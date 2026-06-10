"""
====================================================
ATTITUDE ESTIMATION — All Attitude Estimation Algorithms
====================================================

Purpose:
    Estimates the orientation (attitude) of the camera based on
    identified stars and their known catalog positions.

    Contains:
        - AttitudeEstimator: base class
        - DavenportQAttitude: optimal attitude via eigen-decomposition
        - TriadAttitude: fast 2-star attitude
        - QuestAttitude: fast optimal attitude via Newton-Raphson

Wahba's Problem:
    All algorithms solve Wahba's problem: find the rotation matrix A
    that minimizes:
        J(A) = Σ w_i · |b_i - A·r_i|²
    where b_i are observed star vectors (camera frame) and r_i are
    reference vectors (catalog frame).

Davenport Q-method:
    Converts Wahba's problem to an eigenvalue problem. Builds the
    4×4 K-matrix and finds its largest eigenvalue. The corresponding
    eigenvector is the optimal quaternion. Requires eigen-decomposition.

TRIAD:
    Uses only two stars to build orthonormal frames in both body and
    reference coordinates. The attitude is the rotation between frames.
    Fast but sensitive to noise in either star.

QUEST:
    Uses Newton-Raphson iteration to find the largest eigenvalue of
    the K-matrix without explicit eigen-decomposition. Historically
    important as the de facto standard. Fast and accurate.

Reference:
    Original C++ source: src/attitude-estimators.hpp, src/attitude-estimators.cpp
"""

from __future__ import annotations

import math

import numpy as np

from camera import Camera
from database import Catalog, Stars, StarIdentifiers
from math_utils import (
    IDENTITY_MAT3, Mat3, Vec2, Vec3, Quaternion, Attitude,
    angle_unit,
)


# ─────────────────────────────────────────
# BASE CLASS
# ─────────────────────────────────────────


class AttitudeEstimator:
    """Base class for attitude estimation algorithms.

    Takes identified stars and computes the spacecraft attitude.

    Reference:
        Original C++ source: src/attitude-estimators.hpp:15-25
    """

    def go(self, camera: Camera, stars: Stars, catalog: Catalog,
           star_ids: StarIdentifiers) -> Attitude:
        raise NotImplementedError


# ─────────────────────────────────────────
# DAVENPORT Q-METHOD
# ─────────────────────────────────────────


class DavenportQAttitude(AttitudeEstimator):
    """Davenport Q-method attitude estimation.

    The optimal attitude estimator (in the Wahba sense). Uses all
    identified stars to build the K-matrix and finds its dominant
    eigenvector using numpy's eigen-decomposition.

    Algorithm:
        1. Build attitude profile matrix B = Σ w_i · r_i · b_iᵀ
        2. Compute S = B + Bᵀ, σ = tr(B), Z = [B₂₃-B₃₂, B₃₁-B₁₃, B₁₂-B₂₁]ᵀ
        3. Build 4×4 K-matrix:
           K = [[σ, Zᵀ],
                [Z, S - σI]]
        4. Find largest eigenvalue λ_max and eigenvector q
        5. q is the optimal rotation quaternion

    Complexity: O(N + 4³) where N is the number of identified stars.
    The O(4³) eigen-decomposition is constant-time.

    Reference:
        Original C++ source: src/attitude-estimators.cpp:12-120
    """

    def go(self, camera: Camera, stars: Stars, catalog: Catalog,
           star_ids: StarIdentifiers) -> Attitude:
        if len(star_ids) < 2:
            return Attitude()

        # Build attitude profile matrix B (3x3)
        B = np.zeros((3, 3), dtype=np.float64)

        for s in star_ids:
            b_star = stars[s.star_index]
            b_spatial = camera.camera_to_spatial(Vec2(b_star.x, b_star.y))
            b_vec = np.array([b_spatial.x, b_spatial.y, b_spatial.z])


        for s in star_ids:
            b_star = stars[s.star_index]
            b_spatial = camera.camera_to_spatial(Vec2(b_star.x, b_star.y))
            b_vec = np.array([b_spatial.x, b_spatial.y, b_spatial.z])

            r_star = catalog[s.catalog_index]
            r_vec = np.array([r_star.spatial.x, r_star.spatial.y, r_star.spatial.z])

            B += s.weight * np.outer(r_vec, b_vec)

        # S = B + B^T
        S = B + B.T

        # σ = trace(B)
        sigma = np.trace(B)

        # Z = [B[1,2]-B[2,1], B[2,0]-B[0,2], B[0,1]-B[1,0]]^T
        Z = np.array([
            B[1, 2] - B[2, 1],
            B[2, 0] - B[0, 2],
            B[0, 1] - B[1, 0],
        ])

        # K = 4x4 matrix
        K = np.zeros((4, 4), dtype=np.float64)
        K[0, 0] = sigma
        K[0, 1:4] = Z
        K[1:4, 0] = Z
        K[1:4, 1:4] = S - sigma * np.eye(3)

        # Eigen-decomposition — find largest eigenvalue/vector
        eigenvalues, eigenvectors = np.linalg.eigh(K)
        max_idx = np.argmax(eigenvalues)
        q = eigenvectors[:, max_idx]

        return Attitude(quat=Quaternion(float(q[0]), float(q[1]),
                                        float(q[2]), float(q[3])))


# ─────────────────────────────────────────
# TRIAD
# ─────────────────────────────────────────


def _triad_frame(v1: Vec3, v2: Vec3) -> Mat3:
    """Build an orthonormal coordinate frame from two vectors.

    Frame axes:
        d1 = v1 / |v1|
        d2 = (v1 × v2) / |v1 × v2|
        d3 = d1 × d2

    Returns a matrix whose columns are the frame axes.

    Reference:
        Original C++ source: src/attitude-estimators.cpp:123-132
    """
    d1 = v1.normalize()
    d2 = v1.cross(v2).normalize()
    d3 = d1.cross(d2).normalize()
    return Mat3(
        d1.x, d2.x, d3.x,
        d1.y, d2.y, d3.y,
        d1.z, d2.z, d3.z,
    )


class TriadAttitude(AttitudeEstimator):
    """TRIAD attitude estimation algorithm.

    A fast attitude estimator using only two identified stars.
    Prone to error if either star is noisy. Should not be used
    when more than 2 stars are available and accuracy matters.

    Algorithm:
        1. Select two stars from the identified set.
        2. Build orthonormal frames in both camera (body) and
           catalog (reference) coordinates.
        3. Attitude = body_frame · reference_frameᵀ

    The frame construction ensures orthonormality even with noisy
    observations, at the cost of using only 2 stars' information.

    Reference:
        Original C++ source: src/attitude-estimators.cpp:134-157
    """

    def go(self, camera: Camera, stars: Stars, catalog: Catalog,
           star_ids: StarIdentifiers) -> Attitude:
        if len(star_ids) < 2:
            return Attitude()

        # Pick first star and one near the middle
        a = star_ids[0]
        b = star_ids[len(star_ids) // 2]

        photo_frame = _triad_frame(
            camera.camera_to_spatial(Vec2(stars[a.star_index].x,
                                          stars[a.star_index].y)),
            camera.camera_to_spatial(Vec2(stars[b.star_index].x,
                                          stars[b.star_index].y)),
        )

        catalog_frame = _triad_frame(
            catalog[a.catalog_index].spatial,
            catalog[b.catalog_index].spatial,
        )

        # attitude = photo_frame * catalog_frame^T
        return Attitude(dcm=photo_frame * catalog_frame.transpose())


# ─────────────────────────────────────────
# QUEST
# ─────────────────────────────────────────

_EPSILON = 0.0001  # Newton-Raphson convergence threshold


def _quest_char_poly(x: float, a: float, b: float, c: float,
                     d: float, s: float) -> float:
    """Characteristic polynomial of the QUEST K-matrix.

    p(λ) = (λ² - a)(λ² - b) - cλ + cs - d

    Reference:
        Original C++ source: src/attitude-estimators.cpp:163
    """
    return ((x * x - a) * (x * x - b) - c * x + (c * s) - d)


def _quest_char_poly_prime(x: float, a: float, b: float, c: float) -> float:
    """Derivative of the QUEST characteristic polynomial.

    p'(λ) = 4λ³ - 2(a + b)λ - c

    Reference:
        Original C++ source: src/attitude-estimators.cpp:168
    """
    return 4 * x * x * x - 2 * (a + b) * x - c


def _quest_eigenvalue(guess: float, a: float, b: float, c: float,
                      d: float, s: float) -> float:
    """Find the largest eigenvalue using Newton-Raphson iteration.

    Iterates: λ_{n+1} = λ_n - p(λ_n) / p'(λ_n)

    Reference:
        Original C++ source: src/attitude-estimators.cpp:174-182
    """
    while True:
        h = _quest_char_poly(guess, a, b, c, d, s) / _quest_char_poly_prime(guess, a, b, c)
        guess -= h
        if abs(h) < _EPSILON:
            break
    return guess


class QuestAttitude(AttitudeEstimator):
    """QUEST attitude estimation algorithm.

    A fast alternative to the Davenport Q-method that avoids explicit
    eigen-decomposition. Instead uses Newton-Raphson iteration to find
    the largest eigenvalue, then solves for the optimal quaternion
    analytically.

    Algorithm:
        1. Build attitude profile matrix B (same as DQM).
        2. Compute S, σ, Z (same as DQM).
        3. Compute characteristic polynomial coefficients:
           - δ = det(S), κ = tr(S⁻¹ · δ)
           - a = σ² - κ, b = σ² + |Z|²
           - c = δ + Z·S·Z, d = Z·S²·Z
        4. Find largest eigenvalue via Newton-Raphson.
        5. Solve for optimal quaternion:
           - α = λ² - σ² + κ
           - β = λ - σ
           - γ = (λ + σ)α - δ
           - X = (αI + βS + S²) · Z
           - q = [γ, X] / sqrt(γ² + |X|²)

    Reference:
        Original C++ source: src/attitude-estimators.cpp:184-247

    See also:
        https://ahrs.readthedocs.io/en/latest/filters/quest.html
        https://arc.aiaa.org/doi/pdf/10.2514/1.62549
    """

    def go(self, camera: Camera, stars: Stars, catalog: Catalog,
           star_ids: StarIdentifiers) -> Attitude:
        if len(star_ids) < 2:
            return Attitude()

        # Initial eigenvalue guess = sum of weights
        guess = 0.0

        # Build attitude profile matrix B
        B = np.zeros((3, 3), dtype=np.float64)

        for s in star_ids:
            b_star = stars[s.star_index]
            b_spatial = camera.camera_to_spatial(Vec3(b_star.x, b_star.y, 0))
            b_vec = np.array([b_spatial.x, b_spatial.y, b_spatial.z])

            r_star = catalog[s.catalog_index]
            r_vec = np.array([r_star.spatial.x, r_star.spatial.y, r_star.spatial.z])

            B += s.weight * np.outer(r_vec, b_vec)
            guess += s.weight

        # S = B + B^T
        S = B + B.T
        sigma = np.trace(B)

        Z = np.array([
            B[1, 2] - B[2, 1],
            B[2, 0] - B[0, 2],
            B[0, 1] - B[1, 0],
        ])

        # Characteristic polynomial coefficients
        delta = np.linalg.det(S)
        kappa = np.trace(np.linalg.inv(S) * delta)
        a = sigma * sigma - kappa
        b = sigma * sigma + Z.dot(Z)
        c = delta + Z.dot(S.dot(Z))
        d = Z.dot(S.dot(S.dot(Z)))

        # Newton-Raphson eigenvalue estimation
        eig = _quest_eigenvalue(guess, a, b, c, d, sigma)

        # Solve for optimal quaternion
        alpha = eig * eig - sigma * sigma + kappa
        beta = eig - sigma
        gamma = (eig + sigma) * alpha - delta

        X = (alpha * np.eye(3) + beta * S + S.dot(S)).dot(Z)
        scalar = 1.0 / math.sqrt(gamma * gamma + X.dot(X))
        X = X * scalar
        gamma *= scalar

        return Attitude(quat=Quaternion(gamma, float(X[0]), float(X[1]),
                                        float(X[2])))
