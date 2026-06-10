"""
====================================================
MATH UTILITIES — Vectors, Matrices, Quaternions, Attitude
====================================================

Purpose:
    Provides all mathematical primitives used throughout LOST.
    Ported from C++ attitude-utils.hpp and attitude-utils.cpp.

    Contains:
        - Vec2: 2D vector (pixel coordinates)
        - Vec3: 3D vector (spatial coordinates)
        - Mat3: 3x3 matrix (rotation/DCM)
        - Quaternion: rotation representation
        - EulerAngles: RA/Dec/Roll representation
        - Attitude: combined attitude with dual storage
        - Coordinate transforms (spherical <-> spatial)
        - Unit conversions (rad <-> deg <-> arcsec)
        - Angle computation functions

Mathematical Background:
    All rotations follow the z-y'-x'' Euler angle convention
    (intrinsic rotations). Right ascension is rotation about Z,
    declination is about Y', roll is about X''.

Reference:
    Original C++ source: src/attitude-utils.hpp, src/attitude-utils.cpp
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple


# ─────────────────────────────────────────
# 2D VECTOR
# ─────────────────────────────────────────


class Vec2:
    """A 2D vector with floating-point components.

    Used for pixel coordinates on the camera sensor.
    Origin is top-left corner.
    """

    __slots__ = ('x', 'y')

    def __init__(self, x: float, y: float):
        self.x = float(x)
        self.y = float(y)

    def magnitude_sq(self) -> float:
        return self.x * self.x + self.y * self.y

    def magnitude(self) -> float:
        return math.hypot(self.x, self.y)

    def normalize(self) -> Vec2:
        m = self.magnitude()
        if m == 0:
            return Vec2(0, 0)
        return Vec2(self.x / m, self.y / m)

    def dot(self, other: Vec2) -> float:
        return self.x * other.x + self.y * other.y

    def __add__(self, other: Vec2) -> Vec2:
        return Vec2(self.x + other.x, self.y + other.y)

    def __sub__(self, other: Vec2) -> Vec2:
        return Vec2(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: float) -> Vec2:
        return Vec2(self.x * scalar, self.y * scalar)

    def __rmul__(self, scalar: float) -> Vec2:
        return Vec2(self.x * scalar, self.y * scalar)

    def __neg__(self) -> Vec2:
        return Vec2(-self.x, -self.y)

    def __repr__(self) -> str:
        return f"Vec2({self.x:.6f}, {self.y:.6f})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Vec2):
            return NotImplemented
        return self.x == other.x and self.y == other.y


# ─────────────────────────────────────────
# 3D VECTOR
# ─────────────────────────────────────────


class Vec3:
    """A 3D vector with floating-point components.

    Used for spatial coordinates on the unit sphere.
    For a star at (RA=0, Dec=0), the spatial vector is (1, 0, 0).
    X is the boresight direction of the camera.
    """

    __slots__ = ('x', 'y', 'z')

    def __init__(self, x: float, y: float, z: float):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def magnitude_sq(self) -> float:
        return self.x * self.x + self.y * self.y + self.z * self.z

    def magnitude(self) -> float:
        return math.hypot(math.hypot(self.x, self.y), self.z)

    def normalize(self) -> Vec3:
        m = self.magnitude()
        if m == 0:
            return Vec3(0, 0, 0)
        return Vec3(self.x / m, self.y / m, self.z / m)

    def dot(self, other: Vec3) -> float:
        return math.fma(self.x, other.x,
                        math.fma(self.y, other.y, self.z * other.z))

    def cross(self, other: Vec3) -> Vec3:
        return Vec3(
            self.y * other.z - self.z * other.y,
            -(self.x * other.z - self.z * other.x),
            self.x * other.y - self.y * other.x,
        )

    def outer_product(self, other: Vec3) -> Mat3:
        return Mat3(
            self.x * other.x, self.x * other.y, self.x * other.z,
            self.y * other.x, self.y * other.y, self.y * other.z,
            self.z * other.x, self.z * other.y, self.z * other.z,
        )

    def __add__(self, other: Vec3) -> Vec3:
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Vec3) -> Vec3:
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: float) -> Vec3:
        return Vec3(self.x * scalar, self.y * scalar, self.z * scalar)

    def __rmul__(self, scalar: float) -> Vec3:
        return Vec3(self.x * scalar, self.y * scalar, self.z * scalar)

    def __neg__(self) -> Vec3:
        return Vec3(-self.x, -self.y, -self.z)

    def __repr__(self) -> str:
        return f"Vec3({self.x:.6f}, {self.y:.6f}, {self.z:.6f})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Vec3):
            return NotImplemented
        return self.x == other.x and self.y == other.y and self.z == other.z


# ─────────────────────────────────────────
# 3x3 MATRIX
# ─────────────────────────────────────────


class Mat3:
    """A 3x3 matrix stored as 9 floats in row-major order.

    Used for rotation matrices (Direction Cosine Matrices).
    """

    __slots__ = ('x',)

    def __init__(self,
                 m00: float, m01: float, m02: float,
                 m10: float, m11: float, m12: float,
                 m20: float, m21: float, m22: float):
        self.x = (float(m00), float(m01), float(m02),
                  float(m10), float(m11), float(m12),
                  float(m20), float(m21), float(m22))

    def at(self, i: int, j: int) -> float:
        return self.x[3 * i + j]

    def row(self, i: int) -> Vec3:
        return Vec3(self.at(i, 0), self.at(i, 1), self.at(i, 2))

    def col(self, j: int) -> Vec3:
        return Vec3(self.at(0, j), self.at(1, j), self.at(2, j))

    def transpose(self) -> Mat3:
        return Mat3(
            self.at(0, 0), self.at(1, 0), self.at(2, 0),
            self.at(0, 1), self.at(1, 1), self.at(2, 1),
            self.at(0, 2), self.at(1, 2), self.at(2, 2),
        )

    def trace(self) -> float:
        return self.at(0, 0) + self.at(1, 1) + self.at(2, 2)

    def det(self) -> float:
        return (self.at(0, 0) * (self.at(1, 1) * self.at(2, 2) - self.at(2, 1) * self.at(1, 2))
                - self.at(0, 1) * (self.at(1, 0) * self.at(2, 2) - self.at(2, 0) * self.at(1, 2))
                + self.at(0, 2) * (self.at(1, 0) * self.at(2, 1) - self.at(2, 0) * self.at(1, 1)))

    def inverse(self) -> Mat3:
        d = 1.0 / self.det()
        return Mat3(
            (self.at(1, 1) * self.at(2, 2) - self.at(1, 2) * self.at(2, 1)) * d,
            (self.at(0, 2) * self.at(2, 1) - self.at(0, 1) * self.at(2, 2)) * d,
            (self.at(0, 1) * self.at(1, 2) - self.at(0, 2) * self.at(1, 1)) * d,
            (self.at(1, 2) * self.at(2, 0) - self.at(1, 0) * self.at(2, 2)) * d,
            (self.at(0, 0) * self.at(2, 2) - self.at(0, 2) * self.at(2, 0)) * d,
            (self.at(0, 2) * self.at(1, 0) - self.at(0, 0) * self.at(1, 2)) * d,
            (self.at(1, 0) * self.at(2, 1) - self.at(1, 1) * self.at(2, 0)) * d,
            (self.at(0, 1) * self.at(2, 0) - self.at(0, 0) * self.at(2, 1)) * d,
            (self.at(0, 0) * self.at(1, 1) - self.at(0, 1) * self.at(1, 0)) * d,
        )

    def __add__(self, other: Mat3) -> Mat3:
        return Mat3(
            self.at(0, 0) + other.at(0, 0), self.at(0, 1) + other.at(0, 1), self.at(0, 2) + other.at(0, 2),
            self.at(1, 0) + other.at(1, 0), self.at(1, 1) + other.at(1, 1), self.at(1, 2) + other.at(1, 2),
            self.at(2, 0) + other.at(2, 0), self.at(2, 1) + other.at(2, 1), self.at(2, 2) + other.at(2, 2),
        )

    def __sub__(self, other: Mat3) -> Mat3:
        return Mat3(
            self.at(0, 0) - other.at(0, 0), self.at(0, 1) - other.at(0, 1), self.at(0, 2) - other.at(0, 2),
            self.at(1, 0) - other.at(1, 0), self.at(1, 1) - other.at(1, 1), self.at(1, 2) - other.at(1, 2),
            self.at(2, 0) - other.at(2, 0), self.at(2, 1) - other.at(2, 1), self.at(2, 2) - other.at(2, 2),
        )

    def __mul__(self, other: Mat3 | Vec3 | float) -> Mat3 | Vec3:
        if isinstance(other, Mat3):
            return Mat3(
                self.at(0, 0) * other.at(0, 0) + self.at(0, 1) * other.at(1, 0) + self.at(0, 2) * other.at(2, 0),
                self.at(0, 0) * other.at(0, 1) + self.at(0, 1) * other.at(1, 1) + self.at(0, 2) * other.at(2, 1),
                self.at(0, 0) * other.at(0, 2) + self.at(0, 1) * other.at(1, 2) + self.at(0, 2) * other.at(2, 2),
                self.at(1, 0) * other.at(0, 0) + self.at(1, 1) * other.at(1, 0) + self.at(1, 2) * other.at(2, 0),
                self.at(1, 0) * other.at(0, 1) + self.at(1, 1) * other.at(1, 1) + self.at(1, 2) * other.at(2, 1),
                self.at(1, 0) * other.at(0, 2) + self.at(1, 1) * other.at(1, 2) + self.at(1, 2) * other.at(2, 2),
                self.at(2, 0) * other.at(0, 0) + self.at(2, 1) * other.at(1, 0) + self.at(2, 2) * other.at(2, 0),
                self.at(2, 0) * other.at(0, 1) + self.at(2, 1) * other.at(1, 1) + self.at(2, 2) * other.at(2, 1),
                self.at(2, 0) * other.at(0, 2) + self.at(2, 1) * other.at(1, 2) + self.at(2, 2) * other.at(2, 2),
            )
        elif isinstance(other, Vec3):
            return Vec3(
                other.x * self.at(0, 0) + other.y * self.at(0, 1) + other.z * self.at(0, 2),
                other.x * self.at(1, 0) + other.y * self.at(1, 1) + other.z * self.at(1, 2),
                other.x * self.at(2, 0) + other.y * self.at(2, 1) + other.z * self.at(2, 2),
            )
        elif isinstance(other, (int, float)):
            return Mat3(
                self.at(0, 0) * other, self.at(0, 1) * other, self.at(0, 2) * other,
                self.at(1, 0) * other, self.at(1, 1) * other, self.at(1, 2) * other,
                self.at(2, 0) * other, self.at(2, 1) * other, self.at(2, 2) * other,
            )
        return NotImplemented

    def __rmul__(self, other: Vec3 | float) -> Vec3 | Mat3:
        if isinstance(other, Vec3):
            return Vec3(
                other.x * self.at(0, 0) + other.y * self.at(0, 1) + other.z * self.at(0, 2),
                other.x * self.at(1, 0) + other.y * self.at(1, 1) + other.z * self.at(1, 2),
                other.x * self.at(2, 0) + other.y * self.at(2, 1) + other.z * self.at(2, 2),
            )
        elif isinstance(other, (int, float)):
            return Mat3(
                self.at(0, 0) * other, self.at(0, 1) * other, self.at(0, 2) * other,
                self.at(1, 0) * other, self.at(1, 1) * other, self.at(1, 2) * other,
                self.at(2, 0) * other, self.at(2, 1) * other, self.at(2, 2) * other,
            )
        return NotImplemented

    def __repr__(self) -> str:
        return (f"Mat3({self.at(0, 0):.6f}, {self.at(0, 1):.6f}, {self.at(0, 2):.6f},\n"
                f"     {self.at(1, 0):.6f}, {self.at(1, 1):.6f}, {self.at(1, 2):.6f},\n"
                f"     {self.at(2, 0):.6f}, {self.at(2, 1):.6f}, {self.at(2, 2):.6f})")

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Mat3):
            return NotImplemented
        return self.x == other.x


IDENTITY_MAT3 = Mat3(1, 0, 0, 0, 1, 0, 0, 0, 1)


# ─────────────────────────────────────────
# QUATERNION
# ─────────────────────────────────────────


class Quaternion:
    """A unit quaternion representing a 3D rotation.

    Stored as (real, i, j, k) where real is the scalar part
    and (i, j, k) is the vector part.

    Represents rotation by angle θ about axis (x, y, z):
        q = cos(θ/2) + (x*i + y*j + z*k) * sin(θ/2)

    Reference:
        Original C++ source: src/attitude-utils.hpp:98-125
        src/attitude-utils.cpp:14-133
    """

    __slots__ = ('real', 'i', 'j', 'k')

    def __init__(self, real: float, i: float, j: float, k: float):
        self.real = float(real)
        self.i = float(i)
        self.j = float(j)
        self.k = float(k)

    @classmethod
    def from_axis_angle(cls, axis: Vec3, theta: float) -> Quaternion:
        half = theta / 2.0
        s = math.sin(half)
        return cls(math.cos(half), axis.x * s, axis.y * s, axis.z * s)

    @classmethod
    def from_vector(cls, vec: Vec3) -> Quaternion:
        return cls(0, vec.x, vec.y, vec.z)

    def vector(self) -> Vec3:
        return Vec3(self.i, self.j, self.k)

    def set_vector(self, vec: Vec3) -> None:
        self.i = vec.x
        self.j = vec.y
        self.k = vec.z

    def conjugate(self) -> Quaternion:
        return Quaternion(self.real, -self.i, -self.j, -self.k)

    def __mul__(self, other: Quaternion) -> Quaternion:
        return Quaternion(
            self.real * other.real - self.i * other.i - self.j * other.j - self.k * other.k,
            self.real * other.i + other.real * self.i + self.j * other.k - self.k * other.j,
            self.real * other.j + other.real * self.j + self.k * other.i - self.i * other.k,
            self.real * other.k + other.real * self.k + self.i * other.j - self.j * other.i,
        )

    def rotate(self, vec: Vec3) -> Vec3:
        pure = Quaternion.from_vector(vec)
        result = self * pure * self.conjugate()
        return result.vector()

    def angle(self) -> float:
        if self.real <= -1:
            return 0.0
        if self.real >= 1:
            return 0.0
        return math.acos(self.real) * 2

    def smallest_angle(self) -> float:
        raw = self.angle()
        if raw > math.pi:
            return 2 * math.pi - raw
        return raw

    def set_angle(self, new_angle: float) -> None:
        self.real = math.cos(new_angle / 2.0)
        v = self.vector().normalize() * math.sin(new_angle / 2.0)
        self.set_vector(v)

    def is_unit(self, tolerance: float = 1e-5) -> bool:
        return abs(self.i * self.i + self.j * self.j + self.k * self.k + self.real * self.real - 1) < tolerance

    def canonicalize(self) -> Quaternion:
        if self.real >= 0:
            return self
        return Quaternion(-self.real, -self.i, -self.j, -self.k)

    def to_spherical(self) -> EulerAngles:
        ra = math.atan2(2 * (-self.real * self.k + self.i * self.j),
                        1 - 2 * (self.j * self.j + self.k * self.k))
        if ra < 0:
            ra += 2 * math.pi
        de = -math.asin(2 * (-self.real * self.j - self.i * self.k))
        roll = -math.atan2(2 * (-self.real * self.i + self.j * self.k),
                           1 - 2 * (self.i * self.i + self.j * self.j))
        if roll < 0:
            roll += 2 * math.pi
        return EulerAngles(ra, de, roll)

    def __repr__(self) -> str:
        return f"Quaternion({self.real:.6f}, {self.i:.6f}, {self.j:.6f}, {self.k:.6f})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Quaternion):
            return NotImplemented
        return (self.real == other.real and self.i == other.i
                and self.j == other.j and self.k == other.k)


# ─────────────────────────────────────────
# EULER ANGLES
# ─────────────────────────────────────────


@dataclass
class EulerAngles:
    """Euler angles representing attitude in the z-y'-x'' convention.

    Attributes:
        ra: Right ascension (radians). Rotation about Z axis. [0, 2π)
        de: Declination (radians). Rotation about Y' axis. [-π/2, π/2]
        roll: Roll (radians). Rotation about X'' axis. [0, 2π)
    """

    ra: float
    de: float
    roll: float


# ─────────────────────────────────────────
# ATTITUDE
# ─────────────────────────────────────────


class Attitude:
    """Spacecraft attitude (orientation) with dual representation.

    Stores either a quaternion or direction cosine matrix (DCM).
    Converts automatically between formats as needed.
    May be "unknown" if attitude could not be determined.

    Reference:
        Original C++ source: src/attitude-utils.hpp:138-161
        src/attitude-utils.cpp:348-451
    """

    def __init__(self, quat: Quaternion | None = None, dcm: Mat3 | None = None):
        if quat is not None:
            self._type = 'quaternion'
            self._quaternion = quat
            self._dcm = None
        elif dcm is not None:
            self._type = 'dcm'
            self._dcm = dcm
            self._quaternion = None
        else:
            self._type = 'unknown'
            self._quaternion = None
            self._dcm = None

    def is_known(self) -> bool:
        return self._type != 'unknown'

    def get_quaternion(self) -> Quaternion:
        if self._type == 'quaternion':
            return self._quaternion
        elif self._type == 'dcm':
            return dcm_to_quaternion(self._dcm)
        else:
            raise ValueError("Attitude is unknown")

    def get_dcm(self) -> Mat3:
        if self._type == 'dcm':
            return self._dcm
        elif self._type == 'quaternion':
            return quaternion_to_dcm(self._quaternion)
        else:
            raise ValueError("Attitude is unknown")

    def to_spherical(self) -> EulerAngles:
        if self._type == 'quaternion':
            return self._quaternion.to_spherical()
        elif self._type == 'dcm':
            return dcm_to_quaternion(self._dcm).to_spherical()
        else:
            raise ValueError("Attitude is unknown")

    def rotate(self, vec: Vec3) -> Vec3:
        if self._type == 'dcm':
            return self._dcm * vec
        elif self._type == 'quaternion':
            return self._quaternion.rotate(vec)
        else:
            raise ValueError("Attitude is unknown")

    def __repr__(self) -> str:
        if not self.is_known():
            return "Attitude(unknown)"
        return f"Attitude({self.get_quaternion()})"


# ─────────────────────────────────────────
# COORDINATE TRANSFORMS
# ─────────────────────────────────────────


def spherical_to_spatial(ra: float, de: float) -> Vec3:
    """Convert right ascension and declination to a unit vector.

    A star with RA=0, Dec=0 maps to (1, 0, 0).

    Reference:
        Original C++ source: src/attitude-utils.cpp:136-142
    """
    return Vec3(
        math.cos(ra) * math.cos(de),
        math.sin(ra) * math.cos(de),
        math.sin(de),
    )


def spatial_to_spherical(vec: Vec3) -> Tuple[float, float]:
    """Convert a unit vector to right ascension and declination.

    Returns (ra, de) in radians. RA in [0, 2π), Dec in [-π/2, π/2].

    Reference:
        Original C++ source: src/attitude-utils.cpp:145-150
    """
    ra = math.atan2(vec.y, vec.x)
    if ra < 0:
        ra += 2 * math.pi
    de = math.asin(vec.z)
    return ra, de


def spherical_to_quaternion(ra: float, dec: float, roll: float) -> Quaternion:
    """Create a quaternion from Euler angles (z-y'-x'' convention).

    Applies rotations in order: Z (ra), then Y' (-dec), then X'' (-roll).
    This is an "improper" z-y'-x' Euler rotation.

    Reference:
        Original C++ source: src/attitude-utils.cpp:100-117
    """
    a = Quaternion.from_axis_angle(Vec3(0, 0, 1), ra)
    b = Quaternion.from_axis_angle(Vec3(0, 1, 0), -dec)
    c = Quaternion.from_axis_angle(Vec3(1, 0, 0), -roll)
    result = (a * b * c).conjugate()
    return result


# ─────────────────────────────────────────
# QUATERNION <-> DCM CONVERSIONS
# ─────────────────────────────────────────


def quaternion_to_dcm(quat: Quaternion) -> Mat3:
    """Convert a quaternion to a rotation matrix (DCM).

    Reference:
        Original C++ source: src/attitude-utils.cpp:355-364
    """
    x = quat.rotate(Vec3(1, 0, 0))
    y = quat.rotate(Vec3(0, 1, 0))
    z = quat.rotate(Vec3(0, 0, 1))
    return Mat3(
        x.x, y.x, z.x,
        x.y, y.y, z.y,
        x.z, y.z, z.z,
    )


def dcm_to_quaternion(dcm: Mat3) -> Quaternion:
    """Convert a rotation matrix (DCM) to a quaternion.

    Uses a two-step process:
    1. Find quaternion that aligns X-axis
    2. Find quaternion that aligns Y-axis (with handedness check)

    Reference:
        Original C++ source: src/attitude-utils.cpp:367-393
    """
    old_x = Vec3(1, 0, 0)
    new_x = dcm.col(0)
    x_align_axis = old_x.cross(new_x).normalize()
    x_align_angle = angle_unit(old_x, new_x)
    x_align = Quaternion.from_axis_angle(x_align_axis, x_align_angle)

    old_y = x_align.rotate(Vec3(0, 1, 0))
    new_y = dcm.col(1)
    rotate_clockwise = old_y.cross(new_y).dot(new_x) > 0
    y_angle = angle_unit(old_y, new_y)
    if not rotate_clockwise:
        y_angle = -y_angle
    y_align = Quaternion.from_axis_angle(Vec3(1, 0, 0), y_angle)

    return x_align * y_align


# ─────────────────────────────────────────
# ANGLE COMPUTATION
# ─────────────────────────────────────────


def angle(vec1: Vec3, vec2: Vec3) -> float:
    """Angle between two vectors in radians.
    Normalizes both vectors first.

    Reference:
        Original C++ source: src/attitude-utils.cpp:470-472
    """
    return angle_unit(vec1.normalize(), vec2.normalize())


def angle_unit(vec1: Vec3, vec2: Vec3) -> float:
    """Angle between two unit vectors in radians.
    Slightly faster than angle(); assumes vectors are already unit length.

    Reference:
        Original C++ source: src/attitude-utils.cpp:479-483
    """
    d = vec1.dot(vec2)
    if d >= 1:
        return 0.0
    if d <= -1:
        return math.pi - 0.0000001
    return math.acos(d)


def distance_2d(v1: Vec2, v2: Vec2) -> float:
    return (v1 - v2).magnitude()


def distance_3d(v1: Vec3, v2: Vec3) -> float:
    return (v1 - v2).magnitude()


# ─────────────────────────────────────────
# UNIT CONVERSIONS
# ─────────────────────────────────────────


def rad_to_deg(rad: float) -> float:
    return rad * 180.0 / math.pi


def deg_to_rad(deg: float) -> float:
    return deg / 180.0 * math.pi


def rad_to_arcsec(rad: float) -> float:
    return rad_to_deg(rad) * 3600.0


def arcsec_to_rad(arcsec: float) -> float:
    return deg_to_rad(arcsec / 3600.0)


def decimal_modulo(x: float, mod: float) -> float:
    """Mathematical modulo (not remainder).
    Always returns a value in [0, mod).

    Reference:
        Original C++ source: src/attitude-utils.cpp:168-172
    """
    result = x - mod * math.floor(x / mod)
    return result if result >= 0 else result + mod
