"""
====================================================
STAR IDENTIFICATION — All Star Identification Algorithms
====================================================

Purpose:
    Determines which detected centroids correspond to which catalog stars.

    Contains:
        - StarIdAlgorithm: base class
        - GeometricVotingStarId: voting-based identification
        - PyramidStarId: 4-star pattern matching with probabilistic validation
        - IdentifyRemainingStars: extension after initial identification
        - Internal helper classes for pair-distance matching

    Also contains shared utilities:
        - IRUnidentifiedCentroid: tracks unidentified centroids
        - IdentifyThirdStar: finds candidate catalog stars

Geometric Voting:
    For each star in the image, votes for catalog stars whose pairwise
    distances match. The catalog star with the most votes wins.
    Verification phase uses inter-star distance consistency.

Pyramid:
    Enumerates groups of 4 stars, finds matching 4-star patterns in the
    catalog using pair-distance database queries. Uses a probabilistic
    model to skip likely-false matches. When a unique pyramid is found,
    extends identification to remaining stars.

Reference:
    Original C++ source: src/star-id.hpp, src/star-id.cpp,
    src/star-id-private.hpp
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from camera import Camera
from database import (
    Catalog, CatalogStar, Star, Stars, StarIdentifier, StarIdentifiers,
    PairDistanceKVectorDatabase, MultiDatabase,
)
from math_utils import Vec2, Vec3, angle, angle_unit, deg_to_rad


# ─────────────────────────────────────────
# UTILITY FUNCTIONS
# ─────────────────────────────────────────


def mag_to_brightness(magnitude: int) -> float:
    """Convert catalog magnitude to relative brightness.

    Formula: 10^(-mag / 250)
    Magnitude input is magnitude * 100 (e.g., 150 → 1.50).

    Reference:
        Original C++ source: src/star-utils.cpp:147-149
    """
    return math.pow(10.0, -magnitude / 250.0)


# ─────────────────────────────────────────
# BASE CLASS
# ─────────────────────────────────────────


class StarIdAlgorithm:
    """Base class for star identification algorithms.

    Takes a list of centroids, a catalog, and a camera, and returns
    a list of StarIdentifier objects mapping centroids to catalog stars.

    Reference:
        Original C++ source: src/star-id.hpp:12-23
    """

    def go(self, database: object,
           stars: Stars, catalog: Catalog,
           camera: Camera) -> StarIdentifiers:
        raise NotImplementedError


# ─────────────────────────────────────────
# DUMMY STAR ID (debug)
# ─────────────────────────────────────────


class DummyStarId(StarIdAlgorithm):
    """Returns random star identifications. For debugging.

    Reference:
        Original C++ source: src/star-id.cpp:16-27
    """

    def go(self, database: object,
           stars: Stars, catalog: Catalog,
           camera: Camera) -> StarIdentifiers:
        import random
        rng = random.Random(123456)
        result = []
        for i in range(len(stars)):
            result.append(StarIdentifier(i, rng.randint(0, len(catalog) - 1)))
        return result


# ─────────────────────────────────────────
# GEOMETRIC VOTING STAR ID
# ─────────────────────────────────────────


class GeometricVotingStarId(StarIdAlgorithm):
    """Geometric Voting star identification.

    Algorithm:
        1. For each star i in image:
           a. For each other star j, compute great-circle distance θ_ij.
           b. Query pair-distance DB for catalog pairs within tolerance.
           c. Vote for both catalog stars in each matching pair.
           d. Assign the catalog star with the most votes.
        2. Verification: for each pair of identified stars, check
           inter-star distance consistency. Keep only stars with
           verification votes ≥ 75% of the maximum.

    Complexity: O(N² × M) where N = image stars, M = pairs per query.

    Reference:
        Original C++ source: src/star-id.cpp:29-152
    """

    def __init__(self, tolerance_rad: float):
        self.tolerance = tolerance_rad

    def go(self, database: object,
           stars: Stars, catalog: Catalog,
           camera: Camera) -> StarIdentifiers:

        pair_db = self._load_pair_db(database)
        if pair_db is None:
            print("GeometricVoting: No pair-distance database available.", file=sys.stderr)
            return []

        identified: StarIdentifiers = []

        for i in range(len(stars)):
            votes = [0] * len(catalog)
            i_spatial = camera.camera_to_spatial(Vec2(stars[i].x, stars[i].y)).normalize()

            for j in range(len(stars)):
                if i == j:
                    continue

                j_spatial = camera.camera_to_spatial(Vec2(stars[j].x, stars[j].y)).normalize()
                gc_distance = angle_unit(i_spatial, j_spatial)
                lower = gc_distance - self.tolerance
                upper = gc_distance + self.tolerance

                pairs = pair_db.find_pairs_liberal(lower, upper)
                voted_in_pair = [False] * len(catalog)

                for p in pairs:
                    if not voted_in_pair[p.index1]:
                        votes[p.index1] += 1
                        voted_in_pair[p.index1] = True
                    if not voted_in_pair[p.index2]:
                        votes[p.index2] += 1
                        voted_in_pair[p.index2] = True

            max_votes = max(votes) if votes else 0
            index_of_max = votes.index(max_votes) if max_votes > 0 else 0
            identified.append(StarIdentifier(i, index_of_max))

        # Verification phase
        verification_votes = [0] * len(identified)
        for i in range(len(identified)):
            for j in range(i + 1, len(identified)):
                first = catalog[identified[i].catalog_index]
                second = catalog[identified[j].catalog_index]
                c_dist = angle_unit(first.spatial, second.spatial)

                first_star = stars[identified[i].star_index]
                second_star = stars[identified[j].star_index]
                first_spatial = camera.camera_to_spatial(Vec2(first_star.x, first_star.y))
                second_spatial = camera.camera_to_spatial(Vec2(second_star.x, second_star.y))
                s_dist = angle_unit(first_spatial.normalize(), second_spatial.normalize())

                if abs(s_dist - c_dist) < self.tolerance:
                    verification_votes[i] += 1
                    verification_votes[j] += 1

        max_votes = max(verification_votes) if verification_votes else 0
        threshold = max_votes * 3 // 4

        verified: StarIdentifiers = []
        for i, v in enumerate(verification_votes):
            if v > threshold:
                verified.append(identified[i])

        return verified

    @staticmethod
    def _load_pair_db(database: object) -> Optional[PairDistanceKVectorDatabase]:
        if isinstance(database, PairDistanceKVectorDatabase):
            return database
        return None


# ─────────────────────────────────────────
# PYRAMID STAR ID — INTERNAL HELPERS
# ─────────────────────────────────────────


class _IRUnidentifiedCentroid:
    """Tracks an unidentified centroid during IdentifyRemainingStars.

    Stores the best pair of identified stars that form a triangle
    with this centroid, where one angle is close to 90°.

    The "triangular angle" heuristic: when two identified stars and
    an unidentified star form a triangle, the angle at the unidentified
    star is most discriminating when it's near 90°.

    Reference:
        Original C++ source: src/star-id-private.hpp:17-47
        src/star-id.cpp:268-371
    """

    def __init__(self, star: Star, index: int):
        self.best_angle_from_90 = float('inf')
        self.best_star1: Optional[StarIdentifier] = None
        self.best_star2: Optional[StarIdentifier] = None
        self.index = index
        self.star = star
        self._identified_stars_in_range: List[Tuple[float, StarIdentifier]] = []

    @staticmethod
    def _vertical_angles_to_angle_from_90(v1: float, v2: float) -> float:
        """Compute how far the difference of two angles is from 90°.

        Reference:
            Original C++ source: src/star-id.cpp:268-270
        """
        diff = (v1 - v2) % math.pi
        return abs(diff - math.pi / 2)

    def add_identified_star(self, star_id: StarIdentifier, stars: Stars) -> None:
        """When a centroid within range is identified, register it.

        Checks if this new identified star, combined with any previously
        registered star, forms a triangle angle close to 90°.

        Reference:
            Original C++ source: src/star-id.cpp:276-291
        """
        other_star = stars[star_id.star_index]
        dx = other_star.x - self.star.x
        dy = other_star.y - self.star.y
        angle_from_vertical = math.atan2(dy, dx)

        for other_angle, other_pair in self._identified_stars_in_range:
            cur_angle_from_90 = self._vertical_angles_to_angle_from_90(
                other_angle, angle_from_vertical)
            if cur_angle_from_90 < self.best_angle_from_90:
                self.best_angle_from_90 = cur_angle_from_90
                self.best_star1 = star_id
                self.best_star2 = other_pair

        self._identified_stars_in_range.append((angle_from_vertical, star_id))


class _PairDistanceInvolvingIterator:
    """Iterates over pairs that involve a specific catalog star.

    Given the result of a K-vector query (a list of pairs), this
    iterator returns only those pairs where one of the two stars
    matches the given catalog index.

    Reference:
        Original C++ source: src/star-id.cpp:169-250
    """

    def __init__(self, pairs: List, involving: int):
        self._pairs = pairs
        self._involving = involving
        self._pos = 0

    def __iter__(self):
        return self

    def __next__(self) -> int:
        while self._pos < len(self._pairs):
            p = self._pairs[self._pos]
            self._pos += 1
            if p.index1 == self._involving:
                return p.index2
            if p.index2 == self._involving:
                return p.index1
        raise StopIteration

    def has_value(self) -> bool:
        return self._pos < len(self._pairs)


def _pairs_to_map(pairs: List) -> Dict[int, List[int]]:
    """Convert a list of pairs to a multimap from each star to its partners.

    Result is symmetric: if B is in map[A], then A is also in map[B].

    Reference:
        Original C++ source: src/star-id.cpp:259-266
    """
    result: Dict[int, List[int]] = defaultdict(list)
    for p in pairs:
        result[p.index1].append(p.index2)
        result[p.index2].append(p.index1)
    return result


def _find_unidentified_centroids_in_range(
    centroids: List[_IRUnidentifiedCentroid],
    star: Star, camera: Camera,
    min_distance: float, max_distance: float,
) -> List[int]:
    """Find indices of unidentified centroids within distance bounds of a star.

    Uses cosine comparison to avoid computing acos for every candidate.

    Reference:
        Original C++ source: src/star-id.cpp:299-340
    """
    our_spatial = camera.camera_to_spatial(Vec2(star.x, star.y)).normalize()
    min_cos = math.cos(max_distance)
    max_cos = math.cos(min_distance)

    result = []
    for idx, uc in enumerate(centroids):
        their_spatial = camera.camera_to_spatial(Vec2(uc.star.x, uc.star.y)).normalize()
        angle_cos = our_spatial.dot(their_spatial)
        if min_cos <= angle_cos <= max_cos:
            result.append(idx)
    return result


def _add_to_all_unidentified_centroids(
    star_id: StarIdentifier, stars: Stars,
    above_threshold: List[_IRUnidentifiedCentroid],
    below_threshold: List[_IRUnidentifiedCentroid],
    min_distance: float, max_distance: float,
    angle_from_90_threshold: float,
    camera: Camera,
) -> None:
    """Add an identified star to all nearby unidentified centroids.

    Moves centroids whose best angle drops below the threshold from
    the "above threshold" list to the "below threshold" list.

    Reference:
        Original C++ source: src/star-id.cpp:350-371
    """
    now_below = []
    for idx in _find_unidentified_centroids_in_range(
        above_threshold, stars[star_id.star_index], camera,
        min_distance, max_distance,
    ):
        uc = above_threshold[idx]
        uc.add_identified_star(star_id, stars)
        if uc.best_angle_from_90 <= angle_from_90_threshold:
            below_threshold.append(uc)
            now_below.append(uc.index)

    above_threshold[:] = [uc for uc in above_threshold
                          if uc.index not in now_below]


def _identify_third_star(
    db: PairDistanceKVectorDatabase,
    catalog: Catalog,
    catalog_index1: int, catalog_index2: int,
    distance1: float, distance2: float,
    tolerance: float,
) -> List[int]:
    """Find catalog stars matching distances to two identified stars.

    Given two identified catalog stars and distances from each to an
    unidentified centroid, find all catalog stars that could be the
    unidentified one.

    Checks:
        1. Distance to first star matches within tolerance.
        2. Distance to second star matches within tolerance.
        3. Spectrality: (c1 × c2) · candidate > 0 (correct handedness).

    Reference:
        Original C++ source: src/star-id.cpp:384-427
    """
    query1 = db.find_pairs_exact(
        distance1 - tolerance, distance1 + tolerance)

    spatial1 = catalog[catalog_index1].spatial
    spatial2 = catalog[catalog_index2].spatial
    cross = spatial1.cross(spatial2)

    result = []
    for candidate in _PairDistanceInvolvingIterator(query1, catalog_index1):
        candidate_spatial = catalog[candidate].spatial

        angle2 = angle_unit(candidate_spatial, spatial2)
        if not (distance2 - tolerance <= angle2 <= distance2 + tolerance):
            continue

        spectral_torch = cross.dot(candidate_spatial)
        if spectral_torch <= 0:
            continue

        result.append(candidate)

    return result


def _select_next_unidentified(
    above_threshold: List[_IRUnidentifiedCentroid],
    below_threshold: List[_IRUnidentifiedCentroid],
) -> Optional[_IRUnidentifiedCentroid]:
    """Select the next centroid to attempt identification on.

    Prefers centroids already below the soft threshold (best angle
    from 90 < π/4). Otherwise picks the one with best angle from 90
    in the above-threshold list.

    Reference:
        Original C++ source: src/star-id.cpp:429-451
    """
    # Process by index (brightness) order for deterministic behavior.
    # Using index ensures extra faint stars never jump ahead of brighter
    # ones in the identification queue, keeping results stable regardless
    # of how many dim stars are added via centroid_mag_filter.
    if below_threshold:
        best = min(below_threshold, key=lambda uc: uc.index)
        below_threshold.remove(best)
        return best

    if not above_threshold:
        return None

    # Only select above-threshold centroids that have valid reference pairs.
    # Centroids may have best_star1 = None if ALL identified stars are outside
    # their distance range (e.g. with wider FOV cameras). We must not attempt
    # identification on them since no valid triangle can be formed.
    valid = [uc for uc in above_threshold if uc.best_star1 is not None]
    if not valid:
        return None

    best = min(valid, key=lambda uc: uc.index)
    above_threshold.remove(best)
    return best


_ANGLE_FROM_90_SOFT_THRESHOLD = math.pi / 4
PYRAMID_LIMIT = 10  # Cap pyramid search to brightest N stars for deterministic behavior


def _identify_remaining_stars(
    identifiers: StarIdentifiers,
    stars: Stars,
    db: PairDistanceKVectorDatabase,
    catalog: Catalog,
    camera: Camera,
    tolerance: float,
) -> int:
    """Identify additional stars after initial pyramid match.

    Only processes stars within the first PYRAMID_LIMIT entries of the
    input list. This guarantees that the identification order is invariant
    when extra faint stars are appended at the end (e.g. by lowering
    centroid_mag_filter), preventing instability.

    Algorithm:
        1. Initialize all unidentified centroids.
        2. For each identified star, update nearby unidentified centroids.
        3. While there are candidates:
           a. Pick the centroid with best angle-from-90.
           b. Use pair distances to find unique catalog match.
           c. If exactly one candidate, identify it and update neighbors.

    Reference:
        Original C++ source: src/star-id.cpp:461-569
    """
    all_unidentified = [
        _IRUnidentifiedCentroid(stars[i], i) for i in range(len(stars))
    ]
    identified_indices = {si.star_index for si in identifiers}

    above_threshold = [
        uc for uc in all_unidentified
        if uc.index not in identified_indices
    ]
    below_threshold: List[_IRUnidentifiedCentroid] = []

    for star_id in identifiers:
        _add_to_all_unidentified_centroids(
            star_id, stars,
            above_threshold, below_threshold,
            db.min_distance, db.max_distance,
            _ANGLE_FROM_90_SOFT_THRESHOLD,
            camera,
        )

    num_extra = 0

    while above_threshold or below_threshold:
        next_uc = _select_next_unidentified(above_threshold, below_threshold)
        if next_uc is None:
            break

        unidentified_spatial = camera.camera_to_spatial(
            Vec2(next_uc.star.x, next_uc.star.y))
        spatial1 = camera.camera_to_spatial(
            Vec2(stars[next_uc.best_star1.star_index].x,
                 stars[next_uc.best_star1.star_index].y))
        spatial2 = camera.camera_to_spatial(
            Vec2(stars[next_uc.best_star2.star_index].x,
                 stars[next_uc.best_star2.star_index].y))
        d1 = angle_unit(spatial1.normalize(), unidentified_spatial.normalize())
        d2 = angle_unit(spatial2.normalize(), unidentified_spatial.normalize())
        spectral_torch = spatial1.cross(spatial2).dot(unidentified_spatial)

        if spectral_torch > 0:
            candidates = _identify_third_star(
                db, catalog,
                next_uc.best_star1.catalog_index,
                next_uc.best_star2.catalog_index,
                d1, d2, tolerance)
        else:
            candidates = _identify_third_star(
                db, catalog,
                next_uc.best_star2.catalog_index,
                next_uc.best_star1.catalog_index,
                d2, d1, tolerance)

        if len(candidates) == 1:
            identifiers.append(StarIdentifier(next_uc.index, candidates[0]))

            _add_to_all_unidentified_centroids(
                identifiers[-1], stars,
                above_threshold, below_threshold,
                db.min_distance, db.max_distance,
                _ANGLE_FROM_90_SOFT_THRESHOLD,
                camera,
            )

            num_extra += 1
        else:
            pass  # No candidates or multiple — can't identify

    return num_extra


# ─────────────────────────────────────────
# PYRAMID STAR ID
# ─────────────────────────────────────────


class PyramidStarId(StarIdAlgorithm):
    """Pyramid star identification algorithm.

    The "de facto" star-id algorithm used in many real-world missions.
    Searches for groups of 4 stars (a "pyramid") that uniquely match
    a 4-star pattern in the catalog. Once found, extends identification
    to remaining stars.

    Algorithm:
        1. Enumerate 4-star tuples (i, j, k, r) using interleaved indices.
        2. For each tuple:
           a. Compute 6 inter-star distances.
           b. Compute 3 inner sine angles at each vertex.
           c. Calculate expected mismatches:
              E[m] = N_false⁴ · tol⁵ · sin(θ_ij) / (2π² · sin(min_angle))
           d. If E[m] > max_prob, skip (likely false match).
           e. Query pair-distance DB for three edges (ij, ik, ir).
           f. Find unique (i, j, k, r) catalog match via hash-multimap intersection.
           g. Check spectrality (handedness) of the match.
        3. If unique pyramid found, run IdentifyRemainingStars.

    Reference:
        Original C++ source: src/star-id.cpp:571-770

    The pyramid paper:
        "Pyramid Star Identification Technique" by Mortari et al.
    """

    def __init__(self, tolerance_rad: float,
                 num_false_stars: int = 500,
                 max_mismatch_prob: float = 0.001,
                 cutoff: int = 1000):
        self.tolerance = tolerance_rad
        self.num_false_stars = num_false_stars
        self.max_mismatch_prob = max_mismatch_prob
        self.cutoff = cutoff

    def go(self, database: object,
           stars: Stars, catalog: Catalog,
           camera: Camera) -> StarIdentifiers:

        if database is None or len(stars) < 4:
            print("Not enough stars, or database missing.", file=sys.stderr)
            return []

        pair_db = self._load_pair_db(database, catalog)
        if pair_db is None:
            return []

        expected_mismatches_const = (
            math.pow(self.num_false_stars, 4)
            * math.pow(self.tolerance, 5)
            / 2
            / math.pow(math.pi, 2)
        )

        identified: StarIdentifiers = []

        # Cap the pyramid search to PYRAMID_LIMIT stars so that the iteration
        # pattern (dj, dk, dr, i) is invariant to adding faint stars. The first
        # PYRAMID_LIMIT stars are the same across varying mag_filter settings
        # (they are the stars that pass the brightness filter). The pipeline
        # also caps at PYRAMID_LIMIT before calling this function, so
        # IdentifyRemainingStars never sees stars beyond the limit either,
        # preventing identification-order instability from extra stars.
        sorted_stars = stars
        orig_index = {i: i for i in range(len(stars))}
        num_stars = min(len(stars), PYRAMID_LIMIT)
        across = int(math.sqrt(num_stars)) * 2
        halfway_across = int(math.sqrt(num_stars) / 2)
        total_iterations = 0

        j_max = num_stars - 3
        for j_iter in range(j_max):
            dj = 1 + (j_iter + halfway_across) % j_max

            k_max = num_stars - dj - 2
            for k_iter in range(k_max):
                dk = 1 + (k_iter + across) % k_max

                r_max = num_stars - dj - dk - 1
                for r_iter in range(r_max):
                    dr = 1 + (r_iter + halfway_across) % r_max

                    i_max = num_stars - dj - dk - dr - 1
                    for i_iter in range(i_max + 1):
                        i = (i_iter + i_max // 2) % (i_max + 1)

                        total_iterations += 1
                        if total_iterations > self.cutoff:
                            return identified

                        j = i + dj
                        k = j + dk
                        r = k + dr

                        i_spatial = camera.camera_to_spatial(
                            Vec2(sorted_stars[i].x, sorted_stars[i].y)).normalize()
                        j_spatial = camera.camera_to_spatial(
                            Vec2(sorted_stars[j].x, sorted_stars[j].y)).normalize()
                        k_spatial = camera.camera_to_spatial(
                            Vec2(sorted_stars[k].x, sorted_stars[k].y)).normalize()

                        ij_dist = angle_unit(i_spatial, j_spatial)

                        # Sine inner angles at each vertex
                        # Use angle() not angle_unit() since (j-i) etc. are NOT unit vectors
                        i_sin_inner = math.sin(angle(
                            j_spatial - i_spatial, k_spatial - i_spatial))
                        j_sin_inner = math.sin(angle(
                            i_spatial - j_spatial, k_spatial - j_spatial))
                        k_sin_inner = math.sin(angle(
                            i_spatial - k_spatial, j_spatial - k_spatial))

                        expected_mismatches = (
                            expected_mismatches_const
                            * math.sin(ij_dist)
                            / k_sin_inner
                            / max(i_sin_inner, j_sin_inner, k_sin_inner)
                        )

                        if expected_mismatches > self.max_mismatch_prob:
                            continue

                        r_spatial = camera.camera_to_spatial(
                            Vec2(sorted_stars[r].x, sorted_stars[r].y)).normalize()

                        spectral_torch = i_spatial.cross(j_spatial).dot(k_spatial) > 0

                        ik_dist = angle_unit(i_spatial, k_spatial)
                        ir_dist = angle_unit(i_spatial, r_spatial)
                        jk_dist = angle_unit(j_spatial, k_spatial)
                        jr_dist = angle_unit(j_spatial, r_spatial)
                        kr_dist = angle_unit(k_spatial, r_spatial)

                        # Check all distances are within database bounds
                        for d in [ik_dist, ir_dist, jk_dist, jr_dist, kr_dist]:
                            if (d < pair_db.min_distance + self.tolerance
                                    or d > pair_db.max_distance - self.tolerance):
                                break
                        else:
                            # All distances within bounds — proceed
                            pass
                        # Actually we need to handle the continue properly:
                        distances_ok = True
                        for d in [ik_dist, ir_dist, jk_dist, jr_dist, kr_dist]:
                            if (d < pair_db.min_distance + self.tolerance
                                    or d > pair_db.max_distance - self.tolerance):
                                distances_ok = False
                                break
                        if not distances_ok:
                            continue

                        ij_pairs = pair_db.find_pairs_liberal(
                            ij_dist - self.tolerance, ij_dist + self.tolerance)
                        ik_pairs = pair_db.find_pairs_liberal(
                            ik_dist - self.tolerance, ik_dist + self.tolerance)
                        ir_pairs = pair_db.find_pairs_liberal(
                            ir_dist - self.tolerance, ir_dist + self.tolerance)

                        ik_map = _pairs_to_map(ik_pairs)
                        ir_map = _pairs_to_map(ir_pairs)

                        i_match = j_match = k_match = r_match = -1

                        for p in ij_pairs:
                            # C++ iterates the flat pair array [a,b,a,b,...] element by element,
                            # trying both (a,b) and (b,a) orderings for each pair.
                            for (i_candidate, j_candidate) in [(p.index1, p.index2),
                                                                (p.index2, p.index1)]:
                                i_candidate_spatial = catalog[i_candidate].spatial
                                j_candidate_spatial = catalog[j_candidate].spatial
                                ij_candidate_cross = i_candidate_spatial.cross(j_candidate_spatial)

                                for k_candidate in ik_map.get(i_candidate, []):
                                    k_candidate_spatial = catalog[k_candidate].spatial
                                    candidate_spectral = ij_candidate_cross.dot(k_candidate_spatial) > 0
                                    if candidate_spectral != spectral_torch:
                                        continue

                                    jk_candidate_dist = angle_unit(
                                        j_candidate_spatial, k_candidate_spatial)
                                    if (jk_candidate_dist < jk_dist - self.tolerance
                                            or jk_candidate_dist > jk_dist + self.tolerance):
                                        continue

                                    for r_candidate in ir_map.get(i_candidate, []):
                                        r_candidate_spatial = catalog[r_candidate].spatial
                                        jr_candidate_dist = angle_unit(
                                            j_candidate_spatial, r_candidate_spatial)
                                        if (jr_candidate_dist < jr_dist - self.tolerance
                                                or jr_candidate_dist > jr_dist + self.tolerance):
                                            continue

                                        kr_candidate_dist = angle_unit(
                                            k_candidate_spatial, r_candidate_spatial)
                                        if (kr_candidate_dist < kr_dist - self.tolerance
                                                or kr_candidate_dist > kr_dist + self.tolerance):
                                            continue

                                        if i_match == -1:
                                            i_match = i_candidate
                                            j_match = j_candidate
                                            k_match = k_candidate
                                            r_match = r_candidate
                                        else:
                                            i_match = -2
                                            break

                                    if i_match == -2:
                                        break

                                if i_match == -2:
                                    break

                        if i_match >= 0:
                            identified.append(StarIdentifier(i, i_match))
                            identified.append(StarIdentifier(j, j_match))
                            identified.append(StarIdentifier(k, k_match))
                            identified.append(StarIdentifier(r, r_match))

                            num_extra = _identify_remaining_stars(
                                identified, sorted_stars, pair_db, catalog,
                                camera, self.tolerance)

                            # Remap star indices from sorted order back to original
                            for si in identified:
                                si.star_index = orig_index[si.star_index]

                            return identified

        return identified

    @staticmethod
    def _load_pair_db(database: object, catalog: Catalog) -> Optional[PairDistanceKVectorDatabase]:
        if isinstance(database, PairDistanceKVectorDatabase):
            return database
        return None
