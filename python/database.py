"""
====================================================
DATABASE — Catalog, K-Vector Index, Pair-Distance Database
====================================================

Purpose:
    Provides data structures for loading and querying star catalogs.

    Contains:
        - catalog loading from TSV (BSC format)
        - NarrowCatalog: magnitude/separation/brightest-N filtering
        - CatalogStar, Star, StarIdentifier types
        - KVectorIndex: range-query acceleration structure
        - PairDistanceKVectorDatabase: inter-star distance lookups
        - MagToBrightness conversion

K-Vector Background:
    The K-vector index is a data structure that enables O(1) range queries
    into sorted 1D numerical data. It pre-computes bin boundaries and
    stores cumulative counts, so a query [q_min, q_max] returns the
    index range [bins[bin(q_min)-1], bins[bin(q_max)]].

Reference:
    Original C++ source: src/databases.hpp, src/databases.cpp
                         src/star-utils.hpp, src/star-utils.cpp
"""

from __future__ import annotations

import csv
import math
import os
from typing import List, Optional, Tuple

from math_utils import Vec3, angle_unit, deg_to_rad


# ─────────────────────────────────────────
# TYPES
# ─────────────────────────────────────────


class CatalogStar:
    """A star from the star catalog.

    Attributes:
        spatial: Unit vector position on celestial sphere.
        magnitude: Magnitude * 100 (integer). E.g., magnitude=150 → mag 1.50.
        name: Unique numerical identifier.

    Reference:
        Original C++ source: src/star-utils.hpp:12-43
    """

    __slots__ = ('spatial', 'magnitude', 'name')

    def __init__(self, spatial: Vec3, magnitude: int, name: int):
        self.spatial = spatial
        self.magnitude = magnitude
        self.name = name

    def __repr__(self) -> str:
        return f"CatalogStar(name={self.name}, mag={self.magnitude})"


class Star:
    """A centroid detected in an image.

    Attributes:
        x: X pixel coordinate (top-left origin).
        y: Y pixel coordinate (top-left origin).
        radius_x: Approximate horizontal radius of bright area.
        radius_y: Approximate vertical radius of bright area.
        magnitude: Relative brightness (larger = brighter).

    Reference:
        Original C++ source: src/star-utils.hpp:49-72
    """

    __slots__ = ('x', 'y', 'radius_x', 'radius_y', 'magnitude')

    def __init__(self, x: float, y: float,
                 radius_x: float = 0.0, radius_y: float = 0.0,
                 magnitude: int = 0):
        self.x = float(x)
        self.y = float(y)
        self.radius_x = float(radius_x)
        self.radius_y = float(radius_y)
        self.magnitude = magnitude

    def __repr__(self) -> str:
        return f"Star({self.x:.2f}, {self.y:.2f}, mag={self.magnitude})"


class StarIdentifier:
    """Records that a detected star corresponds to a catalog star.

    Attributes:
        star_index: Index into the Stars array.
        catalog_index: Index into the Catalog array.
        weight: Confidence weight (usually 1.0).

    Reference:
        Original C++ source: src/star-utils.hpp:78-98
    """

    __slots__ = ('star_index', 'catalog_index', 'weight')

    def __init__(self, star_index: int, catalog_index: int, weight: float = 1.0):
        self.star_index = star_index
        self.catalog_index = catalog_index
        self.weight = float(weight)

    def __repr__(self) -> str:
        return f"StarIdentifier(star={self.star_index}, cat={self.catalog_index})"


# Type aliases matching C++ style
Catalog = List[CatalogStar]
Stars = List[Star]
StarIdentifiers = List[StarIdentifier]


# ─────────────────────────────────────────
# CATALOG LOADING
# ─────────────────────────────────────────


def load_catalog(tsv_path: str | None = None) -> Catalog:
    """Load the Bright Star Catalog from a TSV file.

    The TSV format is: raj2000, dej2000, name, weird_char, magnitude_high.magnitude_low
    RA and Dec are in degrees.
    Magnitude is encoded as: magnitudeHigh*100 + (magnitudeHigh < 0 ? -magnitudeLow : magnitudeLow)
    E.g., magnitude 1.50 → magnitudeHigh=1, magnitudeLow=50 → 1*100 + 50 = 150

    If tsv_path is None, checks the LOST_BSC_PATH environment variable,
    then defaults to '../lost/bright-star-catalog.tsv'.

    Reference:
        Original C++ source: src/io.cpp:58-123
    """
    if tsv_path is None:
        tsv_path = os.environ.get('LOST_BSC_PATH', '../lost/bright-star-catalog.tsv')

    catalog: Catalog = []

    with open(tsv_path, 'r') as f:
        reader = csv.reader(f, delimiter='|')
        for row in reader:
            if len(row) < 5:
                continue
            raj2000 = float(row[0])
            dej2000 = float(row[1])
            name = int(row[2])
            mag_str = row[4] if '.' not in row[3] and '.' in row[4] else row[3]
            # The magnitude is stored as two parts: high.low
            # where high is the integer part and low is the fractional part * 100
            mag_parts = row[4].split('.') if '.' in row[4] else row[3].split('.')
            mag_high = int(mag_parts[0])
            mag_low = int(mag_parts[1]) if len(mag_parts) > 1 else 0
            magnitude = mag_high * 100 + (mag_low if mag_high >= 0 else -mag_low)

            catalog.append(CatalogStar(
                spherical_to_spatial(deg_to_rad(raj2000), deg_to_rad(dej2000)),
                magnitude,
                name,
            ))

    # Remove duplicate positions (keep brighter star)
    catalog.sort(key=lambda s: s.spatial.x)
    i = len(catalog) - 1
    while i > 0:
        if (catalog[i].spatial.x - catalog[i - 1].spatial.x) ** 2 + \
           (catalog[i].spatial.y - catalog[i - 1].spatial.y) ** 2 + \
           (catalog[i].spatial.z - catalog[i - 1].spatial.z) ** 2 < (5e-5) ** 2:
            if catalog[i].magnitude > catalog[i - 1].magnitude:
                catalog.pop(i)
            else:
                catalog.pop(i - 1)
            i -= 1
        else:
            i -= 1

    return catalog


def spherical_to_spatial(ra: float, de: float) -> Vec3:
    """Convert right ascension/declination to unit vector.

    Must be defined here to avoid circular imports with math_utils.
    """
    import math
    return Vec3(
        math.cos(ra) * math.cos(de),
        math.sin(ra) * math.cos(de),
        math.sin(de),
    )


# ─────────────────────────────────────────
# CATALOG FILTERING
# ─────────────────────────────────────────


def narrow_catalog(catalog: Catalog,
                   max_magnitude: int,
                   max_stars: int,
                   min_separation: float) -> Catalog:
    """Remove unwanted stars from a catalog.

    Steps:
        1. Remove stars dimmer than max_magnitude.
        2. Remove stars too close to each other (< min_separation radians).
        3. Keep only the max_stars brightest stars.

    Args:
        catalog: Input catalog.
        max_magnitude: Maximum magnitude * 100 (brighter = lower number).
                       Pass 9999 to not filter by magnitude.
        max_stars: Maximum number of stars to keep.
        min_separation: Minimum angular separation in radians.
                       Pass -1 to not filter by separation.

    Reference:
        Original C++ source: src/star-utils.cpp:17-51
    """
    # Step 1: magnitude filter
    result = [s for s in catalog if s.magnitude <= max_magnitude]

    # Step 2: separation filter
    too_close = set()
    for i in range(len(result)):
        for j in range(i + 1, len(result)):
            if angle_unit(result[i].spatial, result[j].spatial) < min_separation:
                too_close.add(i)
                too_close.add(j)

    result = [s for idx, s in enumerate(result) if idx not in too_close]

    # Step 3: brightest-N filter
    if max_stars < len(result):
        result.sort(key=lambda s: s.magnitude)  # lower magnitude = brighter
        result = result[:max_stars]

    return result


def mag_to_brightness(magnitude: int) -> float:
    """Convert catalog magnitude to relative brightness.

    Formula: 10^(-mag / 250)

    The magnitude input is magnitude * 100 (e.g., 150 → 1.50).
    The factor of 250 accounts for this scaling (100 / 0.4 ≈ 250).

    Reference:
        Original C++ source: src/star-utils.cpp:147-149
    """
    return math.pow(10.0, -magnitude / 250.0)


def find_named_star(catalog: Catalog, name: int) -> Optional[CatalogStar]:
    """Find a catalog star by its name.

    Reference:
        Original C++ source: src/star-utils.cpp:54-61
    """
    for star in catalog:
        if star.name == name:
            return star
    return None


# ─────────────────────────────────────────
# K-VECTOR INDEX
# ─────────────────────────────────────────


class KVectorIndex:
    """A data structure enabling constant-time range queries.

    The K-vector index pre-computes cumulative counts for equally-spaced
    bins between min and max. Query [q_min, q_max] returns the index
    range of values that fall (approximately) within the query.

    Reference:
        Original C++ source: src/databases.cpp:64-164
    """

    def __init__(self, values: List[float], min_val: float, max_val: float,
                 num_bins: int):
        num_values = len(values)

        # Generate K-vector bins
        bin_width = (max_val - min_val) / num_bins
        k_vector = [0] * (num_bins + 1)

        last_bin = 0
        for i in range(num_values):
            this_bin = int(math.ceil((values[i] - min_val) / bin_width))
            for bin_idx in range(last_bin, this_bin):
                k_vector[bin_idx] = i
            last_bin = this_bin
        for bin_idx in range(last_bin, num_bins + 1):
            k_vector[bin_idx] = num_values

        self._num_values = num_values
        self._min = min_val
        self._max = max_val
        self._bin_width = bin_width
        self._num_bins = num_bins
        self._bins = k_vector

    @staticmethod
    def deserialize(min_val: float, max_val: float, num_bins: int,
                    bins: List[int], num_values: int) -> KVectorIndex:
        """Create a KVectorIndex from pre-computed values (used when loading)."""
        idx = object.__new__(KVectorIndex)
        idx._num_values = num_values
        idx._min = min_val
        idx._max = max_val
        idx._bin_width = (max_val - min_val) / num_bins if num_bins > 0 else 0
        idx._num_bins = num_bins
        idx._bins = bins
        return idx

    @property
    def num_values(self) -> int:
        return self._num_values

    @property
    def min_val(self) -> float:
        return self._min

    @property
    def max_val(self) -> float:
        return self._max

    @property
    def num_bins(self) -> int:
        return self._num_bins

    def _bin_for(self, query: float) -> int:
        """Find the lowest-indexed bin for a query value.

        Reference:
            Original C++ source: src/databases.cpp:159-164
        """
        result = int(math.ceil((query - self._min) / self._bin_width))
        return result

    def query_liberal(self, min_query: float, max_query: float) -> Tuple[int, int]:
        """Return at least all values in [min_query, max_query].

        May return a few extra entries just outside the range (conservative).

        Returns (lower_index, upper_index) where lower_index is the first
        matching value and upper_index is one-past-the-last matching value.

        Reference:
            Original C++ source: src/databases.cpp:129-156
        """
        if max_query >= self._max:
            max_query = self._max - 0.00001
        if min_query <= self._min:
            min_query = self._min + 0.00001

        if min_query > self._max or max_query < self._min:
            return (0, 0)

        lower_bin = self._bin_for(min_query)
        upper_bin = self._bin_for(max_query)

        lower_idx = self._bins[lower_bin - 1]
        if lower_idx >= self._num_values:
            return (0, 0)

        upper_idx = self._bins[upper_bin]
        return (lower_idx, upper_idx)


# ─────────────────────────────────────────
# PAIR-DISTANCE K-VECTOR DATABASE
# ─────────────────────────────────────────


class KVectorPair:
    """A pair of catalog stars with their angular distance."""

    __slots__ = ('index1', 'index2', 'distance')

    def __init__(self, index1: int, index2: int, distance: float):
        self.index1 = index1
        self.index2 = index2
        self.distance = distance


def catalog_to_pair_distances(catalog: Catalog,
                              min_distance: float,
                              max_distance: float) -> List[KVectorPair]:
    """Compute all inter-star distances within the given range.

    Reference:
        Original C++ source: src/databases.cpp:175-192
    """
    pairs = []
    for i in range(len(catalog)):
        for k in range(i + 1, len(catalog)):
            dist = angle_unit(catalog[i].spatial, catalog[k].spatial)
            if min_distance <= dist <= max_distance:
                pairs.append(KVectorPair(i, k, dist))
    return pairs


class PairDistanceKVectorDatabase:
    """A database of inter-star distances with fast range queries.

    Stores all pairs of catalog stars whose angular separation falls
    within [min_distance, max_distance]. Supports fast lookup of pairs
    by distance.

    Magic value: 0x2536f009

    Reference:
        Original C++ source: src/databases.cpp:220-295
    """

    MAGIC_VALUE = 0x2536F009

    def __init__(self, catalog: Catalog,
                 min_distance: float,
                 max_distance: float,
                 num_bins: int):
        pairs = catalog_to_pair_distances(catalog, min_distance, max_distance)
        pairs.sort(key=lambda p: p.distance)
        distances = [p.distance for p in pairs]

        self._index = KVectorIndex(distances, min_distance, max_distance, num_bins)
        self._pairs = pairs
        self._catalog = catalog

    @property
    def min_distance(self) -> float:
        return self._index.min_val

    @property
    def max_distance(self) -> float:
        return self._index.max_val

    @property
    def num_pairs(self) -> int:
        return self._index.num_values

    def find_pairs_liberal(self, min_query: float, max_query: float) -> List[KVectorPair]:
        """Find all pairs with distance in [min_query, max_query].

        Conservative: may return extra pairs just outside the range.

        Reference:
            Original C++ source: src/databases.cpp:237-246
        """
        lower, upper = self._index.query_liberal(min_query, max_query)
        return self._pairs[lower:upper]

    def find_pairs_exact(self, min_query: float, max_query: float) -> List[KVectorPair]:
        """Find all pairs with distance strictly in [min_query, max_query].

        Uses cosine comparison to eliminate edge pairs that were included
        by the conservative K-vector query.

        Reference:
            Original C++ source: src/databases.cpp:248-279
        """
        max_query_cos = math.cos(min_query)
        min_query_cos = math.cos(max_query)

        lower, upper = self._index.query_liberal(min_query, max_query)

        while lower < upper:
            p = self._pairs[lower]
            cos_val = self._catalog[p.index1].spatial.dot(self._catalog[p.index2].spatial)
            if cos_val < max_query_cos:  # angle > min_query
                break
            lower += 1

        while lower < upper:
            p = self._pairs[upper - 1]
            cos_val = self._catalog[p.index1].spatial.dot(self._catalog[p.index2].spatial)
            if cos_val > min_query_cos:  # angle < max_query
                break
            upper -= 1

        return self._pairs[lower:upper]

    def star_distances(self, star_index: int) -> List[float]:
        """Get all distances from a given star to its paired stars.

        For debugging.

        Reference:
            Original C++ source: src/databases.cpp:287-295
        """
        result = []
        for p in self._pairs:
            if p.index1 == star_index or p.index2 == star_index:
                result.append(p.distance)
        return result


# ─────────────────────────────────────────
# MULTI-DATABASE
# ─────────────────────────────────────────


class MultiDatabase:
    """Container for multiple sub-databases.

    A multi-database maps "magic values" to binary database buffers.
    This is the top-level database format used by LOST.

    Reference:
        Original C++ source: src/databases.cpp:297-364
    """

    CATALOG_MAGIC = 0xF9A283BC

    def __init__(self):
        self._subdatabases: dict[int, bytes] = {}

    def add_subdatabase(self, magic_value: int, data: bytes) -> None:
        self._subdatabases[magic_value] = data

    def get_subdatabase(self, magic_value: int) -> Optional[bytes]:
        return self._subdatabases.get(magic_value)
