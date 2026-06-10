"""
====================================================
CONFIGURATION SYSTEM
====================================================

Purpose:
    Every CLI option from the original LOST pipeline and database
    builder becomes a Python configuration field.

    Uses dataclasses with type hints and default values matching LOST.

    PipelineConfig controls:
        - Camera parameters (focal length, pixel size, FOV)
        - Centroiding algorithm selection and parameters
        - Star ID algorithm selection and parameters
        - Attitude estimation algorithm selection
        - Database path
        - Output/comparison options

    DatabaseConfig controls:
        - Catalog filtering (magnitude, max stars, separation)
        - K-vector database parameters
        - Output settings

Reference:
    Original C++ source: src/pipeline-options.hpp, src/database-options.hpp
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class PipelineConfig:
    """All configuration options for the star-tracking pipeline.

    Default values match the C++ LOST defaults exactly.
    """

    # ── Camera Parameters ──────────────────────────────────────────
    png: str = ""
    focal_length: float = 0.0
    pixel_size: float = -1.0
    fov: float = 20.0

    # ── Centroiding ────────────────────────────────────────────────
    centroid_algo: str = "cog"
    centroid_dummy_stars: int = 5
    centroid_mag_filter: float = -1.0
    centroid_filter_brightest: int = -1

    # ── Database ───────────────────────────────────────────────────
    database_path: str = ""

    # ── Star Identification ────────────────────────────────────────
    star_id_algo: str = "pyramid"
    angular_tolerance: float = 0.04
    false_stars_estimate: int = 500
    max_mismatch_probability: float = 0.001

    # ── Attitude Estimation ────────────────────────────────────────
    attitude_algo: str = "dqm"

    # ── Centroid Magnitude Filter (also in pipeline) ───────────────
    centroid_compare_threshold: float = 2.0
    attitude_compare_threshold: float = 1.0

    # ── Generation (NOT ported, but config preserved) ──────────────
    generate: int = 0
    generate_x_res: int = 1024
    generate_y_res: int = 1024
    generate_centroids_only: bool = False
    generate_zero_mag_photons: float = 20000.0
    generate_saturation_photons: float = 150.0
    generate_spread_stddev: float = 1.0
    generate_shot_noise: bool = True
    generate_dark_current: float = 0.1
    generate_read_noise_stddev: float = 0.05
    generate_ra: float = 88.0
    generate_de: float = 7.0
    generate_roll: float = 0.0
    generate_random_attitudes: bool = False
    generate_blur_ra: float = 0.0
    generate_blur_de: float = 0.0
    generate_blur_roll: float = 0.0
    generate_exposure: float = 0.2
    generate_readout_time: float = 0.0
    generate_oversampling: int = 4
    generate_false_stars: int = 0
    generate_false_min_mag: float = 8.0
    generate_false_max_mag: float = 1.0
    generate_perturb_centroids: float = 0.0
    generate_cutoff_mag: float = 6.0
    generate_seed: int = 394859
    generate_time_based_seed: bool = False

    def __post_init__(self):
        # Validate algorithm choices
        valid_centroid = {"cog", "iwcog", "dummy", ""}
        if self.centroid_algo not in valid_centroid:
            raise ValueError(f"Invalid centroid algorithm: {self.centroid_algo}. "
                             f"Choose from {valid_centroid - {''}}")

        valid_star_id = {"pyramid", "gv", "dummy", ""}
        if self.star_id_algo not in valid_star_id:
            raise ValueError(f"Invalid star ID algorithm: {self.star_id_algo}. "
                             f"Choose from {valid_star_id - {''}}")

        valid_attitude = {"dqm", "triad", "quest", ""}
        if self.attitude_algo not in valid_attitude:
            raise ValueError(f"Invalid attitude algorithm: {self.attitude_algo}. "
                             f"Choose from {valid_attitude - {''}}")


@dataclass
class DatabaseConfig:
    """Configuration for building a star catalog database.

    Reference:
        Original C++ source: src/database-options.hpp
    """

    min_mag: float = 100.0
    max_stars: int = 10000
    min_separation: float = 0.08
    kvector: bool = False
    kvector_min_distance: float = 0.5
    kvector_max_distance: float = 15.0
    kvector_distance_bins: int = 10000
    swap_integer_endianness: bool = False
    swap_decimal_endianness: bool = False
    output: str = "-"


@dataclass
class CatalogStar:
    """A star from the Bright Star Catalog.

    Attributes:
        spatial: Unit vector position on celestial sphere.
        magnitude: Magnitude * 100 (integer). E.g. magnitude=150 → 1.50.
        name: Unique numerical identifier (HD number or HR number).

    Reference:
        Original C++ source: src/star-utils.hpp:12-43
    """
    spatial_x: float
    spatial_y: float
    spatial_z: float
    magnitude: int
    name: int


@dataclass
class Star:
    """A centroid detected in an image.

    Attributes:
        x: X pixel coordinate (top-left origin).
        y: Y pixel coordinate (top-left origin).
        radius_x: Approximate horizontal radius of bright area.
        radius_y: Approximate vertical radius of bright area.
        magnitude: Relative brightness measure (larger = brighter).

    Reference:
        Original C++ source: src/star-utils.hpp:49-72
    """
    x: float
    y: float
    radius_x: float
    radius_y: float
    magnitude: int


@dataclass
class StarIdentifier:
    """Records that a detected star corresponds to a catalog star.

    Attributes:
        star_index: Index into the Stars array.
        catalog_index: Index into the Catalog array.
        weight: Confidence weight (usually 1.0).

    Reference:
        Original C++ source: src/star-utils.hpp:78-98
    """
    star_index: int
    catalog_index: int
    weight: float = 1.0
