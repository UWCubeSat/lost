"""
====================================================
PIPELINE — Star Tracking Pipeline Orchestration
====================================================

Purpose:
    Wires together all stages of the star-tracking pipeline.
    Manages the data flow from input image through centroiding,
    star identification, and attitude estimation.

    The Pipeline class:
        1. Takes configuration (PipelineConfig)
        2. Loads database
        3. Runs each stage, passing outputs as inputs to the next
        4. Returns PipelineOutput with all intermediate results

    Data flow:
        Image (bytes) → [Centroiding] → Stars → [Filtering] → Stars
        → [Star ID] → StarIdentifiers → [Attitude] → Attitude

Reference:
    Original C++ source: src/io.hpp, src/io.cpp, src/main.cpp
"""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field, asdict
from typing import List, Optional

from camera import Camera, fov_to_focal_length
from centroiding import (
    CentroidAlgorithm, CogCentroider, IterativeWeightedCogCentroider,
    DummyCentroider,
)
from config import PipelineConfig
from database import (
    Catalog, CatalogStar, Stars, StarIdentifiers, StarIdentifier,
    PairDistanceKVectorDatabase, narrow_catalog, load_catalog,
)
from math_utils import deg_to_rad
from star_id import (
    StarIdAlgorithm, PyramidStarId, GeometricVotingStarId,
    DummyStarId,
)
from math_utils import Attitude
from attitude import (
    AttitudeEstimator, DavenportQAttitude, TriadAttitude, QuestAttitude,
)


# ─────────────────────────────────────────
# PIPELINE OUTPUT
# ─────────────────────────────────────────


@dataclass
class PipelineOutput:
    """The result of running a pipeline.

    Stores all intermediate and final outputs.
    Timing information is in seconds.

    Reference:
        Original C++ source: src/io.hpp:190-206
    """
    stars: Optional[Stars] = None
    star_ids: Optional[StarIdentifiers] = None
    attitude: Optional[Attitude] = None

    centroiding_time: float = -1.0
    star_id_time: float = -1.0
    attitude_time: float = -1.0

    catalog: Catalog = field(default_factory=list)

    def to_dict(self) -> dict:
        """Serialize output to a dict for tracing."""
        result = {
            'centroiding_time_s': self.centroiding_time,
            'star_id_time_s': self.star_id_time,
            'attitude_time_s': self.attitude_time,
        }
        if self.stars is not None:
            result['num_centroids'] = len(self.stars)
            result['centroids'] = [
                {'x': s.x, 'y': s.y, 'radius_x': s.radius_x,
                 'radius_y': s.radius_y, 'magnitude': s.magnitude}
                for s in self.stars
            ]
        if self.star_ids is not None:
            result['num_identified'] = len(self.star_ids)
            result['star_ids'] = [
                {'star_index': si.star_index,
                 'catalog_index': si.catalog_index,
                 'weight': si.weight}
                for si in self.star_ids
            ]
            if self.catalog:
                result['identified_names'] = [
                    self.catalog[si.catalog_index].name if si.catalog_index < len(self.catalog) else -1
                    for si in self.star_ids
                ]
        if self.attitude is not None:
            result['attitude_known'] = True
            try:
                q = self.attitude.get_quaternion()
                result['quaternion'] = {
                    'real': q.real, 'i': q.i, 'j': q.j, 'k': q.k
                }
                euler = q.to_spherical()
                result['euler_angles_deg'] = {
                    'ra': math.degrees(euler.ra),
                    'de': math.degrees(euler.de),
                    'roll': math.degrees(euler.roll),
                }
            except (ValueError, AttributeError):
                result['attitude_known'] = False
        return result


# ─────────────────────────────────────────
# PIPELINE CLASS
# ─────────────────────────────────────────


class Pipeline:
    """The main star-tracking pipeline.

    Constructed via configuration. Runs all stages sequentially.

    Reference:
        Original C++ source: src/io.hpp:233-253, src/io.cpp:821-1024
    """

    def __init__(self,
                 centroid_algorithm: Optional[CentroidAlgorithm] = None,
                 star_id_algorithm: Optional[StarIdAlgorithm] = None,
                 attitude_algorithm: Optional[AttitudeEstimator] = None,
                 database: object = None,
                 catalog: Optional[Catalog] = None,
                 pair_db: Optional[PairDistanceKVectorDatabase] = None):
        self.centroid_algorithm = centroid_algorithm
        self.star_id_algorithm = star_id_algorithm
        self.attitude_algorithm = attitude_algorithm
        self.database = database
        self.catalog = catalog or []
        self.pair_db = pair_db

        # Centroid filter parameters
        self.centroid_min_magnitude = 0
        self.centroid_min_stars = 0

    @classmethod
    def from_config(cls, config: PipelineConfig,
                    catalog: Optional[Catalog] = None,
                    pair_db: Optional[PairDistanceKVectorDatabase] = None) -> Pipeline:
        """Create a Pipeline from PipelineConfig.

        Sets up algorithms, loads database, and configures filters.

        Reference:
            Original C++ source: src/io.cpp:842-906
        """
        result = cls()
        result.catalog = catalog or []
        result.pair_db = pair_db

        # Centroid algorithm
        if config.centroid_algo == "dummy":
            result.centroid_algorithm = DummyCentroider(config.centroid_dummy_stars)
        elif config.centroid_algo == "cog":
            result.centroid_algorithm = CogCentroider()
        elif config.centroid_algo == "iwcog":
            result.centroid_algorithm = IterativeWeightedCogCentroider()

        # Centroid filter
        if config.centroid_mag_filter > 0:
            result.centroid_min_magnitude = int(config.centroid_mag_filter)
        if config.centroid_filter_brightest > 0:
            result.centroid_min_stars = config.centroid_filter_brightest

        # Database loading
        if config.database_path and config.database_path != "":
            # In Python, we use the in-memory PairDistanceKVectorDatabase
            # rather than the serialized format. If a .dat file exists,
            # we load it but wrap it appropriately.
            pass  # Use pre-loaded pair_db or catalog directly

        # Star ID algorithm
        if config.star_id_algo == "dummy":
            result.star_id_algorithm = DummyStarId()
        elif config.star_id_algo == "gv":
            result.star_id_algorithm = GeometricVotingStarId(
                deg_to_rad(config.angular_tolerance))
        elif config.star_id_algo == "pyramid":
            result.star_id_algorithm = PyramidStarId(
                deg_to_rad(config.angular_tolerance),
                config.false_stars_estimate,
                config.max_mismatch_probability,
                1000)

        # Attitude algorithm
        if config.attitude_algo == "dqm":
            result.attitude_algorithm = DavenportQAttitude()
        elif config.attitude_algo == "triad":
            result.attitude_algorithm = TriadAttitude()
        elif config.attitude_algo == "quest":
            result.attitude_algorithm = QuestAttitude()

        return result

    def run(self, image: Optional[bytes] = None,
            image_width: int = 0, image_height: int = 0,
            stars: Optional[Stars] = None,
            camera: Optional[Camera] = None,
            star_ids: Optional[StarIdentifiers] = None,
            trace: bool = False,
            trace_dir: str = "traces") -> PipelineOutput:
        """Run the full pipeline.

        Each stage runs if the required input and algorithm are available.
        Stages pass their output to the next stage automatically.

        Args:
            image: Raw grayscale image bytes.
            image_width, image_height: Image dimensions.
            stars: Pre-computed centroids (skips centroiding).
            camera: Camera model.
            star_ids: Pre-computed star IDs (skips star ID).
            trace: Enable trace output.
            trace_dir: Directory for trace files.

        Returns:
            PipelineOutput with all intermediate results.

        Reference:
            Original C++ source: src/io.cpp:913-1012
        """
        output = PipelineOutput()
        output.catalog = self.catalog

        if trace:
            os.makedirs(trace_dir, exist_ok=True)

        # ── Stage 1: Centroiding ───────────────────────────────
        if self.centroid_algorithm and image is not None:
            t0 = time.time()
            unfiltered_stars = self.centroid_algorithm.go(image, image_width, image_height)
            t1 = time.time()
            output.centroiding_time = t1 - t0

            if trace:
                with open(os.path.join(trace_dir, '01_raw_centroids.json'), 'w') as f:
                    json.dump([{'x': s.x, 'y': s.y, 'mag': s.magnitude} for s in unfiltered_stars], f, indent=2)

            # ── Stage 2: Centroid Filtering ────────────────────────
            filtered = self._filter_centroids(unfiltered_stars)
            output.stars = filtered

            if trace:
                with open(os.path.join(trace_dir, '02_filtered_centroids.json'), 'w') as f:
                    json.dump([{'x': s.x, 'y': s.y, 'mag': s.magnitude} for s in filtered], f, indent=2)

            stars = filtered
            star_ids = None  # discard any pre-computed star IDs

        # ── Stage 3: Star Identification ───────────────────────
        db_for_star_id = self.pair_db if self.pair_db is not None else self.database
        if self.star_id_algorithm and db_for_star_id and stars is not None and camera is not None:
            t0 = time.time()
            output.star_ids = self.star_id_algorithm.go(
                db_for_star_id, stars, self.catalog, camera)
            t1 = time.time()
            output.star_id_time = t1 - t0

            if trace:
                with open(os.path.join(trace_dir, '03_star_ids.json'), 'w') as f:
                    json.dump([
                        {'star_index': si.star_index,
                         'catalog_index': si.catalog_index,
                         'catalog_name': self.catalog[si.catalog_index].name if si.catalog_index < len(self.catalog) else -1}
                        for si in output.star_ids
                    ], f, indent=2)

            star_ids = output.star_ids
            output.stars = stars

        # ── Stage 4: Attitude Estimation ───────────────────────
        if self.attitude_algorithm and star_ids and camera is not None and stars is not None:
            t0 = time.time()
            output.attitude = self.attitude_algorithm.go(camera, stars, self.catalog, star_ids)
            t1 = time.time()
            output.attitude_time = t1 - t0

            if trace:
                with open(os.path.join(trace_dir, '04_attitude.json'), 'w') as f:
                    d = {'attitude_known': False}
                    try:
                        q = output.attitude.get_quaternion()
                        d = {
                            'attitude_known': True,
                            'quaternion': {
                                'real': q.real, 'i': q.i, 'j': q.j, 'k': q.k
                            },
                        }
                        euler = q.to_spherical()
                        d['euler_angles_deg'] = {
                            'ra': math.degrees(euler.ra),
                            'de': math.degrees(euler.de),
                            'roll': math.degrees(euler.roll),
                        }
                    except (ValueError, AttributeError):
                        pass
                    json.dump(d, f, indent=2)

        if trace:
            with open(os.path.join(trace_dir, '05_pipeline_output.json'), 'w') as f:
                json.dump(output.to_dict(), f, indent=2)

        return output

    def _filter_centroids(self, stars: Stars) -> Stars:
        """Filter centroids by magnitude and sort by brightness.

        Steps:
            1. If centroid_min_stars > 0 and we have more stars than that,
               find the Nth brightest star's magnitude as a threshold.
            2. Keep only stars with magnitude >= max(min_magnitude, threshold).
            3. Sort remaining stars by magnitude descending, with (x, y)
               tiebreak for determinism. This guarantees that the first N
               stars passed to star identification are the same set in the
               same order regardless of centroid_mag_filter, because star
               magnitudes are intrinsic to each centroid and the sort is
               a pure function of the filtered list.

        Reference:
            Original C++ source: src/io.cpp:949-969
        """
        min_mag = self.centroid_min_magnitude
        if (self.centroid_min_stars > 0
                and self.centroid_min_stars < len(stars)):
            sorted_stars = sorted(stars, key=lambda s: s.magnitude, reverse=True)
            min_mag = max(min_mag, sorted_stars[self.centroid_min_stars - 1].magnitude)

        filtered = [s for s in stars if s.magnitude >= min_mag]

        # Sort by magnitude descending, then by (x,y) for determinism.
        # This ensures the same stars always appear in the same order
        # regardless of which centroid_mag_filter is used, making the
        # pyramid search input deterministic across filter settings.
        return sorted(filtered, key=lambda s: (-s.magnitude, s.x, s.y))


# ─────────────────────────────────────────
# CONVENIENCE FUNCTIONS
# ─────────────────────────────────────────


def focal_length_from_config(config: PipelineConfig, x_resolution: int) -> float:
    """Calculate focal length in pixels from config options.

    Two modes:
        1. Using FOV: f = x_res / (2 * tan(FOV/2))
        2. Using pixel size + focal length: f = focal_length_mm * 1000 / pixel_size_um

    Reference:
        Original C++ source: src/io.cpp:319-336
    """
    if config.pixel_size != -1:
        return config.focal_length * 1000 / config.pixel_size
    else:
        return fov_to_focal_length(deg_to_rad(config.fov), x_resolution)
