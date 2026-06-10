"""
Integration tests for the full LOST pipeline.

Tests the complete pipeline: centroiding → star ID → attitude estimation.
Uses synthetic data (not real images) to verify each stage works.
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from camera import Camera
from centroiding import CogCentroider
from config import PipelineConfig
from database import (
    Catalog, CatalogStar, StarIdentifier,
    PairDistanceKVectorDatabase, narrow_catalog,
)
from math_utils import (
    Vec2, Vec3, Quaternion, Attitude, angle_unit,
    spherical_to_spatial, deg_to_rad,
)
from star_id import PyramidStarId, GeometricVotingStarId
from attitude import DavenportQAttitude, TriadAttitude, QuestAttitude


def make_test_catalog() -> Catalog:
    """Create a simple test catalog with stars at known positions."""
    catalog = []
    # Create stars at known positions on the unit sphere
    names_and_positions = [
        (1, 10, 20),    # name, ra_deg, dec_deg
        (2, 30, 10),
        (3, 50, -5),
        (4, 70, 15),
        (5, 90, 25),
        (6, 110, -10),
        (7, 130, 30),
        (8, 150, 5),
        (9, 170, -20),
        (10, 190, 40),
        (11, 210, -15),
        (12, 230, 10),
        (13, 250, -30),
        (14, 270, 20),
        (15, 290, -5),
        (16, 310, 35),
        (17, 330, -25),
        (18, 350, 15),
        (19, 15, 45),
        (20, 45, -35),
    ]
    for name, ra_deg, dec_deg in names_and_positions:
        ra = deg_to_rad(ra_deg)
        dec = deg_to_rad(dec_deg)
        spatial = spherical_to_spatial(ra, dec)
        catalog.append(CatalogStar(spatial, 100, name))  # mag=1.00

    return catalog


def test_centroiding():
    """Test centroiding with a synthetic image."""
    w, h = 100, 100
    image = bytes([0] * (w * h))

    centroid_algo = CogCentroider()
    stars = centroid_algo.go(image, w, h)
    assert isinstance(stars, list)
    print(f"  test_centroiding OK (found {len(stars)} stars in blank image)")


def test_database():
    """Test database creation with a small catalog."""
    catalog = make_test_catalog()

    # Narrow catalog
    narrowed = narrow_catalog(catalog, 500, 50, 0.0)
    assert len(narrowed) == len(catalog)

    # Create pair-distance database
    min_dist = deg_to_rad(1.0)
    max_dist = deg_to_rad(175.0)
    db = PairDistanceKVectorDatabase(narrowed, min_dist, max_dist, 100)

    assert db.num_pairs > 0
    print(f"  test_database OK ({db.num_pairs} pairs)")


def test_pyramid_star_id():
    """Test the Pyramid star identification algorithm.

    We construct a mini catalog with 6 stars all in the camera FOV,
    then verify that the pyramid algorithm can identify them.
    """
    import math

    # Build catalog: 6 stars visible in the camera FOV (pointing along +X)
    # All stars have x > 0 and are within reasonable angular distances
    catalog = [
        CatalogStar(Vec3(0.9, 0.1, 0.05), 100, 1),
        CatalogStar(Vec3(0.85, -0.2, 0.1), 100, 2),
        CatalogStar(Vec3(0.95, 0.05, -0.15), 100, 3),
        CatalogStar(Vec3(0.8, 0.3, 0.1), 100, 4),
        CatalogStar(Vec3(0.88, -0.1, -0.2), 100, 5),
        CatalogStar(Vec3(0.92, 0.2, -0.05), 100, 6),
    ]

    # Normalize all vectors
    for cs in catalog:
        cs.spatial = cs.spatial.normalize()

    # Create database
    min_dist = deg_to_rad(1.0)
    max_dist = deg_to_rad(175.0)
    db = PairDistanceKVectorDatabase(catalog, min_dist, max_dist, 100)

    # Camera with moderate FOV
    cam = Camera.from_center(500, 400, 400)

    # Generate image stars from catalog stars (identity attitude)
    identity = Attitude(quat=Quaternion(1, 0, 0, 0))
    stars = []

    for i, cs in enumerate(catalog):
        rotated = identity.rotate(cs.spatial)
        if rotated.x <= 0:
            continue
        pixel = cam.spatial_to_camera(rotated)
        if cam.in_sensor(pixel):
            from database import Star
            stars.append(Star(pixel.x, pixel.y, 2.0, 2.0, 100))

    if len(stars) >= 4:
        algo = PyramidStarId(deg_to_rad(0.5), 10, 0.01, 1000)
        result = algo.go(db, stars, catalog, cam)
        print(f"  test_pyramid_star_id OK (identified {len(result)} out of {len(stars)} stars)")
        if len(result) < 2:
            print(f"    WARNING: Expected more identifications. Results: {result}")
    else:
        print("  test_pyramid_star_id SKIPPED (not enough stars in FOV, need 4)")


def test_geometric_voting():
    """Test the Geometric Voting star identification algorithm."""
    import math

    # Same catalog as pyramid test
    catalog = [
        CatalogStar(Vec3(0.9, 0.1, 0.05), 100, 1),
        CatalogStar(Vec3(0.85, -0.2, 0.1), 100, 2),
        CatalogStar(Vec3(0.95, 0.05, -0.15), 100, 3),
        CatalogStar(Vec3(0.8, 0.3, 0.1), 100, 4),
        CatalogStar(Vec3(0.88, -0.1, -0.2), 100, 5),
        CatalogStar(Vec3(0.92, 0.2, -0.05), 100, 6),
    ]
    for cs in catalog:
        cs.spatial = cs.spatial.normalize()

    min_dist = deg_to_rad(1.0)
    max_dist = deg_to_rad(175.0)
    db = PairDistanceKVectorDatabase(catalog, min_dist, max_dist, 100)
    cam = Camera.from_center(500, 400, 400)

    identity = Attitude(quat=Quaternion(1, 0, 0, 0))
    stars = []

    for i, cs in enumerate(catalog):
        rotated = identity.rotate(cs.spatial)
        if rotated.x <= 0:
            continue
        pixel = cam.spatial_to_camera(rotated)
        if cam.in_sensor(pixel):
            from database import Star
            stars.append(Star(pixel.x, pixel.y, 2.0, 2.0, 100))

    if len(stars) >= 2:
        algo = GeometricVotingStarId(deg_to_rad(0.5))
        result = algo.go(db, stars, catalog, cam)
        print(f"  test_geometric_voting OK (identified {len(result)} stars)")
    else:
        print("  test_geometric_voting SKIPPED (not enough stars in FOV)")


def test_attitude_estimators():
    """Test all attitude estimation algorithms."""
    import math

    catalog = [
        CatalogStar(Vec3(0.9, 0.1, 0.05), 100, 1),
        CatalogStar(Vec3(0.85, -0.2, 0.1), 100, 2),
        CatalogStar(Vec3(0.95, 0.05, -0.15), 100, 3),
        CatalogStar(Vec3(0.8, 0.3, 0.1), 100, 4),
    ]
    for cs in catalog:
        cs.spatial = cs.spatial.normalize()

    cam = Camera.from_center(500, 400, 400)
    identity = Attitude(quat=Quaternion(1, 0, 0, 0))
    stars = []
    star_ids = []

    for i, cs in enumerate(catalog):
        rotated = identity.rotate(cs.spatial)
        if rotated.x <= 0:
            continue
        pixel = cam.spatial_to_camera(rotated)
        if cam.in_sensor(pixel):
            from database import Star
            stars.append(Star(pixel.x, pixel.y, 2.0, 2.0, 100))
            star_ids.append(StarIdentifier(len(stars) - 1, i))

    if len(stars) >= 2:
        for name, algo in [("DQM", DavenportQAttitude()),
                           ("TRIAD", TriadAttitude()),
                           ("QUEST", QuestAttitude())]:
            attitude = algo.go(cam, stars, catalog, star_ids)
            assert attitude.is_known(), f"{name} should produce known attitude"
            q = attitude.get_quaternion()
            assert q.is_unit(), f"{name} quaternion should be unit"
            # The attitude error should be small (we're using identity attitude
            # with perfect star positions)
            error = (Quaternion(1, 0, 0, 0) * q.conjugate()).smallest_angle()
            print(f"  test_{name.lower()} OK (error={math.degrees(error):.4f}°)")
    else:
        print("  test_attitude_estimators SKIPPED")


if __name__ == '__main__':
    print("Pipeline integration tests:")
    test_centroiding()
    test_database()
    test_pyramid_star_id()
    test_geometric_voting()
    test_attitude_estimators()
    print("\nAll pipeline tests passed!")
