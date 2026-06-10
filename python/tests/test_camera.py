"""
Tests for camera.py
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from camera import Camera, fov_to_focal_length, focal_length_to_fov
from math_utils import Vec2, Vec3, deg_to_rad


def test_camera_construction():
    cam = Camera(1000, 512, 512, 1024, 1024)
    assert cam.x_resolution == 1024
    assert cam.y_resolution == 1024
    assert cam.focal_length == 1000
    print("  test_camera_construction OK")


def test_camera_from_center():
    cam = Camera.from_center(1000, 1024, 1024)
    assert cam.x_resolution == 1024
    assert cam.y_resolution == 1024
    assert cam.focal_length == 1000
    print("  test_camera_from_center OK")


def test_spatial_to_camera():
    # Camera: 1024x1024, focal length 1000px, center at (512, 512)
    cam = Camera.from_center(1000, 1024, 1024)

    # A star on the boresight (x, 0, 0) maps to the center
    v = cam.spatial_to_camera(Vec3(100, 0, 0))
    assert abs(v.x - 512) < 1e-6
    assert abs(v.y - 512) < 1e-6

    # A star to the right maps left (negative y in camera)
    v = cam.spatial_to_camera(Vec3(100, 50, 0))
    assert v.x < 512  # Right in space = Left in image

    # A star above maps down (negative z in camera)
    v = cam.spatial_to_camera(Vec3(100, 0, 50))
    assert v.y < 512  # Up in space = Down in image

    print("  test_spatial_to_camera OK")


def test_camera_to_spatial():
    cam = Camera.from_center(1000, 1024, 1024)

    # Center pixel maps to (1, 0, 0)
    v = cam.camera_to_spatial(Vec2(512, 512))
    assert abs(v.x - 1) < 1e-10
    assert abs(v.y) < 1e-10
    assert abs(v.z) < 1e-10

    # Off-center pixel
    v = cam.camera_to_spatial(Vec2(612, 512))
    # x_pixel = -612 + 512 = -100, y_pixel = -512 + 512 = 0
    # result = (1, -100/1000, 0/1000) = (1, -0.1, 0)
    assert abs(v.x - 1) < 1e-10
    assert abs(v.y - (-0.1)) < 1e-10
    assert abs(v.z - 0) < 1e-10

    print("  test_camera_to_spatial OK")


def test_round_trip():
    cam = Camera.from_center(1000, 1024, 1024)

    # Pick a point, go to spatial, then back to camera
    original = Vec2(400, 600)
    spatial = cam.camera_to_spatial(original)
    back = cam.spatial_to_camera(Vec3(spatial.x, spatial.y, spatial.z))

    assert abs(back.x - original.x) < 1e-6
    assert abs(back.y - original.y) < 1e-6

    print("  test_round_trip OK")


def test_in_sensor():
    cam = Camera.from_center(1000, 1024, 1024)

    assert cam.in_sensor(Vec2(0, 0))
    assert cam.in_sensor(Vec2(1024, 1024))
    assert cam.in_sensor(Vec2(512, 512))
    assert not cam.in_sensor(Vec2(-1, 0))
    assert not cam.in_sensor(Vec2(0, 2000))

    print("  test_in_sensor OK")


def test_fov_conversions():
    fov_deg = 20.0
    fov_rad = deg_to_rad(fov_deg)
    res = 1024

    fl = fov_to_focal_length(fov_rad, res)
    fov_back = focal_length_to_fov(fl, res, 1.0)

    assert abs(rad_to_deg(fov_back) - fov_deg) < 0.01

    print("  test_fov_conversions OK")


def rad_to_deg(r):
    return r * 180 / math.pi


if __name__ == '__main__':
    import math
    print("Camera tests:")
    test_camera_construction()
    test_camera_from_center()
    test_spatial_to_camera()
    test_camera_to_spatial()
    test_round_trip()
    test_in_sensor()
    test_fov_conversions()
    print("\nAll camera tests passed!")
