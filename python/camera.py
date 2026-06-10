"""
====================================================
CAMERA MODEL
====================================================

Purpose:
    Describes the camera's optical properties and provides
    conversions between pixel coordinates and spatial vectors.

    The camera coordinate system has:
        - X pointing away from the sensor (boresight direction)
        - Y and Z spanning the sensor plane
        - Origin at the principal point (center of sensor)

Mathematical Background:
    The pinhole camera model maps a 3D point (x, y, z) in camera
    coordinates to a 2D pixel (u, v) via:
        u = -y * f / x + cx
        v = -z * f / x + cy
    where f is focal length in pixels, (cx, cy) is the principal point.

    The inverse mapping places points on a plane at x=1:
        y' = -(u - cx) / f
        z' = -(v - cy) / f
        returns (1, y', z')

Reference:
    Original C++ source: src/camera.hpp, src/camera.cpp
"""

from __future__ import annotations

import math

from math_utils import Vec2, Vec3, angle_unit


class Camera:
    """Full description of a camera.

    Stores focal length, principal point, and resolution.
    Provides spatial↔pixel coordinate conversions.

    Args:
        focal_length: Focal length in pixels
        x_center: Principal point X coordinate (pixels)
        y_center: Principal point Y coordinate (pixels)
        x_resolution: Sensor width in pixels
        y_resolution: Sensor height in pixels

    Reference:
        Original C++ source: src/camera.hpp:8-58, src/camera.cpp:1-70
    """

    def __init__(self, focal_length: float,
                 x_center: float, y_center: float,
                 x_resolution: int, y_resolution: int):
        self._focal_length = float(focal_length)
        self._x_center = float(x_center)
        self._y_center = float(y_center)
        self._x_resolution = int(x_resolution)
        self._y_resolution = int(y_resolution)

    @classmethod
    def from_center(cls, focal_length: float,
                    x_resolution: int, y_resolution: int) -> Camera:
        """Create a camera with principal point at the sensor center."""
        return cls(focal_length,
                   x_resolution / 2.0, y_resolution / 2.0,
                   x_resolution, y_resolution)

    def spatial_to_camera(self, vector: Vec3) -> Vec2:
        """Convert a 3D spatial point to 2D pixel coordinates.

        Assumes X is the depth direction pointing away from the sensor.
        Any vector (x, 0, 0) maps to the principal point.

        Reference:
            Original C++ source: src/camera.cpp:15-26
        """
        assert vector.x > 0, f"Point behind camera: x={vector.x}"

        focal_factor = self._focal_length / vector.x

        y_pixel = vector.y * focal_factor
        z_pixel = vector.z * focal_factor

        return Vec2(-y_pixel + self._x_center, -z_pixel + self._y_center)

    def camera_to_spatial(self, vector: Vec2) -> Vec3:
        """Convert 2D pixel coordinates to a 3D spatial direction.

        Returns a vector with x-component equal to 1 (placed one unit away
        along the boresight). The returned vector is NOT necessarily
        a unit vector — call .normalize() if you need one.

        WARNING: Other functions rely on x-component being 1.

        Reference:
            Original C++ source: src/camera.cpp:35-48
        """
        x_pixel = -vector.x + self._x_center
        y_pixel = -vector.y + self._y_center

        return Vec3(
            1,
            x_pixel / self._focal_length,
            y_pixel / self._focal_length,
        )

    def in_sensor(self, vector: Vec2) -> bool:
        """Check whether a pixel coordinate is within the sensor bounds.

        Returns True if the point is on or within the sensor edges.

        Reference:
            Original C++ source: src/camera.cpp:51-56
        """
        return (0 <= vector.x <= self._x_resolution
                and 0 <= vector.y <= self._y_resolution)

    @property
    def x_resolution(self) -> int:
        return self._x_resolution

    @property
    def y_resolution(self) -> int:
        return self._y_resolution

    @property
    def focal_length(self) -> float:
        return self._focal_length

    @property
    def fov(self) -> float:
        """Horizontal field of view in radians.

        Reference:
            Original C++ source: src/camera.cpp:66-68
        """
        return focal_length_to_fov(self._focal_length, self._x_resolution, 1.0)

    @focal_length.setter
    def focal_length(self, value: float) -> None:
        self._focal_length = float(value)


# ─────────────────────────────────────────
# FOV <-> FOCAL LENGTH CONVERSIONS
# ─────────────────────────────────────────


def fov_to_focal_length(x_fov: float, x_resolution: int) -> float:
    """Convert horizontal field of view to focal length in pixels.

    Reference:
        Original C++ source: src/camera.cpp:58-60
    """
    return x_resolution / 2.0 / math.tan(x_fov / 2)


def focal_length_to_fov(focal_length: float, x_resolution: int,
                        pixel_size: float) -> float:
    """Convert focal length (with pixel size) to horizontal field of view.

    Reference:
        Original C++ source: src/camera.cpp:62-64
    """
    return math.atan2(x_resolution / 2 * pixel_size, focal_length) * 2
