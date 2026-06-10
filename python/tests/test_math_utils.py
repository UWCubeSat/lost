"""
Tests for math_utils.py — Vectors, Matrices, Quaternions, Attitude
"""

import math
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from math_utils import (
    Vec2, Vec3, Mat3, IDENTITY_MAT3, Quaternion, Attitude,
    EulerAngles, spherical_to_spatial, spatial_to_spherical,
    spherical_to_quaternion, quaternion_to_dcm, dcm_to_quaternion,
    angle, angle_unit, rad_to_deg, deg_to_rad,
    rad_to_arcsec, arcsec_to_rad, decimal_modulo,
)


def test_vec2():
    v = Vec2(3, 4)
    assert v.magnitude() == 5.0
    assert v.magnitude_sq() == 25.0
    n = v.normalize()
    assert abs(n.magnitude() - 1.0) < 1e-10
    assert v.dot(Vec2(1, 0)) == 3.0
    assert Vec2(1, 2) + Vec2(3, 4) == Vec2(4, 6)
    assert Vec2(5, 6) - Vec2(3, 2) == Vec2(2, 4)
    assert Vec2(2, 3) * 2 == Vec2(4, 6)
    assert -Vec2(1, 2) == Vec2(-1, -2)
    print("  test_vec2 OK")


def test_vec3():
    v = Vec3(1, 2, 3)
    assert abs(v.magnitude() - math.sqrt(14)) < 1e-10
    assert abs(v.magnitude_sq() - 14) < 1e-10
    c = v.cross(Vec3(4, 5, 6))
    assert abs(c.x - (-3)) < 1e-10
    assert abs(c.y - 6) < 1e-10
    assert abs(c.z - (-3)) < 1e-10
    d = v.dot(Vec3(1, 0, 0))
    assert abs(d - 1) < 1e-10
    assert Vec3(1, 2, 3) + Vec3(4, 5, 6) == Vec3(5, 7, 9)
    assert Vec3(1, 2, 3) - Vec3(3, 2, 1) == Vec3(-2, 0, 2)
    assert Vec3(2, 3, 4) * 2 == Vec3(4, 6, 8)
    print("  test_vec3 OK")


def test_mat3():
    m = IDENTITY_MAT3
    assert m.det() == 1.0
    assert m.trace() == 3.0
    v = m * Vec3(1, 2, 3)
    assert v == Vec3(1, 2, 3)

    # Rotation matrix about Z by 90 degrees
    c = math.cos(math.pi / 2)
    s = math.sin(math.pi / 2)
    rot_z = Mat3(c, -s, 0, s, c, 0, 0, 0, 1)
    assert abs(rot_z.det() - 1) < 1e-10
    v2 = rot_z * Vec3(1, 0, 0)
    assert abs(v2.x) < 1e-10
    assert abs(v2.y - 1) < 1e-10

    mt = rot_z.transpose()
    assert abs(mt.at(0, 1) - 1) < 1e-10

    m2 = Mat3(4, 3, 2, 1, 0, 5, 2, 1, 3)
    m2i = m2.inverse()
    # m2 * m2i should be identity
    prod = m2 * m2i
    for i in range(3):
        for j in range(3):
            expected = 1.0 if i == j else 0.0
            assert abs(prod.at(i, j) - expected) < 1e-8

    print("  test_mat3 OK")


def test_quaternion():
    # Identity quaternion
    q = Quaternion(1, 0, 0, 0)
    v = q.rotate(Vec3(1, 2, 3))
    assert v == Vec3(1, 2, 3)

    # 90 degrees about Z
    q = Quaternion.from_axis_angle(Vec3(0, 0, 1), math.pi / 2)
    v = q.rotate(Vec3(1, 0, 0))
    assert abs(v.x) < 1e-10
    assert abs(v.y - 1) < 1e-10
    assert abs(v.z) < 1e-10
    assert q.is_unit()

    # Conjugate gives inverse
    q_inv = q.conjugate()
    v2 = q_inv.rotate(v)
    assert abs(v2.x - 1) < 1e-10
    assert abs(v2.y) < 1e-10

    # Canonicalize
    q_neg = Quaternion(-1, 0, 0, 0)
    q_can = q_neg.canonicalize()
    assert q_can.real >= 0

    # Euler angles round-trip
    ra, dec, roll = 1.0, 0.5, 0.3
    q_sph = spherical_to_quaternion(ra, dec, roll)
    euler = q_sph.to_spherical()
    assert abs(euler.ra - ra) < 1e-10
    assert abs(euler.de - dec) < 1e-10

    print("  test_quaternion OK")


def test_attitude():
    q = Quaternion.from_axis_angle(Vec3(0, 0, 1), math.pi / 3)
    a = Attitude(quat=q)
    assert a.is_known()
    assert abs(a.get_quaternion().angle() - math.pi / 3) < 1e-10

    dcm = a.get_dcm()
    a2 = Attitude(dcm=dcm)
    q2 = a2.get_quaternion()
    assert abs(q2.real - q.real) < 1e-10
    assert abs(q2.i - q.i) < 1e-8

    print("  test_attitude OK")


def test_coordinate_transforms():
    ra, de = 2.0, 0.8
    v = spherical_to_spatial(ra, de)
    ra2, de2 = spatial_to_spherical(v)
    assert abs(ra - ra2) < 1e-10
    assert abs(de - de2) < 1e-10

    # Unit sphere check
    assert abs(v.magnitude() - 1) < 1e-10

    print("  test_coordinate_transforms OK")


def test_angle():
    v1 = Vec3(1, 0, 0)
    v2 = Vec3(0, 1, 0)
    assert abs(angle(v1, v2) - math.pi / 2) < 1e-10
    assert abs(angle_unit(v1.normalize(), v2.normalize()) - math.pi / 2) < 1e-10

    # Same vector
    assert angle(v1, v1) == 0.0

    # Opposite
    assert abs(angle(v1, Vec3(-1, 0, 0)) - math.pi) < 1e-6

    print("  test_angle OK")


def test_unit_conversions():
    assert abs(rad_to_deg(math.pi) - 180) < 1e-10
    assert abs(deg_to_rad(180) - math.pi) < 1e-10
    assert abs(rad_to_arcsec(1) - 206264.8) < 1
    assert abs(arcsec_to_rad(206265) - 1) < 0.001

    # Modulo
    assert abs(decimal_modulo(-0.8, 0.6) - 0.4) < 1e-10
    assert abs(decimal_modulo(2.5, 1.0) - 0.5) < 1e-10
    assert abs(decimal_modulo(-0.1, 1.0) - 0.9) < 1e-10

    print("  test_unit_conversions OK")


def test_dcm_quaternion_roundtrip():
    q = Quaternion.from_axis_angle(Vec3(0.6, 0.8, 0), math.radians(45))
    dcm = quaternion_to_dcm(q)
    q2 = dcm_to_quaternion(dcm)
    q2 = q2.canonicalize()
    q = q.canonicalize()
    # Compare
    error = (q * q2.conjugate()).smallest_angle()
    assert error < 1e-8, f"DCM-quaternion round-trip error: {error}"
    print("  test_dcm_quaternion_roundtrip OK")


if __name__ == '__main__':
    print("math_utils tests:")
    test_vec2()
    test_vec3()
    test_mat3()
    test_quaternion()
    test_attitude()
    test_coordinate_transforms()
    test_angle()
    test_unit_conversions()
    test_dcm_quaternion_roundtrip()
    print("\nAll math_utils tests passed!")
