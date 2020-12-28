#include "attitude-utils.hpp"

#include <math.h>
#include <assert.h>

namespace lost {

Quaternion Quaternion::operator*(const Quaternion &other) const {
    return Quaternion(
        real*other.real - i*other.i - j*other.j - k*other.k,
        real*other.i + other.real*i + j*other.k - k*other.j,
        real*other.j + other.real*j - i*other.k + k*other.i,
        real*other.k + other.real*k + i*other.j - other.i*j);
}

Quaternion Quaternion::Conjugate() const {
    return Quaternion(real, -i, -j, -k);
}

Vec3 Quaternion::Vector() const {
    return { i, j, k };
}

Quaternion::Quaternion(const Vec3 &input) {
    real = 0;
    i = input.x;
    j = input.y;
    k = input.z;
}

Quaternion::Quaternion(const Vec3 &input, float theta) {
    real = sin(theta/2);
    // the compiler will optimize it. Right?
    i = input.x * cos(theta/2);
    j = input.y * cos(theta/2);
    k = input.z * cos(theta/2);
}

Vec3 Quaternion::Rotate(const Vec3 &input) const {
    // TODO: optimize
    return ((*this)*Quaternion(input)*Conjugate()).Vector();
}

Quaternion SphericalToQuaternion(float ra, float dec, float roll) {
    assert(roll >= 0.0 && roll <= 2*M_PI);
    assert(ra >= 0 && ra <= 2*M_PI);
    assert(dec >= -M_PI && dec <= M_PI);

    // when we are modifying the coordinate axes, the quaternion multiplication works so that the
    // rotations are applied from left to right. This is the opposite as for modifying vectors.
    return Quaternion({ 0, 0, 1 }, -ra)
        * Quaternion({ 0, 1, 0 }, -dec)
        * Quaternion({ 1, 0, 0 }, -roll);
}

float RadToDeg(float rad) {
    return rad*180.0/M_PI;
}

float DegToRad(float deg) {
    return deg/180.0*M_PI;
}

float RadToArcSec(float rad) {
    return RadToDeg(rad) * 3600.0;
}

float ArcSecToRad(float arcSec) {
    return DegToRad(arcSec / 3600.0);
}

}
