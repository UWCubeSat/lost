#include "attitude-utils.hpp"

#include <math.h>
#include <assert.h>

namespace lost {

Quaternion Quaternion::operator*(const Quaternion &other) const {
    return Quaternion(
        real*other.real - i*other.i - j*other.j - k*other.k,
        real*other.i + other.real*i + j*other.k - k*other.j,
        real*other.j + other.real*j + k*other.i - i*other.k,
        real*other.k + other.real*k + i*other.j - j*other.i);
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
    real = cos(theta/2);
    // the compiler will optimize it. Right?
    i = input.x * sin(theta/2);
    j = input.y * sin(theta/2);
    k = input.z * sin(theta/2);
}

Vec3 Quaternion::Rotate(const Vec3 &input) const {
    // TODO: optimize
    return ((*this)*Quaternion(input)*Conjugate()).Vector();
}

float Quaternion::Angle() const {
    return acos(real)*2;
}

void Quaternion::ToSpherical(float *ra, float *de, float *roll) const {
    *ra = -atan((2*i*j+2*real*k)/(2*real*real+2*i*i-1));
    *de = -atan((2*j*k+2*real*i)/(2*real*real+2*k*k-1));
    *roll = -asin(-2*i*k+2*real*j);
}

Quaternion SphericalToQuaternion(float ra, float dec, float roll) {
    assert(roll >= 0.0 && roll <= 2*M_PI);
    assert(ra >= 0 && ra <= 2*M_PI);
    assert(dec >= -M_PI && dec <= M_PI);

    // when we are modifying the coordinate axes, the quaternion multiplication works so that the
    // rotations are applied from left to right. This is the opposite as for modifying vectors.
    Quaternion a = Quaternion({ 0, 0, 1 }, ra);
    Quaternion b = Quaternion({ 0, 1, 0 }, -dec);
    Quaternion c = Quaternion({ 1, 0, 0 }, -roll);
    return (a*b*c).Conjugate();
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

float GreatCircleDistance(float ra1, float de1, float ra2, float de2) {
    return 2.0*asin(sqrt(pow(sin(abs(de1-de2)/2.0), 2.0)
                         + cos(de1)*cos(de2)*pow(sin(abs(ra1-ra2)/2.0), 2.0)));
}

}
