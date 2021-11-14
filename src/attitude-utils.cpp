#include "attitude-utils.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>

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
    // Working out these equations would be a pain in the ass. Thankfully, this wikipedia page:
    // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Quaternion_to_Euler_angles_conversion
    // uses almost exactly the same euler angle scheme that we do, so we copy their equations almost
    // wholesale. The only differences are that 1, we use de and roll in the opposite directions,
    // and 2, we store the conjugate of the quaternion (double check why?), which means we need to
    // invert the final de and roll terms, as well as negate all the terms involving a mix between
    // the real and imaginary parts.
    *ra = atan2(2*(-real*k+i*j), 1-2*(j*j+k*k));
    if (*ra < 0)
        *ra += 2*M_PI;
    *de = -asin(2*(-real*j-i*k)); // allow de to be positive or negaive, as is convention
    *roll = -atan2(2*(-real*i+j*k), 1-2*(i*i+j*j));
    if (*roll < 0)
        *roll += 2*M_PI;
}

Quaternion SphericalToQuaternion(float ra, float dec, float roll) {
    assert(roll >= 0.0 && roll <= 2*M_PI);
    assert(ra >= 0 && ra <= 2*M_PI);
    assert(dec >= -M_PI && dec <= M_PI);

    // when we are modifying the coordinate axes, the quaternion multiplication works so that the
    // rotations are applied from left to right. This is the opposite as for modifying vectors.

    // It is indeed correct that a positive rotation in our right-handed coordinate frame is in the
    // clockwise direction when looking down/through the axis of rotation. Just like the right hand
    // rule for magnetic field around a current-carrying conductor.
    Quaternion a = Quaternion({ 0, 0, 1 }, ra);
    Quaternion b = Quaternion({ 0, 1, 0 }, -dec);
    Quaternion c = Quaternion({ 1, 0, 0 }, -roll);
    Quaternion result = (a*b*c).Conjugate();
    assert(result.IsUnit(0.00001));
    return result;
}

bool Quaternion::IsUnit(float tolerance) const {
    return abs(i*i+j*j+k*k+real*real - 1) < tolerance;
}

Vec3 SphericalToSpatial(float ra, float de) {
    return {
        cos(ra)*cos(de),
        sin(ra)*cos(de),
        sin(de),
    };
}

void SpatialToSpherical(const Vec3 &vec, float *ra, float *de) {
    *ra = atan2(vec.y, vec.x);
    if (*ra < 0)
        *ra += M_PI*2;
    *de = asin(vec.z);
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

float Vec3::Magnitude() const {
    return sqrt(x*x+y*y+z*z);
}

Vec3 Vec3::Normalize() const {
    float mag = Magnitude();
    return {
        x/mag, y/mag, z/mag,
    };
}

float Vec3::operator*(const Vec3 &other) const {
    return x*other.x + y*other.y + z*other.z;
}

Vec2 Vec2::operator-(const Vec2 &other) const {
    return { x - other.x, y - other.y };
}

Vec3 Vec3::operator-(const Vec3 &other) const {
    return { x - other.x, y - other.y, z - other.z };
}

Vec3 Vec3::crossProduct(const Vec3 &other) const {
    return {
        y*other.z - z*other.y,
        -(x*other.z - z*other.x),
        x*other.y - y*other.x,
    };
}

long SerializeLengthVec3() {
    return sizeof(float)*3;
}

void SerializeVec3(const Vec3 &vec, unsigned char *buffer) {
    float *fBuffer = (float *)buffer;
    *fBuffer++ = vec.x;
    *fBuffer++ = vec.y;
    *fBuffer = vec.z;
}

Vec3 DeserializeVec3(const unsigned char *buffer) {
    Vec3 result;
    const float *fBuffer = (float *)buffer;
    result.x = *fBuffer++;
    result.y = *fBuffer++;
    result.z = *fBuffer;
    return result;
}

float Angle(const Vec3 &vec1, const Vec3 &vec2) {
    return AngleUnit(vec1.Normalize(), vec2.Normalize());
}

float AngleUnit(const Vec3 &vec1, const Vec3 &vec2) {
    float dot = vec1*vec2;
    // TODO: we shouldn't need this nonsense, right? how come acos sometimes gives nan?
    return dot >= 1 ? 0 : dot <= -1 ? M_PI-0.0000001 : acos(dot);
}

float Distance(const Vec2 &v1, const Vec2 &v2) {
    return sqrt(pow(v1.x-v2.x, 2) + pow(v1.y-v2.y, 2));
}

}
