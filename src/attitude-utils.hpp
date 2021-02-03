#ifndef ATTITUDE_UTILS_H
#define ATTITUDE_UTILS_H

#include <memory>

namespace lost {

// At first, I wanted to have two separate Attitude classes, one storing Euler angles and converting
// to Quaterinon, and another storing as Quaternion and converting to Euler. But abstract classes
// make everything more annoying, because you need vectors of pointers...ugh!

struct Vec2 {
    float x;
    float y;
};

struct Vec3 {
    float x;
    float y;
    float z;
};

class Quaternion {
public:
    Quaternion() = default; // I guess this lets you call Attitude on an attitude object?
    Quaternion(const Vec3 &); // Pure quaternion
    Quaternion(const Vec3 &, float);
    Quaternion(float real, float i, float j, float k)
        : real(real), i(i), j(j), k(k) { };

    Quaternion operator*(const Quaternion &other) const;
    Quaternion Conjugate() const;
    Vec3 Vector() const;
    Vec3 Rotate(const Vec3 &) const;
    float Angle() const;
    void ToSpherical(float *ra, float *dec, float *roll) const;

    float real;
    float i;
    float j;
    float k;
};

// Return a quaternion that will reorient the coordinate axes so that the x-axis points at the given
// right ascension and declination, then roll the coordinate axes counterclockwise (i.e., the stars
// will appear to rotate clockwise). This is an "improper" z-y'-x' Euler rotation.
Quaternion SphericalToQuaternion(float ra, float dec, float roll);

float RadToDeg(float);
float DegToRad(float);
float RadToArcSec(float);
float ArcSecToRad(float);

float GreatCircleDistance(float ra1, float de1, float ra2, float de2);

// TODO: quaternion and euler angle conversion, conversion between ascension/declination to rec9tu

}

#endif
