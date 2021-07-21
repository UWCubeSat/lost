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

    float Magnitude() const;
    Vec2 Normalize() const;

    float operator*(const Vec2 &) const;
    Vec2 operator-(const Vec2 &) const;
};

class Vec3 {
public:
    float x;
    float y;
    float z;

    float Magnitude() const;
    Vec3 Normalize() const;

    float operator*(const Vec3 &) const;
    Vec3 operator-(const Vec3 &) const;
};

long SerializeLengthVec3();
void SerializeVec3(const Vec3 &, unsigned char *);
Vec3 DeserializeVec3(const unsigned char *);

float Distance(const Vec2 &, const Vec2 &);
float Distance(const Vec3 &, const Vec3 &);

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

// returns unit vector
Vec3 SphericalToSpatial(float ra, float de);
// angle between two vectors, using dot product and magnitude division
float Angle(const Vec2 &, const Vec2 &);
float Angle(const Vec3 &, const Vec3 &);
// angle between two vectors, /assuming/ that they are already unit length
float AngleUnit(const Vec2 &, const Vec2 &);
float AngleUnit(const Vec3 &, const Vec3 &);

float RadToDeg(float);
float DegToRad(float);
float RadToArcSec(float);
float ArcSecToRad(float);

// TODO: quaternion and euler angle conversion, conversion between ascension/declination to rec9tu

}

#endif
