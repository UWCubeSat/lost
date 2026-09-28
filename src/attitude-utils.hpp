#ifndef ATTITUDE_UTILS_H
#define ATTITUDE_UTILS_H

#include <memory>
#include <vector>

#include "serialize-helpers.hpp"
#include "scalar.hpp"

namespace lost {

// At first, I wanted to have two separate Attitude classes, one storing Euler angles and converting
// to Quaterinon, and another storing as Quaternion and converting to Euler. But abstract classes
// make everything more annoying, because you need vectors of pointers...ugh!

/// A two dimensional vector with floating point components
struct Vec2 {
    scalar x;
    scalar y;

    scalar Magnitude() const;
    scalar MagnitudeSq() const;

    Vec2 Normalize() const;

    scalar operator*(const Vec2 &) const;
    Vec2 operator*(const scalar &) const;
    Vec2 operator-(const Vec2 &) const;
    Vec2 operator+(const Vec2 &) const;
};

class Mat3; // define above so we can use in Vec3 class

/// Three dimensional vector with floating point components
class Vec3 {
public:
    scalar x;
    scalar y;
    scalar z;

    scalar Magnitude() const;
    scalar MagnitudeSq() const;
    Vec3 Normalize() const;

    scalar operator*(const Vec3 &) const;
    Vec3 operator*(const scalar &) const;
    Vec3 operator*(const Mat3 &) const;
    Vec3 operator-(const Vec3 &) const;
    Vec3 CrossProduct(const Vec3 &) const;
    Mat3 OuterProduct(const Vec3 &) const;
};

/// 3x3 matrix with floating point components
class Mat3 {
public:
    scalar x[9];

    scalar At(int i, int j) const;
    Mat3 operator+(const Mat3 &) const;
    Mat3 operator*(const Mat3 &) const;
    Vec3 operator*(const Vec3 &) const;
    Mat3 operator*(const scalar &) const;
    Mat3 Transpose() const;
    Vec3 Column(int) const;
    Vec3 Row(int) const;
    scalar Trace() const;
    scalar Det() const;
    Mat3 Inverse() const;
};

extern const Mat3 kIdentityMat3;

void SerializeVec3(SerializeContext *, const Vec3 &);
Vec3 DeserializeVec3(DeserializeContext *des);

scalar Distance(const Vec2 &, const Vec2 &);
scalar Distance(const Vec3 &, const Vec3 &);

/**
 * A "human-readable" way to represent a 3d rotation or orientation.
 * Euler angles roughly correspond to yaw, pitch, and roll of an airplane, which are easy for humans to understand.
 * There's no one single way to store Euler angles. We use z-y'-x'' angles, according to the notation used on the wikipedia page for euler angles.
 */
class EulerAngles {
public:
    EulerAngles(scalar ra, scalar de, scalar roll)
        : ra(ra), de(de), roll(roll) { };

    /// Right ascension. How far we yaw left. Yaw is performed first.
    scalar ra;
    /// Declination. How far we pitch up (or down if negative). Pitch is performed second, after yaw.
    scalar de;
    /// How far we roll counterclockwise. Roll is performed last (after yaw and pitch).
    scalar roll;
};

/// A quaternion is a common way to represent a 3d rotation.
class Quaternion {
public:
    Quaternion() = default;
    explicit Quaternion(const Vec3 &);
    Quaternion(const Vec3 &, scalar);

    Quaternion(scalar real, scalar i, scalar j, scalar k)
        : real(real), i(i), j(j), k(k) { };

    Quaternion operator*(const Quaternion &other) const;
    Quaternion Conjugate() const;
    Vec3 Vector() const;
    void SetVector(const Vec3 &);
    Vec3 Rotate(const Vec3 &) const;
    scalar Angle() const;
    /// Returns the smallest angle that can be used to represent the rotation represented by the
    /// quaternion. I.e, min(Angle, 2pi-Angle).
    scalar SmallestAngle() const;
    void SetAngle(scalar);
    EulerAngles ToSpherical() const;
    bool IsUnit(scalar tolerance) const;
    Quaternion Canonicalize() const;

    scalar real;
    scalar i;
    scalar j;
    scalar k;
};

//
/**
 * The attitude (orientation) of a spacecraft.
 * The Attitude object stores either a rotation matrix (direction cosine matrix) or a quaternion,
 * and converts automatically to the other format when needed. Importantly, an attitude may also be
 * "unknown", representing for example if there weren't enough stars to determine an attitude. When
 * set to unknown, it is illegal to try and convert it to anything, so check IsKnown first!
 * @note When porting to an embedded device, you may want to get rid of this class and adapt to
 * either quaternions or DCMs exclusively, depending on the natural output format of whatever
 * attitude estimation algorithm you're using.
 */
class Attitude {
public:
    /// constructs unknown attitude:
    Attitude() = default;
    explicit Attitude(const Quaternion &); // NOLINT
    explicit Attitude(const Mat3 &dcm);

    Quaternion GetQuaternion() const;
    Mat3 GetDCM() const;
    EulerAngles ToSpherical() const;
    Vec3 Rotate(const Vec3 &) const;
    bool IsKnown() const;

private:
    enum class AttitudeType {
        UnknownType, /// Doesn't mean "unknown type", but rather "we have no fucking clue where we're pointing"
        QuaternionType,
        DCMType,
    };

    AttitudeType type;
    Quaternion quaternion;
    Mat3 dcm; // direction cosine matrix
};

Mat3 QuaternionToDCM(const Quaternion &);
Quaternion DCMToQuaternion(const Mat3 &);

/// Return a quaternion that will reorient the coordinate axes so that the x-axis points at the given
/// right ascension and declination, then roll the coordinate axes counterclockwise (i.e., the stars
/// will appear to rotate clockwise). This is an "improper" z-y'-x' Euler rotation.
Quaternion SphericalToQuaternion(scalar ra, scalar dec, scalar roll);

/// returns unit vector
Vec3 SphericalToSpatial(scalar ra, scalar de);
void SpatialToSpherical(const Vec3 &, scalar *ra, scalar *de);
/// angle between two vectors, using dot product and magnitude division
scalar Angle(const Vec3 &, const Vec3 &);
/// angle between two vectors, /assuming/ that they are already unit length
scalar AngleUnit(const Vec3 &, const Vec3 &);

scalar RadToDeg(scalar);
scalar DegToRad(scalar);
scalar RadToArcSec(scalar);
scalar ArcSecToRad(scalar);
/// Given a scalar, find it "modulo" another scalar, in the true mathematical sense (not remainder).
/// Always returns something in [0,mod) Eg -0.8 mod 0.6 = 0.4
scalar ScalarModulo(scalar x, scalar mod);

// TODO: quaternion and euler angle conversion, conversion between ascension/declination to rec9tu

}

#endif
