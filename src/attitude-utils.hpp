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
    Vec3 crossProduct(const Vec3 &) const;
};

class Mat3 {
public:
    float x[9];

    float At(int i, int j) const;
    Mat3 operator*(const Mat3 &) const;
    Vec3 operator*(const Vec3 &) const;
    Mat3 Transpose() const;
    Vec3 Column(int) const;
    Vec3 Row(int) const;
};

long SerializeLengthVec3();
void SerializeVec3(const Vec3 &, unsigned char *);
Vec3 DeserializeVec3(const unsigned char *);

float Distance(const Vec2 &, const Vec2 &);
float Distance(const Vec3 &, const Vec3 &);

class EulerAngles {
public:
    EulerAngles(float ra, float de, float roll)
        : ra(ra), de(de), roll(roll) { };

    float ra;
    float de;
    float roll;
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
    EulerAngles ToSpherical() const;
    bool IsUnit(float tolerance) const;
    // A canonical rotation quaternion's first component should be positive. A quaternion and its
    // negative represent the same rotation.
    Quaternion Canonicalize() const;

    float real;
    float i;
    float j;
    float k;
};

// When porting to a embedded device, you'll probalby want to get rid of this class and adapt to
// either quaternions or DCMs exclusively, depending on the natural output format of whatever
// attitude estimation algorithm you're using. This Attitude class stores either quaternion or DCM,
// depending on what's the natural output of the attitude estimation algorithm, then converts to the
// requested format on-demand.
class Attitude {
public:
    Attitude() = default;
    Attitude(const Quaternion &);
    Attitude(const Mat3 &dcm);

    Quaternion GetQuaternion() const;
    Mat3 GetDCM() const;
    EulerAngles ToSpherical() const;
    Vec3 Rotate(const Vec3 &) const;
    
private:
    enum AttitudeType {
        NullType,
        QuaternionType,
        DCMType,
    };

    Quaternion quaternion;
    Mat3 dcm; // direction cosine matrix
    AttitudeType type;
};

Mat3 QuaternionToDCM(const Quaternion &);
Quaternion DCMToQuaternion(const Mat3 &);

// Return a quaternion that will reorient the coordinate axes so that the x-axis points at the given
// right ascension and declination, then roll the coordinate axes counterclockwise (i.e., the stars
// will appear to rotate clockwise). This is an "improper" z-y'-x' Euler rotation.
Quaternion SphericalToQuaternion(float ra, float dec, float roll);

// returns unit vector
Vec3 SphericalToSpatial(float ra, float de);
void SpatialToSpherical(const Vec3 &, float *ra, float *de);
// angle between two vectors, using dot product and magnitude division
float Angle(const Vec3 &, const Vec3 &);
// angle between two vectors, /assuming/ that they are already unit length
float AngleUnit(const Vec3 &, const Vec3 &);

float RadToDeg(float);
float DegToRad(float);
float RadToArcSec(float);
float ArcSecToRad(float);

// TODO: quaternion and euler angle conversion, conversion between ascension/declination to rec9tu

}

#endif
