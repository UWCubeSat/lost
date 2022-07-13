#ifndef ATTITUDE_UTILS_H
#define ATTITUDE_UTILS_H

#include <memory>

namespace lost {

// At first, I wanted to have two separate Attitude classes, one storing Euler angles and converting
// to Quaterinon, and another storing as Quaternion and converting to Euler. But abstract classes
// make everything more annoying, because you need vectors of pointers...ugh!

/**
 * @brief
 * @details
 */
struct Vec2 {
    /// @brief
    float x;

    /// @brief
    float y;

    float Magnitude() const;
    float MagnitudeSq() const;

    /**
     * @brief
     * @return
     */
    Vec2 Normalize() const;

    /**
     * @brief
     * @param vec2
     * @return
     */
    float operator*(const Vec2 &) const;
    Vec2 operator*(const float &) const;
    Vec2 operator-(const Vec2 &) const;
    Vec2 operator+(const Vec2 &) const;
};

class Vec3 {
public:
    float x;
    float y;
    float z;

    float Magnitude() const;
    float MagnitudeSq() const;
    Vec3 Normalize() const;

    float operator*(const Vec3 &) const;
    Vec3 operator*(const float &) const;
    Vec3 operator-(const Vec3 &) const;
    Vec3 crossProduct(const Vec3 &) const;
};

/**
 * @brief
 * @details
 */
class Mat3 {
public:
    /// @brief
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

/**
 * @brief
 * @details
 */
class EulerAngles {
public:
    /**
     * @brief
     * @param ra
     * @param de
     * @param roll
     */
    EulerAngles(float ra, float de, float roll)
        : ra(ra), de(de), roll(roll) { };

    /// @brief
    float ra;

    /// @brief
    float de;

    /// @brief
    float roll;
};

/**
 * @brief
 * @details
 */
class Quaternion {
public:
    /**
     * @brief
     * @note I guess this lets you call Attitude on an attitude object?
     */
    Quaternion() = default;
    Quaternion(const Vec3 &);
    Quaternion(const Vec3 &, float);

    /**
     * @brief
     * @param real
     * @param i
     * @param j
     * @param k
     */
    Quaternion(float real, float i, float j, float k)
        : real(real), i(i), j(j), k(k) { };

    Quaternion operator*(const Quaternion &other) const;
    Quaternion Conjugate() const;
    Vec3 Vector() const;
    void SetVector(const Vec3 &);
    Vec3 Rotate(const Vec3 &) const;
    float Angle() const;
    void SetAngle(float);
    EulerAngles ToSpherical() const;
    bool IsUnit(float tolerance) const;
    Quaternion Canonicalize() const;

    float real;
    float i;
    float j;
    float k;
};

//
/**
 * @brief
 * @details
 * @note When porting to an embedded device, you'll probably want to get rid of this class and adapt to
 * either quaternions or DCMs exclusively, depending on the natural output format of whatever
 * attitude estimation algorithm you're using. This Attitude class stores either quaternion or DCM,
 * depending on what's the natural output of the attitude estimation algorithm, then converts to the
 * requested format on-demand.
 */
class Attitude {
public:
    /// @brief
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
