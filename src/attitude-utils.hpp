#ifndef ATTITUDE_UTILS_H
#define ATTITUDE_UTILS_H

#include <algorithm>
#include <fstream>
#include <memory>
#include <numeric>  // iota
#include <vector>

namespace lost {

// At first, I wanted to have two separate Attitude classes, one storing Euler angles and converting
// to Quaterinon, and another storing as Quaternion and converting to Euler. But abstract classes
// make everything more annoying, because you need vectors of pointers...ugh!

/// A two dimensional vector with floating point components
struct Vec2 {
    float x;
    float y;

    float Magnitude() const;
    float MagnitudeSq() const;

    Vec2 Normalize() const;

    float operator*(const Vec2 &) const;
    Vec2 operator*(const float &) const;
    Vec2 operator-(const Vec2 &) const;
    Vec2 operator+(const Vec2 &) const;

    friend std::ostream &operator<<(std::ostream &output, const Vec2 &vec);
};

class Mat3;  // define above so we can use in Vec3 class

/// Three dimensional vector with floating point components
class Vec3 {
   public:
    float x;
    float y;
    float z;

    Vec3(){};
    Vec3(float x, float y, float z) : x(x), y(y), z(z){};

    float Magnitude() const;
    float MagnitudeSq() const;
    Vec3 Normalize() const;

    bool operator<(const Vec3 &other) const;

    float operator*(const Vec3 &) const;
    Vec3 operator*(const float &) const;
    Vec3 operator*(const Mat3 &) const;
    Vec3 operator-(const Vec3 &) const;
    Vec3 operator+(const Vec3 &) const;
    Vec3 CrossProduct(const Vec3 &) const;
    Mat3 OuterProduct(const Vec3 &) const;

    friend std::ostream &operator<<(std::ostream &output, const Vec3 &vec);
};

/// 3x3 vector with floating point components
class Mat3 {
   public:
    float x[9];

    float At(int i, int j) const;
    Mat3 operator+(const Mat3 &) const;
    Mat3 operator*(const Mat3 &) const;
    Vec3 operator*(const Vec3 &) const;
    Mat3 operator*(const float &) const;
    Mat3 Transpose() const;
    Vec3 Column(int) const;
    Vec3 Row(int) const;
    float Trace() const;
    float Det() const;
    Mat3 Inverse() const;
};

extern const Mat3 kIdentityMat3;

long SerializeLengthVec3();
void SerializeVec3(const Vec3 &, unsigned char *);
Vec3 DeserializeVec3(const unsigned char *);

float Distance(const Vec2 &, const Vec2 &);
float Distance(const Vec3 &, const Vec3 &);

/**
 * A "human-readable" way to represent a 3d rotation or orientation.
 * Euler angles roughly correspond to yaw, pitch, and roll of an airplane, which are easy for humans
 * to understand. There's no one single way to store Euler angles. We use z-y'-x'' angles, according
 * to the notation used on the wikipedia page for euler angles.
 */
class EulerAngles {
   public:
    EulerAngles(float ra, float de, float roll) : ra(ra), de(de), roll(roll){};

    /// Right ascension. How far we yaw left. Yaw is performed first.
    float ra;
    /// Declination. How far we pitch up (or down if negative). Pitch is performed second, after
    /// yaw.
    float de;
    /// How far we roll counterclockwise. Roll is performed last (after yaw and pitch).
    float roll;
};

/// A quaternion is a common way to represent a 3d rotation.
class Quaternion {
   public:
    Quaternion() = default;
    explicit Quaternion(const Vec3 &);
    Quaternion(const Vec3 &, float);

    Quaternion(float real, float i, float j, float k) : real(real), i(i), j(j), k(k){};

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
 * The attitude (orientation) of a spacecraft.
 * The Attitude object stores either a rotation matrix (direction cosine matrix) or a quaternion,
 * and converts automatically to the other format when needed.
 * @note When porting to an embedded device, you'll probably want to get rid of this class and adapt
 * to either quaternions or DCMs exclusively, depending on the natural output format of whatever
 * attitude estimation algorithm you're using.
 */
class Attitude {
   public:
    Attitude() = default;
    explicit Attitude(const Quaternion &);  // NOLINT
    explicit Attitude(const Mat3 &dcm);

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
    Mat3 dcm;  // direction cosine matrix
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
/// Given a float, find it "modulo" another float, in the true mathematical sense (not remainder).
/// Always returns something in [0,mod) Eg -0.8 mod 0.6 = 0.4
float FloatModulo(float x, float mod);

// Argsort function - Tetra
// Sort first vector based on values of second vector (asc)
template <class T, class U>
std::vector<T> ArgsortVector(std::vector<T> arr, std::vector<U> cmp) {
    std::vector<T> res;
    std::vector<int> indices(arr.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&](int a, int b) -> bool { return cmp[a] < cmp[b]; });

    for (int ind : indices) {
        res.push_back(arr[ind]);
    }
    return res;
}

// TODO: quaternion and euler angle conversion, conversion between ascension/declination to rec9tu

}  // namespace lost

#endif
