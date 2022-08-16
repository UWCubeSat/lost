#include "attitude-utils.hpp"

#include <math.h>
#include <assert.h>
#include <iostream>

namespace lost {

/// Multiply two quaternions using the usual definition of quaternion multiplication (effectively composes rotations)
Quaternion Quaternion::operator*(const Quaternion &other) const {
    return Quaternion(
        real*other.real - i*other.i - j*other.j - k*other.k,
        real*other.i + other.real*i + j*other.k - k*other.j,
        real*other.j + other.real*j + k*other.i - i*other.k,
        real*other.k + other.real*k + i*other.j - j*other.i);
}

/// Effectively computes a quaternion representing the inverse rotation of the original.
Quaternion Quaternion::Conjugate() const {
    return Quaternion(real, -i, -j, -k);
}

/// The vector formed by imaginary components of the quaternion. The axis of the represented rotation.
Vec3 Quaternion::Vector() const {
    return { i, j, k };
}

/// Set imaginary components.
void Quaternion::SetVector(const Vec3 &vec) {
    i = vec.x;
    j = vec.y;
    k = vec.z;
}

/// Creates a "pure quaternion" with the given vector for imaginary parts and zero for real part.
Quaternion::Quaternion(const Vec3 &input) {
    real = 0;
    SetVector(input);
}

/// Create a quaternion which represents a rotation of theta around the axis input
Quaternion::Quaternion(const Vec3 &input, float theta) {
    real = cos(theta/2);
    // the compiler will optimize it. Right?
    i = input.x * sin(theta/2);
    j = input.y * sin(theta/2);
    k = input.z * sin(theta/2);
}

/// Rotate a 3d vector according to the rotation represented by the quaternion.
Vec3 Quaternion::Rotate(const Vec3 &input) const {
    // TODO: optimize
    return ((*this)*Quaternion(input)*Conjugate()).Vector();
}

/// How many radians the rotation represented by this quaternion has.
float Quaternion::Angle() const {
    if (real <= -1) {
        return 0; // 180*2=360=0
    }
    // TODO: we shouldn't need this nonsense, right? how come acos sometimes gives nan? (same as in AngleUnit)
    return (real >= 1 ? 0 : acos(real))*2;
}

/// Change the amount of rotation in a quaternion, keeping that rotation around the same axis.
void Quaternion::SetAngle(float newAngle) {
    real = cos(newAngle/2);
    SetVector(Vector().Normalize() * sin(newAngle/2));
}

EulerAngles Quaternion::ToSpherical() const {
    // Working out these equations would be a pain in the ass. Thankfully, this wikipedia page:
    // https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles#Quaternion_to_Euler_angles_conversion
    // uses almost exactly the same euler angle scheme that we do, so we copy their equations almost
    // wholesale. The only differences are that 1, we use de and roll in the opposite directions,
    // and 2, we store the conjugate of the quaternion (double check why?), which means we need to
    // invert the final de and roll terms, as well as negate all the terms involving a mix between
    // the real and imaginary parts.
    float ra = atan2(2*(-real*k+i*j), 1-2*(j*j+k*k));
    if (ra < 0)
        ra += 2*M_PI;
    float de = -asin(2*(-real*j-i*k)); // allow de to be positive or negaive, as is convention
    float roll = -atan2(2*(-real*i+j*k), 1-2*(i*i+j*j));
    if (roll < 0)
        roll += 2*M_PI;

    return EulerAngles(ra, de, roll);
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

/// Whether the quaternion is a unit quaternion. All quaternions representing rotations should be units.
bool Quaternion::IsUnit(float tolerance) const {
    return abs(i*i+j*j+k*k+real*real - 1) < tolerance;
}

/**
 * Ensure that the quaternion's real part is nonnegative. Does not change the rotation represented by the quaternion.
 */
Quaternion Quaternion::Canonicalize() const {
    if (real >= 0) {
        return *this;
    }

    return Quaternion(-real, -i, -j, -k);
}

/// Convert from right ascension & declination to a 3d point on the unit sphere.
Vec3 SphericalToSpatial(float ra, float de) {
    return {
        cos(ra)*cos(de),
        sin(ra)*cos(de),
        sin(de),
    };
}

/// Convert from a 3d point on the unit sphere to right ascension & declination.
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

std::ostream& operator<< (std::ostream& output, const Vec3 &vec){
    output << "x: " << vec.x << "   y: " << vec.y << "  z: " << vec.z;
    return output;
}

/// The square of the magnitude
float Vec3::MagnitudeSq() const {
    return x*x+y*y+z*z;
}

/// The square of the magnitude
float Vec2::MagnitudeSq() const {
    return x*x+y*y;
}

float Vec3::Magnitude() const {
    return sqrt(MagnitudeSq());
}

float Vec2::Magnitude() const {
    return sqrt(MagnitudeSq());
}

/// Create a vector pointing in the same direction with magnitude 1
Vec3 Vec3::Normalize() const {
    float mag = Magnitude();
    return {
        x/mag, y/mag, z/mag,
    };
}

/// Dot product
float Vec3::operator*(const Vec3 &other) const {
    return x*other.x + y*other.y + z*other.z;
}

/// Dot product
Vec2 Vec2::operator*(const float &other) const {
    return { x*other, y*other };
}

/// Vector-Scalar multiplication
Vec3 Vec3::operator*(const float &other) const {
    return { x*other, y*other, z*other };
}

/// Usual vector addition
Vec2 Vec2::operator+(const Vec2 &other) const {
    return {x + other.x, y + other.y };
}

/// Usual vector addition
Vec2 Vec2::operator-(const Vec2 &other) const {
    return { x - other.x, y - other.y };
}

/// Usual vector subtraction
Vec3 Vec3::operator-(const Vec3 &other) const {
    return { x - other.x, y - other.y, z - other.z };
}

Vec3 Vec3::operator+(const Vec3 &other) const{
    return {x + other.x, y + other.y, z + other.z};
}

/// Usual vector cross product
Vec3 Vec3::CrossProduct(const Vec3 &other) const {
    return {
        y*other.z - z*other.y,
        -(x*other.z - z*other.x),
        x*other.y - y*other.x,
    };
}

/// The outer product of two vectors
Mat3 Vec3::OuterProduct(const Vec3 &other) const {
    return {
        x*other.x, x*other.y, x*other.z,
        y*other.x, y*other.y, y*other.z,
        z*other.x, z*other.y, z*other.z
    };
}

/// Vector-matrix multiplication, where the vector is transposed
Vec3 Vec3::operator*(const Mat3 &other) const {
    return {
        x*other.At(0,0) + y*other.At(0,1) + z*other.At(0,2),
        x*other.At(1,0) + y*other.At(1,1) + z*other.At(1,2),
        x*other.At(2,0) + y*other.At(2,1) + z*other.At(2,2),
    };
}

/// Access the i,j-th element of the matrix
float Mat3::At(int i, int j) const {
    return x[3*i+j];
}

/// Get the column at index j
Vec3 Mat3::Column(int j) const {
    return {At(0,j), At(1,j), At(2,j)};
}

/// Get the row at index i
Vec3 Mat3::Row(int i) const {
    return {At(i,0), At(i,1), At(i,2)};
}

/// Normal matrix addition
Mat3 Mat3::operator+(const Mat3 &other) const {
    return {
        At(0,0)+other.At(0,0), At(0,1)+other.At(0,1), At(0,2)+other.At(0,2),
        At(1,0)+other.At(1,0), At(1,1)+other.At(1,1), At(1,2)+other.At(1,2),
        At(2,0)+other.At(2,0), At(2,1)+other.At(2,1), At(2,2)+other.At(2,2)
    };
}

/// Naive matrix multiplication.
Mat3 Mat3::operator*(const Mat3 &other) const {
#define _MATMUL_ENTRY(row, col) At(row,0)*other.At(0,col) + At(row,1)*other.At(1,col) + At(row,2)*other.At(2,col)
    return {
        _MATMUL_ENTRY(0,0), _MATMUL_ENTRY(0,1), _MATMUL_ENTRY(0,2),
        _MATMUL_ENTRY(1,0), _MATMUL_ENTRY(1,1), _MATMUL_ENTRY(1,2),
        _MATMUL_ENTRY(2,0), _MATMUL_ENTRY(2,1), _MATMUL_ENTRY(2,2),
    };
#undef _MATMUL_ENTRY
}

/// Matrix-Vector multiplication
Vec3 Mat3::operator*(const Vec3 &vec) const {
    return {
        vec.x*At(0,0) + vec.y*At(0,1) + vec.z*At(0,2),
        vec.x*At(1,0) + vec.y*At(1,1) + vec.z*At(1,2),
        vec.x*At(2,0) + vec.y*At(2,1) + vec.z*At(2,2),
    };
}

/// Matrix-Scalar multiplication
Mat3 Mat3::operator*(const float &s) const {
    return {
        s*At(0,0), s*At(0,1), s*At(0,2),
        s*At(1,0), s*At(1,1), s*At(1,2),
        s*At(2,0), s*At(2,1), s*At(2,2)
    };
}

/// Transpose of a matrix
Mat3 Mat3::Transpose() const {
    return {
        At(0,0), At(1,0), At(2,0),
        At(0,1), At(1,1), At(2,1),
        At(0,2), At(1,2), At(2,2),
    };
}

/// Trace of a matrix
float Mat3::Trace() const {
    return At(0,0) + At(1,1) + At(2,2);
}

/// Determinant of a matrix
float Mat3::Det() const {
    return (At(0,0) * (At(1,1)*At(2,2) - At(2,1)*At(1,2))) - (At(0,1) * (At(1,0)*At(2,2) - At(2,0)*At(1,2))) + (At(0,2) * (At(1,0)*At(2,1) - At(2,0)*At(1,1)));
}

/// Inverse of a matrix
Mat3 Mat3::Inverse() const {
    // https://ardoris.wordpress.com/2008/07/18/general-formula-for-the-inverse-of-a-3x3-matrix/
    float scalar = 1 / Det();

    Mat3 res = {
        At(1,1)*At(2,2) - At(1,2)*At(2,1), At(0,2)*At(2,1) - At(0,1)*At(2,2), At(0,1)*At(1,2) - At(0,2)*At(1,1),
        At(1,2)*At(2,0) - At(1,0)*At(2,2), At(0,0)*At(2,2) - At(0,2)*At(2,0), At(0,2)*At(1,0) - At(0,0)*At(1,2),
        At(1,0)*At(2,1) - At(1,1)*At(2,0), At(0,1)*At(2,0) - At(0,0)*At(2,1), At(0,0)*At(1,1) - At(0,1)*At(1,0)
    };

    return res * scalar;
}

/// 3x3 identity matrix
 const Mat3 kIdentityMat3 = {1,0,0,
                     0,1,0,
                     0,0,1};

Attitude::Attitude(const Quaternion &quat)
    : quaternion(quat), type(QuaternionType) {}

Attitude::Attitude(const Mat3 &matrix)
    : dcm(matrix), type(DCMType) {}

/// Convert a quaternion to a rotation matrix (Direction Cosine Matrix)
Mat3 QuaternionToDCM(const Quaternion &quat) {
    Vec3 x = quat.Rotate({1, 0, 0});
    Vec3 y = quat.Rotate({0, 1, 0});
    Vec3 z = quat.Rotate({0, 0, 1});
    return {
        x.x, y.x, z.x,
        x.y, y.y, z.y,
        x.z, y.z, z.z,
    };
}

/// Convert a rotation matrix (Direction Cosine Matrix) to a quaternion representing the same rotation.
Quaternion DCMToQuaternion(const Mat3 &dcm) {
    // Make a quaternion that rotates the reference frame X-axis into the dcm's X-axis, just like
    // the DCM itself does
    Vec3 oldXAxis = Vec3({1, 0, 0});
    Vec3 newXAxis = dcm.Column(0); // this is where oldXAxis is mapped to
    assert(abs(newXAxis.Magnitude()-1) < 0.001);
    Vec3 xAlignAxis = oldXAxis.CrossProduct(newXAxis).Normalize();
    float xAlignAngle = AngleUnit(oldXAxis, newXAxis);
    Quaternion xAlign(xAlignAxis, xAlignAngle);

    // Make a quaternion that will rotate the Y-axis into place
    Vec3 oldYAxis = xAlign.Rotate({0, 1, 0});
    Vec3 newYAxis = dcm.Column(1);
    // we still need to take the cross product, because acos returns a value in [0,pi], and thus we
    // need to know which direction to rotate before we rotate. We do this by checking if the cross
    // product of old and new y axes is in the same direction as the new X axis.
    bool rotateClockwise = oldYAxis.CrossProduct(newYAxis) * newXAxis > 0; // * is dot product
    Quaternion yAlign({1, 0, 0}, AngleUnit(oldYAxis, newYAxis) * (rotateClockwise ? 1 : -1));

    // We're done! There's no need to worry about the Z-axis because the handed-ness of the
    // coordinate system is always preserved, which means the Z-axis is uniquely determined as the
    // cross product of the X- and Y-axes. This goes to show that DCMs store redundant information.

    // we want the y alignment to have memory of X, which means we put its multiplication on the
    // right
    return xAlign*yAlign;
}

/// Get the quaternion representing the attitude, converting from whatever format is stored.
Quaternion Attitude::GetQuaternion() const {
    switch (type) {
        case QuaternionType:
            return quaternion;
        case DCMType:
            return DCMToQuaternion(dcm);
        default:
            assert(false);
    }
}

/// Get the rotation matrix (direction cosine matrix) representing the attitude, converting from whatever format is stored.
Mat3 Attitude::GetDCM() const {
    switch (type) {
        case DCMType:
            return dcm;
        case QuaternionType:
            return QuaternionToDCM(quaternion);
        default:
            assert(false);
    }
}

/// Convert a vector from the reference frame to the body frame.
Vec3 Attitude::Rotate(const Vec3 &vec) const {
    switch (type) {
        case DCMType:
            return dcm*vec;
        case QuaternionType:
            return quaternion.Rotate(vec);
        default:
            assert(false);
    }
}

/// Get the euler angles from the attitude, converting from whatever format is stored.
EulerAngles Attitude::ToSpherical() const {
    switch (type) {
        case DCMType:
            return GetQuaternion().ToSpherical();
        case QuaternionType:
            return quaternion.ToSpherical();
        default:
            assert(false);
    }
}

/// The length that a Vec3 will take up when serialized
long SerializeLengthVec3() {
    return sizeof(float)*3;
}

/// Serialize a Vec3 to buffer. Takes up space according to SerializeLengthVec3
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

/// Calculate the inner angle, in radians, between two vectors.
float Angle(const Vec3 &vec1, const Vec3 &vec2) {
    return AngleUnit(vec1.Normalize(), vec2.Normalize());
}

/**
 * Calculate the inner angle, in radians, between two /unit/ vectors.
 * Slightly faster than Angle()
 * @warn If the vectors are not already unit vectors, will return the wrong result!
 */
float AngleUnit(const Vec3 &vec1, const Vec3 &vec2) {
    float dot = vec1*vec2;
    // TODO: we shouldn't need this nonsense, right? how come acos sometimes gives nan?
    return dot >= 1 ? 0 : dot <= -1 ? M_PI-0.0000001 : acos(dot);
}

/// The distance between two vectors, according to the usual distance formula.
float Distance(const Vec2 &v1, const Vec2 &v2) {
    return sqrt(pow(v1.x-v2.x, 2) + pow(v1.y-v2.y, 2));
}


}
