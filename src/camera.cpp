#include "camera.hpp"

#include <math.h>
#include <assert.h>

#include "attitude-utils.hpp"

namespace lost {

/**
 * Converts from a 3D point in space to a 2D point on the camera sensor.
 * Assumes that X is the depth direction and that it points away from the center of the sensor, i.e., any vector (x, 0, 0) will be at (xResolution/2, yResolution/2) on the sensor.
 */
Vec2 Camera::SpatialToCamera(const Vec3 &vector) const {
    // can't handle things behind the camera.
    assert(vector.x > 0);
    // TODO: is there any sort of accuracy problem when vector.y and vector.z are small?

    float focalFactor = focalLength/vector.x;

    float yPixel = vector.y*focalFactor;
    float zPixel = vector.z*focalFactor;

    return { -yPixel + xCenter, -zPixel + yCenter };
}

/**
 * Gives a point in 3d space that could correspond to the given vector, using the same
 * coordinate system described for SpatialToCamera.
 * Not all vectors returned by this function will necessarily have the same magnitude.
 * @return A vector in 3d space corresponding to the given vector, with x-component equal to 1
 * @warning Other functions rely on the fact that returned vectors are placed one unit away (x-component equal to 1). Don't change this behavior!
 */
Vec3 Camera::CameraToSpatial(const Vec2 &vector) const {
    assert(InSensor(vector));

    // isn't it interesting: To convert from center-based to left-corner-based coordinates is the
    // same formula; f(x)=f^{-1}(x) !
    float xPixel = -vector.x + xCenter;
    float yPixel = -vector.y + yCenter;

    return {
        1,
        xPixel / focalLength,
        yPixel / focalLength,
    };
}

Vec3 Camera::CameraToSpatialFov(const Vec2 &vector) const{
    // float fov = DegToRad(this->Fov());
    float fov = this->Fov();
    float centerX = this->XResolution() / 2.0;
    float centerY = this->YResolution() / 2.0;
    float scaleFactor = std::tan(fov / 2) / centerX;

    float jOverI = (centerX - vector.x) * scaleFactor;
    float kOverI = (centerY - vector.y) * scaleFactor;
    float i = 1.0 / std::sqrt(1 + jOverI * jOverI + kOverI * kOverI);
    float j = jOverI * i;
    float k = kOverI * i;
    Vec3 res(i, j, k);

    return res;
}

/// Returns whether a given pixel is actually in the camera's field of view
bool Camera::InSensor(const Vec2 &vector) const {
    // if vector.x == xResolution, then it is at the leftmost point of the pixel that's "hanging
    // off" the edge of the image, so vector is still in the image.
    return vector.x >= 0 && vector.x <= xResolution
        && vector.y >= 0 && vector.y <= yResolution;
}

float FovToFocalLength(float xFov, float xResolution) {
    return xResolution / 2.0f / tan(xFov/2);
}

float FocalLengthToFov(float focalLength, float xResolution, float pixelSize) {
    return atan(xResolution/2 * pixelSize / focalLength) * 2;
}

float Camera::Fov() const {
    return FocalLengthToFov(focalLength, xResolution, 1.0);
}

}
