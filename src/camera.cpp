#include "camera.hpp"

#include <assert.h>
#include <math.h>

#include "attitude-utils.hpp"

namespace lost {

Vec2 Camera::SpatialToCamera(const Vec3 &vector) const {
    // can't handle things behind the camera.
    assert(vector.x > 0);
    // TODO: is there any sort of accuracy problem when vector.y and vector.z
    // are small?

    float focalFactor = focalLength / vector.x;

    float yPixel = vector.y * focalFactor;
    float zPixel = vector.z * focalFactor;

    return {-yPixel + xCenter, -zPixel + yCenter};
}

// we'll just place the points at 1 unit away from the pinhole (x=1).
// Note: Other functions rely on the fact that vectors are placed one unit away.
// Don't change this behavior!
Vec3 Camera::CameraToSpatial(const Vec2 &vector) const {
    assert(InSensor(vector));

    // isn't it interesting: To convert from center-based to left-corner-based
    // coordinates is the same formula; f(x)=f^{-1}(x) !
    float xPixel = -vector.x + xCenter;
    float yPixel = -vector.y + yCenter;

    return {
        1,
        xPixel / focalLength,
        yPixel / focalLength,
    };
}

bool Camera::InSensor(const Vec2 &vector) const {
    // if vector.x == xResolution, then it is at the leftmost point of the pixel
    // that's "hanging off" the edge of the image, so vector is still in the
    // image.
    return vector.x >= 0 && vector.x <= xResolution && vector.y >= 0 &&
           vector.y <= yResolution;
}

float FovToFocalLength(float xFov, float xResolution) {
    return xResolution / 2.0f / tan(xFov / 2);
}

float FocalLengthToFov(float focalLength, float xResolution, float pixelSize) {
    return atan(xResolution / 2 * pixelSize / focalLength) * 2;
}

float Camera::Fov() const {
    return FocalLengthToFov(focalLength, xResolution, 1.0);
}

}  // namespace lost
