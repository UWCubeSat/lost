#include "camera.hpp"

#include <math.h>
#include <assert.h>

#include "attitude-utils.hpp"
#include "scalar.hpp"

namespace lost {

/**
 * Converts from a 3D point in space to a 2D point on the camera sensor.
 * Assumes that X is the depth direction and that it points away from the center of the sensor, i.e., any vector (x, 0, 0) will be at (xResolution/2, yResolution/2) on the sensor.
 */
Vec2 Camera::SpatialToCamera(const Vec3 &vector) const {
    // can't handle things behind the camera.
    assert(vector.x > 0);
    // TODO: is there any sort of accuracy problem when vector.y and vector.z are small?

    scalar focalFactor = focalLength/vector.x;

    scalar yPixel = vector.y*focalFactor;
    scalar zPixel = vector.z*focalFactor;

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
    scalar xPixel = -vector.x + xCenter;
    scalar yPixel = -vector.y + yCenter;

    return {
        1,
        xPixel / focalLength,
        yPixel / focalLength,
    };
}

/// Returns whether a given pixel is actually in the camera's field of view
bool Camera::InSensor(const Vec2 &vector) const {
    // if vector.x == xResolution, then it is at the leftmost point of the pixel that's "hanging
    // off" the edge of the image, so vector is still in the image.
    return vector.x >= 0 && vector.x <= xResolution
        && vector.y >= 0 && vector.y <= yResolution;
}

scalar FovToFocalLength(scalar xFov, scalar xResolution) {
    return xResolution / SCALAR(2.0) / SCALAR_TAN(xFov/2);
}

scalar FocalLengthToFov(scalar focalLength, scalar xResolution, scalar pixelSize) {
    return SCALAR_ATAN(xResolution/2 * pixelSize / focalLength) * 2;
}

scalar Camera::Fov() const {
    return FocalLengthToFov(focalLength, xResolution, 1.0);
}

}
