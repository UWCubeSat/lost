#include "camera.hpp"
#include "attitude-utils.hpp"

#include <math.h>
#include <assert.h>

namespace lost {

Vec2 Camera::SpatialToCamera(const Vec3 &vector) const {
    // can't handle things behind the camera.
    assert(vector.x > 0);
    // TODO: is there any sort of accuracy problem when vector.y and vector.z are small?

    float yTangent = vector.y/vector.x;
    float zTangent = vector.z/vector.x;

    float yPixel = yTangent*xFocalLength;
    float zPixel = zTangent*yFocalLength;

    return { -yPixel + xCenter, -zPixel + yCenter };
}

// we'll just place the points at 1 unit away from the pinhole (x=1)
Vec3 Camera::CameraToSpatial(const Vec2 &vector) const {
    assert(InSensor(vector));

    // isn't it interesting: To convert from center-based to left-corner-based coordinates is the
    // same formula; f(x)=f^{-1}(x) !
    float xPixel = -vector.x + xCenter;
    float yPixel = -vector.y + yCenter;

    return {
        1,
        xPixel / xFocalLength,
        yPixel / yFocalLength,
    };
}

// TODO: reimplement or determine we don't need it
// void Camera::CoordinateAngles(const Vec2 &vector, float *ra, float *de) const {
//     // TODO: off-by-one with xResolution - 1?
//     // TODO: minimize floating point error?

//     // *ra = atan((xResolution/2.0-vector.x)/(xResolution/2.0/tan(xFov/2.0)));
//     // *de = atan((yResolution/2.0-vector.y)/(xResolution/2.0/tan(xFov/2.0)));
// }

bool Camera::InSensor(const Vec2 &vector) const {
    // if vector.x == xResolution, then it is at the leftmost point of the pixel that's "hanging
    // off" the edge of the image, so vector is still in the image.
    return vector.x >= 0 && vector.x <= xResolution
        && vector.y >= 0 && vector.y <= yResolution;
}

int Camera::GetXResolution() const {
    return xResolution;
}

int Camera::GetYResolution() const {
    return yResolution;
}

}
