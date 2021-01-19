#include "camera.hpp"
#include "attitude-utils.hpp"

#include <math.h>
#include <assert.h>

namespace lost {

Vec2 Camera::ConvertCoordinates(const Vec3 &vector) const {
    // can't handle things behind the camera.
    assert(vector.x > 0);
    assert(xFov > 0 && xFov < M_PI);
    // TODO: is there any sort of accuracy problem when vector.y and vector.z are small?
    float yTangent = vector.y/vector.x;
    float zTangent = vector.z/vector.x;
    // these pixels are measured using 0,0 as the /center/

    // TODO: analyze off-by-one errors (is fov for the center of the pixels, in which case we should
    // use xResolution-1??)
    float yPixel = yTangent/tan(xFov/2)*(xResolution-1)/2;
    float zPixel = zTangent/tan(xFov/2)*(xResolution-1)/2;

    // now convert to using 0,0 as the top left, as usual. TODO: analyze off-by-one errors here too.
    // Right now, if we have 5x5, the center becomes (2,2), and if we have 4x4, the center becomes
    // (1.5, 1.5). I think this is right?
    return { -yPixel + (xResolution-1)/2.0f, -zPixel + (yResolution-1)/2.0f };
}

void Camera::CoordinateAngles(const Vec2 &vector, float *ra, float *de) const {
    // TODO: off-by-one with xResolution - 1?
    // TODO: minimize floating point error?
    *ra = atan((xResolution/2.0-vector.x)/(xResolution/2.0/tan(xFov/2.0)));
    *de = atan((yResolution/2.0-vector.y)/(xResolution/2.0/tan(xFov/2.0)));
}

bool Camera::InSensor(const Vec2 &vector) const {
    return vector.x >= 0 && vector.x <= xResolution
        && vector.y >= 0 && vector.y <= yResolution;
}

}
