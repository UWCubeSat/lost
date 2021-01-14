#ifndef CAMERA_H
#define CAMERA_H

#include "attitude-utils.hpp"

namespace lost {

class Camera {
public:
    Camera() = default;
    Camera(float xFov, int xResolution, int yResolution)
        : xFov(xFov), xResolution(xResolution), yResolution(yResolution) { };

    // Converts from a 3D point in space to a 2D point on the camera sensor. Assumes that X is the
    // depth direction and that it points away from the center of the sensor, i.e., any vector (x,
    // 0, 0) will be at (xResolution/2, yResolution/2) on the sensor.
    Vec2 ConvertCoordinates(const Vec3 &vector) const;
    // converts from a 2d point in the camera sensor to right ascension and declination relative to
    // the center of the camera.
    void CoordinateAngles(const Vec2 &vector, float *ra, float *de) const;
    // returns whether a given pixel is actually in the camera's field of view
    bool InSensor(const Vec2 &vector) const;

    float xFov;
    int xResolution;
    int yResolution;
    // TODO: distortion
};

}

#endif
