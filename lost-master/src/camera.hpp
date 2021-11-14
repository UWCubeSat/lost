#ifndef CAMERA_H
#define CAMERA_H

#include "attitude-utils.hpp"

namespace lost {

class Camera {
public:
    // Takes focal lengths in pixels
    Camera(float xFocalLength, float yFocalLength,
           float xCenter, float yCenter,
           int xResolution, int yResolution)
        : xFocalLength(xFocalLength), yFocalLength(yFocalLength),
          xCenter(xCenter), yCenter(yCenter),
          xResolution(xResolution), yResolution(yResolution) { };
    Camera(float xFocalLength, int xResolution, int yResolution)
        : Camera(xFocalLength, xFocalLength,
                 xResolution/(float)2.0, yResolution/(float)2.0,
                 xResolution, yResolution) { };

    // Converts from a 3D point in space to a 2D point on the camera sensor. Assumes that X is the
    // depth direction and that it points away from the center of the sensor, i.e., any vector (x,
    // 0, 0) will be at (xResolution/2, yResolution/2) on the sensor.
    Vec2 SpatialToCamera(const Vec3 &) const;
    // Gives /a/ point in 3d space that could correspond to the given vector, using the same
    // coordinate system described for SpatialToCamera. Not all vectors returned by this function
    // will necessarily have the same magnitude.
    Vec3 CameraToSpatial(const Vec2 &) const;
    // converts from a 2d point in the camera sensor to right ascension and declination relative to
    // the center of the camera.
    // void CoordinateAngles(const Vec2 &vector, float *ra, float *de) const;
    // returns whether a given pixel is actually in the camera's field of view
    bool InSensor(const Vec2 &vector) const;

    int GetXResolution() const;
    int GetYResolution() const;

private:
    // TODO: distortion
    float xFocalLength; float yFocalLength;
    float xCenter; float yCenter;
    int xResolution; int yResolution;
};

float FovToFocalLength(float xFov, float xResolution);

}

#endif
