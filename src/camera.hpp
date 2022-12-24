#ifndef CAMERA_H
#define CAMERA_H

#include "attitude-utils.hpp"

namespace lost {

/// A full description of a camera. Enough information to reconstruct the camera matrix and then some.
class Camera {
public:
    Camera(const Camera &) = default;

    /**
     * @param xCenter,yCenter The "principal point" of the camera. In an ideal camera, just half the resolution, but physical cameras often have a bit of offset.
     */
    Camera(float focalLength,
           float xCenter, float yCenter,
           int xResolution, int yResolution)
        : focalLength(focalLength),
          xCenter(xCenter), yCenter(yCenter),
          xResolution(xResolution), yResolution(yResolution) {};

    Camera(float focalLength, int xResolution, int yResolution)
        : Camera(focalLength,
                 xResolution / (float) 2.0, yResolution / (float) 2.0,
                 xResolution, yResolution) {};

    Vec2 SpatialToCamera(const Vec3 &) const;
    Vec3 CameraToSpatial(const Vec2 &) const;

    Vec3 CameraToSpatialFov(const Vec2 &) const;

    // converts from a 2d point in the camera sensor to right ascension and declination relative to
    // the center of the camera.
    // void CoordinateAngles(const Vec2 &vector, float *ra, float *de) const;

    bool InSensor(const Vec2 &vector) const;

    /// Width of the sensor in pixels
    int XResolution() const { return xResolution; };
    /// Height of the sensor in pixels
    int YResolution() const { return yResolution; };
    /// Focal length in pixels
    float FocalLength() const { return focalLength; };
    /// Horizontal field of view in radians
    float Fov() const; // in radians

    void SetFocalLength(float focalLength) { this->focalLength = focalLength; }

private:
    // TODO: distortion
    float focalLength;
    float xCenter; float yCenter;
    int xResolution; int yResolution;
};

float FovToFocalLength(float xFov, float xResolution);

}

#endif
