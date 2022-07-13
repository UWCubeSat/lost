#ifndef CAMERA_H
#define CAMERA_H

#include "attitude-utils.hpp"

namespace lost {

/**
 * @brief
 * @details
 */
class Camera {
public:
    /**
     * @brief
     * @param
     */
    Camera(const Camera &) = default;

    /**
     * @brief
     * @param focalLength Camera focal lengths in pixels
     * @param xCenter
     * @param yCenter
     * @param xResolution
     * @param yResolution
     */
    Camera(float focalLength,
           float xCenter, float yCenter,
           int xResolution, int yResolution)
        : focalLength(focalLength),
          xCenter(xCenter), yCenter(yCenter),
          xResolution(xResolution), yResolution(yResolution) {};

    /**
     * @brief
     * @param focalLength
     * @param xResolution
     * @param yResolution
     */
    Camera(float focalLength, int xResolution, int yResolution)
        : Camera(focalLength,
                 xResolution / (float) 2.0, yResolution / (float) 2.0,
                 xResolution, yResolution) {};

    Vec2 SpatialToCamera(const Vec3 &) const;
    Vec3 CameraToSpatial(const Vec2 &) const;

    // converts from a 2d point in the camera sensor to right ascension and declination relative to
    // the center of the camera.
    // void CoordinateAngles(const Vec2 &vector, float *ra, float *de) const;

    bool InSensor(const Vec2 &vector) const;

    /**
     * @brief
     * @return
     */
    int XResolution() const { return xResolution; };

    /**
     * @brief
     * @return
     */
    int YResolution() const { return yResolution; };

    /**
     * @brief
     * @return
     */
    float FocalLength() const { return focalLength; };

    /**
     * @brief
     * @return
     */
    float Fov() const;

    /**
     * @brief
     * @param focalLength
     */
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
