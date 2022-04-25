#ifndef CAMERA_H
#define CAMERA_H

#include "attitude-utils.hpp"

namespace lost {

class Camera {
public:
    Camera(const Camera &) = default;
    // Takes focal lengths in pixels
    Camera(float focalLength, float fX, float fY,
           float xCenter, float yCenter,
           float k1, float k2, float k3, float p1, float p2, int iters,
           int xResolution, int yResolution)
        : focalLength(focalLength), focalLengthX(fX), focalLengthY(fY), // just syntax to say this.focalLengthX = fX;
          xCenter(xCenter), yCenter(yCenter),
          k1(k1), k2(k2), k3(k3), p1(p1), p2(p2),
          xResolution(xResolution), yResolution(yResolution) { };
    Camera(float focalLength, int xResolution, int yResolution)
        : Camera(focalLength, focalLength, focalLength,
                 xResolution/(float)2.0, yResolution/(float)2.0,
                 0, 0, 0, 0, 0, 0,
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

    int XResolution() const { return xResolution; };
    int YResolution() const { return yResolution; };
    float FocalLength() const { return focalLength; };

    float FocalLengthX() const {return focalLengthX; };
    float FocalLengthY() const {return focalLengthY; };

    float XCenter() const {return xCenter; };
    float YCenter() const {return yCenter; };

    float K1() const {return k1; };
    float K2() const {return k2; };
    float K3() const {return k3; };
    float P1() const {return p1; };
    float P2() const {return p2; };

    int ITERS() const{return iters; };


    float Fov() const;
    void SetFocalLength(float focalLength) { this->focalLength = focalLength; }

private:
    float focalLength;
    float focalLengthX; float focalLengthY;
    float xCenter; float yCenter;
    float k1, k2, k3, p1, p2;
    float iters;
    int xResolution; int yResolution;
};

float FovToFocalLength(float xFov, float xResolution);

}

#endif
