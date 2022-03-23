#ifndef CAMERA_H
#define CAMERA_H

#include "attitude-utils.hpp"

// Ben's Notes: Questions for Mark or Allen
//      1. Is possible you could explain each argument for the Camera object in terms of where they come from?
//      2. Wait, does my image down below actually represent the 3D space described in the class down below?

namespace lost {

    //Ben's Notes on Focal Length:

        // Imagine that we had an x-axis and y-axis, where the y-axis was the convex lens and the x-axis
        // ran straight through the middle:

        //          (This below is the camera itself)      (y-axis or convex lens of camera)
    //   _____________________________________________o
    //   --------------------------------------------` `
    //                                               ` `
    //                                               ` `
    //                                               ` `
    //                                               ` `
    //                                               ` `
    //               (*)*****************************` `
    // (image sensor) |   *         (Focal Length)   ` `    *
    //                |       *   |------------------`|`          *     (Focal Point, where light ray hits the x-axis)
    // (x-axis)  .....|..........(*).................`.`................(*)................................. (x-axis)
    //              (Focal point,    *               `|`-----------------|    *                     |
    //               where light ray     *           ` ` (Focal length)            *                |
    //               hits the x-axis)        *       ` `                                 *          | (Object you're taking a picture of)
    //                                           *   ` `                                       *    |
    //                                               `*`:******************************************(*) (B)
    //                                               ` `           (Line above is a ray of light coming into
    //                                               ` `             the camera lens from the outside world into the camera
    //                                               ` `               (right to left))
    //                                               ` `
    // ______________________________________________. .
    // -----------------------------------------------o
    //                                               (y-axis or convex lens of camera)

        // The focal length is the distance from the convex lens to the x

        // Common terminology:
        //    - Convex lens: a lens that is thin at the edges and gets thicker towards the center.
        //                   It focuses light that comes in at a straight perpendicualr angle to the lens
        //                   (like the bottom beam of light shown

class Camera { // TODO: Ben | Put the distortCoeffs into some kind of global struct???
public:
    Camera(const Camera &) = default;
    // Takes focal lengths in pixels
    Camera(float focalLength,
           float xCenter, float yCenter,
           int xResolution, int yResolution,
           float k1, float k2, float k3, float k4, float k5, float k6,
           float p1, float p2)
        : focalLength(focalLength),
          xCenter(xCenter), yCenter(yCenter),
          xResolution(xResolution), yResolution(yResolution),
          k1(k1), k2(k2), k3(k3), k4(k4), k5(k5), k6(k6),
          p1(p1), p2(p2) { };
    Camera(float focalLength, int xResolution, int yResolution,
           float k1, float k2, float k3, float k4, float k5, float k6,
           float p1, float p2)
        : Camera(focalLength,
                 xResolution/(float)2.0, yResolution/(float)2.0,
                 xResolution, yResolution,
                 k1, k2, k3, k4, k5, k6,
                 p1, p2) { };

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
    float Fov() const;
    void SetFocalLength(float focalLength) { this->focalLength = focalLength; }
    float K1() const { return this->k1; };
    float K2() const { return this->k2; };
    float K3() const { return this->k3; };
    float K4() const { return this->k4; };
    float K5() const { return this->k5; };
    float K6() const { return this->k6; };
    float P1() const { return this->p1; };
    float P2() const { return this->p2; };

//    DistortCoeffs* distortCoeffs() const { return this.distortCoeffs; } // TODO | Ben Kosa | Take care of this.

private:
    // TODO | Ben Kosa
    //  Find out if the distortion coeffs are integers or floats!
    //
    float focalLength;
    float xCenter; float yCenter;
    int xResolution; int yResolution;
    float k1, k2, k3, k4, k5, k6, p1, p2;
//    struct DistortCoeffs distortCoeffs; // TODO | Ben Kosa | Take care of this.

};

float FovToFocalLength(float xFov, float xResolution);

// Our undistortion function works with up to six radial distortion coefficents (k1 through k6)
// and 2 tangental coeffs (p1 and p2)
// TODO | Ben: Check if this is
struct DistortCoeffs {
    float k1;
    float k2;
    float k3;
    float k4;
    float k5;
    float k6;
//        float p1;
//        float p2;
};

}

#endif
