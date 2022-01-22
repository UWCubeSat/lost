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
        


class Camera { // TODO: Ben | Put the intrinic parameters AND the dist coeffs for the camera here
public:
    Camera(const Camera &) = default;
    // Takes focal lengths in pixels
    Camera(float focalLength,
           float xCenter, float yCenter,
           int xResolution, int yResolution)
        : focalLength(focalLength),
          xCenter(xCenter), yCenter(yCenter),
          xResolution(xResolution), yResolution(yResolution) { };
    Camera(float focalLength, int xResolution, int yResolution)
        : Camera(focalLength,
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

    int XResolution() const { return xResolution; };
    int YResolution() const { return yResolution; };
    float FocalLength() const { return focalLength; };
    float Fov() const;
    void SetFocalLength(float focalLength) { this->focalLength = focalLength; }

private:
    // TODO: distortion
    // Variables to add:
    //      Camera Parameters:
    //      Distortion Coefficient:
    //
    float focalLength;
    float xCenter; float yCenter;
    int xResolution; int yResolution;

    //TODO: Figure out how many distortion coeffs we need (k1 - k6) and if we want to include p1 and p2 (radial distortion coeffs).
    /**
     * <p> Our undistortion function works with up to six radial distortion coefficents. </p>
     * <p> These radial distortion coefficents are found here: TODO: put link to coeff finder here <p>
     * <p> Hardcode these into our distCoeffs enum down below.
     */

    enum distCoeffs {k1 = 0, k2 = 0, k3 = 0, k4 = 0, k5 = 0, k6 = 0, p1 = 0, p2 = 0};

};

float FovToFocalLength(float xFov, float xResolution);

}

#endif
