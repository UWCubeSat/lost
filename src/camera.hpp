#ifndef CAMERA_H
#define CAMERA_H

namespace lost {

class Camera {
public:
    Camera(long xFov, long yFov, long xResolution, long yResolution)
        : xFov(xFov), yFov(yFov), xResolution(xResolution), yResolution(yResolution) { };

    long xFov; // millionths of a degree
    long yFov;   // millionths of a degree
    int xResolution;
    int yResolution;
    // TODO: distortion
};

}

#endif
