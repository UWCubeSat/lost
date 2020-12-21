#ifndef CAMERA_H
#define CAMERA_H

namespace lost {

class Camera {
public:
    long horizontalFov; // millionths of a degree
    long verticalFov;   // millionths of a degree
    // TODO: distortion
};

}

#endif
