#ifndef ATTITUDE_UTILS_H
#define ATTITUDE_UTILS_H

#include "camera.hpp"

class EulerAngles {
public:
    long x;
    long y;
    long z;
};

class QuaternionAngles {
    long real;
    long i;
    long j;
    long k;
};

class Attitude {
public:
    virtual EulerAngles Euler() const;
    virtual QuaternionAngles Quaternion() const;
};

// TODO: quaternion and euler angle conversion, conversion between ascension/declination to rec9tu

#endif
