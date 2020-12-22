#ifndef ATTITUDE_UTILS_H
#define ATTITUDE_UTILS_H

#include "camera.hpp"

namespace lost {

// At first, I wanted to have two separate Attitude classes, one storing Euler angles and converting
// to Quaterinon, and another storing as Quaternion and converting to Euler. But abstract classes
// make everything more annoying, because you need vectors of pointers...ugh!

class EulerAttitude {
public:
    long x;
    long y;
    long z;
};

class Attitude {
public:
    Attitude() = default; // I guess this lets you call Attitude on an attitude object?
    Attitude(EulerAttitude ea) { }; // TODO: implement
    Attitude(long real, long i, long j, long k)
        : real(real), i(i), j(j), k(k) { };

    EulerAttitude Euler() { return { 0 }; }; // TODO

    long real;
    long i;
    long j;
    long k;
};

// TODO: quaternion and euler angle conversion, conversion between ascension/declination to rec9tu

}

#endif
