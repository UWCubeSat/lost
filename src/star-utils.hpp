#ifndef STAR_UTILS_H
#define STAR_UTILS_H

#include <vector>

namespace lost {

// Different centroiding algorithms may return different amounts of data. Because there won't ever
// be /that/ many different centroids in a single image, it's acceptable to waste some space on
// un-used fields. Set unused fields to negative
class Star {
public:
    Star(float x, float y, float radiusX, float radiusY, long magnitude) :
        x(x), y(y), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude), spherePhi(-1.0) { };
    Star(float x, float y, float radiusX) : Star(x, y, radiusX, radiusX, 0) { };
    Star() : Star(0.0, 0.0, 0.0) { };

    float x; // pixels*10^6
    float y; // pixels*10^6
    float radiusX;
    float radiusY;  // if omitted, but x is present, assume circular.
    long  magnitude; // some relative number
    // eccentricity?

    // viewport-relative spherical coordinates. Assumes the center of the viewport is theta = 90deg and phi=0deg
    float sphereTheta;
    float spherePhi; // -1.0 indicates spherical coordinates not set

    int identifiedBscIndex; // if the star has been identified, refers to the index in the BSC,
                            // starting at zero
};

typedef std::vector<Star> Stars;

float StarDistancePixels(Star one, Star two);

}

#endif
