#ifndef CENTROID_H
#define CENTROID_H

#include <iostream>
#include <vector>

namespace lost {

// Different centroiding algorithms may return different amounts of data. Because there won't ever
// be /that/ many different centroids in a single image, it's acceptable to waste some space on
// un-used fields. Set unused fields to negative
class CentroidStar {
public:
    CentroidStar(float x, float y, float radiusX, float radiusY, long magnitude) :
        x(x), y(y), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude) { };
    CentroidStar(float x, float y, float radiusX) : CentroidStar(x, y, radiusX, radiusX, 0) { };

    float x; // pixels*10^6
    float y; // pixels*10^6
    float radiusX;
    float radiusY;  // if omitted, but x is present, assume circular.
    long  magnitude; // some relative number
    // eccentricity?
};

typedef std::vector<CentroidStar> Centroids;

class CentroidAlgorithm {
public:
    virtual Centroids Go(unsigned char *image, int imageWidth, int imageHeight) const = 0;
    virtual ~CentroidAlgorithm() { };
};

class DummyCentroidAlgorithm: public CentroidAlgorithm {
public:
    DummyCentroidAlgorithm(int numStars) : numStars(numStars) { };
    Centroids Go(unsigned char *image, int imageWidth, int imageHeight) const override;
private:
    int numStars;
};

class CenterOfGravityAlgorithm : public CentroidAlgorithm {
    public:
        CenterOfGravityAlgorithm() { };
        Centroids Go(unsigned char *image, int imageWidth, int imageHeight) const override;
};

}

#endif
