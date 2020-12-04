#ifndef CENTROID_H
#define CENTROID_H

#include <iostream>
#include <vector>

namespace lost {

// Different centroiding algorithms may return different amounts of data. Because there won't ever
// be /that/ many different centroids in a single image, it's acceptable to waste some space on
// un-used fields. Set unused fields to negative
class Star {
public:
    Star(float x, float y, float radiusX, float radiusY, long magnitude) :
        x(x), y(y), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude) { };
    Star(float x, float y, float radiusX) : Star(x, y, radiusX, radiusX, 0) { };

    float x; // pixels*10^6
    float y; // pixels*10^6
    float radiusX;
    float radiusY;  // if omitted, but x is present, assume circular.
    long  magnitude; // some relative number
    // eccentricity?
};

class CentroidComparison {
public:
    CentroidComparison() : meanError(0.0f), numExtraStars(0), numMissingStars(0) { };
    void Print();
    float meanError;       // average distance from actual to expected star
    int numExtraStars;    // stars in actual but not expected. Ideally 0
    int numMissingStars;  // stars is expected but not actual. Ideally 0
    // I would add 99th percentile or something, but the really far away stars should really just
    // count in extra_num
};

class CentroidAlgorithm {
public:
    virtual std::vector<Star> Go(unsigned char *image, int imageWidth, int imageHeight) const = 0;
    virtual ~CentroidAlgorithm() { };
};

class DummyCentroidAlgorithm: public CentroidAlgorithm {
public:
    DummyCentroidAlgorithm(int numStars) : numStars(numStars) { };
    std::vector<Star> Go(unsigned char *image, int imageWidth, int imageHeight) const override;
private:
    int numStars;
};

float StarDistancePixels(Star one, Star two);
// compare two lists of centroids to find differences.
CentroidComparison StarCentroidsCompare(float distanceThreshold, // stars further apart than
                                                                      // are different stars.
                                          std::vector<Star> expected,
                                          std::vector<Star> actual);
void CentroidComparisonPrint(CentroidComparison comparison);

}

#endif
