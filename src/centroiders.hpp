#ifndef CENTROID_H
#define CENTROID_H

#include <iostream>
#include <vector>

#include "star-utils.hpp"

namespace lost {

class CentroidAlgorithm {
public:
    virtual Stars Go(unsigned char *image, int imageWidth, int imageHeight) const = 0;
    virtual ~CentroidAlgorithm() { };
};

class DummyCentroidAlgorithm: public CentroidAlgorithm {
public:
    DummyCentroidAlgorithm(int numStars) : numStars(numStars) { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
private:
    int numStars;
};

class CenterOfGravityAlgorithm : public CentroidAlgorithm {
    public:
        CenterOfGravityAlgorithm(int subdivisions) : subdivisions(subdivisions) { };
        Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
    private:
        int subdivisions;
};

class IterativeWeightedCenterOfGravityAlgorithm : public CentroidAlgorithm {
    public:
        IterativeWeightedCenterOfGravityAlgorithm() { };
        Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
};

}

#endif
