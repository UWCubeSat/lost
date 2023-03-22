#ifndef CENTROID_H
#define CENTROID_H

#include <iostream>
#include <vector>

#include "star-utils.hpp"

namespace lost {

/// An algorithm that detects the (x,y) coordinates of bright points in an image, called "centroids"
class CentroidAlgorithm {
public:
    /**
     * Actually perform the centroid detection. This is the "main" function of CentroidAlgorithm.
     * @param image An row-major indexed array of grayscale pixels. 0 is black, 255 is white. The total length is \p imageWidth * \p imageHeight
     * @param imageWidth,imageHeight Image dimensions in pixels.
     */
    virtual Stars Go(unsigned char *image, int imageWidth, int imageHeight) const = 0;

    virtual ~CentroidAlgorithm() { };
};

/// A centroid algorithm for debugging that returns random centroids.
class DummyCentroidAlgorithm: public CentroidAlgorithm {
public:
    explicit DummyCentroidAlgorithm(int numStars) : numStars(numStars) { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
private:
    int numStars;
};

/**
 * A simple, fast, and pretty decent centroid algorithm.
 * Simply finds the weighted average of the coordinates of all bright pixels, the weight being proportional to the brightness of the pixel.
 */
class CenterOfGravityAlgorithm : public CentroidAlgorithm {
public:
    CenterOfGravityAlgorithm() { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
};

/**
 * A more complicated centroid algorithm which doesn't perform much better than CenterOfGravityAlgorithm.
 * Iteratively estimates the center of the centroid. Some papers report that it is slightly more precise than CenterOfGravityAlgorithm, but that has not been our experience.
 */
class IterativeWeightedCenterOfGravityAlgorithm : public CentroidAlgorithm {
    public:
        IterativeWeightedCenterOfGravityAlgorithm() { };
        Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
};

}

#endif
