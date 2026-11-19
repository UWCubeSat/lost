#ifndef CENTROID_H
#define CENTROID_H

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

/**
 * Centroids stars near predicted positions, e.g. where the previous frames say each star will be.
 * For each prediction: estimates the background and noise from the border of a square window around
 * it, seeds at the brightest pixel within gateRadius of the prediction, flood-fills the connected
 * pixels above threshold, and takes their background-subtracted center of gravity. Rejects a star
 * that touches the window or image edge (its centroid would be biased), or that ends up farther than
 * gateRadius from its prediction.
 *
 * Then checks the matches against each other: a star reached from two predictions is ambiguous and
 * dropped, and since the prediction error of a whole frame is almost a common shift (attitude
 * error), stars whose offset from their prediction is more than outlierRadius from the median
 * offset are dropped too.
 */
class WindowedCenterOfGravity {
public:
    /**
     * @param gateRadius Largest expected distance, in pixels, between a star and its prediction.
     * @param maxStarRadius Largest expected star radius in pixels. The window around each
     * prediction extends this much beyond the gate, so that whole stars fit in it.
     * Go finds nothing if either radius is negative or more than kMaxRadius.
     * @param outlierRadius See class description, in pixels.
     */
    WindowedCenterOfGravity(int gateRadius, int maxStarRadius, decimal outlierRadius)
        : gateRadius(gateRadius), maxStarRadius(maxStarRadius), outlierRadius(outlierRadius) { };

    /// Largest gateRadius and maxStarRadius accepted, in pixels.
    static const int kMaxRadius = 1024;

    /**
     * @param predictions Predicted star positions in pixel coordinates (same as Star::position).
     * Predictions outside the image, or not finite, are skipped.
     * @param[out] matchedPredictions For each returned star, the index of its prediction. May be null.
     * @return The stars found, in increasing order of prediction index.
     */
    Stars Go(const unsigned char *image, int imageWidth, int imageHeight,
             const std::vector<Vec2> &predictions, std::vector<int> *matchedPredictions) const;

private:
    int gateRadius;
    int maxStarRadius;
    decimal outlierRadius;
};

}

#endif
