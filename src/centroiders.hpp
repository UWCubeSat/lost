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

/**
 * Centroids stars near predicted positions. For each prediction: takes background and noise from
 * the border of a square window, seeds at the brightest pixel within gateRadius, flood-fills the
 * connected pixels at or above threshold, and takes their background-subtracted center of gravity.
 * Only those pixels count, so the faint wings of a star are left out.
 *
 * Drops a star that touches the window or image edge, that lies farther than gateRadius from its
 * prediction, or whose pixel bounding box overlaps another match's. When at least three stars
 * are left, also drops any whose offset from its prediction is more than outlierRadius from the
 * frame's median offset.
 *
 * Working memory, not measured on the flight computer: one byte per pixel of a window, which is
 * (2 * (gateRadius + maxStarRadius) + 1) pixels square or the image if that is smaller, allocated
 * once per call; up to 8 more bytes per pixel of the largest star, for the flood fill's stack;
 * and a few hundred bytes per prediction. Go also compares every pair of matches. So keep the
 * radii and the number of predictions to what the frame needs.
 */
class WindowedCenterOfGravity {
public:
    /**
     * @param gateRadius Largest expected distance between a star and its prediction, in pixels.
     * @param maxStarRadius Largest expected star radius, in pixels; the window extends this far past the gate.
     * @param outlierRadius Largest accepted distance of a match's offset from the median offset, in pixels.
     * Go finds nothing if it is negative or NaN.
     */
    WindowedCenterOfGravity(int gateRadius, int maxStarRadius, decimal outlierRadius)
        : gateRadius(gateRadius), maxStarRadius(maxStarRadius), outlierRadius(outlierRadius) { };

    /// Largest gateRadius and maxStarRadius accepted, in pixels. Outside 0..kMaxRadius, Go finds nothing.
    static const int kMaxRadius = 1024;
    /// Largest image width and height accepted, in pixels. Go finds nothing in a larger image, or
    /// in one with more pixels than an int can index.
    static const int kMaxImageSide = 1 << 24;

    /**
     * @param predictions Predicted star positions in pixels (as Star::position). Off-image ones are skipped.
     * @param[out] matchedPredictions The prediction index of each returned star. May be null.
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
