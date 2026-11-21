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
 * Centroids stars near predicted positions. For each prediction: background and noise from the
 * border of a square window, seed at the brightest pixel within gateRadius, flood fill of the
 * connected pixels at or above threshold, then their background-subtracted center of gravity.
 * Pixels below threshold do not count, so faint wings are left out.
 *
 * Drops a star that touches the window or image edge, is farther than gateRadius from its
 * prediction, or whose bounding box overlaps another match's. With three or more stars left,
 * also drops any whose offset is more than outlierRadius from the median offset.
 *
 * Memory (not measured on the flight computer): one byte per window pixel, where the window is
 * 2 * (gateRadius + maxStarRadius) + 1 pixels square or the image if smaller; up to 8 bytes per
 * pixel of the largest star; a few hundred bytes per prediction. Every pair of matches is
 * compared. Keep the radii and the number of predictions small.
 */
class WindowedCenterOfGravity {
public:
    /**
     * @param gateRadius Largest expected distance between a star and its prediction, in pixels.
     * @param maxStarRadius Largest expected star radius, in pixels. Pads the window past the gate.
     * @param outlierRadius Largest distance of a match's offset from the median offset, in pixels.
     * Negative or NaN: Go finds nothing.
     */
    WindowedCenterOfGravity(int gateRadius, int maxStarRadius, decimal outlierRadius)
        : gateRadius(gateRadius), maxStarRadius(maxStarRadius), outlierRadius(outlierRadius) { };

    /// Largest gateRadius and maxStarRadius, in pixels. Outside 0..kMaxRadius, Go finds nothing.
    static const int kMaxRadius = 1024;
    /// Largest image width and height, in pixels. Larger: Go finds nothing.
    /// Up to this size a float holds pixel centers exactly, to 1/128 pixel.
    static const int kMaxImageSide = 1 << 16;

    /**
     * @param predictions Predicted positions in pixels, as Star::position. Off-image ones are skipped.
     * @param[out] matchedPredictions Prediction index of each returned star. May be null.
     * @return The stars found, ordered by prediction index.
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
