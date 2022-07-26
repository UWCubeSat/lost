#ifndef CENTROID_H
#define CENTROID_H

#include <iostream>
#include <vector>

#include "star-utils.hpp"

namespace lost {

/**
 * @brief
 * @details
 */
class CentroidAlgorithm {
public:
    /**
     * @brief
     * @param image
     * @param imageWidth
     * @param imageHeight
     * @return
     */
    virtual Stars Go(unsigned char *image, int imageWidth, int imageHeight) const = 0;

    /// @brief
    virtual ~CentroidAlgorithm() { };
};

/**
 * @brief
 * @details
 */
class DummyCentroidAlgorithm: public CentroidAlgorithm {
public:
    /**
     * @brief
     * @param numStars
     */
    explicit DummyCentroidAlgorithm(int numStars) : numStars(numStars) { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
private:
    int numStars;
};

/**
 * @brief
 * @details
 */
class CenterOfGravityAlgorithm : public CentroidAlgorithm {
public:
    /// @brief
    CenterOfGravityAlgorithm() { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
};

/**
 * @brief
 * @details
 */
class IterativeWeightedCenterOfGravityAlgorithm : public CentroidAlgorithm {
    public:
        /// @brief
        IterativeWeightedCenterOfGravityAlgorithm() { };
        Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;
};

}

#endif
