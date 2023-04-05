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

/**
 * @brief Least Squares Gaussian Fit (1D) centroiding algorithm
 *
 * Detect centroids by fitting a 1D Gaussian to the marginals of each window
 * Slightly less accurate than 2D fit
 * Note: increasing window size may introduce noise, making this method less accurate on small stars
 */
class LeastSquaresGaussianFit1D : public CentroidAlgorithm{
public:
    explicit LeastSquaresGaussianFit1D(int nb) : nb(nb), np(nb*2+1) { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;

private:
    const int nb;
    const int np;
};

/**
 * @brief Least Squares Gaussian Fit (2D) centroiding algorithm
 *
 * Detect centroids by fitting a 2D Gaussian to all pixels in each window
 * This is the most accurate centroiding function to date
 * Also more computationally expensive than other methods, scales exponentially with increasing window size
 */
class LeastSquaresGaussianFit2D : public CentroidAlgorithm {
public:
    explicit LeastSquaresGaussianFit2D(int nb) : nb(nb), np(nb*2+1) { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;

private:
    const int nb;
    const int np;
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

class GaussianGrid : public CentroidAlgorithm {
public:
    GaussianGrid(){ };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;

private:
    // Only support 5x5 Gaussian Grid for now
    const int nb = 2;
    const int np = nb*2 + 1;
};

/// Get value of pixel at (x, y) in image with width=w
int Get(int x, int y, const unsigned char *image, int w);

/// Get value of x marginal - keep x fixed, sum pixel values from [y0-nb, y0+nb]
int XMarginal(int x, int y0, int nb, const unsigned char *image, int w);

/// Get value of y marginal - keep y fixed, sum pixel values from [x0-nb, x0+nb]
int YMarginal(int x0, int y, int nb, const unsigned char *image, int w);

/// Compute and subtract noise = mean of all pixels from image
void SubtractNoise(unsigned char *image, int imageWidth, int imageHeight);

/**
 * @brief Given window centered at (x0, y0) in image with width=w, calculate initial guess for sigma
 * parameter of Gaussian function
 * Only used with Least Squares Gaussian Fit algorithms
 */
float FitInitialGuessSigma(int x0, int y0, int maxMag, const int nb, const unsigned char *image,
                           int w);

// DO NOT DELETE
/*
Get initial guess for window centered at (x0, y0) with given nb
Output:
a = max intensity value
(xb, yb) = coordinates of pixel with max intensity
sigma = standard deviation (sigmaX = sigmaY)
*/
// void InitialGuess(int x0, int y0, const int nb, const unsigned char *image, int w, float *a,
//                   float *xb, float *yb, double *sigma);

} // namespace lost

#endif
