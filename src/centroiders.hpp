#ifndef CENTROID_H
#define CENTROID_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_multimap>

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
 * @param Thresholder a class with a constructor that takes a pointer to an image (bytes), int imageWidth, and int imageHeight, then has a method .PixelThreshold(int,int) which returns the threshold for the given x/y coordinate (a pixel must be strictly brighter than the threshold to be considered).
 * @param ConnectedComponentCentroider a class with a zero argument constructor, then two methods: update(int,int,unsigned char), which updates it with information about a new pixel, and merge(const ConnectedComponentThresholder &other), which merges the other into the current, and Vec2 finalize(void), which returns the centroiding result.
 */
template <typename Thresholder, typename ConnectedComponentCentroider>
class ThresholdedConnectedComponentsCentroidAlgorithm : public CentroidAlgorithm {
public:

    virtual Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override {
        Thresholder thresholder = MakeThresholder(image, imageWidth, imageHeight);

        // Weird merging situations should be very rare, so we don't really need a union-find here.

        // owns all the connected components, for garbage collection. We use unique_ptrs so that pointers don't break when the vector is resized.
        std::vector<std::unique_ptr<ConnectedComponentCentroider>> connectedComponentCentroiders;
        std::unordered_map<int, ConnectedComponentCentroider *> idToCC; // this could/should be a vector?
        std::unordered_multimap<ConnectedComponentCentroider *, int> ccToIds;
        
        // At the beginning of iteration for each row, stores the connected component id for each pixel of the preceding row (or -1 if it was below threshold), 
        std::vector<int> rowConnectedComponents(imageWidth, -1);
        for (int y = 0; y < imageHeight; y++) {
            for (int x = 0; x < imageWidth; x++) {
                if (image[y*imageHeight + x] <= thresholder.PixelThreshold(x, y)) {
                    rowConnectedComponents[x] = -1;
                    continue;
                }

                // TODO handle subtracting the threshold. I'm of the opinion we should just pass in the adjusted pixel brightness to the cc centroider and it shouldn't use the image at all.

                int leftCCId = x == 0 ? -1 : rowConnectedComponents[x-1];
                int upCCId = rowConnectedComponents[x]; // already initialized to -1 at start
                int minId = std::min(leftCCId, upCCId);
                int maxId = std::max(leftCCId, upCCId);

                if (maxId < 0) { // neither neighbor is in a CC, create a new one
                    ConnectedComponentCentroider *newCC = new ConnectedComponentCentroider();
                    int newCCId = connectedComponentCentroiders.size();
                    connectedComponentCentroiders.emplace_back(newCC);
                    newCC->update(x, y, image[y*imageHeight + x]);

                    idToCC.emplace(newCCId, newCC);
                    ccToIds.emplace(newCC, newCCId);

                    rowConnectedComponents[x] = newCCId;

                } else if (minId < 0) { // only one neighbor belongs to a CC, assign pixel to that one.
                    int ccId = maxId;
                    rowConnectedComponents[x] = ccId;
                    idToCC[ccId]->update(x, y, image[y*imageHeight + x]);

                } else { // both neighbors in a CC, merge into the one with larger ID
                    ConnectedComponentCentroider *minCC = idToCC[minId];
                    ConnectedComponentCentroider *maxCC = idToCC[maxId];

                    maxCC->merge(minCC);

                    // loop through the ids in the minCC, assign each to the max
                    // using equal_range
                    for (auto it = ccToIds.equal_range(minCC); it.first != it.second; it.first++) {
                        idToCC[it.first->second] = maxCC;
                        ccToIds.emplace(idToCC[maxId], it.first->second);
                    }

                    // we're leaving behind the minCC's old IDs in ccsToIds. But that's okay, because nothing should point to it anymore.

                    rowConnectedComponents[x] = maxId;
                }
            }
        }

        Stars stars;
        // Check which connected components were never merged, and finalize them
        for (auto &cc : connectedComponentCentroiders) {
            // check that it wasn't merged by finding an arbitrary id assigned to the cc, and checking whether that points back to the cc. If the cc was merged, none of its ids point back to it.
            auto idIt = ccToIds.find(cc.get());
            assert(idIt != ccToIds.end());
            if (idToCC[idIt->second] == cc.get()) {
                stars.push_back(cc->finalize());
            }
        }

        return stars;
    }

private:

    /// Make a thresholder. Virtual so that a subclass could override the arguments passed to the thresholder without needing to resort to adding template parameters to the thresholder.
    virtual Thresholder MakeThresholder(unsigned char *image, int imageWidth, int imageHeight) const {
        return Thresholder(image, imageWidth, imageHeight);
    }

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
    explicit LeastSquaresGaussianFit1D(int nb, bool dyn) : nb(nb), dynamic(dyn), np(nb*2+1) { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;

private:
    const int nb;
    bool dynamic;
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
    explicit LeastSquaresGaussianFit2D(int nb, bool dyn) : nb(nb), dynamic(dyn), np(nb*2+1) { };
    Stars Go(unsigned char *image, int imageWidth, int imageHeight) const override;

private:
    const int nb;
    bool dynamic;
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
void SubtractNoise(unsigned char *image, int imageWidth, int imageHeight, float noise);

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
