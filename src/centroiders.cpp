#include "centroiders.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <deque>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
// #include <eigen3/unsupported/Eigen/NumericalDiff>
#include <iostream>
#include <set>  // TODO: remove later, this is just lazy to implement hash for unordered_set
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace lost {

/// A simple thresholder that accepts pixels a certain number of standard deviations above the mean brightness.
class StddevGlobalThresholder {
public:
    StddevGlobalThresholder(unsigned char *image, int imageWidth, int imageHeight, float stddevs) {
        unsigned long totalMag = 0;
        float std = 0;
        long totalPixels = imageHeight * imageWidth;
        for (long i = 0; i < totalPixels; i++) {
            totalMag += image[i];
        }
        float mean = totalMag / totalPixels;

        for (long i = 0; i < totalPixels; i++) {
            std += std::pow(image[i] - mean, 2);
        }
        std = std::sqrt(std / totalPixels);
        threshold = std::min(255, mean + std*stddevs);
    }

    unsigned char PixelThreshold(int x, int y) const {
        return threshold;
    }

private:
    unsigned char threshold;
};

/// A version of StddevGlobalThresholder, but with the threshold fixed at 5 stddevs above the mean (as suggested by a Liebe paper)
class Stddev5GlobalThresholder : public StddevGlobalThresholder {
public:
    Stddev5GlobalThresholder(unsigned char *image, int imageWidth, int imageHeight)
        : StddevGlobalThresholder(image, imageWidth, imageHeight, 5.0) { }
};

/// Each pixel's threshold is a certain number of stddevs above the mean, all calculated according to a local square window around the pixel (truncated near image edges).
class SlidingWindowLocalThresholder {
public:
    SlidingWindowLocalThresholder(unsigned char *image, int imageWidth, int imageHeight, int windowSize, float stddevs)
        : imageWidth(imageWidth), windowSize(windowSize), nwSums(imageWidth*imageHeight) {
        
        assert(windowSize%2 == 1 && windowSize > 1); // only odd window sizes, so we can put an even amount on each side
        int windowRadius = windowSize/2; // rounds down, so the number of pixels to each side not counting the center

        // this sorta uses a lot of memory (4MiB for a 1MiB image). We could just compute for smaller squares at a time and save some memory, but oh well.
        std::vector<unsigned long> nwSums(imageWidth*imageHeight);
        long long sumOfSquares = 0;
        for (int y = 0; y < imageHeight; y++) {
            for (int x = 0; x < imageWidth; x++) {
                long nSum = y == 0 ? 0 : nwSums[(y-1)*imageWidth + x];
                long wSum = x == 0 ? 0 : nwSums[y*imageWidth + (x-1)];
                long nwSum = (x == 0 || y == 0) ? 0 : nwSums[(y-1)*imageWidth + (x-1)];
                nwSums[y*imageWidth + x] = nSum + wSum - nwSum;
                sumOfSquares += image[y*imageWidth + x] * image[y*imageWidth + x];
            }
        }

        // just compute one global stddev
        float mean = (float)nwSums[imageWidth*imageHeight - 1] / (imageWidth*imageHeight);
        float stddev = std::sqrt((float)sumOfSquares / (imageWidth*imageHeight) - mean*mean);

        for (int y = 0; y < imageHeight; y++) {
            for (int x = 0; x < imageWidth; x++) {
                int xMin = std::max(0, x-windowSize);
                int xMax = std::min(imageWidth-1, x+windowSize);
                int yMin = std::max(0, y-windowSize);
                int yMax = std::min(imageHeight-1, y+windowSize);

                thresholds[y*imageWidth + x] = (
                    nwSums[yMax*imageWidth + xMax]
                    - (xMin == 0 ? 0 : nwSums[yMax*imageWidth + xMin-1])
                    - (yMin == 0 ? 0 : nwSums[(yMin-1)*imageWidth + xMax])
                    // if either yMin or xMin is zero, then there was no overlap to be double counted.
                    + (xMin == 0 || yMin == 0 ? 0 : nwSums[(yMin-1)*imageWidth + xMin-1])
                    ) / ((xMax-xMin+1)*(yMax-yMin+1))
                    // and add the user-specified number of stddevs
                    + stddevs*stddev;
            }
        }
    }

    unsigned char PixelThreshold(int x, int y) const {
        return thresholds[y*imageWidth + x];
    }

private:
    int imageWidth;
    std::vector<unsigned char> thresholds;
};

class CenterOfGravityCCCentroider {
public:
    void update(int x, int y, unsigned char intensity) {
        float newMass = mass + intensity;
        center = center*(mass/newMass) + Vec2(x,y)*(intensity/newMass);
        mass = newMass;
    }

    void merge(const CenterOfGravityCCCentroider &other) {
        update(other.center.x(), other.center.y(), other.mass);
    }
private:
    Vec2 center;
    float mass;
};

/**
Prediction of 1D Gaussian model for 1D LSGF method
Note: can be used as f(xi, beta) or f(yi, beta)
*/
static float FitModel(float x, float a, float xb, float sigma);

/// Prediction of 2D Gaussian model for 2D LSGF method
static float FitModel2D(float x, float y, float a, float xb, float yb, float sigmaX, float sigmaY);

static std::vector<unsigned char> CopyImage(const unsigned char* const image, int w, int h);

class Point {
public:
    int x;
    int y;

    bool operator<(const Point &other) const {
        return x<other.x || y<other.y;
    }
};

/**
 * @brief Get a candidate point for each possible star cluster by running a floodfill through image
 *
 * Eliminates possibility of detecting multiple centroids for a single star
 */
std::vector<std::vector<int>> FloodfillPreproc(unsigned char *image, int imageWidth, int imageHeight);

/// A simple, but well tested thresholding algorithm that works well with star images
int BasicThreshold(unsigned char *image, int imageWidth, int imageHeight, float* const noise);

void SubtractNoise(unsigned char *image, int imageWidth, int imageHeight, float noise){
    // float sum = 0;
    // for (long i = 0; i < imageHeight * imageWidth; i++) {
    //     int x = i % imageWidth;
    //     int y = i / imageWidth;
    //     sum += Get(x, y, image, imageWidth);
    // }
    // float noise = sum / (imageWidth * imageHeight);
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        image[i] = std::max(0, int(image[i] - noise));
    }
}

static std::vector<unsigned char> CopyImage(const unsigned char* const image, int w, int h){
    std::vector<unsigned char> res;
    res.insert(res.begin(), &image[0], &image[w*h]);
    return res;
}

struct FloodParams{
    int xMin;
    int xMax;
    int yMin;
    int yMax;
};

std::vector<std::vector<int>> FloodfillPreproc(unsigned char *img, int imageWidth, int imageHeight) {
    std::vector<std::vector<int>> res;
    std::set<Point> checkedPoints;

    std::vector<unsigned char> image = CopyImage(img, imageWidth, imageHeight);

    float noise;
    int cutoff = BasicThreshold(image.data(), imageWidth, imageHeight, &noise);
    cutoff -= noise;
    SubtractNoise(image.data(), imageWidth, imageHeight, noise);

    std::cout << "Cutoff: " << cutoff << std::endl;
    if(cutoff == 0){
        // TODO: scuffed, floodfill will never stop if cutoff=0, so stop it here
        std::cerr << "No stars detected in image, killing process" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (long i = 0; i < imageHeight * imageWidth; i++) {
        int x = i % imageWidth;
        int y = i / imageWidth;
        Point pCurr{x, y};
        if (image[i] >= cutoff && checkedPoints.count(pCurr) == 0) {
            // Floodfill from this point (x, y)
            // Obtain coordinates (x0, y0) of pixel with largest magnitude in the floodfill
            std::vector<Point> pts;
            std::deque<Point> queue;

            FloodParams fp{x, x, y, y};

            int maxMag = -1;
            int x0 = x;
            int y0 = y;

            int floodSize = 0;

            queue.push_back(pCurr);
            while (queue.size() != 0) {
                Point p = queue[0];
                queue.pop_front();

                if (p.x < 0 || p.x >= imageWidth || p.y < 0 || p.y >= imageHeight) continue;
                if (checkedPoints.count(p) != 0) continue;

                int mag = Get(p.x, p.y, image.data(), imageWidth);
                if (mag < cutoff) continue;

                floodSize++;

                fp.xMin = std::min(fp.xMin, p.x);
                fp.xMax = std::max(fp.xMax, p.x);
                fp.yMin = std::min(fp.yMin, p.y);
                fp.yMax = std::max(fp.yMax, p.y);

                // Add this point to pts and checkedPoints
                // We can add to checkedPoints since cutoff is global - ensure no 2 fills collide
                pts.push_back(p);
                checkedPoints.insert(p);

                // Update max pixel value in fill
                if (mag > maxMag) {
                    maxMag = mag;
                    x0 = p.x;
                    y0 = p.y;
                }

                // Add all 8 adjacent points to the queue
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        queue.push_back({p.x + dx, p.y + dy});
                    }
                }
            }

            // Cool, our flood is done
            res.push_back({x0, y0, floodSize, fp.xMax - fp.xMin, fp.yMax - fp.yMin});
        }
    }

    // std::cout << "Floodfill preprocessing: found " << res.size() << " total points" << std::endl;
    return res;
}

// TODO: documentation!
struct CentroidParams {
    float yCoordMagSum;
    float xCoordMagSum;
    long magSum;
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int cutoff;
    bool isValid;
    std::unordered_set<int> checkedIndices;
};

int Get(int x, int y, const unsigned char *image, int w) {
    int ind = y * w + x;
    return image[ind];
}

int XMarginal(int x, int y0, int nb, const unsigned char *image, int w) {
    int sum = 0;
    for (int j = -nb; j <= nb; j++) {
        sum += Get(x, y0 + j, image, w);
    }
    return sum;
}

int YMarginal(int x0, int y, int nb, const unsigned char *image, int w) {
    int sum = 0;
    for (int i = -nb; i <= nb; i++) {
        sum += Get(x0 + i, y, image, w);
    }
    return sum;
}

static float FitModel(float x, float a, float xb, float sigma) {
    return a * exp(-1 * (x - xb) * (x - xb) / (2 * sigma * sigma));
}

static float FitModel2D(float x, float y, float a, float xb, float yb, float sigmaX, float sigmaY) {
    float term1 = exp(-1 * (x - xb) * (x - xb) / (2 * sigmaX * sigmaX));
    float term2 = exp(-1 * (y - yb) * (y - yb) / (2 * sigmaY * sigmaY));
    return a * term1 * term2;
}

// DO NO DELETE
// void InitialGuess(int x0, int y0, const int nb, const unsigned char *image, int w, float *a,
//                   float *xb, float *yb, float *sigma) {
//     // a is set to max intensity value in the window
//     // (xb, yb) = coordinates of pixel with max intensity value
//     int max = -1;
//     for (int i = -nb; i <= nb; i++) {
//         for (int j = -nb; j <= nb; j++) {
//             int pixelValue = Get(x0 + i, y0 + j, image, w);
//             if (pixelValue > max) {
//                 *xb = x0 + i;
//                 *yb = y0 + j;
//                 max = pixelValue;
//             }
//         }
//     }
//     *a = max;
//     // Sigma is a little complicated
//     // sigma = f_whm / (2 * \sqrt{2log(2)})
//     // f_whm \approx sqrt of number of pixels with intensity > 0.5 * a
//     int halfCount = 0;
//     for (int i = -nb; i <= nb; i++) {
//         for (int j = -nb; j <= nb; j++) {
//             if (Get(x0 + i, y0 + j, image, w) > *a / 2) {
//                 halfCount++;
//             }
//         }
//     }
//     float fwhm = std::sqrt(halfCount);
//     *sigma = fwhm / (2 * std::sqrt(2 * std::log(2)));
// }

float FitInitialGuessSigma(int x0, int y0, int maxMag, const int nb, const unsigned char *image,
                           int w) {
    float halfCount = 0;
    for (int i = -nb; i <= nb; i++) {
        for (int j = -nb; j < +nb; j++) {
            if (Get(x0 + i, y0 + j, image, w) > maxMag / 2) {
                halfCount++;
            }
        }
    }
    float fwhm = std::sqrt(halfCount);
    return fwhm / (2 * std::sqrt(2 * std::log(2)));
}

// Generic functor
template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor {
    typedef _Scalar Scalar;
    enum { InputsAtCompileTime = NX, ValuesAtCompileTime = NY };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs;  // Number of parameters in your model
    int m_values;  // Number of data points

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

/// Functor for 1D Least Squares Gaussian Fit
struct LSGF1DFunctor : Functor<float> {
    LSGF1DFunctor(const int nb, const int marg, Eigen::VectorXf X, const unsigned char *image,
                  const int w, const int x0, const int y0)
        : Functor<float>(3, 2 * nb + 1),
          nb(nb),
          marg(marg),
          X(X),
          image(image),
          w(w),
          x0(x0),
          y0(y0) {}

    /*
    Calculate residuals (error = prediction - actual)
    x = parameters (a, xb, sigma)
    */
    int operator()(const Eigen::VectorXf &x, Eigen::VectorXf &fvec) const {
        float a = x(0);
        float xb = x(1);
        float sigma = x(2);
        for (int i = 0; i < X.size(); i++) {
            int marginal;
            if (marg == 0)
                marginal = XMarginal(X(i), y0, nb, image, w);
            else
                marginal = YMarginal(x0, X(i), nb, image, w);
            fvec(i) = FitModel(X(i), a, xb, sigma) - marginal;
        }

        return 0;
    }

    // Calculate Jacobian
    int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const {
        float a = x(0);
        float kb = x(1);
        float sK = x(2);

        int ind = 0;

        for (int i = -nb; i <= nb; i++) {
            int k = (marg == 0) ? (x0 + i) : (y0 + i);

            float expK = exp(-1 * (k - kb) * (k - kb) / (2 * sK * sK));

            fjac(ind, 0) = expK;
            fjac(ind, 1) = a * expK * (k - kb) / (sK * sK);
            fjac(ind, 2) = a * expK * (k - kb) * (k - kb) / std::pow(sK, 3);

            ind++;
        }

        return 0;
    }

    // Window size is (2*nb + 1) x (2*nb + 1)
    const int nb;
    // Flag for which marginal to use (0 = X, 1 = Y)
    const int marg;
    // Data points (set of x or y coordinates)
    Eigen::VectorXf X;
    const unsigned char *image;
    const int w;  // image width in pixels
    // Center coordinates (x0, y0) of this window
    const int x0, y0;
};

/// Functor for 2D Least Squares Gaussian Fit
struct LSGF2DFunctor : Functor<float> {
    // We now have 5 params, beta = (a, xb, yb, sigmaX, sigmaY)
    // Let entire window be (np x np) pixels
    // We have np^2 data points
    LSGF2DFunctor(const int nb, const unsigned char *image, const int w, const int x0, const int y0)
        : Functor<float>(5, (2 * nb + 1) * (2 * nb + 1)),
          nb(nb),
          image(image),
          w(w),
          x0(x0),
          y0(y0) {}

    /*
    x = parameters (a, xb, yb, sigmaX, sigmaY)
    */
    int operator()(const Eigen::VectorXf &x, Eigen::VectorXf &fvec) const {
        float a = x(0);
        float xb = x(1);
        float yb = x(2);
        float sigmaX = x(3);
        float sigmaY = x(4);

        int ind = 0;
        for (int i = -nb; i <= nb; i++) {
            for (int j = -nb; j <= nb; j++) {
                int xi = x0 + i;
                int yi = y0 + j;
                float yPred = FitModel2D(xi, yi, a, xb, yb, sigmaX, sigmaY);
                float yActual = Get(xi, yi, image, w);
                fvec(ind) = yPred - yActual;

                ind++;
            }
        }

        return 0;
    }

    int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const {
        float a = x(0);
        float xb = x(1);
        float yb = x(2);
        float sigmaX = x(3);
        float sigmaY = x(4);

        int ind = 0;

        for (int i = -nb; i <= nb; i++) {
            for (int j = -nb; j <= nb; j++) {
                int xi = x0 + i;
                int yj = y0 + j;

                float expX = exp(-1 * (xi - xb) * (xi - xb) / (2 * sigmaX * sigmaX));
                float expY = exp(-1 * (yj - yb) * (yj - yb) / (2 * sigmaY * sigmaY));

                fjac(ind, 0) = expX * expY;
                fjac(ind, 1) = a * expX * expY * (xi - xb) / (sigmaX * sigmaX);
                fjac(ind, 2) = a * expX * expY * (yj - yb) / (sigmaY * sigmaY);
                fjac(ind, 3) = a * expX * expY * std::pow(xi - xb, 2) / (std::pow(sigmaX, 3));
                fjac(ind, 4) = a * expX * expY * std::pow(yj - yb, 2) / (std::pow(sigmaY, 3));

                ind++;
            }
        }

        return 0;
    }

    // Window size is (2*nb + 1) x (2*nb + 1)
    const int nb;

    const unsigned char *image;
    const int w;  // image width in pixels
    // Center coordinates (x0, y0) of this window
    const int x0, y0;
};

std::vector<Star> LeastSquaresGaussianFit1D::Go(unsigned char *image, int imageWidth,
                                                int imageHeight) const {
    std::vector<Star> result;

    std::cout << "1D Gaussian Fit (" << np << "x" << np << ")" << std::endl;

    std::vector<std::vector<int>> candidatePts = FloodfillPreproc(image, imageWidth, imageHeight);

    for (const std::vector<int> &pt : candidatePts) {
        int x = pt[0];
        int y = pt[1];
        // TODO: just taking largest diameter right now
        int dynamicWinSize = std::max(pt[3], pt[4]);

        int winRadius = (dynamic) ? dynamicWinSize : nb; // /2?
        // int winRadius = nb;

        if (x - winRadius < 0 || x + winRadius >= imageWidth || y - winRadius < 0 ||
            y + winRadius >= imageHeight)
            continue;

        Eigen::VectorXf X(2 * winRadius + 1);
        Eigen::VectorXf Y(2 * winRadius + 1);
        int vInd = 0;
        for (int i = -winRadius; i <= winRadius; i++) {
            X(vInd) = x + i;
            Y(vInd) = y + i;
            vInd++;
        }

        float a = Get(x, y, image, imageWidth);
        float sigma = FitInitialGuessSigma(x, y, a, winRadius, image, imageWidth);

        Eigen::VectorXf betaX(3);
        betaX << a, x, sigma;

        Eigen::VectorXf betaY(3);
        betaY << a, y, sigma;

        LSGF1DFunctor functorX(winRadius, 0, X, image, imageWidth, x, y);
        Eigen::LevenbergMarquardt<LSGF1DFunctor, float> lmX(functorX);

        lmX.parameters.maxfev = 2000;
        lmX.parameters.xtol = 1.0e-10;

        LSGF1DFunctor functorY(winRadius, 1, Y, image, imageWidth, x, y);
        Eigen::LevenbergMarquardt<LSGF1DFunctor, float> lmY(functorY);

        lmY.parameters.maxfev = 2000;
        lmY.parameters.xtol = 1.0e-10;

        lmX.minimize(betaX);
        lmY.minimize(betaY);

        a = betaX(0);

        float xb = betaX(1);
        float yb = betaY(1);

        sigma = betaX(2);

        result.push_back(Star(xb + 0.5, yb + 0.5, sigma, sigma, pt[2]));
    }

    return result;
}

std::vector<Star> LeastSquaresGaussianFit2D::Go(unsigned char *image, int imageWidth,
                                                int imageHeight) const {
    std::vector<Star> result;

    std::cout << "2D Gaussian Fit (" << np << "x" << np << ")" << std::endl;

    std::vector<std::vector<int>> candidatePts = FloodfillPreproc(image, imageWidth, imageHeight);

    for (const std::vector<int> &pt : candidatePts) {
        int x = pt[0];
        int y = pt[1];

        // TODO: just taking largest diameter right now
        int dynamicWinSize = std::max(pt[3], pt[4]);

        int winRadius = (dynamic) ? dynamicWinSize : nb; // /2?
        // int winRadius = nb;

        if (x - winRadius < 0 || x + winRadius >= imageWidth || y - winRadius < 0 ||
            y + winRadius >= imageHeight)
            continue;

        float a = Get(x, y, image, imageWidth);
        float sigma = FitInitialGuessSigma(x, y, a, winRadius, image, imageWidth);

        Eigen::VectorXf beta(5);
        beta << a, x, y, sigma, sigma;

        LSGF2DFunctor functor(winRadius, image, imageWidth, x, y);
        Eigen::LevenbergMarquardt<LSGF2DFunctor, float> lm(functor);

        lm.parameters.maxfev = 2000;
        lm.parameters.xtol = 1.0e-10;

        lm.minimize(beta);

        float aRes = beta(0);
        float xRes = beta(1);
        float yRes = beta(2);
        float sigmaX = beta(3);
        float sigmaY = beta(4);

        result.push_back(Star(xRes + 0.5, yRes + 0.5, sigmaX, sigmaY, pt[2]));
    }

    return result;
}

static void GetGGridCoeffs(const std::vector<int> &w, std::vector<int> *const a,
                           std::vector<int> *const b) {
    (*a)[0] = 2 * w[2] * w[1] * w[0] + 6 * w[2] * w[0] * w[3] +
              12 * (w[3] * w[4] * w[0] + w[3] * w[1] * w[0]) + 32 * w[2] * w[4] * w[0] +
              36 * w[4] * w[1] * w[0];
    (*a)[1] = 2 * w[2] * w[3] * w[1] + 6 * w[3] * w[4] * w[1] - 4 * w[2] * w[1] * w[0] +
              12 * w[2] * w[4] * w[1] - 18 * w[3] * w[1] * w[0] - 48 * w[4] * w[1] * w[0];
    (*a)[2] = 2 * (w[2] * w[3] * w[4] - w[2] * w[1] * w[0]) - 4 * w[1] * w[2] * w[3] -
              18 * (w[0] * w[2] * w[3] + w[1] * w[2] * w[4]) - 64 * w[0] * w[2] * w[4];
    (*a)[3] = 2 * w[2] * w[1] * w[3] + 6 * w[3] * w[1] * w[0] - 4 * w[2] * w[3] * w[4] +
              12 * w[2] * w[3] * w[0] - 18 * w[3] * w[4] * w[1] - 48 * w[3] * w[4] * w[0];
    (*a)[4] = 2 * w[2] * w[3] * w[4] + 6 * w[2] * w[4] * w[1] +
              12 * (w[3] * w[4] * w[1] + w[4] * w[1] * w[0]) + 32 * w[2] * w[4] * w[0] +
              36 * w[3] * w[4] * w[0];

    (*b)[0] = -w[2] * w[1] * w[0] + 3 * w[2] * w[3] * w[0] +
              18 * (w[3] * w[4] * w[0] + w[4] * w[1] * w[0]) + 32 * w[2] * w[4] * w[0];
    (*b)[1] = w[2] * w[3] * w[1] + 4 * w[2] * w[1] * w[0] +
              9 * (w[3] * w[4] * w[1] + w[3] * w[1] * w[0]) + 12 * w[2] * w[4] * w[1];
    (*b)[2] = 3 * (w[2] * w[3] * w[4] - w[2] * w[1] * w[0]) +
              9 * (w[2] * w[3] * w[0] - w[2] * w[1] * w[4]);
    (*b)[3] = -w[2] * w[3] * w[1] - 4 * w[2] * w[3] * w[4] -
              9 * (w[3] * w[4] * w[1] + w[3] * w[0] * w[1]) - 12 * w[2] * w[3] * w[0];
    (*b)[4] = w[2] * w[3] * w[4] - 3 * w[2] * w[1] * w[4] -
              18 * (w[3] * w[0] * w[4] + w[1] * w[0] * w[4]) - 32 * w[2] * w[4] * w[0];
}

Stars GaussianGrid::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    std::vector<Star> result;

    std::cout << "Gaussian Grid (" << np << "x" << np << ")" << std::endl;

    std::vector<std::vector<int>> candidatePts = FloodfillPreproc(image, imageWidth, imageHeight);

    for (const std::vector<int> &pt : candidatePts) {
        int x = pt[0];
        int y = pt[1];
        if (x - nb < 0 || x + nb >= imageWidth || y - nb < 0 || y + nb >= imageHeight) continue;

        float nom = 0;
        float denom = 0;

        for (int j = -nb; j <= nb; j++) {
            std::vector<int> w;
            for (int k = -nb; k <= nb; k++) {
                w.push_back(std::pow(Get(x + k, y + j, image, imageWidth), 1));
            }
            std::vector<int> a(np);
            std::vector<int> b(np);
            GetGGridCoeffs(w, &a, &b);
            for (int i = -nb; i <= nb; i++) {
                float v = Get(x + i, y + j, image, imageWidth);
                nom += b[i + nb] * ((v == 0) ? -1e6 : std::log(v));
                denom += a[i + nb] * ((v == 0) ? -1e6 : std::log(v));
            }
        }

        float xRes = x + nom / denom;

        nom = 0;
        denom = 0;
        for (int i = -nb; i <= nb; i++) {
            std::vector<int> w;
            for (int k = -nb; k <= nb; k++) {
                w.push_back(Get(x + i, y + k, image, imageWidth));
            }
            std::vector<int> a(np);
            std::vector<int> b(np);
            GetGGridCoeffs(w, &a, &b);
            for (int j = -nb; j <= nb; j++) {
                float v = Get(x + i, y + j, image, imageWidth);
                nom += b[j + nb] * ((v == 0) ? -1e6 : std::log(v));
                denom += a[j + nb] * ((v == 0) ? -1e6 : std::log(v));
            }
        }

        float yRes = y + nom / denom;

        result.push_back(Star(xRes + 0.5, yRes + 0.5, nb, nb, pt[2]));
    }

    return result;
}

// DUMMY
std::vector<Star> DummyCentroidAlgorithm::Go(unsigned char *, int imageWidth,
                                             int imageHeight) const {
    std::vector<Star> result;

    unsigned int randomSeed = 123456;
    for (int i = 0; i < numStars; i++) {
        result.push_back(
            Star(rand_r(&randomSeed) % imageWidth, rand_r(&randomSeed) % imageHeight, 10.0));
    }

    return result;
}

// a poorly designed thresholding algorithm
int BadThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    // loop through entire array, find sum of magnitudes
    long totalMag = 0;
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    return (((totalMag / (imageHeight * imageWidth)) + 1) * 15) / 10;
}

// a more sophisticated thresholding algorithm, not tailored to star images
int OtsusThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    // code here, duh
    long total = imageWidth * imageHeight;
    // float top = 255;
    float sumB = 0;
    float sum1 = 0;
    float wB = 0;
    float maximum = 0;
    int level = 0;
    // make the histogram (array length 256)
    int histogram[256];

    memset(histogram, 0, sizeof(int) * 256);

    for (long i = 0; i < total; i++) {
        histogram[image[i]]++;
    }
    for (int i = 0; i < 256; i++) {
        sum1 += i * histogram[i];
    }
    for (int i = 0; i < 256; i++) {
        float wF = total - wB;
        // std::cout << "wF\n" << wB << "\n";
        // std::cout << "wB\n" << wF << "\n";
        if (wB > 0 && wF > 0) {
            float mF = (sum1 - sumB) / wF;
            float val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
            // std::cout << val << "\n";
            if (val >= maximum) {
                level = i;
                maximum = val;
            }
        }
        wB = wB + histogram[i];
        sumB = sumB + i * histogram[i];
    }
    return level;
}

int BasicThreshold(unsigned char *image, int imageWidth, int imageHeight, float* const noise) {
    unsigned long totalMag = 0;
    float std = 0;
    long totalPixels = imageHeight * imageWidth;
    for (long i = 0; i < totalPixels; i++) {
        totalMag += image[i];
    }
    float mean = totalMag / totalPixels;
    *noise = mean;

    for (long i = 0; i < totalPixels; i++) {
        std += std::pow(image[i] - mean, 2);
    }
    std = std::sqrt(std / totalPixels);
    return mean + (std * 5);
}

// basic thresholding, but do it faster (trade off of some accuracy?)
int BasicThresholdOnePass(unsigned char *image, int imageWidth, int imageHeight) {
    unsigned long totalMag = 0;
    float std = 0;
    float sq_totalMag = 0;
    long totalPixels = imageHeight * imageWidth;
    for (long i = 0; i < totalPixels; i++) {
        totalMag += image[i];
        sq_totalMag += image[i] * image[i];
    }
    float mean = totalMag / totalPixels;
    float variance = (sq_totalMag / totalPixels) - (mean * mean);
    std = std::sqrt(variance);
    return mean + (std * 5);
}

// recursive helper here
void CogHelper(CentroidParams *p, long i, unsigned char *image, int imageWidth, int imageHeight) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p->cutoff &&
        p->checkedIndices.count(i) == 0) {
        // check if pixel is on the edge of the image, if it is, we dont want to centroid this star
        if (i % imageWidth == 0 || i % imageWidth == imageWidth - 1 || i / imageWidth == 0 ||
            i / imageWidth == imageHeight - 1) {
            p->isValid = false;
        }
        p->checkedIndices.insert(i);
        if (i % imageWidth > p->xMax) {
            p->xMax = i % imageWidth;
        } else if (i % imageWidth < p->xMin) {
            p->xMin = i % imageWidth;
        }
        if (i / imageWidth > p->yMax) {
            p->yMax = i / imageWidth;
        } else if (i / imageWidth < p->yMin) {
            p->yMin = i / imageWidth;
        }
        p->magSum += image[i];
        p->xCoordMagSum += ((i % imageWidth)) * image[i];
        p->yCoordMagSum += ((i / imageWidth)) * image[i];
        if (i % imageWidth != imageWidth - 1) {
            CogHelper(p, i + 1, image, imageWidth, imageHeight);
        }
        if (i % imageWidth != 0) {
            CogHelper(p, i - 1, image, imageWidth, imageHeight);
        }
        CogHelper(p, i + imageWidth, image, imageWidth, imageHeight);
        CogHelper(p, i - imageWidth, image, imageWidth, imageHeight);
    }
}

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *img, int imageWidth,
                                               int imageHeight) const {
    CentroidParams p;

    std::vector<Star> result;
    std::vector<unsigned char> image = CopyImage(img, imageWidth, imageHeight);

    float noise;
    p.cutoff = BasicThreshold(image.data(), imageWidth, imageHeight, &noise);
    p.cutoff -= noise;
    SubtractNoise(image.data(), imageWidth, imageHeight, noise);

    for (long i = 0; i < imageHeight * imageWidth; i++) {
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
            // iterate over pixels that are part of the star
            int xDiameter = 0;  // radius of current star
            int yDiameter = 0;
            p.yCoordMagSum = 0;  // y coordinate of current star
            p.xCoordMagSum = 0;  // x coordinate of current star
            p.magSum = 0;        // sum of magnitudes of current star

            p.xMax = i % imageWidth;
            p.xMin = i % imageWidth;
            p.yMax = i / imageWidth;
            p.yMin = i / imageWidth;
            p.isValid = true;

            int sizeBefore = p.checkedIndices.size();

            CogHelper(&p, i, image.data(), imageWidth, imageHeight);
            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            // use the sums to finish CoG equation and add stars to the result
            float xCoord = (p.xCoordMagSum / (p.magSum * 1.0));
            float yCoord = (p.yCoordMagSum / (p.magSum * 1.0));

            if (p.isValid) {
                result.push_back(Star(xCoord + 0.5f, yCoord + 0.5f, ((float)(xDiameter)) / 2.0f,
                                      ((float)(yDiameter)) / 2.0f,
                                      p.checkedIndices.size() - sizeBefore));
            }
        }
    }
    return result;
}

// Determines how accurate and how much iteration is done by the IWCoG algorithm,
// smaller means more accurate and more iterations.
float iWCoGMinChange = 0.0002;

struct IWCoGParams {
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int cutoff;
    int maxIntensity;
    int guess;
    bool isValid;
    std::unordered_set<int> checkedIndices;
};

void IWCoGHelper(IWCoGParams *p, long i, unsigned char *image, int imageWidth, int imageHeight,
                 std::vector<int> *starIndices) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p->cutoff &&
        p->checkedIndices.count(i) == 0) {
        // check if pixel is on the edge of the image, if it is, we dont want to centroid this star
        if (i % imageWidth == 0 || i % imageWidth == imageWidth - 1 || i / imageWidth == 0 ||
            i / imageWidth == imageHeight - 1) {
            p->isValid = false;
        }
        p->checkedIndices.insert(i);
        starIndices->push_back(i);
        if (image[i] > p->maxIntensity) {
            p->maxIntensity = image[i];
            p->guess = i;
        }
        if (i % imageWidth > p->xMax) {
            p->xMax = i % imageWidth;
        } else if (i % imageWidth < p->xMin) {
            p->xMin = i % imageWidth;
        }
        if (i / imageWidth > p->yMax) {
            p->yMax = i / imageWidth;
        } else if (i / imageWidth < p->yMin) {
            p->yMin = i / imageWidth;
        }
        if (i % imageWidth != imageWidth - 1) {
            IWCoGHelper(p, i + 1, image, imageWidth, imageHeight, starIndices);
        }
        if (i % imageWidth != 0) {
            IWCoGHelper(p, i - 1, image, imageWidth, imageHeight, starIndices);
        }
        IWCoGHelper(p, i + imageWidth, image, imageWidth, imageHeight, starIndices);
        IWCoGHelper(p, i - imageWidth, image, imageWidth, imageHeight, starIndices);
    }
}

Stars IterativeWeightedCenterOfGravityAlgorithm::Go(unsigned char *img, int imageWidth,
                                                    int imageHeight) const {
    IWCoGParams p;
    std::vector<Star> result;

    std::vector<unsigned char> image = CopyImage(img, imageWidth, imageHeight);
    float noise;
    p.cutoff = BasicThreshold(image.data(), imageWidth, imageHeight, &noise);
    p.cutoff -= noise;
    SubtractNoise(image.data(), imageWidth, imageHeight, noise);

    for (long i = 0; i < imageHeight * imageWidth; i++) {
        // check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
            // TODO: store longs --Mark
            std::vector<int> starIndices;  // indices of the current star
            p.maxIntensity = 0;
            int xDiameter = 0;
            int yDiameter = 0;
            float yWeightedCoordMagSum = 0;
            float xWeightedCoordMagSum = 0;
            float weightedMagSum = 0;
            float fwhm;  // fwhm variable
            float standardDeviation;
            float w;  // weight value

            p.xMax = i % imageWidth;
            p.xMin = i % imageWidth;
            p.yMax = i / imageWidth;
            p.yMin = i / imageWidth;
            p.isValid = true;

            IWCoGHelper(&p, i, image.data(), imageWidth, imageHeight, &starIndices);

            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            // calculate fwhm
            float count = 0;
            for (int j = 0; j < (int)starIndices.size(); j++) {
                if (image[starIndices.at(j)] > p.maxIntensity / 2) {
                    count++;
                }
            }
            fwhm = sqrt(count);
            standardDeviation = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
            float modifiedStdDev = 2.0 * pow(standardDeviation, 2);
            // TODO: Why are these floats? --Mark
            float guessXCoord = (float)(p.guess % imageWidth);
            float guessYCoord = (float)(p.guess / imageWidth);
            // how much our new centroid estimate changes w each iteration
            float change = INFINITY;
            int stop = 0;
            // while we see some large enough change in estimated, maybe make it a global variable
            while (change > iWCoGMinChange && stop < 100000) {
                // traverse through star indices, calculate W at each coordinate, add to final
                // coordinate sums
                yWeightedCoordMagSum = 0;
                xWeightedCoordMagSum = 0;
                weightedMagSum = 0;
                stop++;
                for (long j = 0; j < (long)starIndices.size(); j++) {
                    // calculate w
                    float currXCoord = (float)(starIndices.at(j) % imageWidth);
                    float currYCoord = (float)(starIndices.at(j) / imageWidth);
                    w = p.maxIntensity *
                        exp(-1.0 * ((pow(currXCoord - guessXCoord, 2) / modifiedStdDev) +
                                    (pow(currYCoord - guessYCoord, 2) / modifiedStdDev)));

                    xWeightedCoordMagSum += w * currXCoord * ((float)image[starIndices.at(j)]);
                    yWeightedCoordMagSum += w * currYCoord * ((float)image[starIndices.at(j)]);
                    weightedMagSum += w * ((float)image[starIndices.at(j)]);
                }
                float xTemp = xWeightedCoordMagSum / weightedMagSum;
                float yTemp = yWeightedCoordMagSum / weightedMagSum;

                change = abs(guessXCoord - xTemp) + abs(guessYCoord - yTemp);

                guessXCoord = xTemp;
                guessYCoord = yTemp;
            }
            if (p.isValid) {
                result.push_back(Star(guessXCoord + 0.5f, guessYCoord + 0.5f,
                                      ((float)(xDiameter)) / 2.0f, ((float)(yDiameter)) / 2.0f,
                                      starIndices.size()));
            }
        }
    }
    return result;
}

}  // namespace lost
