#include "centroiders.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <deque>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
#include <eigen3/unsupported/Eigen/NumericalDiff>
#include <iostream>
#include <set>  // TODO: remove later, this is just lazy to implement hash for unordered_set
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace lost {

typedef std::vector<int> Point;

std::vector<Point> FloodfillPreproc(unsigned char *image, int imageWidth, int imageHeight);

int BasicThreshold(unsigned char *image, int imageWidth, int imageHeight);

std::vector<Point> FloodfillPreproc(unsigned char *image, int imageWidth, int imageHeight) {
    std::vector<Point> res;
    std::set<Point> checkedPoints;

    int cutoff = BasicThreshold(image, imageWidth, imageHeight);
    std::cout << "cutoff: " << cutoff << std::endl;

    for (long i = 0; i < imageHeight * imageWidth; i++) {
        int x = i % imageWidth;
        int y = i / imageWidth;
        Point pCurr{x, y};
        if (image[i] >= cutoff && checkedPoints.count(pCurr) == 0) {
            // Floodfill from this point (x, y)
            // Obtain coordinates (x0, y0) of pixel with largest magnitude in the floodfill
            std::vector<Point> pts;
            std::deque<Point> queue;

            int maxMag = -1;
            int x0 = x;
            int y0 = y;

            queue.push_back(pCurr);
            while (queue.size() != 0) {
                Point p = queue[0];
                queue.pop_front();

                // TODO: make Point a struct probably
                int px = p[0];
                int py = p[1];
                if (px < 0 || px >= imageWidth || py < 0 || py >= imageHeight) continue;
                if (checkedPoints.count(p) != 0) continue;

                int mag = Get(px, py, image, imageWidth);
                if (mag < cutoff) continue;

                // Add this point to pts and checkedPoints
                // We can add to checkedPoints since cutoff is global - ensure no 2 fills collide
                pts.push_back(p);
                checkedPoints.insert(p);

                // Update max pixel value in fill
                if (mag > maxMag) {
                    maxMag = mag;
                    x0 = px;
                    y0 = py;
                }

                // Add all 8 adjacent points to the queue
                for (int dx = -1; dx <= 1; dx++) {
                    for (int dy = -1; dy <= 1; dy++) {
                        queue.push_back({px + dx, py + dy});
                    }
                }
            }

            // Cool, our flood is done
            res.push_back({x0, y0});
        }
    }

    std::cout << "Floodfill preprocessing: found " << res.size() << " total points" << std::endl;
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

/*
Get pixel value at (x, y) where x and y are 0-based
Note: top left is (0, 0)
*/
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

/*
Can be used as f(xi, beta) or f(yi, beta)
*/
float FitModel(float x, float a, float xb, float sigma) {
    return a * exp(-1 * (x - xb) * (x - xb) / (2 * sigma * sigma));
}

///
float FitModel2D(float x, float y, float a, float xb, float yb, float sigmaX, float sigmaY){
    float term1 = exp(-1 * (x - xb) * (x - xb) / (2 * sigmaX * sigmaX));
    float term2 = exp(-1 * (y - yb) * (y - yb) / (2 * sigmaY * sigmaY));
    return a * term1 * term2;
}


void InitialGuess(int x0, int y0, const int nb, const unsigned char *image, int w, float *a,
                  float *xb, float *yb, double *sigma) {
    // a is set to max intensity value in the window
    // (xb, yb) = coordinates of pixel with max intensity value
    int max = -1;
    for (int i = -nb; i <= nb; i++) {
        for (int j = -nb; j <= nb; j++) {
            int pixelValue = Get(x0 + i, y0 + j, image, w);
            if (pixelValue > max) {
                *xb = x0 + i;
                *yb = y0 + j;
                max = pixelValue;
            }
        }
    }
    *a = max;
    // Sigma is a little complicated
    // sigma = f_whm / (2 * \sqrt{2log(2)})
    // f_whm \approx sqrt of number of pixels with intensity > 0.5 * a
    int halfCount = 0;
    for (int i = -nb; i <= nb; i++) {
        for (int j = -nb; j <= nb; j++) {
            if (Get(x0 + i, y0 + j, image, w) > *a / 2) {
                halfCount++;
            }
        }
    }
    double fwhm = std::sqrt(halfCount);
    *sigma = fwhm / (2 * std::sqrt(2 * std::log(2)));
}

float InitialGuess2(int x0, int y0, int maxMag, const int nb, const unsigned char *image, int w){
    float halfCount = 0;
    for(int i = -nb; i <= nb; i++){
        for(int j = -nb; j <+ nb; j++){
            if(Get(x0+i, y0+j, image, w) > maxMag / 2){
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

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct LSGFFunctor : Functor<double> {
    // First param = number of parameters in your model
    // Second param = number of data points you want to test over
    // for us, = 2 * nb + 1
    LSGFFunctor(const int nb, const int marg, Eigen::VectorXd X, const unsigned char *image,
                const int w, const int x0, const int y0)
        : Functor<double>(3, 2 * nb + 1),
          nb(nb),
          marg(marg),
          X(X),
          image(image),
          w(w),
          x0(x0),
          y0(y0) {}

    /*
    x = parameters (a, xb, sigma)
    fvec = residual (one for every data point, of size = 2*nb+1)
    */
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {
        double a = x(0);
        double xb = x(1);
        double sigma = x(2);
        for (int i = 0; i < X.size(); i++) {
            int marginal;
            if (marg == 0)
                marginal = XMarginal(X(i), y0, nb, image, w);
            else
                marginal = YMarginal(x0, X(i), nb, image, w);
            fvec(i) = marginal - FitModel(X(i), a, xb, sigma);
        }

        return 0;
    }

    // Window size is (2*nb + 1) x (2*nb + 1)
    const int nb;
    // Flag for which marginal to use (0 = X, 1 = Y)
    const int marg;
    // Data points (set of x or y coordinates)
    Eigen::VectorXd X;
    const unsigned char *image;
    const int w;  // image width in pixels
    // Center coordinates (x0, y0) of this window
    const int x0, y0;
};

/// Functor for 2D Least-squares Gaussian Fit algo
struct LSGF2DFunctor : Functor<double> {
    // We now have 5 params, beta = (a, xb, yb, sigmaX, sigmaY)
    // Let entire window be (np x np) pixels
    // We have np^2 data points
    LSGF2DFunctor(const int nb, const unsigned char *image,
                const int w, const int x0, const int y0)
        : Functor<double>(5, (2 * nb + 1) * (2 * nb + 1)),
          nb(nb),
          image(image),
          w(w),
          x0(x0),
          y0(y0) {}

    /*
    x = parameters (a, xb, yb, sigmaX, sigmaY)
    */
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {
        double a = x(0);
        double xb = x(1);
        double yb = x(2);
        double sigmaX = x(3);
        double sigmaY = x(4);

        int ind = 0;
        for(int i = -nb; i <= nb; i++){
            for(int j = -nb; j <= nb; j++){
                int xi = x0 + i;
                int yi = y0 + j;
                // TODO: other way around, yPred is our modelPred, yActual is our pixel intensity
                float yPred = FitModel2D(xi, yi, a, xb, yb, sigmaX, sigmaY);
                float yActual = Get(xi, yi, image, w);
                fvec(ind) = yPred - yActual;

                ind++;
            }
        }

        return 0;
    }

    int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const{
        double a = x(0);
        double xb = x(1);
        double yb = x(2);
        double sigmaX = x(3);
        double sigmaY = x(4);

        for (int i = -nb; i <= nb; i++) {
            for (int j = -nb; j <= nb; j++) {
                int xi = x0 + i;
                int yj = y0 + j;

                float yActual = Get(xi, yj, image, w);
                float yPred = FitModel2D(xi, yj, a, xb, yb, sigmaX, sigmaY);
                float ei = yActual - yPred;

                float expX = exp(-1 * (xi - xb) * (xi - xb) / (2 * sigmaX * sigmaX));
                float expY = exp(-1 * (yj - yb) * (yj - yb) / (2 * sigmaY * sigmaY));

                int row = i + nb;
                fjac(row, 0) = -expX * expY;
                fjac(row, 1) = -a * expY * (xi - xb) / (sigmaX * sigmaX) * expX;
                fjac(row, 2) = -a * expX * (yj - yb) / (sigmaY * sigmaY) * expY;
                fjac(row, 3) = -a * expX * expY * std::pow(xi - xb, 2) / (std::pow(sigmaX, 3));
                fjac(row, 4) = -a * expX * expY * std::pow(yj - yb, 2) / (std::pow(sigmaY, 3));
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

    std::cout << "1D Gaussian Fit" << std::endl;

    const int nb = 2;

    // std::set<Point> checkedPoints;

    std::vector<Point> candidatePts = FloodfillPreproc(image, imageWidth, imageHeight);

    for (const Point &pt : candidatePts) {
        int x = pt[0];
        int y = pt[1];
        if (x - nb < 0 || x + nb >= imageWidth || y - nb < 0 || y + nb >= imageHeight) continue;

        Eigen::VectorXd X(2 * nb + 1);
        Eigen::VectorXd Y(2 * nb + 1);
        int vInd = 0;
        for (int i = -nb; i <= nb; i++) {
            X(vInd) = x + i;
            Y(vInd) = y + i;
            vInd++;
        }

        float a;
        float xb, yb;
        double sigma;
        InitialGuess(x, y, nb, image, imageWidth, &a, &xb, &yb, &sigma);

        Eigen::VectorXd betaX(3);
        betaX << a, xb, sigma;

        Eigen::VectorXd betaY(3);
        betaY << a, yb, sigma;

        LSGFFunctor functorX(nb, 0, X, image, imageWidth, x, y);
        Eigen::NumericalDiff<LSGFFunctor> numDiffX(functorX);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LSGFFunctor>, double> lmX(numDiffX);

        lmX.parameters.maxfev = 2000;
        lmX.parameters.xtol = 1.0e-10;

        LSGFFunctor functorY(nb, 1, Y, image, imageWidth, x, y);
        Eigen::NumericalDiff<LSGFFunctor> numDiffY(functorY);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LSGFFunctor>, double> lmY(numDiffY);

        lmY.parameters.maxfev = 2000;
        lmY.parameters.xtol = 1.0e-10;

        lmX.minimize(betaX);
        lmY.minimize(betaY);

        a = betaX(0);

        xb = betaX(1);
        yb = betaY(1);

        sigma = betaX(2);

        // std::cout << "Original: " << x << ", " << y << std::endl;
        // std::cout << "final: " << xb << ", " << yb << std::endl;
        result.push_back(Star(xb, yb, 0));
        // result.push_back(Star(x, y, 0));
    }

    std::cout << "Number of centroids: " << result.size() << std::endl;

    return result;
}

std::vector<Star> LeastSquaresGaussianFit2D::Go(unsigned char *image, int imageWidth,
                                                int imageHeight) const {
    std::vector<Star> result;

    std::cout << "2D GAUSSIAN FIT" << std::endl;

    const int nb = 2;
    const int np = 2 * nb + 1;

    // std::set<Point> checkedPoints;

    std::vector<Point> candidatePts = FloodfillPreproc(image, imageWidth, imageHeight);

    // TODO: remove, testing
    int same = 0;

    for (const Point &pt : candidatePts) {
        int x = pt[0];
        int y = pt[1];
        if (x - nb < 0 || x + nb >= imageWidth || y - nb < 0 || y + nb >= imageHeight) continue;

        float a = Get(x, y, image, imageWidth);
        double sigma = InitialGuess2(x, y, a, nb, image, imageWidth);
        // InitialGuess(x, y, nb, image, imageWidth, &a, &xb, &yb, &sigmaX);

        Eigen::VectorXd beta(5);
        beta << a, x, y, sigma, sigma;

        // LSGF2DFunctor functor(nb, 0, X, image, imageWidth, x, y);
        LSGF2DFunctor functor(nb, image, imageWidth, x, y);
        // Eigen::NumericalDiff<LSGF2DFunctor> numDiff(functor);
        // Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LSGF2DFunctor>, double> lm(numDiff);
        Eigen::LevenbergMarquardt<LSGF2DFunctor, double> lm(functor);

        lm.parameters.maxfev = 2000;
        lm.parameters.xtol = 1.0e-10;

        lm.minimize(beta);

        a = beta(0);
        float xRes = beta(1);
        float yRes = beta(2);
        float sigmaX = beta(3);
        float sigmaY = beta(4);

        // std::cout << "Original: " << x << ", " << y << std::endl;
        std::cout << xRes << ", " << yRes << std::endl;
        // if(x == xRes && y == yRes) same++;
        result.push_back(Star(xRes+0.5, yRes+0.5, 0));


        // result.push_back(Star(x, y, nb));
    }

    std::cout << "Number of centroids: " << result.size() << std::endl;
    std::cout << "same: " << same << std::endl;

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

// a simple, but well tested thresholding algorithm that works well with star images
int BasicThreshold(unsigned char *image, int imageWidth, int imageHeight) {
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

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth,
                                               int imageHeight) const {
    CentroidParams p;

    std::vector<Star> result;

    p.cutoff = BasicThreshold(image, imageWidth, imageHeight);
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

            CogHelper(&p, i, image, imageWidth, imageHeight);
            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            // use the sums to finish CoG equation and add stars to the result
            float xCoord = (p.xCoordMagSum / (p.magSum * 1.0));
            float yCoord = (p.yCoordMagSum / (p.magSum * 1.0));

            if (p.isValid) {
                result.push_back(Star(xCoord + 0.5f, yCoord + 0.5f, ((float)(xDiameter)) / 2.0f,
                                      ((float)(yDiameter)) / 2.0f,
                                      p.checkedIndices.size() - sizeBefore));

                // std::cout << result.back().position.x << ", " << result.back().position.y << std::endl;
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

Stars IterativeWeightedCenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth,
                                                    int imageHeight) const {
    IWCoGParams p;
    std::vector<Star> result;
    p.cutoff = BasicThreshold(image, imageWidth, imageHeight);
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

            IWCoGHelper(&p, i, image, imageWidth, imageHeight, &starIndices);

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
