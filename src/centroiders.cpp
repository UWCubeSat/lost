#include "centroiders.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
#include <eigen3/unsupported/Eigen/NumericalDiff>
#include <iostream>
#include <set>  // TODO: remove later, this is just lazy to implement hash for unordered_set
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace lost {

int BasicThreshold(unsigned char *image, int imageWidth, int imageHeight);

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

struct LSGF1DFunctor : Functor<double> {
 public:
  LSGF1DFunctor(const Eigen::VectorXd X, int x0, int y0, int w, int h, const unsigned char *image)
      : Functor<double>(3, 3), X(X), x0(x0), y0(y0), w(w), h(h), image(image) {}

  int operator()(const Eigen::VectorXd &params, Eigen::VectorXd &fvec) const {
    double a = params(0);
    int xb = params(1);
    double sigma = params(2);
    for (int i = 0; i < X.size(); i++) {
      fvec(i) = XMarginal(X(i)) - Model(X(i), a, xb, sigma);
    }
    return 0;
  }

 private:
  const Eigen::VectorXd X;
  // Eigen::VectorXd Y; // Don't need, calculate Y_i = V_{m, x_i}
  const int nb = 2;
  const int x0, y0;  // coordinates of centroid window
  const int w, h;    // w = number of pixels in width, h = number of pixels in height
  const unsigned char *image;

  // Get value of pixel at universal coordinates (x, y)
  // TODO: is this correct representation?
  // TODO: will need to deal with case where on boundary of acceptable nb
  unsigned char Get(int x, int y) const {
    // int ind = i * w + j;
    int ind = x + y * w;
    return image[ind];
  }

  // Assuming 0, 0 is the center
  int XMarginal(int x) const {
    int res = 0;
    for (int j = -nb; j <= nb; j++) {
      res += Get(x0 + x, y0 + j);
    }
    return res;
  }

  double Model(int xi, double a, int xb, double sigma) const {
    return a * exp(-1 * (xi - xb) * (xi - xb) / (2 * sigma * sigma));
  }
};

typedef std::vector<int> Point;

// TODO: duplicate, factor out better
unsigned char Get(int x, int y, unsigned char *image, int w, int h) {
  int ind = x + y * w;
  return image[ind];
}

void InitialGuess(int x, int y, unsigned char *image, int w, int h, double *a, int *xb,
                 double *sigma) {
  // a is set to max intensity value in the window
  // (xb, yb) = coordinates of pixel with max intensity value
  const int nb = 2;
  int max = -1;
  for (int i = -nb; i <= nb; i++) {
    for (int j = -nb; j <= nb; j++) {
      int pixelValue = Get(x + i, y + j, image, w, h);
      if (pixelValue > max) {
        *xb = x + i;
        max = pixelValue;
      }
    }
  }
  *a = max;
  // Sigma is kinda complicated
  // sigma = f_whm / (2 * \sqrt{2log(2)})
  // f_whm \approx sqrt of number of pixels with intensity > 0.5 * a
  int halfCount = 0;
  for (int i = -nb; i <= nb; i++) {
    for (int j = -nb; j <= nb; j++) {
      if(Get(x + i, y + j, image, w, h) > *a / 2){
        halfCount++;
      }
    }
  }
  double fwhm = std::sqrt(halfCount);
  *sigma = fwhm / (2 * std::sqrt(2 * std::log(2)));
}

std::vector<Star> LeastSquaresGaussianFit1D::Go(unsigned char *image, int imageWidth,
                                                int imageHeight) const {
  std::vector<Star> result;

  const int cb = 2;
  std::set<Point> checkedPoints;

  int cutoff = BasicThreshold(image, imageWidth, imageHeight);
  for (int i = 0; i < imageHeight * imageWidth; i++) {
    int x = i / imageWidth;
    int y = i % imageHeight;
    Point p{x, y};
    if (x - cb < 0 || x + cb >= imageWidth || y - cb < 0 || y + cb >= imageHeight) continue;
    // TODO: check that entire window is not in checkedPoints? - no, just check the point itself
    if (image[i] >= cutoff && checkedPoints.count(p) == 0) {
      checkedPoints.insert(p);

      std::vector<int> xs;
      std::vector<int> ys;
      Eigen::VectorXd X(2 * cb + 1);
      for (int i = -cb; i <= cb; i++) {
        xs.push_back(x + cb);
        X << (x + cb);
        ys.push_back(y + cb);
      }

      // LSGF1DFunctor(const Eigen::VectorXd X, int x0, int y0, int w, int h,
      //               const unsigned char *image)
      //     : X(X), x0(x0), y0(y0), w(w), h(h), image(image) {}

      Eigen::VectorXd beta(3);  // initial guess
      double a, sigma;
      int xb;
      InitialGuess(x, y, image, imageWidth, imageHeight, &a, &xb, &sigma);
      // TODO: fill beta with described values
      beta << a, xb, sigma;

      LSGF1DFunctor functor(X, x, y, imageWidth, imageHeight, image);
      Eigen::NumericalDiff<LSGF1DFunctor> numDiff(functor);
      Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LSGF1DFunctor>, double> lm(numDiff);

      lm.parameters.maxfev = 1000;
      lm.parameters.xtol = 1.0e-10;

      lm.minimize(beta);

      a = beta(0);
      xb = beta(1);
      sigma = beta(2);

      std::cout << a << ", " << xb << ", " << sigma << std::endl;
    }
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
                              ((float)(yDiameter)) / 2.0f, p.checkedIndices.size() - sizeBefore));
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
        // traverse through star indices, calculate W at each coordinate, add to final coordinate
        // sums
        yWeightedCoordMagSum = 0;
        xWeightedCoordMagSum = 0;
        weightedMagSum = 0;
        stop++;
        for (long j = 0; j < (long)starIndices.size(); j++) {
          // calculate w
          float currXCoord = (float)(starIndices.at(j) % imageWidth);
          float currYCoord = (float)(starIndices.at(j) / imageWidth);
          w = p.maxIntensity * exp(-1.0 * ((pow(currXCoord - guessXCoord, 2) / modifiedStdDev) +
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
        result.push_back(Star(guessXCoord + 0.5f, guessYCoord + 0.5f, ((float)(xDiameter)) / 2.0f,
                              ((float)(yDiameter)) / 2.0f, starIndices.size()));
      }
    }
  }
  return result;
}

}  // namespace lost
