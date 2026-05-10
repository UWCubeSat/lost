#include "centroiders.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "decimal.hpp"

namespace lost {

// DUMMY

Stars DummyCentroidAlgorithm::Go(unsigned char *, int imageWidth, int imageHeight) const {
    Stars result;

    unsigned int randomSeed = 123456;
    for (int i = 0; i < numStars; i++) {
        result.push_back(Star(rand_r(&randomSeed) % imageWidth, rand_r(&randomSeed) % imageHeight, DECIMAL(10.0)));
    }

    return result;
}

// a poorly designed thresholding algorithm
int BadThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    //loop through entire array, find sum of magnitudes
    long totalMag = 0;
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    return (((totalMag/(imageHeight * imageWidth)) + 1) * 15) / 10;
}

// a more sophisticated thresholding algorithm, not tailored to star images
int OtsusThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    // code here, duh
    long total = imageWidth * imageHeight;
    //decimal top = 255;
    decimal sumB = 0;
    decimal sum1 = 0;
    decimal wB = 0;
    decimal maximum = 0;
    int level = 0;
    // make the histogram (array length 256)
    int histogram[256];

    memset(histogram, 0, sizeof(int)*256);

    for (long i = 0; i < total; i++) {
        histogram[image[i]]++;
    }
    for (int i = 0; i < 256; i ++) {
        sum1 += i * histogram[i];
    }
    for (int i = 0; i < 256; i ++) {
        decimal wF = total - wB;
        //std::cout << "wF\n" << wB << "\n";
        //std::cout << "wB\n" << wF << "\n";
        if (wB > 0 && wF > 0) {
            decimal mF = (sum1 - sumB) / wF;
            decimal val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
            //std::cout << val << "\n";
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
    decimal std = 0;
    long totalPixels = imageHeight * imageWidth;
    for (long i = 0; i < totalPixels; i++) {
        totalMag += image[i];
    }
    decimal mean = totalMag / totalPixels;
    for (long i = 0; i < totalPixels; i++) {
        std += DECIMAL_POW(image[i] - mean, 2);
    }
    std = DECIMAL_SQRT(std / totalPixels);
    return mean + (std * 5);
}

// basic thresholding, but do it faster (trade off of some accuracy?)
int BasicThresholdOnePass(unsigned char *image, int imageWidth, int imageHeight) {
    unsigned long totalMag = 0;
    decimal std = 0;
    decimal sq_totalMag = 0;
    long totalPixels = imageHeight * imageWidth;
    for (long i = 0; i < totalPixels; i++) {
        totalMag += image[i];
        sq_totalMag += image[i] * image[i];
    }
    decimal mean = totalMag / totalPixels;
    decimal variance = (sq_totalMag / totalPixels) - (mean * mean);
    std = DECIMAL_SQRT(variance);
    return mean + (std * 5);
}

struct CentroidParams {
    decimal yCoordMagSum;
    decimal xCoordMagSum;
    long magSum;
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int cutoff;
    bool isValid;
    size_t checkedCount;
};

template <typename VisitPixelFunc>
void TraverseConnectedPixels(long startIndex,
                             unsigned char *image,
                             int imageWidth,
                             int imageHeight,
                             int cutoff,
                             vector<unsigned char, LOST_ETL_MAX_IMAGE_PIXELS> *checked,
                             size_t *checkedCount,
                             bool *isValid,
                             VisitPixelFunc visitPixel) {
    vector<int, LOST_ETL_MAX_IMAGE_PIXELS> stack;
    stack.push_back(static_cast<int>(startIndex));
    long imageSize = static_cast<long>(imageWidth) * imageHeight;

    while (!stack.empty()) {
        int i = stack.back();
        stack.pop_back();

        if (i < 0 || i >= imageSize) {
            continue;
        }
        if (image[i] < cutoff || (*checked)[i] != 0) {
            continue;
        }

        int x = i % imageWidth;
        int y = i / imageWidth;
        // Reject components touching the image border as before.
        if (x == 0 || x == imageWidth - 1 || y == 0 || y == imageHeight - 1) {
            *isValid = false;
        }

        (*checked)[i] = 1;
        (*checkedCount)++;
        visitPixel(i, x, y);

        if (x < imageWidth - 1) {
            stack.push_back(i + 1);
        }
        if (x > 0) {
            stack.push_back(i - 1);
        }
        if (y < imageHeight - 1) {
            stack.push_back(i + imageWidth);
        }
        if (y > 0) {
            stack.push_back(i - imageWidth);
        }
    }
}

Stars CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    CentroidParams p;

    Stars result;
    EtlRuntimeBoundCheck((size_t)imageWidth * (size_t)imageHeight,
                         LOST_ETL_MAX_IMAGE_PIXELS,
                         "centroid image pixel count");
    vector<unsigned char, LOST_ETL_MAX_IMAGE_PIXELS> checked(imageHeight * imageWidth, 0);
    p.checkedCount = 0;

    p.cutoff = BasicThreshold(image, imageWidth, imageHeight);
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        if (image[i] >= p.cutoff && checked[i] == 0) {

            //iterate over pixels that are part of the star
            int xDiameter = 0; //radius of current star
            int yDiameter = 0;
            p.yCoordMagSum = 0; //y coordinate of current star
            p.xCoordMagSum = 0; //x coordinate of current star
            p.magSum = 0; //sum of magnitudes of current star

            p.xMax = i % imageWidth;
            p.xMin = i % imageWidth;
            p.yMax = i / imageWidth;
            p.yMin = i / imageWidth;
            p.isValid = true;

            size_t sizeBefore = p.checkedCount;

            TraverseConnectedPixels(
                i,
                image,
                imageWidth,
                imageHeight,
                p.cutoff,
                &checked,
                &p.checkedCount,
                &p.isValid,
                [&](long index, int x, int y) {
                    if (x > p.xMax) {
                        p.xMax = x;
                    } else if (x < p.xMin) {
                        p.xMin = x;
                    }
                    if (y > p.yMax) {
                        p.yMax = y;
                    } else if (y < p.yMin) {
                        p.yMin = y;
                    }
                    p.magSum += image[index];
                    p.xCoordMagSum += x * image[index];
                    p.yCoordMagSum += y * image[index];
                });
            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            //use the sums to finish CoG equation and add stars to the result
            decimal xCoord = (p.xCoordMagSum / (p.magSum * DECIMAL(1.0)));
            decimal yCoord = (p.yCoordMagSum / (p.magSum * DECIMAL(1.0)));

            if (p.isValid) {
                result.push_back(Star(xCoord + DECIMAL(0.5), yCoord + DECIMAL(0.5), (xDiameter)/DECIMAL(2.0), (yDiameter)/DECIMAL(2.0), p.checkedCount - sizeBefore));
            }
        }
    }
    return result;
}

//Determines how accurate and how much iteration is done by the IWCoG algorithm,
//smaller means more accurate and more iterations.
decimal iWCoGMinChange = DECIMAL(0.0002);

struct IWCoGParams {
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int cutoff;
    int maxIntensity;
    int guess;
    bool isValid;
    size_t checkedCount;
};

Stars IterativeWeightedCenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    IWCoGParams p;
    Stars result;
    EtlRuntimeBoundCheck((size_t)imageWidth * (size_t)imageHeight,
                         LOST_ETL_MAX_IMAGE_PIXELS,
                         "iwcog image pixel count");
    vector<unsigned char, LOST_ETL_MAX_IMAGE_PIXELS> checked(imageHeight * imageWidth, 0);
    p.checkedCount = 0;
    p.cutoff = BasicThreshold(image, imageWidth, imageHeight);
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && checked[i] == 0) {
            // TODO: store longs --Mark
            vector<int, LOST_ETL_MAX_STARS> starIndices; //indices of the current star
            p.maxIntensity = 0;
            int xDiameter = 0;
            int yDiameter = 0;
            decimal yWeightedCoordMagSum = 0;
            decimal xWeightedCoordMagSum = 0;
            decimal weightedMagSum = 0;
            decimal fwhm; //fwhm variable
            decimal standardDeviation;
            decimal w; //weight value

            p.xMax = i % imageWidth;
            p.xMin = i % imageWidth;
            p.yMax = i / imageWidth;
            p.yMin = i / imageWidth;
            p.isValid = true;


            TraverseConnectedPixels(
                i,
                image,
                imageWidth,
                imageHeight,
                p.cutoff,
                &checked,
                &p.checkedCount,
                &p.isValid,
                [&](long index, int x, int y) {
                    starIndices.push_back(static_cast<int>(index));
                    if (image[index] > p.maxIntensity) {
                        p.maxIntensity = image[index];
                        p.guess = index;
                    }
                    if (x > p.xMax) {
                        p.xMax = x;
                    } else if (x < p.xMin) {
                        p.xMin = x;
                    }
                    if (y > p.yMax) {
                        p.yMax = y;
                    } else if (y < p.yMin) {
                        p.yMin = y;
                    }
                });

            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            //calculate fwhm
            decimal count = 0;
            for (int j = 0; j < (int) starIndices.size(); j++) {
                if (image[starIndices.at(j)] > p.maxIntensity / 2) {
                    count++;
                }
            }
            fwhm = DECIMAL_SQRT(count);
            standardDeviation = fwhm / (DECIMAL(2.0) * DECIMAL_SQRT(DECIMAL(2.0) * DECIMAL_LOG(2.0)));
            decimal modifiedStdDev = DECIMAL(2.0) * DECIMAL_POW(standardDeviation, 2);
            // TODO: Why are these decimals? --Mark
            decimal guessXCoord = (p.guess % imageWidth);
            decimal guessYCoord = (p.guess / imageWidth);
            //how much our new centroid estimate changes w each iteration
            decimal change = INFINITY;
            int stop = 0;
            //while we see some large enough change in estimated, maybe make it a global variable
            while (change > iWCoGMinChange && stop < 100000) {
            //traverse through star indices, calculate W at each coordinate, add to final coordinate sums
                yWeightedCoordMagSum = 0;
                xWeightedCoordMagSum = 0;
                weightedMagSum = 0;
                stop++;
                for (long j = 0; j < (long)starIndices.size(); j++) {
                    //calculate w
                    decimal currXCoord = starIndices.at(j) % imageWidth;
                    decimal currYCoord = starIndices.at(j) / imageWidth;
                    w = p.maxIntensity * DECIMAL_EXP(DECIMAL(-1.0) * ((DECIMAL_POW(currXCoord - guessXCoord, 2) / modifiedStdDev) + (DECIMAL_POW(currYCoord - guessYCoord, 2) / modifiedStdDev)));

                    xWeightedCoordMagSum += w * currXCoord * DECIMAL(image[starIndices.at(j)]);
                    yWeightedCoordMagSum += w * currYCoord * DECIMAL(image[starIndices.at(j)]);
                    weightedMagSum += w * DECIMAL(image[starIndices.at(j)]);
                }
                decimal xTemp = xWeightedCoordMagSum / weightedMagSum;
                decimal yTemp = yWeightedCoordMagSum / weightedMagSum;

                change = abs(guessXCoord - xTemp) + abs(guessYCoord - yTemp);

                guessXCoord = xTemp;
                guessYCoord = yTemp;
            }
            if (p.isValid) {
                result.push_back(Star(guessXCoord + DECIMAL(0.5), guessYCoord + DECIMAL(0.5), xDiameter/DECIMAL(2.0), yDiameter/DECIMAL(2.0), starIndices.size()));
            }
        }
    }
    return result;
}

}
