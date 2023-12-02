#include "centroiders.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "decimal.hpp"

namespace lost {

// DUMMY

std::vector<Star> DummyCentroidAlgorithm::Go(unsigned char *, int imageWidth, int imageHeight) const {
    std::vector<Star> result;

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
    std::unordered_set<int> checkedIndices;
};

//recursive helper here
void CogHelper(CentroidParams *p, long i, unsigned char *image, int imageWidth, int imageHeight) {

    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p->cutoff && p->checkedIndices.count(i) == 0) {
        //check if pixel is on the edge of the image, if it is, we dont want to centroid this star
        if (i % imageWidth == 0 || i % imageWidth == imageWidth - 1 || i / imageWidth == 0 || i / imageWidth == imageHeight - 1) {
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

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    CentroidParams p;

    std::vector<Star> result;

    p.cutoff = BasicThreshold(image, imageWidth, imageHeight);
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {

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

            int sizeBefore = p.checkedIndices.size();

            CogHelper(&p, i, image, imageWidth, imageHeight);
            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            //use the sums to finish CoG equation and add stars to the result
            decimal xCoord = (p.xCoordMagSum / (p.magSum * DECIMAL(1.0)));
            decimal yCoord = (p.yCoordMagSum / (p.magSum * DECIMAL(1.0)));

            if (p.isValid) {
                result.push_back(Star(xCoord + DECIMAL(0.5), yCoord + DECIMAL(0.5), (xDiameter)/DECIMAL(2.0), (yDiameter)/DECIMAL(2.0), p.checkedIndices.size() - sizeBefore));
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
    std::unordered_set<int> checkedIndices;
};

void IWCoGHelper(IWCoGParams *p, long i, unsigned char *image, int imageWidth, int imageHeight, std::vector<int> *starIndices) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p->cutoff && p->checkedIndices.count(i) == 0) {
        //check if pixel is on the edge of the image, if it is, we dont want to centroid this star
        if (i % imageWidth == 0 || i % imageWidth == imageWidth - 1 || i / imageWidth == 0 || i / imageWidth == imageHeight - 1) {
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

Stars IterativeWeightedCenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    IWCoGParams p;
    std::vector<Star> result;
    p.cutoff = BasicThreshold(image, imageWidth, imageHeight);
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
            // TODO: store longs --Mark
            std::vector<int> starIndices; //indices of the current star
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


            IWCoGHelper(&p, i, image, imageWidth, imageHeight, &starIndices);

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
