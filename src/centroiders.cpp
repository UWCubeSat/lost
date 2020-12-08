#include "centroiders.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <iostream>

#include <unordered_map>

#include <unordered_set>

namespace lost {

// DUMMY

std::vector<Star> DummyCentroidAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    std::vector<Star> result;

    for (int i = 0; i < numStars; i++) {
        result.push_back(Star(rand() % imageWidth, rand() % imageHeight, 10.0));
    }

    return result;
}

//recursive helper here
void cogHelper(float &yCoordMagSum,  
    float &xCoordMagSum, 
    int &magSum, 
    int &xMin, 
    int &xMax, 
    int &yMin, 
    int &yMax, 
    int &cutoff, 
    std::unordered_set<int> &checkedIndeces, int i, unsigned char *image, int imageWidth, int imageHeight) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= cutoff && checkedIndeces.count(i) == 0) {
        //std::cout << i << "\n";
        checkedIndeces.insert(i);
        if (i % imageWidth > xMax) {
            xMax = i % imageWidth;
        } else if (i % imageWidth < xMin) {
            xMin = i % imageWidth;
        }
        if (i / imageWidth > yMax) {
            yMax = i / imageWidth;
        } else if (i / imageWidth < yMin) {
            yMin = i / imageWidth;
        }
        magSum += image[i];
        xCoordMagSum += ((i % imageWidth) + 1) * image[i];
        yCoordMagSum += ((i / imageWidth) + 1) * image[i];
        if((i + 1) % imageWidth != 0) {
            cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i + 1, image, imageWidth, imageHeight);
        }
        if ((i - 1) % imageWidth != (imageWidth - 1)) {
            cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i - 1, image, imageWidth, imageHeight);
        }
        cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i + imageWidth, image, imageWidth, imageHeight);
        cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i - imageWidth, image, imageWidth, imageHeight);
    }
}

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    float yCoordMagSum = 0; 
    float xCoordMagSum = 0;
    int magSum = 0;
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int cutoff;
    std::unordered_set<int> checkedIndeces;
    
    std::vector<Star> result;
    //loop through entire array, find sum of magnitudes
    int totalMag = 0;
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    // cutoff might need a new equation
    cutoff = (((totalMag/(imageHeight * imageWidth)) + 1) * 15) / 10;
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= cutoff && checkedIndeces.count(i) == 0) {
            checkedIndeces.insert(i);
            //iterate over pixels that are part of the star

            int xDiameter = 0; //radius of current star
            int yDiameter = 0;
            yCoordMagSum = 0; //y coordinate of current star
            xCoordMagSum = 0; //x coordinate of current star
            magSum = 0; //sum of magnitudes of current star

            //computes indices to skip after done w current star
            int j = i;
            while(j < imageWidth * imageHeight && (j + 1) % imageWidth != 0 && image[j] >= cutoff) {
                j++;
            }

            magSum += image[i];
            xMax = i % imageWidth;
            xMin = i % imageWidth;
            yMax = i / imageWidth;
            yMin = i / imageWidth;
            xCoordMagSum += ((i % imageWidth) + 1) * image[i];
            yCoordMagSum += ((i / imageWidth) + 1) * image[i];
            if((i + 1) % imageWidth != 0) {
                cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i + 1, image, imageWidth, imageHeight);
            }
            if ((i - 1) % imageWidth != (imageWidth - 1)) {
                cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i - 1, image, imageWidth, imageHeight);
            }
            cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i + imageWidth, image, imageWidth, imageHeight);
            cogHelper(yCoordMagSum, xCoordMagSum, magSum, xMin, xMax, yMin, yMax, cutoff, checkedIndeces, i - imageWidth, image, imageWidth, imageHeight);
            xDiameter = (xMax - xMin) + 1;
            yDiameter = (yMax - yMin) + 1;
            //use the sums to finish CoG equation and add stars to the result
            float xCoord = (xCoordMagSum / (magSum * 1.0));      
            float yCoord = (yCoordMagSum / (magSum * 1.0));
            result.push_back(Star(xCoord, yCoord, ((double)(xDiameter * 1.0))/2.0, ((double)(yDiameter * 1.0))/2.0, 0));
            i = j - 1;
        }
    }
    return result;
}

// OTHER CENTROID RELATED FUNCTIONS

float StarDistancePixels(Star one, Star two) {
    float distX = one.x - two.x;
    float distY = one.y - two.y;
    return sqrt(distX*distX + distY*distY);
}

// helper for StarCentroidsCompare
static std::vector<int> FindClosestCentroids(float threshold,
                                              std::vector<Star> one,
                                              std::vector<Star> two) {
    std::vector<int> result;

    for (int i = 0; i < (int)one.size(); i++) {
        float closestDistance = INFINITY;
        int   closestIndex = -1;

        for (int k = 0; k < (int)two.size(); k++) {
            float currDistance = StarDistancePixels(one[i], two[k]);
            if (currDistance < threshold && currDistance < closestDistance) {
                closestDistance = currDistance;
                closestIndex = k;
            }
        }
        
        result.push_back(closestIndex);
    }

    return result;
}

CentroidComparison StarCentroidsCompare(float threshold,
                                    std::vector<Star> expected,
                                    std::vector<Star> actual) {

    CentroidComparison result;
    // maps from indexes in each list to the closest centroid from other list
    std::vector<int> expectedToActual, actualToExpected;

    expectedToActual = FindClosestCentroids(threshold, expected, actual);
    actualToExpected = FindClosestCentroids(threshold, actual, expected);

    // any expected stars whose closest actual star does not refer to them are missing
    for (int i = 0; i < (int)expectedToActual.size(); i++) {
        if (expectedToActual[i] == -1 ||
            actualToExpected[expectedToActual[i]] != i) {
            result.numMissingStars++;
        } else {
            result.meanError += StarDistancePixels(expected[i], actual[expectedToActual[i]]);
        }
    }
    result.meanError /= (expected.size() - result.numMissingStars);


    // any actual star whose closest expected star does not refer to them is extra
    for (int i = 0; i < (int)actual.size(); i++) {
        if (actualToExpected[i] == -1 ||
            expectedToActual[actualToExpected[i]] != i) {
            result.numExtraStars++;
        }
    }

    return result;
}

void CentroidComparison::Print() {
    printf("Extra stars: %d\nMissing stars: %d\nMean error: %f\n",
           numExtraStars, numMissingStars, meanError);
}

}
