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

//I dont think this follows style guidelines but its the only way I could
//think about doing recursion
float yCoordMagSum = 0; 
float xCoordMagSum = 0;
int magSum = 0;
int xMin;
int xMax;
int cutoff;
std::unordered_set<int> checkedIndeces;

//recursive helper here

void cogHelper(int i, unsigned char *image, int imageWidth) {
    if (image[i] >= cutoff && checkedIndeces.count(i) == 0) {
        checkedIndeces.insert(i);
        if (i % imageWidth > xMax) {
            xMax = i % imageWidth;
        } else if (i % imageWidth < xMin) {
            xMin = i % imageWidth;
        }
        magSum += image[i];
        xCoordMagSum += ((i % imageWidth) + 1) * image[i];
        yCoordMagSum += ((i / imageWidth) + 1) * image[i];
        cogHelper(i + 1, image, imageWidth);
        cogHelper(i - 1, image, imageWidth);
        cogHelper(i + imageWidth, image, imageWidth);
    }
}

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    std::vector<Star> result;
    //loop through entire array, find sum of magnitudes
    int totalMag = 0;
    
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    cutoff = totalMag/(imageHeight * imageWidth) + 1;

    for (int i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= cutoff && checkedIndeces.count(i) == 0) {
            checkedIndeces.insert(i);
            //iterate over pixels that are part of the star

            int radius = 0; //radius of current star
            yCoordMagSum = 0; //y coordinate of current star
            xCoordMagSum = 0; //x coordinate of current star
            magSum = 0; //sum of magnitudes of current star

            magSum += image[i];
            xMax = i % imageWidth;
            xMin = i % imageWidth;
            xCoordMagSum += ((i % imageWidth) + 1) * image[i];
            yCoordMagSum += ((i / imageWidth) + 1) * image[i];
            cogHelper(i + 1, image, imageWidth);
            cogHelper(i - 1, image, imageWidth);
            cogHelper(i + imageWidth, image, imageWidth);
            radius = (xMax - xMin) + 1;

            //use the sums to finish CoG equation and add stars to the result
            float xCoord = (xCoordMagSum / (magSum * 1.0));      
            float yCoord = (yCoordMagSum / (magSum * 1.0));
            result.push_back(Star(xCoord, yCoord, ((double)(radius * 1.0))/2));
            i += radius;
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
