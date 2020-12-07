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

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    std::vector<Star> result;
    //loop through entire array, find sum of magnitudes
    int totalMag;
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    int cutoff = totalMag/(imageHeight * imageWidth) + 1;
    cutoff = 150;

    //std::unordered_map<int, std::unordered_set<int>> indicesAlreadyChecked;
    std::unordered_set<int> checkedIndeces;

    for (int i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= cutoff && checkedIndeces.count(i) == 0) {
            checkedIndeces.insert(i);
            //iterate over pixels that are part of the star

            int radius = 0; //radius of current star
            int magSum; //current magnitude sum of star pixels
            float yCoordMagSum; //sum of magnitudes * their coordinate
            float xCoordMagSum;

            int j = i;
            
            //find radius
            while (image[j] >= cutoff) {
                radius++;
                j++;
            }

            //iterate over star pixels within the radius to make a square
            for (int k = 0; k < radius; k++) {
                for (int l = 0; l < radius; l++) {
                    magSum += image[i + k + (imageWidth * l)];
                    xCoordMagSum += (((i + k + (imageWidth * l)) % imageWidth) + 1) * image[i + k + (imageWidth * l)];
                    yCoordMagSum += (((i + k + (imageWidth * l)) / imageHeight) + 1) * image[i + k + (imageWidth * l)];
                }
            }

            //use the sums to finish CoG equation and add stars to the result
            float xCoord = (xCoordMagSum / (magSum * 1.0)) - 1.0;
            std::cout << xCoord;
            std::cout << "\n";
            
            float yCoord = (yCoordMagSum / (magSum * 1.0)) - 1.0;
            std::cout << yCoord;
            std::cout << "\n";

            //Star *currentStar = new Star(xCoord, yCoord, (double)(radius * 1.0));
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
