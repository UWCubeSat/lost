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

struct CentroidParams {
    float yCoordMagSum; 
    float xCoordMagSum;
    int magSum;
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int cutoff;
    std::unordered_set<int> checkedIndeces;
};

//recursive helper here
void CogHelper(CentroidParams &p, int i, unsigned char *image, int imageWidth, int imageHeight) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p.cutoff && p.checkedIndeces.count(i) == 0) {
        p.checkedIndeces.insert(i);
        if (i % imageWidth > p.xMax) {
            p.xMax = i % imageWidth;
        } else if (i % imageWidth < p.xMin) {
            p.xMin = i % imageWidth;
        }
        if (i / imageWidth > p.yMax) {
            p.yMax = i / imageWidth;
        } else if (i / imageWidth < p.yMin) {
            p.yMin = i / imageWidth;
        }
        p.magSum += image[i];
        p.xCoordMagSum += ((i % imageWidth) + 1) * image[i];
        p.yCoordMagSum += ((i / imageWidth) + 1) * image[i];
        if((i + 1) % imageWidth != 0) {
            CogHelper(p, i + 1, image, imageWidth, imageHeight);
        }
        if ((i - 1) % imageWidth != (imageWidth - 1)) {
            CogHelper(p, i - 1, image, imageWidth, imageHeight);
        }
        CogHelper(p, i + imageWidth, image, imageWidth, imageHeight);
        CogHelper(p, i - imageWidth, image, imageWidth, imageHeight);
    }
}

int DetermineCutoff(unsigned char *image, int imageWidth, int imageHeight) {
    //loop through entire array, find sum of magnitudes
    int totalMag = 0;
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    return (((totalMag/(imageHeight * imageWidth)) + 1) * 15) / 10;
}

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    CentroidParams p;
    std::unordered_set<int> checkedIndeces;
    
    std::vector<Star> result;
    
    p.cutoff = DetermineCutoff(image, imageWidth, imageHeight);
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && p.checkedIndeces.count(i) == 0) {
            checkedIndeces.insert(i);
            //iterate over pixels that are part of the star

            int xDiameter = 0; //radius of current star
            int yDiameter = 0;
            p.yCoordMagSum = 0; //y coordinate of current star
            p.xCoordMagSum = 0; //x coordinate of current star
            p.magSum = 0; //sum of magnitudes of current star

            //computes indices to skip after done w current star
            int j = i;
            while(j < imageWidth * imageHeight && (j + 1) % imageWidth != 0 && image[j] >= p.cutoff) {
                j++;
            }

            p.magSum += image[i];
            p.xMax = i % imageWidth;
            p.xMin = i % imageWidth;
            p.yMax = i / imageWidth;
            p.yMin = i / imageWidth;
            p.xCoordMagSum += ((i % imageWidth)) * image[i];
            p.yCoordMagSum += ((i / imageWidth)) * image[i];

            CogHelper(p, i, image, imageWidth, imageHeight);
            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;
            //use the sums to finish CoG equation and add stars to the result
            float xCoord = (p.xCoordMagSum / (p.magSum * 1.0));      
            float yCoord = (p.yCoordMagSum / (p.magSum * 1.0));
            result.push_back(Star(xCoord, yCoord, ((float)(xDiameter * 1.0))/2.0, ((double)(yDiameter * 1.0))/2.0, 0));
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
