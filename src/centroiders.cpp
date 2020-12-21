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
    std::unordered_set<int> checkedIndices;
};

//recursive helper here
void CogHelper(CentroidParams &p, int i, unsigned char *image, int imageWidth, int imageHeight) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
        p.checkedIndices.insert(i);
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
        p.xCoordMagSum += ((i % imageWidth)) * image[i];
        p.yCoordMagSum += ((i / imageWidth)) * image[i];
        if(i % imageWidth != imageWidth - 1) {
            CogHelper(p, i + 1, image, imageWidth, imageHeight);
        }
        if (i % imageWidth != 0) {
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
    std::unordered_set<int> checkedIndices;
    
    std::vector<Star> result;
    
    p.cutoff = DetermineCutoff(image, imageWidth, imageHeight);
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
            checkedIndices.insert(i);
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
            result.push_back(Star(xCoord, yCoord, ((float)(xDiameter * 1.0))/2.0, ((float)(yDiameter * 1.0))/2.0, 0));
            i = j - 1;
        }
    }
    return result;
}

}
