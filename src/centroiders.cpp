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
    
    std::vector<Star> result;
    
    p.cutoff = DetermineCutoff(image, imageWidth, imageHeight);
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
            //iterate over pixels that are part of the star
            int xDiameter = 0; //radius of current star
            int yDiameter = 0;
            p.yCoordMagSum = 0; //y coordinate of current star
            p.xCoordMagSum = 0; //x coordinate of current star
            p.magSum = 0; //sum of magnitudes of current star

            //computes indices to skip after done w current star

            p.xMax = i % imageWidth;
            p.xMin = i % imageWidth;
            p.yMax = i / imageWidth;
            p.yMin = i / imageWidth;

            CogHelper(p, i, image, imageWidth, imageHeight);
            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;
            //use the sums to finish CoG equation and add stars to the result
            float xCoord = (p.xCoordMagSum / (p.magSum * 1.0));      
            float yCoord = (p.yCoordMagSum / (p.magSum * 1.0));
            result.push_back(Star(xCoord, yCoord, ((float)(xDiameter * 1.0))/2.0, ((float)(yDiameter * 1.0))/2.0, 0));
        }
    }
    return result;
}

//Determines how accurate and how much iteration is done by the IWCoG algorithm,
//smaller means more accurate and more iterations.
float iWCoGMinChange = 0.00001;

struct IWCoGParams {
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int cutoff;
    int maxIntensity;
    int guess;
    std::unordered_set<int> checkedIndices;
};

void IWCoGHelper(IWCoGParams &p, int i, unsigned char *image, int imageWidth, int imageHeight, std::vector<int> &starIndices) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
        p.checkedIndices.insert(i);
        starIndices.push_back(i);
        if (image[i] > p.maxIntensity) {
            p.maxIntensity = image[i];
            p.guess = i;
        }
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
        if(i % imageWidth != imageWidth - 1) {
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
    p.cutoff = DetermineCutoff(image, imageWidth, imageHeight);
    for (int i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
            std::vector<int> starIndices; //indices of the current star
            p.maxIntensity = 0;
            int xDiameter = 0; 
            int yDiameter = 0;
            float yWeightedCoordMagSum = 0; 
            float xWeightedCoordMagSum = 0; 
            float weightedMagSum = 0; 
            float fwhm; //fwhm variable
            float standardDeviation;
            float w; //weight value

            p.xMax = i % imageWidth;
            p.xMin = i % imageWidth;
            p.yMax = i / imageWidth;
            p.yMin = i / imageWidth;

            IWCoGHelper(p, i, image, imageWidth, imageHeight, starIndices);

            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            //calculate fwhm
            int count = 0;
            for (int j = 0; j < (int) starIndices.size(); j++) {
                if (image[starIndices.at(j)] > p.maxIntensity / 2) {
                    count++;
                }
            }
            fwhm = sqrt((float) count);
            standardDeviation = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
            float guessXCoord = (float) (p.guess % imageWidth);
            float guessYCoord = (float) (p.guess / imageWidth);
            //how much our new centroid estimate changes w each iteration
            float change = INFINITY;
            //while we see some large enough change in estimated, maybe make it a global variable
            while (change > 0.0001) {
            //traverse through star indices, calculate W at each coordinate, add to final coordinate sums
                for (int j = 0; j < (int) starIndices.size(); j++) {
                    //calculate w
                    float currXCoord = (float) (starIndices.at(j) % imageWidth);
                    float currYCoord = (float) (starIndices.at(j) / imageWidth);
                    float modifiedStdDev = 2.0 * pow(standardDeviation, 2);
                    w = p.maxIntensity * exp(-1.0 * ((pow(currXCoord - guessXCoord, 2) / modifiedStdDev) + (pow(currYCoord - guessYCoord, 2) / modifiedStdDev)));

                    xWeightedCoordMagSum += w * currXCoord * ((float) image[starIndices.at(j)]);
                    yWeightedCoordMagSum += w * currYCoord * ((float) image[starIndices.at(j)]);
                    weightedMagSum += w * ((float) image[starIndices.at(j)]);

                    float xTemp = xWeightedCoordMagSum / weightedMagSum;
                    float yTemp = yWeightedCoordMagSum / weightedMagSum;

                    change = abs(guessXCoord - xTemp) + abs(guessYCoord - yTemp);

                    guessXCoord = xTemp;
                    guessYCoord = yTemp;
                }
            }
            result.push_back(Star(guessXCoord, guessYCoord, ((float)(xDiameter * 1.0))/2.0, ((float)(yDiameter * 1.0))/2.0, 0));
        }
    }
    return result;
}

}
