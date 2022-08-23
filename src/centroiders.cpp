#include "centroiders.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <iostream>

#include <unordered_map>

#include <unordered_set>
#include <string.h>
#include <cassert>

namespace lost {

// DUMMYS
std::vector<Star> DummyCentroidAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    std::vector<Star> result;

    for (int i = 0; i < numStars; i++) {
        result.push_back(Star(rand() % imageWidth, rand() % imageHeight, 10.0));
    }

    return result;
}

// UNUSED: a poorly designed thresholding algorithm
int BadThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    //loop through entire array, find sum of magnitudes
    long totalMag = 0;
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        totalMag += image[i];
    }
    return (((totalMag/(imageHeight * imageWidth)) + 1) * 15) / 10;
}

// UNUSED: a more sophisticated thresholding algorithm, not tailored to star images
int OtsusThreshold(unsigned char *image, int imageWidth, int imageHeight) {
    // code here, duh
    long total = imageWidth * imageHeight;
    //float top = 255;
    float sumB = 0;
    float sum1 = 0;
    float wB = 0;
    float maximum = 0;
    int level = 0;
    // make the histogram (array length 256)
    int histogram[256];

    memset(histogram, 0, sizeof(int)*256);

    for (long i = 0; i < total; i++) {
        histogram[image[i]] ++;
    }
    for (int i = 0; i < 256; i ++) {
        sum1 += i * histogram[i];
    }
    for (int i = 0; i < 256; i ++) {
        float wF = total - wB;
        //std::cout << "wF\n" << wB << "\n";
        //std::cout << "wB\n" << wF << "\n";
        if (wB > 0 && wF > 0) {
            float mF = (sum1 - sumB) / wF;
            float val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
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

// UNUSED: a simple, but well tested thresholding algorithm that works well with star images.
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

// UNUSED: Accepts an array of brightnesses and image dimensions, and returns the 
// threshold brightness for the entire image, which is defined as 5 standard deviations
// above the mean.
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
    return std::round(mean + (std * 5));
}

//                          *****START OF COG ALGORITHM*****

// This algorithm takes an image and finds all "centroids" or potential stars in it. It divides
// the image into a square grid, calculates the threshold brightness of each subsection, or "box",
// and uses that information to find the centroids. A threshold is defined as 5 standard deviations
// above the mean brightness of pixels.

// Commonly used terms:
    // box: A subsection in the subdivision scheme
    // subdivision scheme: The square grid that contains the boxes the image is divided into. Each
        // box can be assigned an index in this scheme, which follows row-major order.
    // subdivisions: The number of slices the image is sliced into on each side. The square of
        // this quantity indicates the number of boxes.
    // leftover: When doing this division, we often cannot give each box the same number of
        // rows/columns to each box. This quantity tell us how many rows/columns in the subdivision
        // scheme will have 1 extra row/column, and these always start with the first rows/columns.
    // div: The minimum rows/columns each box will have.
    // For leftover and div, "horizontal" refers to rows and "vertical" refers to columns

// By default, this algorithm ensures that regardless of the subdivisions entered, the number of
// pixels in each subdivision is imageWidth*imageHeight <= size <= 100 to ensure that enough
// pixles will be sampled in each subdivision

// Accepts an index in the subdivision scheme, and a horizontal/vertical leftover and div
// This method uses these quantities to find the first row/column index in the actual image
// that "i" references.
int StartOfSubdivision(int i, int leftover, int div) {
    if(i < leftover) {
        return i * (div + 1);
    } else {
        return leftover * (div + 1) + (i - leftover) * div;
    }
}

// Refer to definition for details
long FindSubdivision(long i, int imageWidth, int imageHeight, int subdivisions);


// Accepts a subdivsion, or "box" number, the image's brightness array and associated dimensions, 
// and the subdivisions.
// Returns the threshold for the corresponding box.
int LocalBasicThreshold(int box, unsigned char *image, int imageWidth, int imageHeight, 
        int subdivisions) {

    int horizontalDiv = imageHeight / subdivisions;
    int verticalDiv = imageWidth / subdivisions;
    int horizontalLeftover = imageHeight % subdivisions;
    int verticalLeftover = imageWidth % subdivisions;
    int row = box / subdivisions; // Finds the current row in the subdivision scheme
    int col = box % subdivisions; // Finds the current column in the subdivision scheme
    double average = 0;
    double squareSum = 0;
    long count = 0;

    // Runs through "box" in row-major order
    for(long i = StartOfSubdivision(row, horizontalLeftover, horizontalDiv); 
            i < StartOfSubdivision(row + 1, horizontalLeftover, horizontalDiv); i++) {
        for(long j = StartOfSubdivision(col, verticalLeftover, verticalDiv); 
                j < StartOfSubdivision(col + 1, verticalLeftover, verticalDiv); j++) {
            average += image[i * imageWidth + j];
            squareSum += image[i * imageWidth + j] * image[i * imageWidth + j];
            count++;

            // Checks for Error
            assert(FindSubdivision(i*imageWidth+j, imageWidth, imageHeight, subdivisions) == box);
        }
    }

    average /= count;
    return average + (5 * std::sqrt((squareSum - count * average * average) / (count - 1)));
}

// Accepts the image's brightness array and dimensions, and the subdivisions, and returns
// a vector with the threshold for each "box" in row-major order
std::vector<int> LocalThresholding(unsigned char *image, int imageWidth, int imageHeight, 
        int subdivisions) {
    
    std::vector<int> thresholds;
    for(int i = 0; i < subdivisions * subdivisions; i++) {
        thresholds.push_back(LocalBasicThreshold(i, image, imageWidth, imageHeight, 
                subdivisions));
    }

    return thresholds;
}

// Used in function CenterOfGravityAlgorithm:Go, and is useful for keeping track of
// star dimensions and characteristics, as well as the database of thresholds among
// other quantities relating to the image
struct CentroidParams {
    float yCoordMagSum; 
    float xCoordMagSum;
    long magSum;
    int xMin;
    int xMax;
    int yMin;
    int yMax;
    int magnitude;
    bool isValid;
};

// Accepts the row/column of an index on the image, the imageWidth/imageHeight, the subdivisions
// and the horizontal/vertical div and leftover
// Returns the index of the row/column in the subdivision scheme for the corresponding "i".
// NOTE: This function is the inverse of the StartOfSubdivision function in a way that
// RowOrColumn(StartOfSubdivision(x)) = x for any given row/column index "x" in the subdivision 
// scheme, but StartOfSubdivision(RowOrColumn(x)) only equals x if x = StartOfSubdivision(y),
// for any given row/column index "y" in the subdivision scheme 
long RowOrColumn(long i, int leftover, int div) {
    if(i < (div + 1) * leftover) {
        return i / (div + 1);
    } else {
        return leftover + (i - (div + 1) * leftover) / (div);
    }
}

// Accepts an index in the image, its dimensions and the subdivisions, and returns the 
// corresponding "box" number that the index "i" is in
long FindSubdivision(long i, int imageWidth, int imageHeight, int subdivisions) {
    long row = RowOrColumn(i / imageWidth, imageHeight % subdivisions, imageHeight / subdivisions);
    long col = RowOrColumn(i % imageWidth, imageWidth % subdivisions, imageWidth / subdivisions);
    return  row * subdivisions + col;
}


// Accepts an array of the image's brightnesses, and the image's dimensions, and finds all
// stars in the image and returns the stars as an array
std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const {
    
    // Program will use divisions to represent the subdivisions
    int divisions = subdivisions;
    int min = 0;
    if(divisions < 1) {
        divisions = 1;
    }
    if(imageWidth > imageHeight) {
        min = imageWidth;
    } else {
        min = imageHeight;
    }
    if(min / subdivisions < 10) {
        divisions = min / 10;
    }
    // Make an array-map of stars
    std::vector<int> localCutoff = LocalThresholding(image, imageWidth, imageHeight, divisions);
    std::unordered_map<int, int> equivalencies;
    
    // Does this ensure 0/null based values?
    unsigned char *stars = (unsigned char *) std::malloc(imageWidth * imageHeight * sizeof(unsigned char));
    stars[0] = (image[0] >= localCutoff.at(0)) ? 1 : 0;
    int L = stars[0];

    for(long i = 1; i < imageHeight * imageWidth; i++) {
        if(image[i] >= localCutoff.at(FindSubdivision(i, imageWidth, imageHeight, divisions))) {
            int up = i - imageWidth;
            int left = i - 1;
            bool leftEq = stars[left] != 0;
            bool upEq = stars[up] != 0;
            if(i / imageWidth == 0) {
                if(leftEq) {
                    stars[i] = stars[left];
                } else {
                    stars[i] = ++L;
                }
            } else if(i % imageWidth == 0) {
                if(upEq) {
                    stars[i] = stars[up];
                } else {
                    stars[i] = ++L;
                }
            } else {
                if(leftEq && upEq && stars[left] == stars[up]) {
                    stars[i] = stars[left];
                } else if(leftEq && upEq && stars[left] != stars[up]) {
                    stars[i] = std::min(stars[left], stars[up]);
                    equivalencies.insert(std::pair<int, int>(int(std::max(stars[left], stars[up])), int(stars[i])));
                    assert(equivalencies.find(stars[i]) == equivalencies.end());
                } else if(leftEq) {
                    stars[i] = stars[left];
                } else if(upEq) {
                    stars[i] = stars[up];
                } else {
                    stars[i] = ++L;
                }
            }
        }
    }
     // Get statistics of each star
    std::unordered_map<int, CentroidParams> params; // Star # to param
   for(int i = 0; i < imageWidth * imageHeight; i++) {
        if(stars[i] != 0) {
            int starNumber = equivalencies.find(stars[i]) == equivalencies.end() ? 
                    stars[i] : equivalencies.find(stars[i]) -> second;
            if(params.find(starNumber) == params.end()) {
                CentroidParams p;
                p.magnitude = 0;
                p.magSum = 0;
                p.xCoordMagSum = 0;
                p.yCoordMagSum = 0;
                p.xMax = 0;
                p.xMin = imageWidth;
                p.yMax = 0;
                p.yMin = imageHeight;
                p.isValid = true;
                params.insert(std::pair<int, CentroidParams>(starNumber, p));
            }
            CentroidParams *centroid_ptr = &params.at(starNumber);
            int row = i / imageWidth;
            int col = i % imageWidth;
            centroid_ptr -> xCoordMagSum += col * image[i];
            centroid_ptr -> yCoordMagSum += row * image[i];
            centroid_ptr -> magnitude++;
            centroid_ptr -> magSum += image[i];
            if(col > centroid_ptr -> xMax) {
                centroid_ptr -> xMax = col;
            }
            if(col < centroid_ptr -> xMin) {
                centroid_ptr -> xMin = col;
            }
            if(row > centroid_ptr -> yMax) {
                centroid_ptr -> yMax = row;
            }
            if(row < centroid_ptr -> yMin) {
                centroid_ptr -> yMin = row;
            }
            if (i % imageWidth == 0 || i % imageWidth == imageWidth - 1 || i / imageWidth == 0 || 
                i / imageWidth == imageHeight - 1) {
                centroid_ptr -> isValid = false;
            }

        }
   }

    // Now we actually make stars
    std::vector<Star> result;
    for(auto itr = params.begin(); itr != params.end(); itr++) {
        if(itr -> second.isValid && itr -> first != 0) {
            CentroidParams p = itr -> second;
            float xRadius = p.xMax - p.xMin + 1;
            float yRadius = p.yMax - p.yMin + 1;
            float xCoord = (p.xCoordMagSum / (p.magSum * 1.0));
            float yCoord = (p.yCoordMagSum / (p.magSum * 1.0));
            result.push_back(Star(xCoord + 0.5f, yCoord + 0.5f, xRadius / 2.0f, yRadius / 2.0f, p.magnitude));
        }
    }

    return result;
}

//                          *****END OF COG ALGORITHM*****

//Determines how accurate and how much iteration is done by the IWCoG algorithm,
//smaller means more accurate and more iterations.
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

void IWCoGHelper(IWCoGParams &p, long i, unsigned char *image, int imageWidth, int imageHeight, std::vector<int> &starIndices) {
    if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
        //check if pixel is on the edge of the image, if it is, we dont want to centroid this star
        if (i % imageWidth == 0 || i % imageWidth == imageWidth - 1 || i / imageWidth == 0 
                || i / imageWidth == imageHeight - 1) {
            p.isValid = false;
        }
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
    p.cutoff = BasicThreshold(image, imageWidth, imageHeight);
    for (long i = 0; i < imageHeight * imageWidth; i++) {
        //check if pixel is part of a "star" and has not been iterated over
        if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0) {
            // TODO: store longs --Mark
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
            p.isValid = true;


            IWCoGHelper(p, i, image, imageWidth, imageHeight, starIndices);

            xDiameter = (p.xMax - p.xMin) + 1;
            yDiameter = (p.yMax - p.yMin) + 1;

            //calculate fwhm
            float count = 0;
            for (int j = 0; j < (int) starIndices.size(); j++) {
                if (image[starIndices.at(j)] > p.maxIntensity / 2) {
                    count++;
                }
            }
            fwhm = sqrt(count);
            standardDeviation = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
            float modifiedStdDev = 2.0 * pow(standardDeviation, 2);
            // TODO: Why are these floats? --Mark
            float guessXCoord = (float) (p.guess % imageWidth);
            float guessYCoord = (float) (p.guess / imageWidth);
            //how much our new centroid estimate changes w each iteration
            float change = INFINITY;
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
                    float currXCoord = (float) (starIndices.at(j) % imageWidth);
                    float currYCoord = (float) (starIndices.at(j) / imageWidth);
                    w = p.maxIntensity * exp(-1.0 * ((pow(currXCoord - guessXCoord, 2) / modifiedStdDev) + (pow(currYCoord - guessYCoord, 2) / modifiedStdDev)));

                    xWeightedCoordMagSum += w * currXCoord * ((float) image[starIndices.at(j)]);
                    yWeightedCoordMagSum += w * currYCoord * ((float) image[starIndices.at(j)]);
                    weightedMagSum += w * ((float) image[starIndices.at(j)]);
                }
                float xTemp = xWeightedCoordMagSum / weightedMagSum;
                float yTemp = yWeightedCoordMagSum / weightedMagSum;

                change = abs(guessXCoord - xTemp) + abs(guessYCoord - yTemp);

                guessXCoord = xTemp;
                guessYCoord = yTemp;
            }
            if (p.isValid) {
                result.push_back(Star(guessXCoord + 0.5f, guessYCoord + 0.5f, ((float)(xDiameter))/2.0f, ((float)(yDiameter))/2.0f, starIndices.size()));
            }
        }
    }
    return result;
}

}
