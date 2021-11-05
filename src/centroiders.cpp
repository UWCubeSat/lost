#include "centroiders.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <iostream>

#include <unordered_map>

#include <unordered_set>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>

#include <eigen3/unsupported/Eigen/NonLinearOptimization>
#include <eigen3/unsupported/Eigen/NumericalDiff>

namespace lost
{

    // DUMMY

    std::vector<Star> DummyCentroidAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const
    {
        std::vector<Star> result;

        for (int i = 0; i < numStars; i++)
        {
            result.push_back(Star(rand() % imageWidth, rand() % imageHeight, 10.0));
        }

        return result;
    }

    struct CentroidParams
    {
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
    void CogHelper(CentroidParams &p, int i, unsigned char *image, int imageWidth, int imageHeight)
    {
        if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p.cutoff && p.checkedIndices.count(i) == 0)
        {
            p.checkedIndices.insert(i);
            if (i % imageWidth > p.xMax)
            {
                p.xMax = i % imageWidth;
            }
            else if (i % imageWidth < p.xMin)
            {
                p.xMin = i % imageWidth;
            }
            if (i / imageWidth > p.yMax)
            {
                p.yMax = i / imageWidth;
            }
            else if (i / imageWidth < p.yMin)
            {
                p.yMin = i / imageWidth;
            }
            p.magSum += image[i];
            p.xCoordMagSum += ((i % imageWidth)) * image[i];
            p.yCoordMagSum += ((i / imageWidth)) * image[i];
            if (i % imageWidth != imageWidth - 1)
            {
                CogHelper(p, i + 1, image, imageWidth, imageHeight);
            }
            if (i % imageWidth != 0)
            {
                CogHelper(p, i - 1, image, imageWidth, imageHeight);
            }
            CogHelper(p, i + imageWidth, image, imageWidth, imageHeight);
            CogHelper(p, i - imageWidth, image, imageWidth, imageHeight);
        }
    }

    int DetermineCutoff(unsigned char *image, int imageWidth, int imageHeight)
    {
        //loop through entire array, find sum of magnitudes
        int totalMag = 0;
        for (int i = 0; i < imageHeight * imageWidth; i++)
        {
            totalMag += image[i];
        }
        return (((totalMag / (imageHeight * imageWidth)) + 1) * 15) / 10;
    }

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

    std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const
    {
        CentroidParams p;

        std::vector<Star> result;

        p.cutoff = DetermineCutoff(image, imageWidth, imageHeight);
        for (int i = 0; i < imageHeight * imageWidth; i++)
        {
            //check if pixel is part of a "star" and has not been iterated over
            if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0)
            {
                //iterate over pixels that are part of the star
                int xDiameter = 0; //radius of current star
                int yDiameter = 0;
                p.yCoordMagSum = 0; //y coordinate of current star
                p.xCoordMagSum = 0; //x coordinate of current star
                p.magSum = 0;       //sum of magnitudes of current star

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
                result.push_back(Star(xCoord, yCoord, ((float)(xDiameter * 1.0)) / 2.0, ((float)(yDiameter * 1.0)) / 2.0, 0));
            }
        }
        return result;
    }

    //Determines how accurate and how much iteration is done by the IWCoG algorithm,
    //smaller means more accurate and more iterations.
    float iWCoGMinChange = 0.0002;

    struct IWCoGParams
    {
        int xMin;
        int xMax;
        int yMin;
        int yMax;
        int cutoff;
        int maxIntensity;
        int guess;
        std::unordered_set<int> checkedIndices;
    };

    void IWCoGHelper(IWCoGParams &p, int i, unsigned char *image, int imageWidth, int imageHeight, std::vector<int> &starIndices)
    {
        if (i >= 0 && i < imageWidth * imageHeight && image[i] >= p.cutoff && p.checkedIndices.count(i) == 0)
        {
            p.checkedIndices.insert(i);
            starIndices.push_back(i);
            if (image[i] > p.maxIntensity)
            {
                p.maxIntensity = image[i];
                p.guess = i;
            }
            if (i % imageWidth > p.xMax)
            {
                p.xMax = i % imageWidth;
            }
            else if (i % imageWidth < p.xMin)
            {
                p.xMin = i % imageWidth;
            }
            if (i / imageWidth > p.yMax)
            {
                p.yMax = i / imageWidth;
            }
            else if (i / imageWidth < p.yMin)
            {
                p.yMin = i / imageWidth;
            }
            if (i % imageWidth != imageWidth - 1)
            {
                IWCoGHelper(p, i + 1, image, imageWidth, imageHeight, starIndices);
            }
            if (i % imageWidth != 0)
            {
                IWCoGHelper(p, i - 1, image, imageWidth, imageHeight, starIndices);
            }
            IWCoGHelper(p, i + imageWidth, image, imageWidth, imageHeight, starIndices);
            IWCoGHelper(p, i - imageWidth, image, imageWidth, imageHeight, starIndices);
        }
    }

    Stars IterativeWeightedCenterOfGravityAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const
    {
        IWCoGParams p;
        std::vector<Star> result;
        p.cutoff = DetermineCutoff(image, imageWidth, imageHeight);
        for (int i = 0; i < imageHeight * imageWidth; i++)
        {
            //check if pixel is part of a "star" and has not been iterated over
            if (image[i] >= p.cutoff && p.checkedIndices.count(i) == 0)
            {
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
                float count = 0;
                for (int j = 0; j < (int)starIndices.size(); j++)
                {
                    if (image[starIndices.at(j)] > p.maxIntensity / 2)
                    {
                        count++;
                    }
                }
                fwhm = sqrt(count);
                standardDeviation = fwhm / (2.0 * sqrt(2.0 * log(2.0)));
                float modifiedStdDev = 2.0 * pow(standardDeviation, 2);
                float guessXCoord = (float)(p.guess % imageWidth);
                float guessYCoord = (float)(p.guess / imageWidth);
                //how much our new centroid estimate changes w each iteration
                float change = INFINITY;
                int stop = 0;
                //while we see some large enough change in estimated, maybe make it a global variable
                while (change > iWCoGMinChange && stop < 100000)
                {
                    //traverse through star indices, calculate W at each coordinate, add to final coordinate sums
                    yWeightedCoordMagSum = 0;
                    xWeightedCoordMagSum = 0;
                    weightedMagSum = 0;
                    stop++;
                    for (int j = 0; j < (int)starIndices.size(); j++)
                    {
                        //calculate w
                        float currXCoord = (float)(starIndices.at(j) % imageWidth);
                        float currYCoord = (float)(starIndices.at(j) / imageWidth);
                        w = p.maxIntensity * exp(-1.0 * ((pow(currXCoord - guessXCoord, 2) / modifiedStdDev) + (pow(currYCoord - guessYCoord, 2) / modifiedStdDev)));

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
                result.push_back(Star(guessXCoord, guessYCoord, ((float)(xDiameter * 1.0)) / 2.0, ((float)(yDiameter * 1.0)) / 2.0, 0));
            }
        }
        return result;
    }

    void Gauss1DHelper(int i, int cutoff, unsigned char *image, int imageWidth, int imageHeight, std::vector<int> &starIndices, std::unordered_set<int> &checkedIndices)
    {
        if (i >= 0 && i < imageWidth * imageHeight && image[i] >= cutoff && checkedIndices.count(i) == 0)
        {
            checkedIndices.insert(i);
            starIndices.push_back(i);
            if (i % imageWidth != imageWidth - 1)
            {
                Gauss1DHelper(i + 1, cutoff, image, imageWidth, imageHeight, starIndices, checkedIndices);
            }
            if (i % imageWidth != 0)
            {
                Gauss1DHelper(i - 1, cutoff, image, imageWidth, imageHeight, starIndices, checkedIndices);
            }
            Gauss1DHelper(i + imageWidth, cutoff, image, imageWidth, imageHeight, starIndices, checkedIndices);
            Gauss1DHelper(i - imageWidth, cutoff, image, imageWidth, imageHeight, starIndices, checkedIndices);
        }
    }

    struct LMFunctor
    {
        Eigen::MatrixXf measuredValues;
        // computes m errors, one per data point, for the given params in x
        int operator()(const Eigen::VectorXf &x, Eigen::VectorXf &fvec) const
        {
            float aParam = x(0);
            float y_bParam = x(1);
            float sdParam = x(2);

            for (int j = 0; j < values(); j++)
            {
                float yCoord = measuredValues(j, 0);
                float magSum = measuredValues(j, 1);
                fvec(j) = magSum - (aParam * (exp((-1.0 * pow(yCoord - y_bParam, 2.0)) / (2.0 * pow(sdParam, 2.0)))));
            }
            return 0;
        }
        // Compute jacobian of teh errors
        int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const
        {
            float epsilon;
            epsilon = 1e-5f;

            for (int i = 0; i < x.size(); i++)
            {
                Eigen::VectorXf xPlus(x);
                xPlus(i) += epsilon;
                Eigen::VectorXf xMinus(x);
                xMinus(i) -= epsilon;

                Eigen::VectorXf fvecPlus(values());
                operator()(xPlus, fvecPlus);

                Eigen::VectorXf fvecMinus(values());
                operator()(xMinus, fvecMinus);

                Eigen::VectorXf fvecDiff(values());
                fvecDiff = (fvecPlus - fvecMinus) / (2.0f * epsilon);

                fjac.block(0, i, values(), 1) = fvecDiff;
            }
            return 0;
        }
        // num data points
        int m;

        // gets m
        int values() const { return m; }

        // num parameters
        int n;

        int inputs() const { return n; }
    };

    Stars GaussianFit1DAlgorithm::Go(unsigned char *image, int imageWidth, int imageHeight) const
    {
        std::vector<Star> result;
        std::unordered_set<int> checkedIndices;
        int cutoff = BasicThreshold(image, imageWidth, imageHeight);
        for (int i = 0; i < imageHeight * imageWidth; i++)
        {
            //check if pixel is part of a "star" and has not been iterated over
            if (image[i] >= cutoff && checkedIndices.count(i) == 0)
            {
                std::vector<int> starIndices; //vector of star coordinates
                Gauss1DHelper(i, cutoff, image, imageWidth, imageHeight, starIndices, checkedIndices);

                int maxXInd = 0;
                int minXInd = 10000;
                int maxYInd = 0;
                int minYInd = 10000;
                for (int j = 0; j < (int)starIndices.size(); j++)
                {
                    int temp = starIndices.at(j);
                    if (temp % imageWidth > maxXInd)
                    {
                        maxXInd = temp % imageWidth;
                    }
                    if (temp % imageWidth < minXInd)
                    {
                        minXInd = temp % imageWidth;
                    }
                    if (temp / imageWidth > maxYInd)
                    {
                        maxYInd = temp / imageWidth;
                    }
                    if (temp / imageWidth < minYInd)
                    {
                        minYInd = temp / imageWidth;
                    }
                }

                int xRange = 1 + maxXInd - minXInd;
                int yRange = 1 + maxYInd - minYInd;
                float xMid = (maxXInd + minXInd) / 2.0;
                float yMid = (maxYInd + minYInd) / 2.0;


                Eigen::MatrixXf xMeasuredValues(xRange, 2);
                Eigen::MatrixXf yMeasuredValues(yRange, 2);


                float fwhm = 0.0;

                int start = (minYInd * imageWidth) + minXInd;
                int end = (maxYInd * imageWidth) + maxXInd;

                for (int j = start; j <= end; j++) {
                    if (image[j] >= 128) {
                        fwhm++;
                    }
                    xMeasuredValues((j % imageWidth) - minXInd, 0) = (float) (j % imageWidth);
                    xMeasuredValues((j % imageWidth) - minXInd, 1) += (float) image[j];

                    yMeasuredValues((j / imageWidth) - minYInd, 0) = (float) (j / imageWidth);
                    yMeasuredValues((j / imageWidth) - minYInd, 1) += (float) image[j];

                    if (j % imageWidth >= maxXInd) {
                        j = j + imageWidth - xRange + 1;
                    }
                }

                std::cout << "xMeasuredValues:\n" << xMeasuredValues << "\n";
                std::cout << "yMeasuredValues:\n" << yMeasuredValues << "\n";


                float initStd = fwhm / (2.0 * sqrt(2.0 * log(2.0)));

                // eigen stuff here:

                Eigen::VectorXf x(3);

                x(0) = 255.0;
                x(1) = xMid;
                x(2) = initStd;

                int n = 3; // num parameters
                int xm = xRange;

                LMFunctor xFunctor;
                xFunctor.measuredValues = xMeasuredValues;
                xFunctor.m = xm;
                xFunctor.n = n;

                Eigen::LevenbergMarquardt<LMFunctor, float> xlm(xFunctor);
                xlm.minimize(x);
                float xCoord = x(1);

                int ym = yRange; // row data points

                x(0) = 255.0;
                x(1) = yMid;
                x(2) = initStd;

                LMFunctor yFunctor;
                yFunctor.measuredValues = yMeasuredValues;
                yFunctor.m = ym;
                yFunctor.n = n;

                Eigen::LevenbergMarquardt<LMFunctor, float> ylm(yFunctor);
                ylm.minimize(x);

                float yCoord = x(1);

                result.push_back(Star(xCoord, yCoord, (float)xRange / 2.0, (float)yRange / 2, 0));
            }
        }
        return result;
    }

}