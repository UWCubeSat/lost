#include "centroiders.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
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
int BasicThreshold(unsigned char *image, int imageWidth, int imageHeight,
                   std::vector<unsigned char> *rowMax = nullptr) {
    unsigned long totalMag = 0;
    decimal std = 0;
    long totalPixels = imageHeight * imageWidth;
    if (rowMax != nullptr) rowMax->assign(imageHeight, 0);
    for (int y = 0; y < imageHeight; y++) {
        const unsigned char *row = image + (long)y * imageWidth;
        unsigned char brightest = 0;
        for (int x = 0; x < imageWidth; x++) {
            totalMag += row[x];
            brightest = std::max(brightest, row[x]);
        }
        if (rowMax != nullptr) (*rowMax)[y] = brightest;
    }
    decimal mean = DECIMAL(totalMag) / totalPixels;
    for (long i = 0; i < totalPixels; i++) {
        std += DECIMAL_POW(image[i] - mean, 2);
    }
    std = DECIMAL_SQRT(std / totalPixels);
    return mean + (std * 5);
}

// Is q + k <= mean + 5 * std, where mean = q + r/n?
static bool WithinThreshold(uint64_t k, uint64_t n, uint64_t r, uint64_t scaledSquares) {
    uint64_t added = scaledSquares + 2 * k * r;
    uint64_t subtracted = k * k * n;
    if (added < subtracted) return false;
    uint64_t rest = added - subtracted;
    if (rest >= 26 * n) return true;
    return 26 * r * r <= n * rest;
}

// floor(mean + 5 * standard deviation) of all pixels.
// exact for images up to 840,000,000 pixels
// larger images use floating point.
// If rowMax is given, it is filled with the brightest pixel of each row, like BasicThreshold.
int BasicThresholdOnePass(unsigned char *image, int imageWidth, int imageHeight,
                          std::vector<unsigned char> *rowMax = nullptr) {
    const uint64_t kMaxExactPixels = 840000000;
    const int kChunk = 65536;

    if (rowMax != nullptr) rowMax->assign(std::max(imageHeight, 0), 0);
    if (imageWidth <= 0 || imageHeight <= 0) return 0;
    uint64_t n = (uint64_t)imageWidth * imageHeight;

    uint64_t sum = 0, sumSquares = 0;
    for (int y = 0; y < imageHeight; y++) {
        const unsigned char *row = image + (size_t)y * imageWidth;
        unsigned char brightest = 0;
        for (int x0 = 0; x0 < imageWidth;) {
            int x1 = x0 + std::min(kChunk, imageWidth - x0);
            uint32_t chunkSum = 0, chunkSquares = 0;
            for (int x = x0; x < x1; x++) {
                uint32_t pixel = row[x];
                chunkSum += pixel;
                chunkSquares += pixel * pixel;
                brightest = std::max(brightest, row[x]);
            }
            sum += chunkSum;
            sumSquares += chunkSquares;
            x0 = x1;
        }
        if (rowMax != nullptr) (*rowMax)[y] = brightest;
    }

    if (n > kMaxExactPixels) {
        double mean = (double)sum / n;
        return (int)(mean + 5 * std::sqrt(std::max((double)sumSquares / n - mean * mean, 0.0)));
    }

    uint64_t q = sum / n, r = sum % n;
    uint64_t scaledSquares = 25 * (sumSquares + q * q * n - 2 * q * sum);

    uint64_t k = 0;
    while (WithinThreshold(k + 1, n, r, scaledSquares)) k++;
    return (int)(q + k);
}

struct CentroidParams {
    decimal yCoordMagSum = 0;
    decimal xCoordMagSum = 0;
    long magSum = 0;
    int xMin;
    int xMax;
    int yMin;
    int yMax;

    CentroidParams(int x, int y) : xMin(x), xMax(x), yMin(y), yMax(y) {}
};

std::vector<Star> CenterOfGravityAlgorithm::Go(unsigned char* image, int imageWidth, int imageHeight) const {
    std::vector<Star> result;

    std::vector<unsigned char> rowMax;
    int cutoff = BasicThreshold(image, imageWidth, imageHeight, &rowMax);
    const int offsets[] = {1, -1, imageWidth, -imageWidth};
    const long n = imageWidth * imageHeight;
    std::vector<bool> visited(n, false);

    for (int row = 0; row < imageHeight; row++) {
        if (rowMax[row] < cutoff) continue;
        const long rowEnd = (long)(row + 1) * imageWidth;
        for (long i = (long)row * imageWidth; i < rowEnd; i++) {
            if (image[i] >= cutoff && !visited[i]) {
                visited[i] = true;

                // If any pixels are on the edge of the image then we don't centroid
                // this star, in case part of the star is missing
                bool isValid = true;
                int centroidSize = 0;
                std::vector<long> queue = {i};
                CentroidParams p(i % imageWidth, i / imageWidth);

                while (!queue.empty()) {
                    long ind = queue.back();
                    queue.pop_back();
                    centroidSize++;
                    long x = ind % imageWidth, y = ind / imageWidth;
                    if (x == 0 || x == imageWidth - 1 || y == 0 || y == imageHeight - 1) {
                        isValid = false;
                    }

                    p.magSum += image[ind];
                    p.xCoordMagSum += x * image[ind];
                    p.yCoordMagSum += y * image[ind];

                    if (x > p.xMax) p.xMax = x;
                    if (x < p.xMin) p.xMin = x;
                    if (y > p.yMax) p.yMax = y;
                    if (y < p.yMin) p.yMin = y;

                    for (int offset : offsets) {
                        long ni = ind + offset;
                        if (ni >= 0 && ni < n && std::abs((ni % imageWidth)-x) <= 1 && image[ni] >= cutoff && !visited[ni]) {
                            queue.push_back(ni);
                            visited[ni] = true;
                        }
                    }
                }

                int xDiameter = (p.xMax - p.xMin) + 1;
                int yDiameter = (p.yMax - p.yMin) + 1;

                if (isValid && p.magSum > 0) {
                    decimal xCoord = p.xCoordMagSum / DECIMAL(p.magSum);
                    decimal yCoord = p.yCoordMagSum / DECIMAL(p.magSum);
                    result.push_back(Star(xCoord + DECIMAL(0.5), yCoord + DECIMAL(0.5), (xDiameter)/DECIMAL(2.0), (yDiameter)/DECIMAL(2.0), centroidSize));
                }
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

// WINDOWED CENTER OF GRAVITY

struct WindowedMatch {
    int prediction;
    int xMin, xMax, yMin, yMax;  ///< bounding box of the star's pixels
    Star star;
    Vec2 offset;  ///< star position minus prediction
};

/// Threshold, in noise standard deviations above the background.
static const int kThresholdSigmas = 5;
/// Added to the MAD, in gray levels, before converting it to a standard deviation. Gray levels
/// are whole numbers and the border is a small sample, so the MAD often comes out too small.
/// Also keeps the threshold at least 3.7 gray levels above a noiseless background.
static const decimal kMadAllowance = DECIMAL(0.5);

/**
 * Twice the median of `count` samples in a histogram: the sum of the two middle values (the
 * middle value twice if count is odd).
 */
static int TwiceMedian(const int *histogram, int bins, int count) {
    int lower = -1;
    int seen = 0;
    for (int value = 0; value < bins; value++) {
        seen += histogram[value];
        if (lower < 0 && seen > (count - 1) / 2) lower = value;
        if (seen > count / 2) return lower + value;
    }
    return 0;  // not reached if the histogram holds `count` samples
}

/**
 * Median and median absolute deviation (MAD) of a histogram of `count` gray levels. An even
 * count averages the two middle values. Returned as twice the median and four times the MAD,
 * which are whole numbers.
 */
static void HistogramMedianAndMad(const int histogram[256], int count, int *twiceMedian, int *fourTimesMad) {
    *twiceMedian = TwiceMedian(histogram, 256, count);
    // histogram of |2 * level - twiceMedian|, twice each pixel's distance from the median
    int deviations[511] = {};
    for (int level = 0; level < 256; level++) {
        deviations[std::abs(2 * level - *twiceMedian)] += histogram[level];
    }
    *fourTimesMad = TwiceMedian(deviations, 511, count);
}

Stars WindowedCenterOfGravity::Go(const unsigned char *image, int imageWidth, int imageHeight,
                                  const std::vector<Vec2> &predictions,
                                  std::vector<int> *matchedPredictions) const {
    if (matchedPredictions != nullptr) {
        matchedPredictions->clear();
    }
    // unusable settings; the last check refuses a NaN outlierRadius too
    if (image == nullptr || gateRadius < 0 || gateRadius > kMaxRadius || maxStarRadius < 0 || maxStarRadius > kMaxRadius
        || !(outlierRadius >= 0)) {
        return Stars();
    }
    // pixel indices and window coordinates are ints; pixel centers must be exact in a float
    if (imageWidth <= 0 || imageHeight <= 0 || imageWidth > kMaxImageSide || imageHeight > kMaxImageSide
        || imageWidth > std::numeric_limits<int>::max() / imageHeight) {
        return Stars();
    }

    const int halfWindow = gateRadius + maxStarRadius;
    const int windowSide = 2 * halfWindow + 1;
    // flood fill marks, one per window pixel; a window is never larger than the image
    const int visitedWidth = std::min(windowSide, imageWidth);
    std::vector<unsigned char> visited(visitedWidth * std::min(windowSide, imageHeight));
    std::vector<int> stack;
    std::vector<WindowedMatch> matches;

    for (int p = 0; p < (int)predictions.size(); p++) {
        const Vec2 &prediction = predictions[p];
        // skip predictions off the image; non-finite ones are skipped too
        if (!(prediction.x >= 0 && prediction.x <= imageWidth && prediction.y >= 0 && prediction.y <= imageHeight)) {
            continue;
        }
        // pixel (x, y) covers [x, x+1) * [y, y+1): this is the pixel containing the prediction
        const int windowX0 = (int)DECIMAL_FLOOR(prediction.x) - halfWindow;
        const int windowY0 = (int)DECIMAL_FLOOR(prediction.y) - halfWindow;
        const int windowX1 = windowX0 + windowSide - 1;
        const int windowY1 = windowY0 + windowSide - 1;
        const int xMin = std::max(0, windowX0), xMax = std::min(imageWidth - 1, windowX1);
        const int yMin = std::max(0, windowY0), yMax = std::min(imageHeight - 1, windowY1);
        if (xMin >= xMax || yMin >= yMax) {
            continue;
        }

        // background and noise from the window border
        // (a star crossing it barely moves the median and MAD)
        int histogram[256] = {};
        int count = 0;
        for (int x = xMin; x <= xMax; x++) {
            histogram[image[yMin * imageWidth + x]]++;
            histogram[image[yMax * imageWidth + x]]++;
            count += 2;
        }
        for (int y = yMin + 1; y < yMax; y++) {
            histogram[image[y * imageWidth + xMin]]++;
            histogram[image[y * imageWidth + xMax]]++;
            count += 2;
        }
        int twiceBackground, fourTimesMad;
        HistogramMedianAndMad(histogram, count, &twiceBackground, &fourTimesMad);
        // 1.4826 * MAD estimates the standard deviation of Gaussian noise
        const decimal sigma = DECIMAL(1.4826) * (fourTimesMad / DECIMAL(4.0) + kMadAllowance);
        // star pixels are at least this bright
        const int threshold = (int)DECIMAL_CEIL(twiceBackground / DECIMAL(2.0) + kThresholdSigmas * sigma);

        // seed: brightest pixel whose center is within the gate
        int seed = -1;
        int seedValue = -1;
        for (int y = std::max(yMin, windowY0 + maxStarRadius); y <= std::min(yMax, windowY1 - maxStarRadius); y++) {
            for (int x = std::max(xMin, windowX0 + maxStarRadius); x <= std::min(xMax, windowX1 - maxStarRadius); x++) {
                decimal dx = x + DECIMAL(0.5) - prediction.x;
                decimal dy = y + DECIMAL(0.5) - prediction.y;
                if (dx * dx + dy * dy <= gateRadius * gateRadius && image[y * imageWidth + x] > seedValue) {
                    seedValue = image[y * imageWidth + x];
                    seed = y * imageWidth + x;
                }
            }
        }
        if (seedValue < threshold) {
            continue;
        }

        // flood fill inside the window
        std::fill(visited.begin(), visited.end(), 0);
        stack.assign(1, seed);
        visited[(seed / imageWidth - yMin) * visitedWidth + (seed % imageWidth - xMin)] = 1;
        int64_t weightSum = 0, xWeightSum = 0, yWeightSum = 0;
        int pixels = 0;
        int starXMin = seed % imageWidth, starXMax = starXMin;
        int starYMin = seed / imageWidth, starYMax = starYMin;
        bool touchesEdge = false;
        while (!stack.empty()) {
            const int index = stack.back();
            stack.pop_back();
            const int x = index % imageWidth;
            const int y = index / imageWidth;
            if (x == xMin || x == xMax || y == yMin || y == yMax) {
                touchesEdge = true;  // possibly clipped by the window or the image
            }
            // twice the level above the background; the factor cancels in the centroid
            const int weight = 2 * image[index] - twiceBackground;
            weightSum += weight;
            xWeightSum += (int64_t)weight * x;
            yWeightSum += (int64_t)weight * y;
            pixels++;
            starXMin = std::min(starXMin, x);
            starXMax = std::max(starXMax, x);
            starYMin = std::min(starYMin, y);
            starYMax = std::max(starYMax, y);

            const int neighborX[4] = {x + 1, x - 1, x, x};
            const int neighborY[4] = {y, y, y + 1, y - 1};
            for (int n = 0; n < 4; n++) {
                const int nx = neighborX[n], ny = neighborY[n];
                if (nx < xMin || nx > xMax || ny < yMin || ny > yMax) {
                    continue;
                }
                unsigned char &seen = visited[(ny - yMin) * visitedWidth + (nx - xMin)];
                if (!seen && image[ny * imageWidth + nx] >= threshold) {
                    seen = 1;
                    stack.push_back(ny * imageWidth + nx);
                }
            }
        }
        if (touchesEdge) {
            continue;
        }

        const Vec2 position = {
            (decimal)xWeightSum / weightSum + DECIMAL(0.5),
            (decimal)yWeightSum / weightSum + DECIMAL(0.5),
        };
        const Vec2 offset = position - prediction;
        if (offset.MagnitudeSq() > gateRadius * gateRadius) {
            continue;
        }
        matches.push_back({p, starXMin, starXMax, starYMin, starYMax,
                           Star(position.x, position.y,
                                (starXMax - starXMin + 1) / DECIMAL(2.0), (starYMax - starYMin + 1) / DECIMAL(2.0),
                                pixels),
                           offset});
    }

    // A star reached from more than one prediction is ambiguous: drop all its matches.
    // Thresholds differ per prediction, so the same star can come out as different pixels. But a
    // star is everything connected to its seed above its threshold, so two that share or touch
    // at a pixel always have overlapping bounding boxes.
    std::vector<bool> ambiguous(matches.size(), false);
    for (size_t a = 0; a < matches.size(); a++) {
        for (size_t b = a + 1; b < matches.size(); b++) {
            if (matches[a].xMin <= matches[b].xMax && matches[b].xMin <= matches[a].xMax
                && matches[a].yMin <= matches[b].yMax && matches[b].yMin <= matches[a].yMax) {
                ambiguous[a] = ambiguous[b] = true;
            }
        }
    }
    std::vector<WindowedMatch> unique;
    for (size_t i = 0; i < matches.size(); i++) {
        if (!ambiguous[i]) unique.push_back(matches[i]);
    }

    // drop matches that disagree with the frame's common shift
    if (unique.size() >= 3) {
        std::vector<decimal> xOffsets, yOffsets;
        for (const WindowedMatch &match : unique) {
            xOffsets.push_back(match.offset.x);
            yOffsets.push_back(match.offset.y);
        }
        std::nth_element(xOffsets.begin(), xOffsets.begin() + xOffsets.size() / 2, xOffsets.end());
        std::nth_element(yOffsets.begin(), yOffsets.begin() + yOffsets.size() / 2, yOffsets.end());
        const Vec2 medianOffset = {xOffsets[xOffsets.size() / 2], yOffsets[yOffsets.size() / 2]};
        unique.erase(std::remove_if(unique.begin(), unique.end(), [&](const WindowedMatch &match) {
            return (match.offset - medianOffset).MagnitudeSq() > outlierRadius * outlierRadius;
        }), unique.end());
    }

    std::sort(unique.begin(), unique.end(),
              [](const WindowedMatch &a, const WindowedMatch &b) { return a.prediction < b.prediction; });
    Stars result;
    for (const WindowedMatch &match : unique) {
        result.push_back(match.star);
        if (matchedPredictions != nullptr) {
            matchedPredictions->push_back(match.prediction);
        }
    }
    return result;
}

}
