#include "star-id.hpp"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <algorithm>
#include <random>
#include <set>
#include <vector>

#include "attitude-utils.hpp"
#include "databases.hpp"

namespace lost {

TetraStarIdAlgorithm::TetraStarIdAlgorithm() {}

/**
 * Calculates exponent base for logarithmic binning
 * TODO: what is log binning? what do the params mean
 * @param errorSlope
 * @param errorOffset
 * @return
 */
float TetraStarIdAlgorithm::getLogBinningBase(float errorSlope,
                                              float errorOffset) const {
    if (errorOffset <= 0) {
        cerr << "Error: errorOffset must be positive" << endl;
        exit(EXIT_FAILURE);
    }
    // Calculate base of logarithmic binning function
    // fmax(a, b) returns the larger of 2 floating point arguments
    float base = (1 + errorSlope) / fmax(1 - errorSlope, 0);
    return base;
}

/**
 * TODO
 * @param input
 * @param errorSlope
 * @param errorOffset
 * @return
 */
int TetraStarIdAlgorithm::logBin(float input, float errorSlope,
                                 float errorOffset) const {
    if (isinf(errorSlope) || isinf(errorOffset)) {
        return 0;
    }
    // get base of log binning function
    float base = getLogBinningBase(errorSlope, errorOffset);
    int bin;
    // TODO: linear binning?
    if (base <= 1 + errorOffset * tetra::binSizeRatio / 10.0) {
        bin = input / (2 * (errorSlope + errorOffset) * tetra::binSizeRatio);
    } else {
        bin = (log(input * errorSlope / errorOffset + 1) / log(base)) /
              tetra::binSizeRatio;
    }
    return bin;
}

/**
 * Calculates min possible input value given its log bin TODO: how does it work
 * @param bin
 * @param errorSlope
 * @param errorOffset
 * @return
 */
float TetraStarIdAlgorithm::logUnbin(int bin, float errorSlope,
                                     float errorOffset) const {
    if (isinf(errorSlope) || isinf(errorOffset)) {
        return 0;
    }
    float base = getLogBinningBase(errorSlope, errorOffset);
    float minInput;
    if (base <= 1 + errorOffset * tetra::binSizeRatio / 10.0) {
        minInput = bin * 2 * (errorSlope + errorOffset) * tetra::binSizeRatio;
    } else {
        // pow(base, pwr) returns base^pwr
        minInput = (pow(base, bin * tetra::binSizeRatio) - 1) * errorOffset /
                   errorSlope;
    }
    return minInput;
}

/**
 * Bin largest edge length TODO: how?
 * @param largestEdgeLen
 * @param errorRatio
 * @return
 */
int TetraStarIdAlgorithm::binLargestEdge(int largestEdgeLen,
                                         int errorRatio) const {
    float leRatio = largestEdgeLen / ((1 << 16) - 1.0);
    leRatio +=
        errorRatio * (leRatio * tetra::leErrorSlope + tetra::leErrorOffset);
    return logBin(leRatio, tetra::leErrorSlope, tetra::leErrorOffset);
}

/**
 * Return min possible LE ratio within the bin
 * TODO: what is this
 * @param leBin
 * @return
 */
float TetraStarIdAlgorithm::unbinLargestEdge(int leBin) const {
    return logUnbin(leBin, tetra::leErrorSlope, tetra::leErrorOffset);
}

/**
 * Returns y coordinate bin TODO
 * @param y
 * @param leBin
 * @param errorRatio
 * @return
 */
int TetraStarIdAlgorithm::binYCoord(int y, int leBin, int errorRatio) const {
    float minLERatio = unbinLargestEdge(leBin);
    float errorConst = tetra::leErrorOffset / (2 - tetra::maxScaleFactor);
    float errorSlope = errorConst / fmax(minLERatio - errorConst, 0);
    float errorOffset = errorSlope;

    float yRatio = y / ((1 << 14) - 1.0);
    // copysign(a, b) returns a floating point value with magnitude of a & sign
    // of b fabs(x) returns absolute value of a floating point x
    yRatio +=
        errorRatio * copysign(fabs(yRatio) * errorSlope + errorOffset, yRatio);
    int bin = logBin(fabs(yRatio), errorSlope, errorOffset);
    if (yRatio < 0) {
        // ~ is bitwise NOT, inverts all bits in binary number
        bin = ~bin;
    }
    return bin;
}

/**
 * TODO
 * @param yBin
 * @param leBin
 * @return
 */
float TetraStarIdAlgorithm::unbinYCoord(int yBin, int leBin) const {
    float minLERatio = unbinLargestEdge(leBin);
    float errorConst = tetra::leErrorOffset / (2 - tetra::maxScaleFactor);
    float errorSlope = errorConst / fmax(minLERatio - errorConst, 0);
    float errorOffset = errorSlope;

    float maxYRatio =
        logUnbin(yBin >= 0 ? yBin + 1 : (~yBin) + 1, errorSlope, errorOffset);
    return maxYRatio;
}

/**
 * Bin x coordinate using y coordinate and LE bins
 * TODO: math?
 * @param x
 * @param leBin
 * @param yBin
 * @param errorRatio
 * @return
 */
int TetraStarIdAlgorithm::binXCoord(int x, int leBin, int yBin,
                                    int errorRatio) const {
    float minLERatio = unbinLargestEdge(leBin);
    float maxYRatio = unbinYCoord(yBin, leBin);
    float errorConst = tetra::leErrorOffset / (2 - tetra::maxScaleFactor);

    float errorSlope = errorConst / fmax(minLERatio - errorConst, 0);
    float errorOffset =
        errorSlope * (1 + 2 * sqrt((1.0 / 4) + maxYRatio * maxYRatio)) / 2;
    float xRatio = x / ((1 << 14) - 1.0);
    xRatio +=
        errorRatio * copysign(fabs(xRatio) * errorSlope + errorOffset, xRatio);

    int bin = logBin(fabs(xRatio), errorSlope, errorOffset);
    if (xRatio < 0) {
        bin = ~bin;
    }
    return bin;
}

/**
 * Hash function
 * TODO: uint64_t = unsigned long long- replace with some other hash function?
 * @param oldHash
 * @param key
 * @return
 */
uint64_t TetraStarIdAlgorithm::hashInt(uint64_t oldHash, uint64_t key) const {
    key = key * 11400714819323198549ULL;
    return oldHash ^ (oldHash >> 13) ^ key;
}

/**
 * Hash function, takes in a Pattern and produces a corresponding catalog
 * position Hash is based on Pattern's bins NOTE: this function just computes an
 * index, doesn't put anything in the actual catalog
 * @param patt
 * @return
 */
// TODO: pass Pattern by reference?
uint64_t TetraStarIdAlgorithm::hashPattern(tetra::Pattern patt) const {
    // initialize hash value = largest edge bin
    int leBin = binLargestEdge(patt.leLength, 0);
    uint64_t hash = hashInt(0, leBin);

    // Update hash using each Feature's x and y bins
    for (int i = 0; i < tetra::numPattStars - 2; i++) {
        int yBin = binYCoord(patt.features[i].y, leBin, 0);
        hash = hashInt(hash, yBin + (1 << 31));
        int xBin = binXCoord(patt.features[i].x, leBin, yBin, 0);
        hash = hashInt(hash, xBin + (1 << 31));
    }
    // Could result in a collision
    return hash % tetra::numCatalogPatts;
}

/**
 * Checks that 2 Patterns (newPattern and catPattern) have the same bin pairs
 * (LE, x, y) in case of a collision Used in isMatch() as a kind of pre-check
 * @param newPattern Pattern newly created from image
 * @param catPattern Pattern stored in catlog
 * @return true if all bin pairs are the same; false otherwise
 */
bool TetraStarIdAlgorithm::checkSameBins(tetra::Pattern newPattern,
                                         tetra::Pattern catPattern) const {
    // check that both Patterns have same largest edge bin
    int newLEBin = binLargestEdge(newPattern.leLength, 0);
    int catLEBin =
        binLargestEdge(catPattern.leLength, 2 * catPattern.leBinOffset - 1);
    // cout << "LEBINs: " << newLEBin << " " << catLEBin << endl;
    if (newLEBin != catLEBin) {
        return false;
    }

    // check each Feature, confirm they have same x and y bins
    for (int i = 0; i < tetra::numPattStars - 2; i++) {
        tetra::Feature newFeature = newPattern.features[i];
        tetra::Feature catFeature = catPattern.features[i];
        int newYBin = binYCoord(newFeature.y, newLEBin, 0);
        int catYBin =
            binYCoord(catFeature.y, catLEBin, 2 * catFeature.yBinOffset - 1);

        int newXBin = binXCoord(newFeature.x, newLEBin, newYBin, 0);
        int catXBin = binXCoord(catFeature.x, catLEBin, catYBin,
                                2 * catFeature.xBinOffset - 1);

        // cout << "BINS: " << newYBin << " " << catYBin << " " << newXBin << "
        // "
        //      << catXBin << endl;
        if ((newYBin != catYBin) || (newXBin != catXBin)) {
            // cout << "BINS NOT MATCHED" << endl;
            return false;
        }
    }
    return true;
}

/**
 * TODO: define "match"- coordinates nearly the same?
 * TODO: math
 * @param newPattern
 * @param catPattern
 * @return
 */
bool TetraStarIdAlgorithm::isMatch(tetra::Pattern newPattern,
                                   tetra::Pattern catPattern) const {
    // cout << "Matching..." << endl;
    float newLERatio = newPattern.leLength / ((1 << 16) - 1.0);
    float catLERatio = catPattern.leLength / ((1 << 16) - 1.0);
    float maxLEError = catLERatio * tetra::leErrorSlope + tetra::leErrorOffset;
    if (fabs(newLERatio - catLERatio) > maxLEError) {
        return false;
    }

    float coordErrorConst = tetra::leErrorOffset / (2 - tetra::maxScaleFactor);
    float coordErrorSlope =
        coordErrorConst / fmax(newLERatio - coordErrorConst, 0);
    float coordErrorOffsetY = coordErrorSlope;

    for (int i = 0; i < tetra::numPattStars - 2; i++) {
        float newYCoord = newPattern.features[i].y / ((1 << 14) - 1.0);
        float catYCoord = catPattern.features[i].y / ((1 << 14) - 1.0);
        float maxYError = fabs(catYCoord) * coordErrorSlope + coordErrorOffsetY;
        if (fabs(newYCoord - catYCoord) > maxYError) {
            return false;
        }
    }

    int catLEBin =
        binLargestEdge(catPattern.leLength, 2 * catPattern.leBinOffset - 1);
    for (int i = 0; i < tetra::numPattStars - 2; i++) {
        int catYBin = binYCoord(catPattern.features[i].y, catLEBin,
                                2 * catPattern.features[i].yBinOffset - 1);
        float maxYRatio = unbinYCoord(catYBin, catLEBin);
        float coordErrorOffsetX =
            coordErrorSlope *
            (1 + 2 * sqrt((1.0 / 4) + maxYRatio * maxYRatio)) / 2;

        float newXCoord = newPattern.features[i].x / ((1 << 14) - 1.0);
        float catXCoord = catPattern.features[i].x / ((1 << 14) - 1.0);
        float maxXError = fabs(catXCoord) * coordErrorSlope + coordErrorOffsetX;
        if (fabs(newXCoord - catXCoord) > maxXError) {
            return false;
        }
    }
    // cout << "isMatch() FAILED" << endl;
    return true;
}

/**
 * Update where we're looking in the Pattern cache, possibly update
 * contents of cache itself
 * @param pattCatalog
 * @param catCache
 * @param offset Where the cache starts in the CATALOG
 * @param cacheOffset Where we are in the CACHE
 * @param probeStep How far to step in the cache to find next possible match
 * @return
 */
bool TetraStarIdAlgorithm::incrementOffset(
    FILE *pattCatalog, tetra::Pattern catCache[tetra::pattCacheSize],
    uint64_t *offset, int *cacheOffset, int *probeStep) const {
    if (((*probeStep) * (*probeStep + 1)) / 2 > tetra::maxProbeDepth) {
        return false;  // probe went outside probe bounds
    }
    // Update cache offset
    *cacheOffset += *probeStep;
    // If we probe outside our cache, update cache to next probe offset
    if (*cacheOffset >= tetra::pattCacheSize) {
        // Update offset (i.e., where we are in the CATALOG)
        *offset += *cacheOffset;
        // Reset our cache offset to 0
        *cacheOffset = 0;
        // Moves file pointer {offset} bytes from {origin}
        // SEEK_SET = 0 = beginning of file
        fseek(pattCatalog, *offset * sizeof(tetra::Pattern), SEEK_SET);
        // Reads an array of size = {pattCacheSize = 16} Patterns from catalog,
        // stores them in our cache
        fread(catCache, sizeof(tetra::Pattern), tetra::pattCacheSize,
              pattCatalog);
    }
    *probeStep += 1;
    return true;  // probe stayed within probe bounds
}

/**
 * Looks through catalog to see if there is a catalog Pattern
 * matching the Pattern we constructed from our image
 * Note: we use a Pattern cache which stores a section of the catalog to look
 * for matches
 * @param imgPattern Our constructed Pattern that we want to find a match for
 * @param catalogPattern Output, pattern in catalog that matches ours
 * @param pattCatalog Catalog storing all precomputed Patterns
 * @return true if a UNIQUE match is found; else return false if multiple or no
 * matches found
 */
bool TetraStarIdAlgorithm::getMatchingPattern(tetra::Pattern imgPattern,
                                              tetra::Pattern *catalogPattern,
                                              FILE *pattCatalog) const {
    // TODO: static?
    // Cache of Patterns from catalog
    // TODO: make this a vector instead?
    static tetra::Pattern catCache[tetra::pattCacheSize];
    // Explore cache from beginning
    int cacheOffset = 0;
    // probeStep grows linearly (0, 1, 3, 6, 10, 15, ... ) which results in
    // quadratic probing
    int probeStep = 1;
    // Track whether a matching catalog Pattern has been found yet
    bool foundMatch = false;
    // Initialize beginning of cache in catalog to hash of our imgPattern
    uint64_t offset = hashPattern(imgPattern);
    // Start pointer to catalog at offset * sizeof(Pattern)
    // TODO: replace this stuff with database object
    fseek(pattCatalog, offset * sizeof(tetra::Pattern), SEEK_SET);
    // Fill in our cache
    fread(catCache, sizeof(tetra::Pattern), tetra::pattCacheSize, pattCatalog);
    // Perform quadratic probing through catalog
    // TODO: while hasPattern(p) => { return p.leLength > 0 }- make own function

    // cout << "Offset: " << offset << endl;
    // cout << "CAT LEN: " << catCache[cacheOffset].leLength
    //      << endl;  // if 0, no pattern
    // cout << "CAT id check: " << catCache[cacheOffset].leStarID1 << endl;

    while (catCache[cacheOffset].leLength > 0) {
        tetra::Pattern cp = catCache[cacheOffset];
        // Matching Patterns will have same bins, so check that first

        // cout << "CAT: " << cp.leStarID1 << " END" << endl;
        // cout << "CAT 2: " << cp.leStarID2 << " END" << endl;
        // cout << "OURS: " << imgPattern.leStarID1 << endl;
        // cout << "OURS 2: " << imgPattern.leStarID2
        //      << endl;  // NOTE: indices, not IDs

        // This can actually succeed
        // cout << "CHECK BINS: " << checkSameBins(imgPattern, cp) << endl;
        if (checkSameBins(imgPattern, cp)) {
            if (isMatch(imgPattern, cp)) {
                if (foundMatch) {
                    // found multiple matching Patterns in catalog, so return
                    // false
                    // cout << "ERROR: Multiple matches found" << endl;
                    return false;
                }
                // if a unique matching catalog Pattern was found, update output
                *catalogPattern = cp;
                foundMatch = true;
            }
            if (cp.isLast) {
                // cout << "LAST" << endl; // TODO: possible
                // cout << "islast" << endl;
                break;
            }
        }
        // Catalog Pattern we're looking at right now isn't a match, so
        // go to next location in cache/catalog
        cout << "increment" << endl;
        if (!incrementOffset(pattCatalog, catCache, &offset, &cacheOffset,
                             &probeStep)) {
            // This is fine?
            // cout << "increment fail" << endl;
            return false;
        }
    }
    if (foundMatch) {
        return true;
    }
    return false;
}

// TODO: used in compareBins, set in identifyStars()
// is there a better way to do this?
int leBin = 0;

/**
 * Compare 2 Features based on their bins
 * @param p
 * @param q
 * @return
 */
int TetraStarIdAlgorithm::compareBins(const tetra::Feature &p,
                                      const tetra::Feature &q) const {
    int pYBin = binYCoord(p.y, leBin, 0);
    int qYBin = binYCoord(q.y, leBin, 0);
    int pXBin = binXCoord(p.x, leBin, pYBin, 0);
    int qXBin = binXCoord(q.x, leBin, qYBin, 0);

    if (pXBin != qXBin) {
        return pXBin - qXBin;
    }
    return pYBin - qYBin;
}

/**
 *
 * @param imageStars TODO: centroided stars given in terms of Camera vectors
 * @param imageStarsInds TODO: which stars to choose as part of our Pattern
 * @param patternCatalog
 * @param matches Output
 * @return
 */
bool TetraStarIdAlgorithm::identifyStars(
    const vector<Vec3> &imageStars, int imageStarInds[tetra::numPattStars],
    FILE *patternCatalog, int matches[tetra::numPattStars][2]) const {
    // TODO: imageStarInds is always [3, 2, 1, 0] or [4, 2, 1, 0]?
    // TODO: try DIY, just pass in 4 random indices in [0, maxStars)

    // cout << "HERE" << endl; // ok

    // This is the Pattern we construct
    tetra::Pattern newPattern;
    // Catalog Pattern that uniquely matches our Pattern
    // catPattern defined by getMatchingPattern() call
    tetra::Pattern catPattern;

    // Iterate over all pairs of stairs to find and build largest edge
    float largestEdgeLength = 0.0;
    for (int i = 0; i < tetra::numPattStars; i++) {
        for (int j = i + 1; j < tetra::numPattStars; j++) {
            const Vec3 star1 = imageStars[imageStarInds[i]];
            const Vec3 star2 = imageStars[imageStarInds[j]];
            float newEdgeLength = Distance(star1, star2);
            if (newEdgeLength > largestEdgeLength) {
                largestEdgeLength = newEdgeLength;

                newPattern.leStarID1 = imageStarInds[i];
                newPattern.leStarID2 = imageStarInds[j];
            }
        }
    }

    newPattern.leLength =
        (largestEdgeLength / tetra::maxLELength) * ((1 << 16) - 1);

    // Calculate vector along x-axis of Pattern coordinate system
    Vec3 xAxis =
        imageStars[newPattern.leStarID2] - imageStars[newPattern.leStarID1];

    // Calculate vector along y-axis of Pattern coordinate system
    Vec3 yAxis = imageStars[newPattern.leStarID2].crossProduct(
        imageStars[newPattern.leStarID1]);

    xAxis = xAxis.Normalize();
    yAxis = yAxis.Normalize();

    // Initialize the Pattern Features
    int featureInd = 0;
    for (int i = 0; i < tetra::numPattStars; i++) {
        // Skip the largest edge stars
        int imageStarInd = imageStarInds[i];
        if (imageStarInd != newPattern.leStarID1 &&
            imageStarInd != newPattern.leStarID2) {
            newPattern.features[featureInd].starInd = imageStarInd;
            // Calculate normalized x, y coordinates
            // Uses vector projection
            float x = (xAxis * imageStars[imageStarInd]) / largestEdgeLength;
            float y = (yAxis * imageStars[imageStarInd]) / largestEdgeLength;
            // TODO: why are we doing this conversion
            newPattern.features[featureInd].x = x * ((1 << 14) - 1);
            newPattern.features[featureInd].y = y * ((1 << 14) - 1);

            // TODO:
            tetra::Feature *f = &newPattern.features[featureInd];
            if (f->x == 0) {
                f->x = 1;
            }
            if (f->y == 0) {
                f->y = 1;
            }
            featureInd++;
        }
    }

    int patternRotation;
    leBin = binLargestEdge(newPattern.leLength, 0);

    // Sort Pattern's Features- by x bin, then y bin
    // TODO: error
    sort(newPattern.features, newPattern.features + tetra::numPattStars - 2,
         [this](const tetra::Feature &a, const tetra::Feature &b) -> int {
             return this->compareBins(a, b);
         });

    // qsort(newPattern.features, tetra::numPattStars - 2,

    tetra::Feature firstFeature = newPattern.features[0];
    // cout << "FIRST FEATURE 1: " << firstFeature.x << endl;
    // TODO: does this not update?s
    firstFeature.x *= -1;
    firstFeature.y *= -1;
    // cout << "FIRST FEATURE 1-1: " << firstFeature.x << endl;

    patternRotation = compareBins(
        firstFeature, (newPattern.features[tetra::numPattStars - 3]));

    if (patternRotation >= 0) {
        for (int i = 0; i < tetra::numPattStars - 2; i++) {
            newPattern.features[i].x *= -1;
            newPattern.features[i].y *= -1;
        }
        for (int i = 0; i < (tetra::numPattStars - 2) / 2; i++) {
            tetra::Feature temp = newPattern.features[i];
            newPattern.features[i] =
                newPattern.features[tetra::numPattStars - 3 - i];
            newPattern.features[tetra::numPattStars - 3 - i] = temp;
        }
        int temp = newPattern.leStarID1;
        newPattern.leStarID1 = newPattern.leStarID2;
        newPattern.leStarID2 = temp;
    }

    // Find matching catalog pattern
    if (!getMatchingPattern(newPattern, &catPattern, patternCatalog)) {
        // TODO: error, failing here
        // cout << "Match fail" << endl;
        return false;
    }

    matches[0][0] = newPattern.leStarID1;
    matches[1][0] = newPattern.leStarID2;
    matches[0][1] = catPattern.leStarID1;
    matches[1][1] = catPattern.leStarID2;
    cout << matches[0][1] << ", " << matches[1][1] << endl;
    for (int i = 0; i < tetra::numPattStars - 2; i++) {
        matches[i + 2][0] = newPattern.features[i].starInd;
        matches[i + 2][1] = catPattern.features[i].starInd;
    }

    return true;
}

/*
 * Algo:
 * 0. Given list of centroids, camera, catalog
 * 1. Pick 4 random stars from list of centroids
 *  TODO: imageStars[maxStars], but what if number of centroids < maxStars:
 * choose index within appropriate range if number of centroids < numPattStars,
 * then report too few stars, skip
 * 2. Convert each centroid in stars to a Vec3, put in a vector
 * 3. Create matches array, this is our result
 * 4. Pass to identifyStars() to fill matches
 */
StarIdentifiers TetraStarIdAlgorithm::Go(const unsigned char *database,
                                         const Stars &stars,
                                         const Catalog &catalog,
                                         const Camera &camera) const {
    StarIdentifiers res;

    if ((int)stars.size() < tetra::numPattStars) {
        cerr << "Too few stars in image" << endl;
        return res;  // TODO: throw an error?
    }

    // TESTING # 1 ////
    cout << "STAR: " << sizeof(tetra::Star) << endl;
    // cout << "FEATURE SIZE: " << sizeof(tetra::Feature) << endl;

    // TESTING: remove ////////////////////////////////////////////
    bool testOn = false;
    if (testOn) {
        StarIdentifier si(18, 270);
        // correct = 399
        res.push_back(si);
        return res;
    }
    //////////////////////////////////////////////////////////////

    // TODO: change later
    FILE *patternCatalog;
    // lost/src/pattern_catalog
    patternCatalog = fopen("pattern_catalog", "rb");
    if (!patternCatalog) {
        perror("ERROR");  // No such file or directory
        // cerr << "Error: could not open pattern catalog" << endl;
        exit(EXIT_FAILURE);
    }

    // TODO: remove, testing
    // float positions[24] = {502.883148,  -210.082062, -406.207733,
    // -277.607758,
    //                        118.289345,  -249.106827, -83.203239,  103.738808,
    //                        -34.876026,  27.294966,   191.579407,  302.748718,
    //                        -341.711548, -307.523010, 27.620153,   212.193848,
    //                        -502.522675, -210.006012, 506.891174,  423.762390,
    //                        152.898392,  480.741241,  -289.067657,
    //                        304.056427};

    vector<Vec3> imageStars;
    for (int i = 0; i < (int)stars.size(); i++) {
        // for (int i = 0; i < 12; i++) {
        // Vec2 testPos;
        // testPos.x = positions[i * 2];
        // testPos.y = positions[i * 2 + 1];

        Vec3 spatialStarCoord =
            camera.CameraToSpatial(stars[i].position).Normalize();
        // Vec3 spatialStarCoord = camera.CameraToSpatial(testPos).Normalize();
        // cout << "Vec2 coordinates: " << testPos.x << ", " << testPos.y <<
        // endl;

        // cout << "Vec2 coordinates: " << stars[i].position.x << ", "
        //      << stars[i].position.y << endl;
        // cout << "Vec3 star: " << spatialStarCoord.x << ", "
        //      << spatialStarCoord.y << ", " << spatialStarCoord.z << endl;
        imageStars.push_back(spatialStarCoord);
    }
    int matches[tetra::numPattStars][2] = {{0}};
    int imageStarInds[tetra::numPattStars] = {0};

    random_device rd;
    uniform_int_distribution<unsigned> u(0, (int)stars.size() - 1);
    // default_random_engine e(time(0)); // BUG, always generates same first
    // number this is horrible, could cause star-id to ALWAYS FAIL
    set<int> chosen;

    // TODO: repeat this step as long as we failed to identify?

    // TODO: uncomment this block for actual usage
    // for (int i = 0; i < tetra::numPattStars; i++) {
    //     int r = u(rd);
    //     while (chosen.count(r) != 0) {
    //         r = u(rd);
    //     }
    //     // cout << r << endl;
    //     chosen.insert(r);

    //     imageStarInds[i] = r;

    //     cout << imageStarInds[i] << endl;
    // }
    imageStarInds[0] = 9;
    imageStarInds[1] = 18;
    imageStarInds[2] = 3;
    imageStarInds[3] = 17;

    // TODO: fill in
    // TODO: error
    if (!identifyStars(imageStars, imageStarInds, patternCatalog, matches)) {
        cerr << "ERROR: Failed to identify stars" << endl;
    };

    for (int i = 0; i < tetra::numPattStars; i++) {
        // TODO: why is this printing not-HIP numbers?
        cout << matches[i][0] << ", " << matches[i][1] << endl;
        StarIdentifier si(matches[i][0], matches[i][1]);
        // TODO: bug, matches[i][0] and matches[i][1] is always 0
        res.push_back(si);
    }
    return res;
}

StarIdentifiers DummyStarIdAlgorithm::Go(const unsigned char *database,
                                         const Stars &stars,
                                         const Catalog &catalog,
                                         const Camera &camera) const {
    StarIdentifiers result;

    for (int i = 0; i < (int)stars.size(); i++) {
        result.push_back(StarIdentifier(i, rand() % catalog.size()));
    }

    return result;
}

StarIdentifiers GeometricVotingStarIdAlgorithm::Go(
    const unsigned char *database, const Stars &stars, const Catalog &catalog,
    const Camera &camera) const {
    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);
    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(
        PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL) {
        return identified;
    }
    PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);

    for (int i = 0; i < (int)stars.size(); i++) {
        std::vector<int16_t> votes(catalog.size(), 0);
        Vec3 iSpatial = camera.CameraToSpatial(stars[i].position).Normalize();
        for (int j = 0; j < (int)stars.size(); j++) {
            if (i != j) {
                // TODO: find a faster way to do this:
                std::vector<bool> votedInPair(catalog.size(), false);
                Vec3 jSpatial =
                    camera.CameraToSpatial(stars[j].position).Normalize();
                float greatCircleDistance = AngleUnit(iSpatial, jSpatial);
                // give a greater range for min-max Query for bigger radius
                // (GreatCircleDistance)
                float lowerBoundRange = greatCircleDistance - tolerance;
                float upperBoundRange = greatCircleDistance + tolerance;
                const int16_t *upperBoundSearch;
                const int16_t *lowerBoundSearch =
                    vectorDatabase.FindPairsLiberal(
                        lowerBoundRange, upperBoundRange, &upperBoundSearch);
                // loop from lowerBoundSearch till numReturnedPairs, add one
                // vote to each star in the pairs in the datastructure
                for (const int16_t *k = lowerBoundSearch; k != upperBoundSearch;
                     k++) {
                    if ((k - lowerBoundSearch) % 2 == 0) {
                        float actualAngle = AngleUnit(
                            catalog[*k].spatial, catalog[*(k + 1)].spatial);
                        assert(actualAngle <=
                               greatCircleDistance + tolerance * 2);
                        assert(actualAngle >=
                               greatCircleDistance - tolerance * 2);
                    }
                    if (!votedInPair[*k] || true) {
                        // if (i == 542 && *k == 9085) {
                        //     printf("INC, distance %f from query %f to %f\n",
                        //     greatCircleDistance,
                        //         lowerBoundRange, upperBoundRange);
                        // }
                        votes[*k]++;
                        votedInPair[*k] = true;
                    }
                }
                // US voting system
            }
        }
        // Find star w most votes
        int16_t maxVotes = votes[0];
        int indexOfMax = 0;
        for (int v = 1; v < (int)votes.size(); v++) {
            if (votes[v] > maxVotes) {
                maxVotes = votes[v];
                indexOfMax = v;
            }
        }
        // if (i == 542) {
        //     for (float dist : vectorDatabase.StarDistances(9085, catalog)) {
        //         printf("Actual 9085 distance: %f\n", dist);
        //     }
        //     puts("Debug star.");
        //     for (int i = 0; i < (int)votes.size(); i++) {
        //         if (votes[i] > maxVotes/2) {
        //             printf("Star %4d received %d votes.\n", catalog[i].name,
        //             votes[i]);
        //         }
        //     }
        //     printf("Debug star: Actually voted for %d with %d votes\n",
        //            catalog[indexOfMax].name, maxVotes);
        // }
        // printf("Max votes: %d\n", maxVotes);
        // starIndex = i, catalog index = indexOfMax
        StarIdentifier newStar(i, indexOfMax);
        // Set identified[i] to value of catalog index of star w most votesr
        identified.push_back(newStar);
    }
    // optimizations? N^2
    // https://www.researchgate.net/publication/3007679_Geometric_voting_algorithm_for_star_trackers
    //
    //  Do we have a metric for localization uncertainty? Star brighntess?
    // loop i from 1 through n
    std::vector<int16_t> verificationVotes(identified.size(), 0);
    for (int i = 0; i < (int)identified.size(); i++) {
        // loop j from i+1 through n
        for (int j = i + 1; j < (int)identified.size(); j++) {
            // Calculate distance between catalog stars
            CatalogStar first = catalog[identified[i].catalogIndex];
            CatalogStar second = catalog[identified[j].catalogIndex];
            float cDist = AngleUnit(first.spatial, second.spatial);

            Star firstIdentified = stars[identified[i].starIndex];
            Star secondIdentified = stars[identified[j].starIndex];
            Vec3 firstSpatial =
                camera.CameraToSpatial(firstIdentified.position);
            Vec3 secondSpatial =
                camera.CameraToSpatial(secondIdentified.position);
            float sDist = Angle(firstSpatial, secondSpatial);

            // if sDist is in the range of (distance between stars in the image
            // +- R) add a vote for the match
            if (abs(sDist - cDist) < tolerance) {
                verificationVotes[i]++;
                verificationVotes[j]++;
            }
        }
    }
    // Find star w most votes
    int maxVotes = verificationVotes.size() > 0 ? verificationVotes[0] : 0;
    for (int v = 1; v < (int)verificationVotes.size(); v++) {
        if (verificationVotes[v] > maxVotes) {
            maxVotes = verificationVotes[v];
        }
    }

    // If the stars are within a certain range of the maximal number of votes,
    // we consider it correct.
    // maximal votes = maxVotes
    StarIdentifiers verified;
    int thresholdVotes = maxVotes * 3 / 4;
    printf("Verification threshold: %d\n", thresholdVotes);
    for (int i = 0; i < (int)verificationVotes.size(); i++) {
        if (verificationVotes[i] > thresholdVotes) {
            verified.push_back(identified[i]);
        }
    }

    return verified;
}

/**
 * Strategies:
 *
 * 1. For each star, enumerate all stars which have the same combination of
 * distances to some other stars, getting down to a hopefully small (<10) list
 * of candidates for each star, then do a quad-nested loop to correlate them.
 *
 * 2. Loop through all possible stars in the catalog for star i. Then look at
 * edge ij, using this to select possible j-th stars. If ever there is not a
 * possible j-th star, continue the i-loop. When a possible ij combination is
 * found, loop through k stars according to ik. IF none are found, continue the
 * outer i loop. If some are found, check jk for each one. For each possible ijk
 * triangle,
 */

class PairDistanceInvolvingIterator {
   public:
    // unqualified constructor makes a "past-the-end" iterator
    PairDistanceInvolvingIterator() : pairs(NULL), end(NULL){};

    PairDistanceInvolvingIterator(const int16_t *pairs, const int16_t *end,
                                  int16_t involving)
        : pairs(pairs), end(end), involving(involving) {
        assert((end - pairs) % 2 == 0);
        forwardUntilInvolving();
    };

    // PairDistanceInvolvingIterator operator++() {
    //     PairDistanceInvolvingIterator result(*this);
    //     ++(*this);
    //     return result;
    // }

    PairDistanceInvolvingIterator &operator++() {
        assert(hasValue());
        pairs += 2;
        forwardUntilInvolving();
        return *this;
    }

    int16_t operator*() const { return curValue; }

    bool hasValue() { return pairs != end; }

    // bool operator==(const PairDistanceInvolvingIterator &other) const {
    //     return ()other.pairs == pairs;
    // }

    // bool operator!=(const PairDistanceInvolvingIterator &other) const {
    //     return !(*this == other);
    // }
   private:
    const int16_t *pairs;
    const int16_t *end;
    int16_t involving;
    int16_t curValue;

    // like postfix++, except it's a no-op if already on a valid spot.
    void forwardUntilInvolving() {
        while (pairs != end) {
            if (pairs[0] == involving) {
                curValue = pairs[1];
                return;
            }
            if (pairs[1] == involving) {
                curValue = pairs[0];
                return;
            }
            pairs += 2;
        }
    }
};

void PyramidIdentifyRemainingStars(StarIdentifiers *identifiers,
                                   const Stars &stars, const Catalog &catalog,
                                   const PairDistanceKVectorDatabase &db,
                                   const Camera &camera, float tolerance) {
    assert(identifiers->size() == 4);
    StarIdentifiers pyramidIdentifiers =
        *identifiers;  // copy with only the pyramid's high confidence stars
    Vec3 pyramidActualSpatials[4];
    for (int l = 0; l < 4; l++) {
        pyramidActualSpatials[l] =
            camera
                .CameraToSpatial(
                    stars[pyramidIdentifiers[l].starIndex].position)
                .Normalize();
    }

    for (int p = 0; p < (int)stars.size(); p++) {
        // ensure this star isn't in the pyramid
        bool pInPyramid = false;
        for (const StarIdentifier &id : pyramidIdentifiers) {
            if (id.starIndex == p) {
                pInPyramid = true;
                break;
            }
        }
        if (pInPyramid) {
            continue;
        }

        Vec3 pSpatial = camera.CameraToSpatial(stars[p].position).Normalize();
        float ipDist = AngleUnit(pyramidActualSpatials[0], pSpatial);
        const int16_t *ipEnd;
        const int16_t *ipPairs =
            db.FindPairsLiberal(ipDist - tolerance, ipDist + tolerance, &ipEnd);
        PairDistanceInvolvingIterator pIterator(
            ipPairs, ipEnd, pyramidIdentifiers[0].catalogIndex);

        std::vector<int16_t>
            pCandidates;  // collect them all in the loop, at the
                          // end only identify the star if unique
        while (pIterator.hasValue()) {
            bool ok = true;
            for (int l = 1; l < 4; l++) {
                float actualDist =
                    AngleUnit(pSpatial, pyramidActualSpatials[l]);
                float expectedDist = AngleUnit(
                    catalog[*pIterator].spatial,
                    catalog[pyramidIdentifiers[l].catalogIndex].spatial);
                if (actualDist < expectedDist - tolerance ||
                    actualDist > expectedDist + tolerance) {
                    ok = false;
                }
            }
            if (ok) {
                pCandidates.push_back(*pIterator);
            }
            ++pIterator;
        }

        if (pCandidates.size() == 1) {
            identifiers->push_back(StarIdentifier(p, pCandidates[0]));
        }
        if (pCandidates.size() > 1) {
            std::cerr << "duplicate other star??" << std::endl;
            for (int16_t c : pCandidates) {
                std::cerr << catalog[c].name << std::endl;
            }
        }
    }
}

StarIdentifiers PyramidStarIdAlgorithm::Go(const unsigned char *database,
                                           const Stars &stars,
                                           const Catalog &catalog,
                                           const Camera &camera) const {
    StarIdentifiers identified;
    MultiDatabase multiDatabase(database);

    // cout << sizeof(tetra::Pattern) << endl;

    const unsigned char *databaseBuffer = multiDatabase.SubDatabasePointer(
        PairDistanceKVectorDatabase::kMagicValue);
    if (databaseBuffer == NULL || stars.size() < 4) {
        std::cerr << "Not enough stars, or database missing." << std::endl;
        return identified;
    }
    PairDistanceKVectorDatabase vectorDatabase(databaseBuffer);

    // smallest normal single-precision float is around 10^-38 so we should be
    // all good. See Analytic_Star_Pattern_Probability on the HSL wiki for
    // details.
    float expectedMismatchesConstant =
        pow(numFalseStars, 4) * pow(tolerance, 5) / 2 / pow(M_PI, 2);

    // this iteration technique is described in the Pyramid paper. Briefly: i
    // will always be the lowest index, then dj and dk are how many indexes
    // ahead the j-th star is from the i-th, and k-th from the j-th. In
    // addition, we here add some other numbers so that the pyramids are not
    // weird lines in wide FOV images. TODO: Select the starting points to
    // ensure that the first pyramids are all within measurement tolerance.
    int numStars = (int)stars.size();
    // the idea is that the square root is about across the FOV horizontally
    int across = floor(sqrt(numStars)) * 2;
    int halfwayAcross = floor(sqrt(numStars) / 2);
    long totalIterations = 0;

    int jMax = numStars - 3;
    for (int jIter = 0; jIter < jMax; jIter++) {
        int dj = 1 + (jIter + halfwayAcross) % jMax;

        int kMax = numStars - dj - 2;
        for (int kIter = 0; kIter < kMax; kIter++) {
            int dk = 1 + (kIter + across) % kMax;

            int rMax = numStars - dj - dk - 1;
            for (int rIter = 0; rIter < rMax; rIter++) {
                int dr = 1 + (rIter + halfwayAcross) % rMax;

                int iMax = numStars - dj - dk - dr - 1;
                for (int iIter = 0; iIter <= iMax; iIter++) {
                    int i = (iIter + iMax / 2) %
                            (iMax + 1);  // start near the center of the photo

                    // identification failure due to cutoff
                    if (++totalIterations > cutoff) {
                        std::cerr << "Cutoff reached." << std::endl;
                        return identified;
                    }

                    int j = i + dj;
                    int k = j + dk;
                    int r = k + dr;

                    assert(i != j && j != k && k != r && i != k && i != r &&
                           j != r);

                    // TODO: move this out of the loop?
                    Vec3 iSpatial =
                        camera.CameraToSpatial(stars[i].position).Normalize();
                    Vec3 jSpatial =
                        camera.CameraToSpatial(stars[j].position).Normalize();
                    Vec3 kSpatial =
                        camera.CameraToSpatial(stars[k].position).Normalize();

                    float ijDist = AngleUnit(iSpatial, jSpatial);

                    float iSinInner =
                        sin(Angle(jSpatial - iSpatial, kSpatial - iSpatial));
                    float jSinInner =
                        sin(Angle(iSpatial - jSpatial, kSpatial - jSpatial));
                    float kSinInner =
                        sin(Angle(iSpatial - kSpatial, jSpatial - kSpatial));

                    // if we made it this far, all 6 angles are confirmed! Now
                    // check that this match would not often occur due to
                    // chance. See Analytic_Star_Pattern_Probability on the HSL
                    // wiki for details
                    float expectedMismatches =
                        expectedMismatchesConstant * sin(ijDist) / kSinInner /
                        std::max(std::max(iSinInner, jSinInner), kSinInner);

                    if (expectedMismatches > maxMismatchProbability) {
                        std::cout << "skip: mismatch prob." << std::endl;
                        continue;
                    }

                    Vec3 rSpatial =
                        camera.CameraToSpatial(stars[r].position).Normalize();

                    // sign of determinant, to detect flipped patterns
                    bool spectralTorch =
                        iSpatial.crossProduct(jSpatial) * kSpatial > 0;

                    float ikDist = AngleUnit(iSpatial, kSpatial);
                    float irDist = AngleUnit(iSpatial, rSpatial);
                    float jkDist = AngleUnit(jSpatial, kSpatial);
                    float jrDist = AngleUnit(jSpatial, rSpatial);
                    float krDist = AngleUnit(
                        kSpatial, rSpatial);  // TODO: we don't really need to
                                              // check krDist, if k has been
                                              // verified by i and j it's fine.

                    // we check the distances with the extra tolerance
                    // requirement to ensure that there isn't some pyramid
                    // that's just outside the database's bounds, but within
                    // measurement tolerance of the observed pyramid, since that
                    // would possibly cause a non-unique pyramid to be
                    // identified as unique.
#define _CHECK_DISTANCE(_dist)                              \
    if (_dist < vectorDatabase.MinDistance() + tolerance || \
        _dist > vectorDatabase.MaxDistance() - tolerance) { \
        continue;                                           \
    }
                    _CHECK_DISTANCE(ikDist);
                    _CHECK_DISTANCE(irDist);
                    _CHECK_DISTANCE(jkDist);
                    _CHECK_DISTANCE(jrDist);
                    _CHECK_DISTANCE(krDist);
#undef _CHECK_DISTANCE

                    const int16_t *ijEnd, *ikEnd, *irEnd;
                    const int16_t *const ijQuery =
                        vectorDatabase.FindPairsLiberal(
                            ijDist - tolerance, ijDist + tolerance, &ijEnd);
                    const int16_t *const ikQuery =
                        vectorDatabase.FindPairsLiberal(
                            ikDist - tolerance, ikDist + tolerance, &ikEnd);
                    const int16_t *const irQuery =
                        vectorDatabase.FindPairsLiberal(
                            irDist - tolerance, irDist + tolerance, &irEnd);

                    int iMatch = -1, jMatch = -1, kMatch = -1, rMatch = -1;
                    std::vector<bool> iSeen(catalog.size(), false);
                    for (const int16_t *iCandidateQuery = ijQuery;
                         iCandidateQuery != ijEnd; iCandidateQuery++) {
                        int iCandidate = *iCandidateQuery;
                        if (iSeen[iCandidate]) {
                            continue;
                        }
                        iSeen[iCandidate] = true;

                        const Vec3 &iCandidateSpatial =
                            catalog[iCandidate].spatial;

                        // TODO: caching these iterator results into vectors can
                        // improve performance, but at the cost of memory. It
                        // would be best to put some kind of guarantee on the
                        // memory usage, and then switch to using the iterator
                        // without caching if that memory limit is exceeded.
                        PairDistanceInvolvingIterator jIterator(ijQuery, ijEnd,
                                                                iCandidate);
                        PairDistanceInvolvingIterator kIterator(ikQuery, ikEnd,
                                                                iCandidate);
                        PairDistanceInvolvingIterator rIterator(irQuery, irEnd,
                                                                iCandidate);
                        std::vector<int16_t> jCandidates;
                        std::vector<int16_t> kCandidates;
                        std::vector<int16_t> rCandidates;
                        while (jIterator.hasValue()) {
                            jCandidates.push_back(*jIterator);
                            ++jIterator;
                        }
                        while (kIterator.hasValue()) {
                            kCandidates.push_back(*kIterator);
                            ++kIterator;
                        }
                        while (rIterator.hasValue()) {
                            rCandidates.push_back(*rIterator);
                            ++rIterator;
                        }
                        // TODO: break fast if any of the iterators are empty,
                        // if it's any significant performance improvement.
                        for (int16_t jCandidate : jCandidates) {
                            const Vec3 &jCandidateSpatial =
                                catalog[jCandidate].spatial;
                            Vec3 ijCandidateCross =
                                iCandidateSpatial.crossProduct(
                                    jCandidateSpatial);

                            for (int16_t kCandidate : kCandidates) {
                                Vec3 kCandidateSpatial =
                                    catalog[kCandidate].spatial;
                                bool candidateSpectralTorch =
                                    ijCandidateCross * kCandidateSpatial > 0;
                                // checking the spectral-ity early to fail fast
                                if (candidateSpectralTorch != spectralTorch) {
                                    continue;
                                }

                                // small optimization: We can calculate jk
                                // before iterating through r, so we will!
                                float jkCandidateDist = AngleUnit(
                                    jCandidateSpatial, kCandidateSpatial);
                                if (jkCandidateDist < jkDist - tolerance ||
                                    jkCandidateDist > jkDist + tolerance) {
                                    continue;
                                }

                                // TODO: if there are no jr matches, there's no
                                // reason to continue iterating through all the
                                // other k-s. Possibly enumarete all r matches,
                                // according to ir, before this loop
                                for (int16_t rCandidate : rCandidates) {
                                    const Vec3 &rCandidateSpatial =
                                        catalog[rCandidate].spatial;
                                    float jrCandidateDist = AngleUnit(
                                        jCandidateSpatial, rCandidateSpatial);
                                    float krCandidateDist;
                                    if (jrCandidateDist < jrDist - tolerance ||
                                        jrCandidateDist > jrDist + tolerance) {
                                        continue;
                                    }
                                    krCandidateDist = AngleUnit(
                                        kCandidateSpatial, rCandidateSpatial);
                                    if (krCandidateDist < krDist - tolerance ||
                                        krCandidateDist > krDist + tolerance) {
                                        continue;
                                    }

                                    // we have a match!

                                    if (iMatch == -1) {
                                        iMatch = iCandidate;
                                        jMatch = jCandidate;
                                        kMatch = kCandidate;
                                        rMatch = rCandidate;
                                    } else {
                                        // uh-oh, stinky!
                                        // TODO: test duplicate detection, it's
                                        // hard to cause it in the real
                                        // catalog...
                                        std::cerr
                                            << "Pyramid not unique, skipping..."
                                            << std::endl;
                                        goto sensorContinue;
                                    }
                                }
                            }
                        }
                    }

                    if (iMatch != -1) {
                        printf(
                            "Matched unique pyramid!\nExpected mismatches: "
                            "%e\n",
                            expectedMismatches);
                        identified.push_back(StarIdentifier(i, iMatch));
                        identified.push_back(StarIdentifier(j, jMatch));
                        identified.push_back(StarIdentifier(k, kMatch));
                        identified.push_back(StarIdentifier(r, rMatch));

                        PyramidIdentifyRemainingStars(&identified, stars,
                                                      catalog, vectorDatabase,
                                                      camera, tolerance);
                        printf("Identified an additional %d stars\n",
                               (int)identified.size() - 4);

                        return identified;
                    }

                sensorContinue:;
                }
            }
        }
    }

    std::cerr << "Tried all pyramids; none matched." << std::endl;
    return identified;
}

}  // namespace lost
