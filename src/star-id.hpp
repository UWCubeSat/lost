#ifndef STAR_ID_H
#define STAR_ID_H

#include <cmath>  // TODO: added in here
#include <vector>

#include "attitude-utils.hpp"
#include "camera.hpp"
#include "centroiders.hpp"
#include "star-utils.hpp"

using namespace std;

namespace lost {

// TODO: ok? nested namespace for Tetra constants, structs
namespace tetra {

/////////////// Constants for GenerateTetraCatalog.cpp //////

// Minimum star brightness for inclusion in catalog
// TODO: units? "lower mag values are brighter?"
const float minStarMag = 6.2;
// Minimum angle between stars (in radians)
// Can be used to remove double stars
const float minStarSeparation = 0.004;
// Load factor of catalog- how much of catalog is used
const int loadFactor = 0.6;
// We generate the catalog in cached pieces
const int maxCatCache = 100000000;
// TODO: why is different than maxProbeDepth?
const int maxCatProbeDepth = 50000;
// Number of entries in HIPPARCHOS catalog
const int hipCatalogSize = 9110;
// Current year
const int currYear = 2022;

/////////////////////////////////////////////////////////////

// Number of stars in a Pattern- recommended >= 4
// TODO: have these be passed through constructor for TetraStarIdAlgorithm?
const int numPattStars = 4;
const int numCatalogPatts = 770708495;
const int pattCacheSize = 16;
const int maxProbeDepth = 4278;

const float binSizeRatio = 3.0;

const int maxStars = 12;  // TODO: max number of stars to consider?
const float maxFov = 0.247;
// Errors:
const float maxFovError = 0.01;  // TODO: put in constructor
const float maxCentroidError = .00069054;

const float maxScaleFactor =
    fmax(tan(maxFov * (1 + maxFovError) / 2.0f) / tan(maxFov / 2.0f),
         1 - tan(maxFov * (1 - maxFovError) / 2.0f) / tan(maxFov / 2.0f));

const float maxLELength = 2 * sin(maxFov * (1 + maxFovError) / 2.0f);
const float leErrorSlope = maxScaleFactor - 1;
const float leErrorOffset = 2 * maxCentroidError / (2 - maxScaleFactor);

//////////////////// structs for Tetra ////////////////////////////////
/**
 * Feature represents a star that is NOT part of the largest edge in a Pattern
 * TODO: size of these structs matter for the catalog- could need to regenerate
 */
struct Feature {
    // normalized x, y coordinates w.r.t our constructed coordinate system
    int x : 15;
    int xBinOffset : 1;  // TODO: what are these
    int y : 15;
    int yBinOffset : 1;
    int starInd : 15;  // TODO: index?
    int pad : 1;
};

/**
 * Pattern represents the star pattern and associated coordinate system
 * This is what is stored in a given catalog position
 */
struct Pattern {
    // Features represent the stars that are NOT part of the largest edge
    // stored in order of increasing x bin, then y bin
    Feature features[numPattStars - 2];
    // Length of largest edge in Pattern TODO: might not be entirely accurate
    uint16_t leLength;
    int leBinOffset : 1;
    int leStarID1 : 15;  // ID of first star that forms the largest edge
    int leStarID2 : 15;  // ID of second star that forms the largest edge
    // Indicate whether this catalog position contains the last matching
    // Pattern in the catalog. If so, no need to probe further
    int isLast : 1;
};

////// Used in GenerateTetraCatalog.cpp only //////

// sizeof(Star) = 20
struct Star {
    Vec3 vec;
    float mag;
    int starID;  // HIP number, NOT HR
};

//////////////////////////////////////////////////////////////////

}  // namespace tetra

class StarIdAlgorithm {
   public:
    virtual StarIdentifiers Go(const unsigned char *database, const Stars &,
                               const Catalog &, const Camera &) const = 0;
    virtual ~StarIdAlgorithm(){};
};

class TetraStarIdAlgorithm : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &,
                       const Catalog &, const Camera &) const;
    // Constructor
    TetraStarIdAlgorithm();

   private:
    float getLogBinningBase(float errorSlope, float errorOffset) const;
    int logBin(float input, float errorSlope, float errorOffset) const;
    float logUnbin(int bin, float errorSlope, float errorOffset) const;
    int binLargestEdge(int largestEdgeLen, int errorRatio) const;
    float unbinLargestEdge(int bin) const;
    int binYCoord(int y, int leBin, int errorRatio) const;
    float unbinYCoord(int yBin, int leBin) const;
    int binXCoord(int x, int leBin, int yBin, int errorRatio) const;

    uint64_t hashInt(uint64_t oldHash, uint64_t key) const;
    uint64_t hashPattern(tetra::Pattern patt) const;

    bool checkSameBins(tetra::Pattern newPattern,
                       tetra::Pattern catPattern) const;
    bool isMatch(tetra::Pattern newPattern, tetra::Pattern catPattern) const;

    bool incrementOffset(FILE *pattCatalog,
                         tetra::Pattern catCache[tetra::pattCacheSize],
                         uint64_t *offset, int *cacheOffset,
                         int *probeStep) const;
    bool getMatchingPattern(tetra::Pattern imgPattern,
                            tetra::Pattern *catalogPattern,
                            FILE *pattCatalog) const;

    int compareBins(const tetra::Feature &f1, const tetra::Feature &f2) const;
    bool identifyStars(const vector<Vec3> &imageStars,
                       int imageStarInds[tetra::numPattStars],
                       FILE *patternCatalog,
                       int matches[tetra::numPattStars][2]) const;
};

class DummyStarIdAlgorithm final : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &,
                       const Catalog &, const Camera &) const;
};

class GeometricVotingStarIdAlgorithm : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &,
                       const Catalog &, const Camera &) const;
    GeometricVotingStarIdAlgorithm(float tolerance) : tolerance(tolerance){};

   private:
    float tolerance;
};

class PyramidStarIdAlgorithm final : public StarIdAlgorithm {
   public:
    StarIdentifiers Go(const unsigned char *database, const Stars &,
                       const Catalog &, const Camera &) const;
    /**
     * @param tolerance Angular tolerance in distances (measurement error)
     * @param numFalseStars an estimate of the number of false stars in the
     * whole celestial sphere (not just the field of view). Eg, if you estimate
     * 10 dead pixels in a 40 degree FOV, you'd want to multiply that up to a
     * hundred-something numFalseStars.
     * @param maxMismatchProbability The maximum allowable probability for any
     * star to be mis-id'd.
     * @param cutoff Maximum number of pyramids to iterate through.
     */
    PyramidStarIdAlgorithm(float tolerance, int numFalseStars,
                           float maxMismatchProbability, long cutoff)
        : tolerance(tolerance),
          numFalseStars(numFalseStars),
          maxMismatchProbability(maxMismatchProbability),
          cutoff(cutoff){};

   private:
    float tolerance;
    int numFalseStars;
    float maxMismatchProbability;
    long cutoff;
};

}  // namespace lost

#endif
