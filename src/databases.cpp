#include "databases.hpp"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

// TODO: check which ones are useless
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "attitude-utils.hpp"
#include "star-utils.hpp"

namespace lost {

struct KVectorPair {
    int16_t index1;
    int16_t index2;
    float distance;
};

bool CompareKVectorPairs(const KVectorPair &p1, const KVectorPair &p2) {
    return p1.distance < p2.distance;
}

// just the index part of the kvector, doesn't store the sorted list it refers to. This makes it
// flexible whether used to store star distance angles, an angle of a triangle, etc.
/**
 K-vector index layout. No magic value because its never used on its own.
 | size          | name       | description                                                 |
 |---------------+------------+-------------------------------------------------------------|
 | 4             | numEntries |                                                             |
 | sizeof float  | min        | minimum value contained in the database                     |
 | sizeof float  | max        | max value contained in index                                |
 | 4             | numBins    |                                                             |
 | 4*(numBins+1) | bins       | The `i'th bin (starting from zero) stores how many pairs of |
 |               |            | stars have a distance lesst han or equal to:                |
 |               |            | min+i*(max-min)/numBins                                     |
 */

/// The number of bytes that a kvector index will take up whe serialized
long SerializeLengthKVectorIndex(long numBins) {
    return 4 + sizeof(float) + sizeof(float) + 4 + 4 * (numBins + 1);
}

// apparently there's no easy way to accept an iterator argument. Hate java all you want, but at
// least it can do that!
// https://stackoverflow.com/questions/5054087/declare-a-function-accepting-generic-iterator or
// rather, the correct way is to use a template and just use the ++ and * operators willy-nilly,
// which will throw an error if the operators are not implemented.

/**
 * Serialize a KVector index to disk.
 * Use SerializeLengthKVectorIndex to determine how long the buffer should be.
 * @param values The actual entries the kvector should be referring to, sorted in ascending order.
 * @pre values must be sorted in ascending order!
 * @param min,max Guaranteed bounds on the entries of values
 * @todo Consider replacing min and max parameters by just calculating the min and max of values?
 * @param numBins the number of "bins" the KVector should use. A higher number makes query results
 * "tighter" but takes up more disk space. Usually should be set somewhat smaller than (max-min)
 * divided by the "width" of the typical query.
 * @param buffer[out] index is written here.
 */
void SerializeKVectorIndex(const std::vector<float> &values, float min, float max, long numBins,
                           unsigned char *buffer) {
    std::vector<int32_t> kVector(numBins + 1);  // numBins = length, all elements zero
    float binWidth = (max - min) / numBins;

    // generate the k-vector part
    // Idea: When we find the first star that's across any bin boundary, we want to update all the
    // newly sealed bins
    long lastBin = 0;  // first bin the last star belonged to
    for (int32_t i = 0; i < (int32_t)values.size(); i++) {
        if (i > 0) {
            assert(values[i] >= values[i - 1]);
        }
        assert(values[i] >= min);
        assert(values[i] <= max);
        long thisBin = (long)ceil((values[i] - min) / binWidth);  // first bin we belong to
        assert(thisBin >= 0);
        assert(thisBin <=
               numBins);  // thisBin == numBins is acceptable since kvector length == numBins + 1
        // if thisBin > lastBin, then no more stars can be added to those bins between thisBin and
        // lastBin, so set them.
        for (long bin = lastBin; bin < thisBin; bin++) {
            kVector[bin] = i;  // our index is the number of pairs with shorter distance
        }
        lastBin = thisBin;
    }
    for (long bin = lastBin; bin <= numBins; bin++) {
        kVector[bin] = values.size();
    }

    // verify kvector
    int32_t lastBinVal = -1;
    for (const int32_t &bin : kVector) {
        assert(bin >= lastBinVal);
        lastBinVal = bin;
    }

    unsigned char *bufferStart = buffer;
    // metadata fields
    *(int32_t *)buffer = values.size();
    buffer += sizeof(int32_t);
    *(float *)buffer = min;
    buffer += sizeof(float);
    *(float *)buffer = max;
    buffer += sizeof(float);
    *(int32_t *)buffer = numBins;
    buffer += sizeof(int32_t);

    // kvector index field
    // you could probably do this with memcpy instead, but the explicit loop is necessary for endian
    // concerns? TODO endianness
    for (const int32_t &bin : kVector) {
        *(int32_t *)buffer = bin;
        buffer += sizeof(int32_t);
    }

    // verify length
    assert(buffer - bufferStart == SerializeLengthKVectorIndex(numBins));
}

/// Construct from serialized buffer.
KVectorIndex::KVectorIndex(const unsigned char *buffer) {
    numValues = *(int32_t *)buffer;
    buffer += sizeof(int32_t);
    min = *(float *)buffer;
    buffer += sizeof(float);
    max = *(float *)buffer;
    buffer += sizeof(float);
    numBins = *(int32_t *)buffer;
    buffer += sizeof(int32_t);

    assert(min >= 0.0f);
    assert(max > min);
    binWidth = (max - min) / numBins;

    bins = (const int32_t *)buffer;
}

/**
 * Finds all the entries in the given range, and possibly a few just outside the range on the ends.
 * @param upperIndex[out] Is set to the index of the last returned value +1.
 * @return the index (starting from zero) of the first value matching the query
 */
long KVectorIndex::QueryLiberal(float minQueryDistance, float maxQueryDistance,
                                long *upperIndex) const {
    assert(maxQueryDistance > minQueryDistance);
    if (maxQueryDistance >= max) {
        maxQueryDistance = max - 0.00001;  // TODO: better way to avoid hitting the bottom bin
    }
    if (minQueryDistance <= min) {
        minQueryDistance = min + 0.00001;
    }
    if (minQueryDistance > max || maxQueryDistance < min) {
        *upperIndex = 0;
        return 0;
    }
    long lowerBin = BinFor(minQueryDistance);
    long upperBin = BinFor(maxQueryDistance);
    assert(upperBin >= lowerBin);
    assert(upperBin <= numBins);
    // bins[lowerBin-1]=number of pairs <= r < query distance, so it is the index of the
    // first possible item that might be equal to the query distance
    int lowerIndex = bins[lowerBin - 1];
    if (lowerIndex >= numValues) {
        // all pairs have distance less than queried. Return value is irrelevant as long as
        // numReturned=0
        return 0;
    }
    // bins[upperBin]=number of pairs <= r >= query distance
    *upperIndex = bins[upperBin];
    return lowerIndex;
}

/// return the lowest-indexed bin that contains the number of pairs with distance <= dist
long KVectorIndex::BinFor(float query) const {
    long result = (long)ceil((query - min) / binWidth);
    assert(result >= 0);
    assert(result <= numBins);
    return result;
}

/**
 pair K-vector database layout. The kvector appears before the bulk pair data because it contains
 the number of pairs, which is necessary to read the bulk pair data.

     | size (bytes)             | name         | description |
     |--------------------------+--------------+-------------------------------------------------------------|
     | sizeof kvectorIndex      | kVectorIndex | Serialized KVector index | |
 2*sizeof(int16)*numPairs | pairs        | Bulk pair data |
 */
std::vector<KVectorPair> CatalogToPairDistances(const Catalog &catalog, float minDistance,
                                                float maxDistance) {
    std::vector<KVectorPair> result;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t k = i + 1; k < (int16_t)catalog.size(); k++) {
            KVectorPair pair = {i, k, AngleUnit(catalog[i].spatial, catalog[k].spatial)};
            assert(isfinite(pair.distance));
            assert(pair.distance >= 0);
            assert(pair.distance <= M_PI);

            if (pair.distance >= minDistance && pair.distance <= maxDistance) {
                // we'll sort later
                result.push_back(pair);
            }
        }
    }
    return result;
}

long SerializeLengthPairDistanceKVector(long numPairs, long numBins) {
    return SerializeLengthKVectorIndex(numBins) + 2 * sizeof(int16_t) * numPairs;
}

/// Number of bytes that a serialized KVectorDatabase will take up
long SerializeLengthPairDistanceKVector(const Catalog &catalog, float minDistance,
                                        float maxDistance, long numBins) {
    return SerializeLengthPairDistanceKVector(
        CatalogToPairDistances(catalog, minDistance, maxDistance).size(), numBins);
}

/**
 * Serialize a pair-distance KVector into buffer.
 * Use SerializeLengthPairDistanceKVector to determine how large the buffer needs to be. See command
 * line documentation for other options.
 */
void SerializePairDistanceKVector(const Catalog &catalog, float minDistance, float maxDistance,
                                  long numBins, unsigned char *buffer) {
    std::vector<int32_t> kVector(numBins + 1);  // numBins = length, all elements zero
    std::vector<KVectorPair> pairs = CatalogToPairDistances(catalog, minDistance, maxDistance);

    // sort pairs in increasing order.
    std::sort(pairs.begin(), pairs.end(), CompareKVectorPairs);
    std::vector<float> distances;

    for (const KVectorPair &pair : pairs) {
        distances.push_back(pair.distance);
    }

    unsigned char *bufferStart = buffer;

    // index field
    SerializeKVectorIndex(distances, minDistance, maxDistance, numBins, buffer);
    buffer += SerializeLengthKVectorIndex(numBins);

    // bulk pairs field
    for (const KVectorPair &pair : pairs) {
        *(int16_t *)buffer = pair.index1;
        buffer += sizeof(int16_t);
        *(int16_t *)buffer = pair.index2;
        buffer += sizeof(int16_t);
    }

    // verify length
    assert(buffer - bufferStart == SerializeLengthPairDistanceKVector(pairs.size(), numBins));
}

std::pair<std::vector<uint16_t>, std::vector<uint16_t>> TetraPreparePattCat(const Catalog &catalog,
                                                                            const float maxFovDeg) {
    // Would not recommend changing these parameters!
    // Currently optimal to generate many patterns with smaller FOV
    // If you make the maxFOV of Tetra database larger, you need to scale the following 2 parameters
    // up note that this will cause number of patterns to grow exponentially
    const int pattStarsPerFOV = 10;
    const int verificationStarsPerFOV = 20;

    // To eliminate double stars, specify that star must be > 0.05 degrees apart
    const float starMinSep = 0.05;

    const float maxFOV = DegToRad(maxFovDeg);

    int numEntries = catalog.size();
    int keepForPattCount = 1;
    std::vector<bool> keepForPatterns(numEntries, false);
    std::vector<bool> keepForVerifying(numEntries, false);

    // NOTE: pattern set will always be a subset of verification set
    // We technically don't even need the verification set unless we plan
    // to calculate probability of mismatch - (which we probably should)

    // Definitely keep the first star
    keepForPatterns[0] = true;
    keepForVerifying[0] = true;

    for (int i = 1; i < numEntries; i++) {
        // Spatial vector representing new star
        Vec3 vec = catalog[i].spatial;

        bool anglesForPattOK = true;
        int numPattStarsInFov = 0;

        // We should test each new star's angular distance to all stars
        // we've already selected to be kept for pattern construction
        // Stop early if:
        // a) double star: angle < min separation allowed
        // b) Number of stars in region maxFov/2 >= pattStarsPerFOV
        for (int j = 0; j < i; j++) {
            if (keepForPatterns[j]) {
                float angle = Angle(vec, catalog[j].spatial);
                if(angle < DegToRad(starMinSep)){
                    anglesForPattOK = false;
                    break;
                }
                // If angle between new star i and old star j is less than maxFov/2, OK
                if(angle < maxFOV/2){
                    numPattStarsInFov++;
                    if (numPattStarsInFov >= pattStarsPerFOV) {
                        anglesForPattOK = false;
                        break;
                    }
                }
            }
        }

        if (anglesForPattOK) {
            keepForPatterns[i] = true;
            keepForVerifying[i] = true;
            keepForPattCount++;
            continue;
        }

        bool anglesForVerifOK = true;
        int numVerStarsInFov = 0;

        // Same thing here, we should test each new star's angular distance to all
        // stars we've already selected to be kept for verification
        // std::vector<float> angsVerifying;
        for (int j = 0; j < i; j++) {
            if (keepForVerifying[j]) {
                float angle = Angle(vec, catalog[j].spatial);
                if(angle < DegToRad(starMinSep)){
                    anglesForVerifOK = false;
                    break;
                }
                if(angle < maxFOV/2){
                    numVerStarsInFov++;
                    if (numVerStarsInFov >= verificationStarsPerFOV) {
                        anglesForVerifOK = false;
                        break;
                    }
                }
            }
        }

        if (anglesForVerifOK) {
            keepForVerifying[i] = true;
        }
    }

    // List of indices, track which stars in our star catalog are OK for Tetra to use
    std::vector<uint16_t> finalCatIndices;
    // List of indices, indexing into finalCatIndices
    // Track which stars in the final star table should be used for pattern construction
    // later in Tetra's database generation step
    std::vector<uint16_t> pattStarIndices;

    // finalCat is the final version of the star table
    for (int i = 0; i < (int)keepForVerifying.size(); i++) {
        if (keepForVerifying[i]) {
            finalCatIndices.push_back(i);
        }
    }

    // Pretty clever way of finding which stars in the final star table
    // should be used for pattern construction later in Tetra's database generation step
    uint16_t cumulativeSum = -1;
    for (int i = 0; i < (int)keepForVerifying.size(); i++) {
        if (keepForVerifying[i]) {
            cumulativeSum++;
        }
        if (keepForPatterns[i]) {
            pattStarIndices.push_back(cumulativeSum);
        }
    }
    return std::pair<std::vector<uint16_t>, std::vector<uint16_t>>{finalCatIndices, pattStarIndices};
}

TetraDatabase::TetraDatabase(const unsigned char *buffer) : buffer_(buffer) {
    maxAngle_ = *(float*)buffer;
    buffer += sizeof(float);
    catalogSize_ = *(int32_t *)buffer;
    std::cout << "Tetra database, size= " << catalogSize_ << std::endl;
}

float TetraDatabase::MaxAngle() const { return maxAngle_; }

int TetraDatabase::Size() const { return catalogSize_; }

std::vector<int> TetraDatabase::GetPattern(int index) const {
    std::vector<int> res;
    const unsigned char *p = buffer_ + headerSize;
    for (int i = 0; i < 4; i++) {
        res.push_back(*((uint16_t *)p + 4 * index + i));
    }

    return res;
}

uint16_t TetraDatabase::GetTrueCatInd(int tetraInd) const {
    // TODO: don't harcode this 4
    const unsigned char *p = buffer_ + headerSize + Size() * 4 * sizeof(uint16_t);
    return *((uint16_t *)p + tetraInd);
}

/// Create the database from a serialized buffer.
PairDistanceKVectorDatabase::PairDistanceKVectorDatabase(const unsigned char *buffer)
    : index(KVectorIndex(buffer)) {
    // TODO: errors? (not even sure what i meant by this comment anymore)
    buffer += SerializeLengthKVectorIndex(index.NumBins());
    pairs = (const int16_t *)buffer;
}

/// Return the value in the range [low,high] which is closest to num
float Clamp(float num, float low, float high) { return num < low ? low : num > high ? high : num; }

/**
 * Return at least all the star pairs whose inter-star distance is between min and max
 * @param end[out] Is set to an "off-the-end" pointer, one past the last pair being returned by the
 * query.
 * @return A pointer to the start of the matched pairs. Each pair is stored as simply two 16-bit
 * integers, each of which is a catalog index. (you must increment the pointer twice to get to the
 * next pair).
 */
const int16_t *PairDistanceKVectorDatabase::FindPairsLiberal(float minQueryDistance,
                                                             float maxQueryDistance,
                                                             const int16_t **end) const {
    assert(maxQueryDistance <= M_PI);

    long upperIndex = -1;
    long lowerIndex = index.QueryLiberal(minQueryDistance, maxQueryDistance, &upperIndex);
    *end = &pairs[upperIndex * 2];
    return &pairs[lowerIndex * 2];
}

const int16_t *PairDistanceKVectorDatabase::FindPairsExact(const Catalog &catalog,
                                                           float minQueryDistance,
                                                           float maxQueryDistance,
                                                           const int16_t **end) const {
    // Instead of computing the angle for every pair in the database, we pre-compute the /cosines/
    // of the min and max query distances so that we can compare against dot products directly! As
    // angle increases, cosine decreases, up to M_PI (and queries larger than that don't really make
    // sense anyway)
    assert(maxQueryDistance <= M_PI);

    float maxQueryCos = cos(minQueryDistance);
    float minQueryCos = cos(maxQueryDistance);

    long liberalUpperIndex;
    long liberalLowerIndex =
        index.QueryLiberal(minQueryDistance, maxQueryDistance, &liberalUpperIndex);
    // now we need to find the first and last index that actually matches the query
    // step the lower index forward
    // There's no good reason to be using >= and <= for the comparison against max/min, but the
    // tests fail otherwise (because they use angle, with its acos, instead of forward cos like us).
    // It's an insignificant difference.
    while (liberalLowerIndex < liberalUpperIndex &&
           catalog[pairs[liberalLowerIndex * 2]].spatial *
                   catalog[pairs[liberalLowerIndex * 2 + 1]].spatial >=
               maxQueryCos) {
        liberalLowerIndex++;
    }

    // step the upper index backward
    while (liberalLowerIndex < liberalUpperIndex
           // the liberalUpperIndex is past the end of the logically returned range, so we need to
           // subtract 1
           && catalog[pairs[(liberalUpperIndex - 1) * 2]].spatial *
                      catalog[pairs[(liberalUpperIndex - 1) * 2 + 1]].spatial <=
                  minQueryCos) {
        liberalUpperIndex--;
    }

    *end = &pairs[liberalUpperIndex * 2];
    return &pairs[liberalLowerIndex * 2];
}

/// Number of star pairs stored in the database
long PairDistanceKVectorDatabase::NumPairs() const { return index.NumValues(); }

/// Return the distances from the given star to each star it's paired with in the database (for
/// debugging).
std::vector<float> PairDistanceKVectorDatabase::StarDistances(int16_t star,
                                                              const Catalog &catalog) const {
    std::vector<float> result;
    for (int i = 0; i < NumPairs(); i++) {
        if (pairs[i * 2] == star || pairs[i * 2 + 1] == star) {
            result.push_back(
                AngleUnit(catalog[pairs[i * 2]].spatial, catalog[pairs[i * 2 + 1]].spatial));
        }
    }
    return result;
}

///////////////////// Tetra database //////////////////////

static long SerializeTetraHelper(const Catalog &catalog, float maxFovDeg, unsigned char *buffer,
                            const std::vector<uint16_t> &pattStarIndices,
                            const std::vector<uint16_t> &catIndices, bool ser) {
    const float maxFov = DegToRad(maxFovDeg);

    // TODO: pattBins here and numPattBins in TetraStarIDAlgorithm::numPattBins should be the same
    // Otherwise we break things
    const uint16_t pattBins = 50;
    const int tempBins = 4;

    Catalog tetraCatalog;
    for (uint16_t ind : catIndices) {
        tetraCatalog.push_back(catalog[ind]);
    }

    std::map<Vec3, std::vector<uint16_t>> tempCoarseSkyMap;

    for (const uint16_t starID : pattStarIndices) {
        Vec3 v(tetraCatalog[starID].spatial);
        Vec3 hash{floor((v.x + 1) * tempBins), floor((v.y + 1) * tempBins),
                  floor((v.z + 1) * tempBins)};
        tempCoarseSkyMap[hash].push_back(starID);
    }

    // Return a list of stars that are nearby
    auto tempGetNearbyStars = [&tempCoarseSkyMap, &tetraCatalog](const Vec3 &vec, float radius) {
        std::vector<float> components{vec.x, vec.y, vec.z};

        std::vector<std::vector<int>> hcSpace;
        for (float x : components) {
            std::vector<int> range;
            int lo = int((x + 1 - radius) * tempBins);
            lo = std::max(lo, 0);
            int hi = int((x + 1 + radius) * tempBins);
            hi = std::min(hi + 1, 2 * tempBins);
            range.push_back(lo);
            range.push_back(hi);
            hcSpace.push_back(range);
        }

        // hashcode space has 3 ranges, one for each of [x, y, z]

        std::vector<uint16_t> nearbyStarIDs;

        for (int a = hcSpace[0][0]; a < hcSpace[0][1]; a++) {
            for (int b = hcSpace[1][0]; b < hcSpace[1][1]; b++) {
                for (int c = hcSpace[2][0]; c < hcSpace[2][1]; c++) {
                    Vec3 code{static_cast<float>(a), static_cast<float>(b), static_cast<float>(c)};

                    // For each star j in partition with key=code,
                    // see if our star and j have angle < radius. If so, they are nearby
                    for (uint16_t starID : tempCoarseSkyMap[code]) {
                        float dotProd = vec * tetraCatalog[starID].spatial;
                        if (dotProd > std::cos(radius)) {
                            nearbyStarIDs.push_back(starID);
                        }
                    }
                }
            }
        }

        return nearbyStarIDs;
    };

    const int pattSize = 4;
    typedef std::array<uint16_t, pattSize> Pattern;
    std::vector<Pattern> pattList;
    Pattern patt{0, 0, 0, 0};

    // Construct all possible patterns
    for (int ind = 0; ind < (int)pattStarIndices.size(); ind++) {
        uint16_t firstStarID = pattStarIndices[ind];
        patt[0] = firstStarID;

        Vec3 v(tetraCatalog[firstStarID].spatial);
        Vec3 hashCode{floor((v.x + 1) * tempBins), floor((v.y + 1) * tempBins),
                      floor((v.z + 1) * tempBins)};

        // Remove star i from its sky map partition
        auto removeIt = std::find(tempCoarseSkyMap[hashCode].begin(),
                                  tempCoarseSkyMap[hashCode].end(), firstStarID);
        if (removeIt == tempCoarseSkyMap[hashCode].end()) {
            std::cerr << "Fatal error with hashing" << std::endl;
            exit(EXIT_FAILURE);
        }
        tempCoarseSkyMap[hashCode].erase(removeIt);

        auto nearbyStars = tempGetNearbyStars(v, maxFov);
        const int n = nearbyStars.size();

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                for (int k = j + 1; k < n; k++) {
                    patt[1] = nearbyStars[i];
                    patt[2] = nearbyStars[j];
                    patt[3] = nearbyStars[k];

                    bool pattFits = true;
                    for (int pair1 = 0; pair1 < pattSize; pair1++) {
                        for (int pair2 = pair1 + 1; pair2 < pattSize; pair2++) {
                            float dotProd = tetraCatalog[patt[pair1]].spatial *
                                            tetraCatalog[patt[pair2]].spatial;
                            if (dotProd <= std::cos(maxFov)) {
                                pattFits = false;
                                break;
                            }
                        }
                    }

                    if (pattFits) {
                        pattList.push_back(patt);
                    }
                }
            }
        }
    }

    std::cout << "Found " << pattList.size() << " patterns" << std::endl;

    // Load factor of 0.5
    long long pattCatalogLength = 2 * (int)pattList.size();
    std::vector<Pattern> pattCatalog(pattCatalogLength);

    for (Pattern patt : pattList) {
        std::vector<float> pattEdgeLengths;
        for (int i = 0; i < pattSize; i++) {
            CatalogStar star1 = tetraCatalog[patt[i]];
            for (int j = i + 1; j < pattSize; j++) {
                // calculate distance between vectors
                CatalogStar star2 = tetraCatalog[patt[j]];
                float edgeLen = (star2.spatial - star1.spatial).Magnitude();
                pattEdgeLengths.push_back(edgeLen);
            }
        }
        std::sort(pattEdgeLengths.begin(), pattEdgeLengths.end());

        float pattLargestEdge = pattEdgeLengths[(int)(pattEdgeLengths.size() - 1)];
        std::vector<float> pattEdgeRatios;
        for (int i = 0; i < (int)pattEdgeLengths.size() - 1;
             i++) {  // size()-1, since we ignore the largest edge
            pattEdgeRatios.push_back(pattEdgeLengths[i] / pattLargestEdge);
        }

        std::vector<int> key;
        for (float edgeRatio : pattEdgeRatios) {
            key.push_back(int(edgeRatio * pattBins));
        }

        int hashIndex = KeyToIndex(key, pattBins, pattCatalogLength);
        long long offset = 0;
        // Quadratic probing to find next available bucket for element with key=hashIndex
        while (true) {
            int index = int(hashIndex + std::pow(offset, 2)) % pattCatalogLength;
            offset++;
            if (pattCatalog[index][0] == 0 && pattCatalog[index][1] == 0) {
                pattCatalog[index] = patt;
                break;
            }
        }
    }

    // Done with everything, write to buffer
    // See TetraDatabase for layout of db

    if (ser) {
        *((float *)buffer) = maxFovDeg;
        buffer += sizeof(float);
        *((int32_t *)buffer) = (int)pattCatalog.size();
        buffer += sizeof(int32_t);

        for (Pattern patt : pattCatalog) {
            for (int i = 0; i < pattSize; i++) {
                *((uint16_t *)buffer) = patt[i];
                buffer += sizeof(uint16_t);
            }
        }

        for (uint16_t ind : catIndices) {
            *((uint16_t *)buffer) = ind;
            buffer += sizeof(uint16_t);
        }

        return -1;
    } else {
        return sizeof(float) + sizeof(int32_t) + sizeof(uint16_t) * pattCatalog.size() * pattSize +
               catIndices.size() * sizeof(uint16_t);
    }
}

long SerializeLengthTetraDatabase(const Catalog &cat, float maxFov,
                                  const std::vector<uint16_t> &pattStarIndices,
                                  const std::vector<uint16_t> &catIndices) {
    return SerializeTetraHelper(cat, maxFov, nullptr, pattStarIndices, catIndices, false);
}

void SerializeTetraDatabase(const Catalog &cat, float maxFov, unsigned char *buffer,
                            const std::vector<uint16_t> &pattStarIndices,
                            const std::vector<uint16_t> &catIndices) {
    SerializeTetraHelper(cat, maxFov, buffer, pattStarIndices, catIndices, true);
}

/**
   MultiDatabase memory layout:

   | size           | name              | description                                             |
   |----------------+-------------------+---------------------------------------------------------|
   | 8*maxDatabases | table of contents | each 8-byte entry is the 4-byte magic value followed by |
   |                |                   | a 4-byte index into the bulk where that db begins       |
   | Large          | databases         | the database contents                                   |
 */

/**
 * @brief return a pointer to the start of the database type indicated by the magic value, if such
 * a sub-database is present in the database
 * @param magicValue
 * @return Returns a pointer to the start of the database type indicated by the magic value, null if
 * not found
 */
const unsigned char *MultiDatabase::SubDatabasePointer(int32_t magicValue) const {
    long databaseIndex = -1;
    int32_t *toc = (int32_t *)buffer;
    for (int i = 0; i < kMultiDatabaseMaxDatabases; i++) {
        int32_t curMagicValue = *toc;
        toc++;
        if (curMagicValue == magicValue) {
            databaseIndex = *toc;
            break;
        }
        toc++;
    }
    // the database was not found
    if (databaseIndex < 0) {
        return NULL;
    }

    return buffer + kMultiDatabaseTocLength + databaseIndex;
}

/**
 * Add a database to a MultiDatabase
 * @param magicValue A value unique to this type of database which is used to extract it out of the
 * database later.
 * @param length The number of bytes to allocate for this database.
 * @return Pointer to the start of the space allocated for said database. Return null if full (too
 * many databases).
 */
unsigned char *MultiDatabaseBuilder::AddSubDatabase(int32_t magicValue, long length) {
    // find unused spot in toc and take it!
    int32_t *toc = (int32_t *)buffer;
    bool foundSpot = false;
    for (int i = 0; i < kMultiDatabaseMaxDatabases; i++) {
        if (*toc == 0) {
            *toc = magicValue;
            toc++;
            *toc = bulkLength;
            foundSpot = true;
            break;
        }
        // skip the entry
        toc += 2;
    }

    // database is full
    if (!foundSpot) {
        return NULL;
    }

    buffer = (unsigned char *)realloc(buffer, kMultiDatabaseTocLength + bulkLength + length);
    // just past the end of the last database
    unsigned char *result = buffer + kMultiDatabaseTocLength + bulkLength;
    bulkLength += length;
    return result;
}

MultiDatabaseBuilder::~MultiDatabaseBuilder() { free(buffer); }

}  // namespace lost

// TODO: after creating the database, print more statistics, such as average number of pairs per
// star, stars per bin, space distribution between array and index.
