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
#include <utility>  // pair
#include <vector>

#include "attitude-utils.hpp"
#include "serialize-helpers.hpp"
#include "star-utils.hpp"

namespace lost {

const int32_t PairDistanceKVectorDatabase::kMagicValue = 0x2536f009;
const int32_t TetraDatabase::kMagicValue = 0x26683787;

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
 * @param numBins the number of "bins" the KVector should use. A higher number makes query results "tighter" but takes up more disk space. Usually should be set somewhat smaller than (max-min) divided by the "width" of the typical query.
 * @param buffer[out] index is written here.
 */
void SerializeKVectorIndex(SerializeContext *ser, const std::vector<float> &values, float min, float max, long numBins) {
    std::vector<int32_t> kVector(numBins+1); // We store sums before and after each bin
    float binWidth = (max - min) / numBins;

    // generate the k-vector part
    // Idea: When we find the first star that's across any bin boundary, we want to update all the newly sealed bins
    long lastBin = 0; // first bin the last star belonged to
    for (int32_t i = 0; i < (int32_t)values.size(); i++) {
        if (i > 0) {
            assert(values[i] >= values[i-1]);
        }
        assert(values[i] >= min);
        assert(values[i] <= max);
        long thisBin = (long)ceil((values[i] - min) / binWidth); // first bin we belong to
        assert(thisBin >= 0);
        assert(thisBin <= numBins); // thisBin == numBins is acceptable since kvector length == numBins + 1
        // if thisBin > lastBin, then no more stars can be added to those bins between thisBin and lastBin, so set them.
        for (long bin = lastBin; bin < thisBin; bin++) {
            kVector[bin] = i; // our index is the number of pairs with shorter distance
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

    // metadata fields
    SerializePrimitive<int32_t>(ser, values.size());
    SerializePrimitive<float>(ser, min);
    SerializePrimitive<float>(ser, max);
    SerializePrimitive<int32_t>(ser, numBins);

    // kvector index field
    for (const int32_t &bin : kVector) {
        SerializePrimitive<int32_t>(ser, bin);
    }
}

/// Construct from serialized buffer.
KVectorIndex::KVectorIndex(DeserializeContext *des) {

    numValues = DeserializePrimitive<int32_t>(des);
    min = DeserializePrimitive<float>(des);
    max = DeserializePrimitive<float>(des);
    numBins = DeserializePrimitive<int32_t>(des);

    assert(min >= 0.0f);
    assert(max > min);
    binWidth = (max - min) / numBins;

    bins = DeserializeArray<int32_t>(des, numBins+1);
}

/**
 * Finds all the entries in the given range, and possibly a few just outside the range on the ends.
 * @param upperIndex[out] Is set to the index of the last returned value +1.
 * @return the index (starting from zero) of the first value matching the query
 */
long KVectorIndex::QueryLiberal(float minQueryDistance, float maxQueryDistance, long *upperIndex) const {
    assert(maxQueryDistance > minQueryDistance);
    if (maxQueryDistance >= max) {
        maxQueryDistance = max - 0.00001; // TODO: better way to avoid hitting the bottom bin
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
    int lowerIndex = bins[lowerBin-1];
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
 pair K-vector database layout. The kvector appears before the bulk pair data because it contains the
 number of pairs, which is necessary to read the bulk pair data.

     | size (bytes)             | name         | description                                                 |
     |--------------------------+--------------+-------------------------------------------------------------|
     | sizeof kvectorIndex      | kVectorIndex | Serialized KVector index                                    |
     | 2*sizeof(int16)*numPairs | pairs        | Bulk pair data                                              |
 */
std::vector<KVectorPair> CatalogToPairDistances(const Catalog &catalog, float minDistance, float maxDistance) {
    std::vector<KVectorPair> result;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t k = i+1; k < (int16_t)catalog.size(); k++) {

            KVectorPair pair = { i, k, AngleUnit(catalog[i].spatial, catalog[k].spatial) };
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

/**
 * Serialize a pair-distance KVector into buffer.
 * Use SerializeLengthPairDistanceKVector to determine how large the buffer needs to be. See command line documentation for other options.
 */
void SerializePairDistanceKVector(SerializeContext *ser, const Catalog &catalog, float minDistance, float maxDistance, long numBins) {
    std::vector<int32_t> kVector(numBins+1); // numBins = length, all elements zero
    std::vector<KVectorPair> pairs = CatalogToPairDistances(catalog, minDistance, maxDistance);

    // sort pairs in increasing order.
    std::sort(pairs.begin(), pairs.end(), CompareKVectorPairs);
    std::vector<float> distances;

    for (const KVectorPair &pair : pairs) {
        distances.push_back(pair.distance);
    }

    // index field
    SerializeKVectorIndex(ser, distances, minDistance, maxDistance, numBins);

    // bulk pairs field
    for (const KVectorPair &pair : pairs) {
        SerializePrimitive<int16_t>(ser, pair.index1);
        SerializePrimitive<int16_t>(ser, pair.index2);
    }
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
                if (angle < DegToRad(starMinSep)) {
                    anglesForPattOK = false;
                    break;
                }
                // If angle between new star i and old star j is less than maxFov/2, OK
                if (angle < maxFOV / 2) {
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
                if (angle < DegToRad(starMinSep)) {
                    anglesForVerifOK = false;
                    break;
                }
                if (angle < maxFOV / 2) {
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

    // Find which stars in the final star table
    // should be used for pattern construction later in Tetra's database generation step
    // pattStarIndices will be a double-index:
    //  our "final star table" is a list of indices into the original (unchanged) catalog
    //  pattStarIndices is a list of indices into the final star table to tell Tetra which ones
    //  to use for pattern construction
    int cumulativeSum = -1;
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

TetraDatabase::TetraDatabase(DeserializeContext *des) {
    maxAngle_ = DeserializePrimitive<float>(des);
    pattCatSize_ = DeserializePrimitive<uint64_t>(des);
    tetraStarCatSize_ = DeserializePrimitive<uint64_t>(des);
    pattCats_ = DeserializeArray<uint16_t>(des, 4 * pattCatSize_);
    starCatInds_ = DeserializeArray<uint16_t>(des, tetraStarCatSize_);
}

TetraPatt TetraDatabase::GetPattern(int index) const {
    TetraPatt res;
    for (int i = 0; i < 4; i++) {
        res.push_back(pattCats_[4*index + i]);
    }
    return res;
}

uint16_t TetraDatabase::GetTrueCatInd(int tetraInd) const {
    return starCatInds_[tetraInd];
}

std::vector<TetraPatt> TetraDatabase::GetPatternMatches(int index) const {
    std::vector<TetraPatt> res;
    for (int c = 0;; c++) {
        int i = (index + c * c) % PattCatSize();

        TetraPatt tableRow = GetPattern(i);

        if (tableRow[0] == 0 && tableRow[1] == 0) {
            break;
        } else {
            res.push_back(tableRow);
        }
    }
    return res;
}

/// Create the database from a serialized buffer.
PairDistanceKVectorDatabase::PairDistanceKVectorDatabase(DeserializeContext *des)
    : index(KVectorIndex(des)) {

    pairs = DeserializeArray<int16_t>(des, 2*index.NumValues());
}

/// Return the value in the range [low,high] which is closest to num
float Clamp(float num, float low, float high) {
    return num < low ? low : num > high ? high : num;
}

/**
 * Return at least all the star pairs whose inter-star distance is between min and max
 * @param end[out] Is set to an "off-the-end" pointer, one past the last pair being returned by the query.
 * @return A pointer to the start of the matched pairs. Each pair is stored as simply two 16-bit integers, each of which is a catalog index. (you must increment the pointer twice to get to the next pair).
 */
const int16_t *PairDistanceKVectorDatabase::FindPairsLiberal(
    float minQueryDistance, float maxQueryDistance, const int16_t **end) const {

    assert(maxQueryDistance <= M_PI);

    long upperIndex = -1;
    long lowerIndex = index.QueryLiberal(minQueryDistance, maxQueryDistance, &upperIndex);
    *end = &pairs[upperIndex * 2];
    return &pairs[lowerIndex * 2];
}

const int16_t *PairDistanceKVectorDatabase::FindPairsExact(const Catalog &catalog,
                                                           float minQueryDistance, float maxQueryDistance, const int16_t **end) const {

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

void SerializeTetraDatabase(SerializeContext *ser, const Catalog &catalog, float maxFovDeg,
                            const std::vector<uint16_t> &pattStarIndices,
                            const std::vector<uint16_t> &catIndices) {
    const float maxFovRad = DegToRad(maxFovDeg);

    // TODO: these are hardcoded values
    // pattBins here and numPattBins in TetraStarIDAlgorithm::numPattBins must be the same
    const int pattBins = 50;
    const int tempBins = 4;
    const int pattSize = 4;

    Catalog tetraCatalog;
    for (int ind : catIndices){
        tetraCatalog.push_back(catalog[ind]);
    }

    auto spatialHash = [](const Vec3 &vec){
        std::hash<float> hasher;
        return hasher(vec.x) ^ hasher(vec.y) ^ hasher(vec.z);
    };
    auto spatialEquals = [](const Vec3 &vec1, const Vec3 &vec2){
        return vec1.x==vec2.x && vec1.y==vec2.y && vec1.z==vec2.z;
    };
    std::unordered_map<Vec3, std::vector<uint16_t>, decltype(spatialHash), decltype(spatialEquals)>
        tempCoarseSkyMap(8, spatialHash, spatialEquals);

    for (int starID : pattStarIndices){
        Vec3 v{tetraCatalog[starID].spatial};
        Vec3 hash{
            floor((v.x+1) * tempBins),
            floor((v.y+1) * tempBins),
            floor((v.z+1) * tempBins)
        };
        tempCoarseSkyMap[hash].push_back(starID);
    }

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
        // Hashcode space has 3 ranges, one for each of [x, y, z]
        std::vector<uint16_t> nearbyStarIDs; // TODO: typedef this
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
    using Pattern = std::array<uint16_t, pattSize>;
    std::vector<Pattern> pattList;
    Pattern patt{0, 0, 0, 0};

    // Construct all possible patterns
    for(int firstStarInd : pattStarIndices){
        patt[0] = firstStarInd;

        Vec3 v{tetraCatalog[firstStarInd].spatial};
        Vec3 hashCode{floor((v.x + 1) * tempBins), floor((v.y + 1) * tempBins),
                      floor((v.z + 1) * tempBins)};

        // Remove star=firstStarInd from its sky map partition
        auto removeIt = std::find(tempCoarseSkyMap[hashCode].begin(),
                                  tempCoarseSkyMap[hashCode].end(), firstStarInd);
        assert(removeIt != tempCoarseSkyMap[hashCode].end());
        tempCoarseSkyMap[hashCode].erase(removeIt);

        auto nearbyStars = tempGetNearbyStars(v, maxFovRad);
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
                            if (dotProd <= std::cos(maxFovRad)) {
                                pattFits = false;
                                break;
                            }
                        }
                        if (!pattFits) break;
                    }

                    if (pattFits) {
                        pattList.push_back(patt);
                    }
                }
            }
        }
    }

    std::cerr << "Tetra found " << pattList.size() << " patterns" << std::endl;
    // Ensure load factor just < 0.5
    long long pattCatalogLen = 2 * (int)pattList.size() + 1;
    std::vector<Pattern> pattCatalog(pattCatalogLen);
    for (Pattern patt : pattList) {
        std::vector<double> pattEdgeLengths;
        for (int i = 0; i < pattSize; i++) {
            CatalogStar star1 = tetraCatalog[patt[i]];
            for (int j = i + 1; j < pattSize; j++) {
                // calculate distance between vectors
                CatalogStar star2 = tetraCatalog[patt[j]];
                double edgeLen = (star2.spatial - star1.spatial).Magnitude();
                pattEdgeLengths.push_back(edgeLen);
            }
        }
        std::sort(pattEdgeLengths.begin(), pattEdgeLengths.end());
        double pattLargestEdge = pattEdgeLengths.back();
        std::vector<double> pattEdgeRatios;
        // Skip last edge since ratio is just 1
        for (int i = 0; i < pattEdgeLengths.size() - 1; i++) {
            pattEdgeRatios.push_back(pattEdgeLengths[i] / pattLargestEdge);
        }

        std::vector<int> key;
        for (double edgeRatio : pattEdgeRatios) {
            key.push_back(int(edgeRatio * pattBins));
        }

        int hashIndex = KeyToIndex(key, pattBins, pattCatalogLen);
        long long offset = 0;
        // Quadratic probing to find next available bucket for element with key=hashIndex
        while (true) {
            int index = int(hashIndex + std::pow(offset, 2)) % pattCatalogLen;
            offset++;
            if (pattCatalog[index][0] == 0 && pattCatalog[index][1] == 0) {
                pattCatalog[index] = patt;
                break;
            }
        }
    }
    SerializePrimitive<float>(ser, maxFovDeg);
    SerializePrimitive<uint64_t>(ser, pattCatalog.size());
    SerializePrimitive<uint64_t>(ser, catIndices.size());
    for (Pattern patt : pattCatalog) {
        for (int i = 0; i < pattSize; i++) {
            SerializePrimitive<uint16_t>(ser, patt[i]);
        }
    }
    for (uint16_t ind : catIndices) {
        SerializePrimitive<uint16_t>(ser, ind);
    }
}

/**
   MultiDatabase memory layout:

   | size | name           | description                                 |
   |------+----------------+---------------------------------------------|
   |    4 | magicValue     | unique database identifier                  |
   |    4 | databaseLength | length in bytes (32-bit unsigned)           |
   |    n | database       | the entire database. 8-byte aligned         |
   |  ... | ...            | More databases (each has value, length, db) |
   |    4 | caboose        | 4 null bytes indicate the end               |
 */

/**
 * @brief return a pointer to the start of the database type indicated by the magic value, if such
 * a sub-database is present in the database
 * @param magicValue
 * @return Returns a pointer to the start of the database type indicated by the magic value, null if not found
 */
const unsigned char *MultiDatabase::SubDatabasePointer(int32_t magicValue) const {
    DeserializeContext desValue(buffer);
    DeserializeContext *des = &desValue; // just for naming consistency with how we use `des` elsewhere

    assert(magicValue != 0);
    while (true) {
        int32_t curMagicValue = DeserializePrimitive<int32_t>(des);
        if (curMagicValue == 0) {
            return nullptr;
        }
        uint32_t dbLength = DeserializePrimitive<uint32_t>(des);
        assert(dbLength > 0);
        DeserializePadding<uint64_t>(des); // align to an 8-byte boundary
        const unsigned char *curSubDatabasePointer = DeserializeArray<unsigned char>(des, dbLength);
        if (curMagicValue == magicValue) {
            return curSubDatabasePointer;
        }
    }
    // shouldn't ever make it here. Compiler should remove this assertion as unreachable.
    assert(false);
}

void SerializeMultiDatabase(SerializeContext *ser,
                            const MultiDatabaseDescriptor &dbs) {
    for (const MultiDatabaseEntry &multiDbEntry : dbs) {
        SerializePrimitive<int32_t>(ser, multiDbEntry.magicValue);
        SerializePrimitive<uint32_t>(ser, multiDbEntry.bytes.size());
        SerializePadding<uint64_t>(ser);
        std::copy(multiDbEntry.bytes.cbegin(), multiDbEntry.bytes.cend(), std::back_inserter(ser->buffer));
    }
    SerializePrimitive<int32_t>(ser, 0); // caboose
}

}

// TODO: after creating the database, print more statistics, such as average number of pairs per
// star, stars per bin, space distribution between array and index.
