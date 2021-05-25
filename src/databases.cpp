#include "databases.hpp"
#include "attitude-utils.hpp"
#include "star-utils.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <math.h>

#define kKVectorMagicNumber 0x4253f009

namespace lost {

struct KVectorPair {
    int16_t index1;
    int16_t index2;
    float distance;
};

bool CompareKVectorPairs(const KVectorPair &p1, const KVectorPair &p2) {
    return p1.distance < p2.distance;
}

/**
 K-vector database layout
     | size (bytes) | name        | description                                                 |
     |--------------+-------------+-------------------------------------------------------------|
     |            4 | magic       | Magic number kKvectorMagicNumber (ensure database validity) |
     |            4 | numPairs    | Number of star pairs that are in the min/max distance range |
     |            4 | minDistance | Floating point, radians                                     |
     |            4 | maxDistance | Floating point, radians                                     |
     |            4 | numBins     | Number of distance bins                                     |
     | 2*2*numPairs | pairs       | Pairs, sorted by distance between the stars in the pair.    |
     |              |             | Simply stores the BSC index of the first star immediately   |
     |              |             | followed by the BSC index of the other star.                |
     |4*(numBins+1) | bins        | The `i'th bin (starting from zero) stores how many pairs of |
     |              |             | stars have a distance less than or equal to:                |
     |              |             | minDistance+i*(maxDistance-minDistance)/numBins             |
 */
unsigned char *BuildKVectorDatabase(const Catalog &catalog, long *length,
                           float minDistance, float maxDistance, long numBins) {
    std::vector<int32_t> kVector(numBins+1); // numBins = length, all elements zero
    std::vector<KVectorPair> pairs;
    float binWidth = (maxDistance - minDistance) / numBins;

    long dummyLength;
    if (length == NULL) {
        length = &dummyLength;
    }

    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t k = i+1; k < (int16_t)catalog.size(); k++) {

            KVectorPair pair = { i, k, AngleUnit(catalog[i].spatial, catalog[k].spatial) };
            assert(isfinite(pair.distance));
            assert(pair.distance >= 0);
            assert(pair.distance <= M_PI);

            if (pair.distance >= minDistance && pair.distance <= maxDistance) {
                // we'll sort later
                pairs.push_back(pair);
            }
        }
    }

    // sort pairs in increasing order.
    std::sort(pairs.begin(), pairs.end(), CompareKVectorPairs);

    // generate the k-vector part
    // Idea: When we find the first star that's across any bin boundary, we want to update all the newly sealed bins
    long lastBin = 0; // first bin the last star belonged to
    for (int32_t i = 0; i < (int32_t)pairs.size(); i++) {
        long thisBin = (long)ceil((pairs[i].distance - minDistance) / binWidth); // first bin we belong to
        assert(thisBin >= 0);
        assert(thisBin <= numBins); // thisBin == numBins is acceptable since kvector length == numBins + 1
        // if thisBin > lastBin, then no more stars can be added to those bins between thisBin and lastBin, so set them.
        for (long bin = lastBin; bin < thisBin; bin++) {
            kVector[bin]=i; // our index is the number of pairs with shorter distance
        }
        lastBin = thisBin;
    }
    for (long bin = lastBin; bin <= numBins; bin++) {
        kVector[bin] = pairs.size();
    }

    // verify kvector
    int32_t lastBinVal = -1;
    for (const int32_t &bin : kVector) {
        assert(bin >= lastBinVal);
        lastBinVal = bin;
    }

    // TODO: determine the correct length, copy the correct vectors and numbers into result
    *length = 4+4+4+4+4 + 2*sizeof(int16_t)*pairs.size() + sizeof(int32_t)*(numBins+1);
    
    unsigned char *result = new unsigned char[*length];
    int32_t *resultMagicValue = (int32_t *)result;
    int32_t *resultNumPairs = (int32_t *)resultMagicValue + 1;
    float *resultMinDistance = (float *)(resultNumPairs + 1);
    float *resultMaxDistance = (float *)(resultMinDistance + 1);
    int32_t *resultNumBins = (int32_t *)(resultMaxDistance + 1);
    int16_t *resultPairs = (int16_t *)(resultNumBins + 1);
    int32_t *resultKVector = (int32_t *)(resultPairs + pairs.size()*2);

    *resultNumPairs = pairs.size();
    *resultMinDistance = minDistance;
    *resultMaxDistance = maxDistance;
    *resultNumBins = numBins;

    for (const KVectorPair &pair : pairs) {
        resultPairs[0] = pair.index1;
        resultPairs[1] = pair.index2;
        resultPairs += 2;
    }
    assert((void *)resultPairs == (void *)resultKVector);

    // you could probably do this with memcpy instead, but the explicit loop is necessary for endian
    // concerns? TODO endianness
    for (const int32_t &bin : kVector) {
        *resultKVector = bin;
        resultKVector++;
    }
    assert((unsigned char *)resultKVector - result == *length);

    *resultMagicValue = kKVectorMagicNumber;
    return result;
}

KVectorDatabase::KVectorDatabase(const unsigned char *databaseBytes) {
    // TODO: errors?
    int32_t *bytesMagicValue = (int32_t *)databaseBytes;
    int32_t *bytesNumPairs = (int32_t *)bytesMagicValue + 1;
    float *bytesMinDistance = (float *)(bytesNumPairs + 1);
    float *bytesMaxDistance = (float *)(bytesMinDistance + 1);
    int32_t *bytesNumBins = (int32_t *)(bytesMaxDistance + 1);

    assert(*bytesMagicValue == kKVectorMagicNumber);
    numPairs = *bytesNumPairs;
    minDistance = *bytesMinDistance;
    assert(minDistance > 0.0f);
    maxDistance = *bytesMaxDistance;
    assert(maxDistance > minDistance);
    numBins = *bytesNumBins;
    binWidth = (maxDistance - minDistance) / numBins;

    pairs = (int16_t *)(bytesNumBins + 1);
    bins = (int32_t *)(pairs + 2*numPairs);
}

float Clamp(float num, float low, float high) {
    return num < low ? low : num > high ? high : num;
}

int16_t *KVectorDatabase::FindPossibleStarPairsApprox(
    float minQueryDistance, float maxQueryDistance, long *numReturnedPairs) const {

    assert(maxQueryDistance > minQueryDistance);
    if (maxQueryDistance >= maxDistance) {
        maxQueryDistance = maxDistance - 0.00001; // TODO: better way to avoid hitting the bottom bin
    }
    if (minQueryDistance <= minDistance) {
        minQueryDistance = minDistance + 0.00001;
    }
    if (minQueryDistance > maxDistance || maxQueryDistance < minDistance) {
        *numReturnedPairs = 0;
        return NULL;
    }
    long lowerBin = BinForDistance(minQueryDistance);
    long upperBin = BinForDistance(maxQueryDistance);
    assert(upperBin >= lowerBin);
    assert(upperBin <= numBins);
    // bins[lowerBin-1]=number of pairs <= r < query distance, so it is the index of the
    // first possible item that might be equal to the query distance
    int lowerPair = bins[lowerBin-1];
    if (lowerPair >= numPairs) {
        // all pairs have distance less than queried
        return NULL;
    }
    // bins[upperBin]=number of pairs <= r >= query distance
    int upperPair = bins[upperBin] - 1;
    *numReturnedPairs = upperPair - lowerPair + 1;
    assert(*numReturnedPairs >= 0);
    return &pairs[lowerPair * 2];
}

long KVectorDatabase::BinForDistance(float distance) const {
    long result = (long)ceil((distance - minDistance) / binWidth);
    assert(result >= 0);
    assert(result <= numBins);
    return result;
}

long KVectorDatabase::NumPairs() const {
    return numPairs;
}

std::vector<float> KVectorDatabase::StarDistances(int16_t star, const Catalog &catalog) const {
    std::vector<float> result;
    for (int i = 0; i < numPairs; i++) {
        if (pairs[i*2] == star || pairs[i*2+1] == star) {
            result.push_back(AngleUnit(catalog[pairs[i*2]].spatial, catalog[pairs[i*2+1]].spatial));
        }
    }
    return result;
}

}

// TODO: after creating the database, print more statistics, such as average number of pairs per
// star, stars per bin, space distribution between array and index.
