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
     |              |             | minDistance+(maxDistance-minDistance)*i/numBins             |
 */
unsigned char *BuildKVectorDatabase(const Catalog &catalog, long *length,
                           float minDistance, float maxDistance, long numBins) {
    std::vector<int32_t> kVector(numBins+1); // numBins = length, all elements zero
    std::vector<KVectorPair> pairs;
    float binWidth = (maxDistance - minDistance) / numBins;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t k = i+1; k < (int16_t)catalog.size(); k++) {

            KVectorPair pair = { i, k, GreatCircleDistance(catalog[i].raj2000, catalog[i].dej2000,
                                                           catalog[k].raj2000, catalog[k].dej2000) };
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
    long lastBin = 0;
    for (int32_t i = 0; i < (int32_t)pairs.size(); i++) {
        long thisBin = (long)floor((pairs[i].distance - minDistance) / binWidth);
        if (thisBin == numBins) {
            // rounding error? TODO
            thisBin--;
        }
        assert(thisBin >= 0);
        assert(thisBin < numBins);
        for (long bin = lastBin; bin <= thisBin; bin++) {
            kVector[bin+1]=i;
        }
        lastBin = thisBin;
    }
    for (long bin = lastBin + 1; bin < numBins; bin++) {
        kVector[bin+1] = kVector[lastBin+1];
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

    pairs = (int16_t *)(bytesNumBins + 1);
    bins = (int32_t *)(pairs + 2*numPairs);
}

float Clamp(float num, float low, float high) {
    return num < low ? low : num > high ? high : num;
}

int16_t *KVectorDatabase::FindPossibleStarPairsApprox(
    float minQueryDistance, float maxQueryDistance, long *numReturnedPairs) const {

    assert(maxQueryDistance > minQueryDistance);
    if (minQueryDistance < minDistance || minQueryDistance > maxDistance ||
        maxQueryDistance < minDistance || maxQueryDistance > maxDistance) {
        *numReturnedPairs = 0;
        return NULL;
    }
    //tbr v
    long lowerBin = BinForDistance(minQueryDistance);
    long upperBin = BinForDistance(maxQueryDistance);
    assert(upperBin >= lowerBin);
    assert(upperBin < numBins);
    int lowerPair = bins[lowerBin];
    assert(lowerPair < numPairs);
    int upperPair = bins[upperBin+1];
    *numReturnedPairs = upperPair - lowerPair;
    assert(*numReturnedPairs >= 0);
    return &pairs[lowerPair * 2];
}

int16_t *KVectorDatabase::FindPossibleStarPairsExact(
    float minDistance, float maxDistance, const Catalog &catalog, int *numReturnedPairs) const {

    
}

void KVectorDatabase::BinBounds(int bin, float *min, float *max) const {
    assert(bin >= 0 && bin < numBins);

    if (min != NULL) {
        *min = (maxDistance-minDistance)/numBins*bin;
    }

    if (max != NULL) {
        *max = (maxDistance-minDistance)/numBins*(bin+1);
    }
}

long KVectorDatabase::BinForDistance(float distance) const {
    float binWidth = (maxDistance - minDistance) / numBins;
    long result = (long)floor((distance - minDistance) / binWidth);
    assert (result <= numBins);
    if (result == numBins) {
        return numBins - 1;
    } else {
        return result;
    }
}

}
