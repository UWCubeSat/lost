#include "databases.hpp"
#include "attitude-utils.hpp"
#include "star-utils.hpp"

#include <vector>
#include <algorithm>
#include <assert.h>
#include <math.h>

#define kKVectorMagicNumber 0x4253f009

namespace lost {

static float GreatCircleDistance(const CatalogStar &one, const CatalogStar &two) {
    return 2.0*asin(sqrt(pow(sin(abs(one.dej2000-two.dej2000)/2.0), 2.0)
                         + cos(one.dej2000)*cos(two.dej2000)*pow(sin(abs(one.raj2000-two.raj2000)/2.0), 2.0)));
}

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
     |    4*numBins | bins        | The `i'th bin (starting from zero) stores how many pairs of |
     |              |             | stars have a distance less than or equal to:                |
     |              |             | minDistance+(maxDistance-minDistance)*i/numBins             |
 */
unsigned char *BuildKVectorDatabase(const Catalog &catalog, long *length,
                           float minDistance, float maxDistance, long numBins) {
    std::vector<int32_t> kVector(numBins); // numBins = length, all elements zero
    std::vector<KVectorPair> pairs;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t k = i+1; k < (int16_t)catalog.size(); k++) {

            KVectorPair pair = { i, k, GreatCircleDistance(catalog[i], catalog[k]) };
            assert(isfinite(pair.distance));
            assert(pair.distance >= 0);
            assert(pair.distance <= M_PI);

            if (pair.distance >= minDistance && pair.distance <= maxDistance) {
                // we'll sort later
                pairs.push_back(pair);

                for (int j = 0;
                     j < numBins && pair.distance <= minDistance+(maxDistance-minDistance)*j/numBins;
                     j++) {
                    kVector[j]++;
                }
            }
        }
    }

    // sort pairs in increasing order.
    std::sort(pairs.begin(), pairs.end(), CompareKVectorPairs);

    // TODO: determine the correct length, copy the correct vectors and numbers into result
    *length = 4+4+4+4+4 + 2*sizeof(int16_t)*pairs.size() + sizeof(int32_t)*numBins;
    
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

KVectorDatabase::KVectorDatabase(unsigned char *databaseBytes) {
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

int16_t *KVectorDatabase::FindPossibleStarPairsApprox(float minQueryDistance, float maxQueryDistance, int *numReturnedPairs) const {
    //minDistance is going the starting point
    float binWidth = (maxDistance - minDistance) / numBins;
    //tbr v
    int lowerBound = bins[(int16_t)floor((minQueryDistance - minDistance) / binWidth)];
    int upperBound = bins[(int16_t)floor((maxQueryDistance - maxDistance) / binWidth)];
    *numReturnedPairs = upperBound - lowerBound;
    return &pairs[lowerBound * 2];
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

int KVectorDatabase::BinForDistance(float distance) const {
    // TODO: Be careful about rounding!
}

}
