#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <inttypes.h>
#include <vector>

#include "star-utils.hpp"

namespace lost {

unsigned char *BuildKVectorDatabase(const Catalog &catalog, long *length,
                                    float minDistance, float maxDistance, long numBins);

class KVectorDatabase {
    std::vector<int> StarsWithDistance(long minDistance, long maxDistance);
private:
    long numPairs;
    long minDistance;
    long maxDistance;
    long numBins;
    // TODO: endianness
    int16_t *pairs;
    int32_t *bins;
};

}

#endif
