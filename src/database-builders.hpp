#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <inttypes.h>
#include <vector>

#include "catalog-generic.hpp"

namespace lost {

/**
 K-vector database layout
     | size (bytes) | name        | description                                                 |
     |--------------+-------------+-------------------------------------------------------------|
     |            4 | numPairs    | Number of star pairs that are in the min/max distance range |
     |            4 | minDistance | Millionths of a degree                                      |
     |            4 | maxDistance | Millionths of a degree                                      |
     |            4 | numBins     | Number of distance bins                                     |
     | 2*2*numPairs | pairs       | Pairs, sorted by distance between the stars in the pair.    |
     |              |             | Simply stores the BSC index of the first star immediately   |
     |              |             | followed by the BSC index of the other star.                |
     |    4*numBins | bins        | The `i'th bin (starting from zero) stores how many pairs of |
     |              |             | stars have a distance less than or equal to:                |
     |              |             | minDistance+(maxDistance-minDistance)*i/numBins             |
 */

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
