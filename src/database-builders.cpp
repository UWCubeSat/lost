#include "database-builders.hpp"
#include "star-utils.hpp"

#include <vector>

namespace lost {


/**
 K-vector database layout
     | size (bytes) | name        | description                                                 |
     |--------------+-------------+-------------------------------------------------------------|
     |            4 | magic       | Magic number 0x4253f009 (ensure database validity)          |
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
                           float minDistance, float maxDistance, long numBins) {
    *length = 0;

    std::vector<long> kVector(numBins); // numBins = length, all elements zero
    std::vector<std::pair<int, int>> pairs;
    for (const CatalogStar &one : catalog) {
        for (const CatalogStar &two : catalog) {
            // TODO: great circle (haversine) distance between the two stars, check if within
            // minDistance-maxDistance, if so add to `pairs` and increment the appropriate bins of
            // kVector.
        }
    }

    // TODO: determine the correct length, copy the correct vectors and numbers into result
    *length = 6969;
    unsigned char *result = new unsigned char[*length];

    return result;
}

}
