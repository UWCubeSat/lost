#include "database-builders.hpp"
#include "star-utils.hpp"

#include <vector>

namespace lost {

unsigned char *BuildKVectorDatabase(const Catalog &catalog, long *length,
                           float minDistance, float maxDistance, long numBins) {
    long numPairs = 0;
    for (const CatalogStar &one : catalog) {
        for (const CatalogStar &two : catalog) {
            // TODO: great circle distance between the two stars, check if within
            // minDistance-maxDistance, if so add to vector.
        }
    }

    *length = 6969;
    return new unsigned char[*length];
}

}
