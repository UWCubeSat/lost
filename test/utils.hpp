#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "databases.hpp"

namespace lost {

unsigned char *BuildPairDistanceKVectorDatabase(const Catalog &catalog, long *length, float minDistance, float maxDistance, long numBins);
/// simple O(n^2) check
bool AreStarIdentifiersEquivalent(const StarIdentifiers &, const StarIdentifiers &);

}

#endif
