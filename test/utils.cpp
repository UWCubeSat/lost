#include <algorithm>

#include "databases.hpp"

namespace lost {

unsigned char *BuildPairDistanceKVectorDatabase(
    const Catalog &catalog, long *length, float minDistance, float maxDistance, long numBins) {

    long dummyLength;
    if (length == NULL)
        length = &dummyLength;

    *length = SerializeLengthPairDistanceKVector(catalog, minDistance, maxDistance, numBins);
    unsigned char *result = new unsigned char[*length];
    SerializePairDistanceKVector(catalog, minDistance, maxDistance, numBins, result);
    return result;
}

bool AreStarIdentifiersEquivalent(const StarIdentifiers &ids1, const StarIdentifiers &ids2) {
    if (ids1.size() != ids2.size()) {
        return false;
    }

    for (const StarIdentifier &id1 : ids1) {
        if (std::find(ids2.cbegin(), ids2.cend(), id1) == ids2.cend()) {
            return false;
        }
    }
    return true;
}

}
