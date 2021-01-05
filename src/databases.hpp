#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <inttypes.h>
#include <vector>

#include "star-utils.hpp"

namespace lost {

unsigned char *BuildKVectorDatabase(const Catalog &catalog, long *length,
                                    float minDistance, float maxDistance, long numBins);

class KVectorDatabase {
public:
    KVectorDatabase(unsigned char *databaseBytes);

    // the "Approx" functions assume star distances are uniformly distributed through the bins. The
    // "Exact" functions actually check the distances between the returned stars to make sure it's correct.

    std::vector<int16_t> FindPossibleStarsApprox(
        float minDistance, float maxDistance) const;
    std::vector<int16_t> FindPossibleStarsExact(
        float minDistance, float maxDistance, const Catalog &) const;
    int16_t *FindPossibleStarPairsApprox(
        float minDistance, float maxDistance, int *numReturnedPairs) const;
    int16_t *FindPossibleStarPairsExact(
        float minDistance, float maxDistance, const Catalog &, int *numReturnedPairs) const;

    int NumStars() const;
private:
    void BinBounds(int bin, float *min, float *max) const;
    int BinForDistance(float dist) const;

    long numPairs;
    float minDistance;
    float maxDistance;
    long numBins;
    // TODO: endianness
    int16_t *pairs;
    int32_t *bins;
};

}

#endif
