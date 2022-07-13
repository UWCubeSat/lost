#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <stdlib.h>
#include <inttypes.h>
#include <vector>

#include "star-utils.hpp"

namespace lost {

const int32_t kCatalogMagicValue = 0xF9A283BC;

/**
 * @brief
 * @details
 * @note Not an instantiable database on its own -- used in other databases
 * @todo QueryConservative, and QueryTrapezoidal which interpolates linearly between endpoints
 */
class KVectorIndex {
public:
    KVectorIndex(const unsigned char *);

    long QueryLiberal(float minQueryDistance, float maxQueryDistance, long *upperIndex) const;

    /**
     * @brief
     * @return
     */
    long NumValues() const { return numValues; };

    /**
     * @brief
     * @return
     */
    long NumBins() const { return numBins; };

    /**
     * @brief
     * @return
     */
    float Max() const { return max; };

    /**
     * @brief
     * @return
     */
    float Min() const { return min; };
private:
    // return the lowest-indexed bin that contains the number of pairs with distance <= dist
    long BinFor(float dist) const;

    long numValues;
    float min;
    float max;
    float binWidth;
    long numBins;
    const int32_t *bins;
};

long SerializeLengthPairDistanceKVector(const Catalog &, float minDistance, float maxDistance, long numBins);
void SerializePairDistanceKVector(const Catalog &, float minDistance, float maxDistance, long numBins, unsigned char *buffer);

/**
 * @brief Stores angular distance between pairs of stars
 * @details
 * @warning Sensitive to uncalibrated camera parameters
 * @todo Trapezoidal interpolation
 */
class PairDistanceKVectorDatabase {
public:
    PairDistanceKVectorDatabase(const unsigned char *databaseBytes);

    const int16_t *FindPairsLiberal(float min, float max, const int16_t **end) const;

    std::vector<float> StarDistances(int16_t star, const Catalog &) const;

    /**
     * @brief
     * @return
     */
    float MaxDistance() const { return index.Max(); };

    /**
     * @brief
     * @return
     */
    float MinDistance() const { return index.Min(); };

    /**
     * @brief
     * @return
     */
    long NumPairs() const;

    /// @brief
    const static int32_t kMagicValue = 0x2536f009;
private:
    KVectorIndex index;
    // TODO: endianness
    const int16_t *pairs;
};

// stores "inner angles" between star triples. Unsensitive to first-order error in basic camera
// parameters (eg, wrong FOV or principal point), can be sensitive to second-order errors (eg,
// camera distortion, which may cause the effective FOV or principal point to be different in
// different parts of the image). Used for Mortari's Non-Dimensional Star-ID
class TripleInnerKVectorDatabase {
public:
    TripleInnerKVectorDatabase(const unsigned char *databaseBytes);

    // return at least all the triples with inner angle in the given range. The numReturnedTriples*3
    // ints from the returned pointer are valid to read.
    void FindTriplesLiberal(float min, float max, long **begin, long **end) const;
    // TODO: trapezoidal interpolation
private:
    KVectorIndex index;
    int16_t *triples;
};

// maximum number of databases in a MultiDatabase
const int kMultiDatabaseMaxDatabases = 64;
const long kMultiDatabaseTocLength = 8*kMultiDatabaseMaxDatabases;

// ,
/**
 * @brief Represents a database that contains multiple databases
 * @details This is almost always what will be used in the real world,
 * since you'll want to store at least the catalog plus one specific database.
 */
class MultiDatabase {
public:
    /**
     * @brief
     * @param buffer
     */
    MultiDatabase(const unsigned char *buffer) : buffer(buffer) { };
    const unsigned char *SubDatabasePointer(int32_t magicValue) const;
private:
    const unsigned char *buffer;
};

/**
 * @brief
 * @details
 */
class MultiDatabaseBuilder {
public:
    /**
     * @brief
     * @note the () after new ensures it's zero-initialized
     */
    MultiDatabaseBuilder()
        : buffer((unsigned char *)calloc(1, kMultiDatabaseTocLength)), bulkLength(0) { };
    ~MultiDatabaseBuilder();
    unsigned char *AddSubDatabase(int32_t magicValue, long length);

    /**
     * @brief
     * @return
     */
    unsigned char *Buffer() { return buffer; };

    /**
     * @brief
     * @return
     */
    long BufferLength() { return kMultiDatabaseTocLength+bulkLength; };
private:
    // Throughout LOST, most dynamic memory is managed with `new` and `delete` to make it easier to
    // use unique pointers. Here, however, we use realloc, so C-style memory management.
    unsigned char *buffer;
    // how many bytes are presently allocated for databases (excluding map)
    long bulkLength;
};

}

#endif
