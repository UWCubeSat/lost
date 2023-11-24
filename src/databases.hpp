#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <stdlib.h>
#include <inttypes.h>
#include <vector>

#include "star-utils.hpp"
#include "serialize-helpers.hpp"

namespace lost {

const int32_t kCatalogMagicValue = 0xF9A283BC;

/**
 * A data structure enabling constant-time range queries into fixed numerical data.
 * 
 * @note Not an instantiable database on its own -- used in other databases
 */
// TODO: QueryConservative, QueryExact, QueryTrapezoidal?
class KVectorIndex {
public:
    explicit KVectorIndex(DeserializeContext *des);

    long QueryLiberal(decimal minQueryDistance, decimal maxQueryDistance, long *upperIndex) const;

    /// The number of data points in the data referred to by the kvector
    long NumValues() const { return numValues; };
    long NumBins() const { return numBins; };
    /// Upper bound on elements
    decimal Max() const { return max; };
    // Lower bound on elements
    decimal Min() const { return min; };
private:
    long BinFor(decimal dist) const;

    long numValues;
    decimal min;
    decimal max;
    decimal binWidth;
    long numBins;
    const int32_t *bins;
};

void SerializePairDistanceKVector(SerializeContext *, const Catalog &, decimal minDistance, decimal maxDistance, long numBins);

/**
 * A database storing distances between pairs of stars.
 * Supports fast range queries to find all pairs of stars separated by approximately a certain distance.
 * @warning Sensitive to uncalibrated camera parameters
 */
class PairDistanceKVectorDatabase {
public:
    explicit PairDistanceKVectorDatabase(DeserializeContext *des);

    const int16_t *FindPairsLiberal(decimal min, decimal max, const int16_t **end) const;
    const int16_t *FindPairsExact(const Catalog &, decimal min, decimal max, const int16_t **end) const;
    std::vector<decimal> StarDistances(int16_t star, const Catalog &) const;

    /// Upper bound on stored star pair distances
    decimal MaxDistance() const { return index.Max(); };
    /// Lower bound on stored star pair distances
    decimal MinDistance() const { return index.Min(); };
    /// Exact number of stored pairs
    long NumPairs() const;

    /// Magic value to use when storing inside a MultiDatabase
    static const int32_t kMagicValue; // 0x2536f009
    // apparently you're supposed to not actually put the value of the static variables here, but
    // rather in a cpp implementation file.
private:
    KVectorIndex index;
    // TODO: endianness
    const int16_t *pairs;
};

// /**
//  * @brief Stores "inner angles" between star triples
//  * @details Unsensitive to first-order error in basic camera
//  * parameters (eg, wrong FOV or principal point), can be sensitive to second-order errors (eg,
//  * camera distortion, which may cause the effective FOV or principal point to be different in
//  * different parts of the image). Used for Mortari's Non-Dimensional Star-ID
//  */
// class TripleInnerKVectorDatabase {
// public:
//     explicit TripleInnerKVectorDatabase(const unsigned char *databaseBytes);

//     void FindTriplesLiberal(decimal min, decimal max, long **begin, long **end) const;
// private:
//     KVectorIndex index;
//     int16_t *triples;
// };

/**
 * A database that contains multiple databases
 * This is almost always the database that is actually passed to star-id algorithms in the real world, since you'll want to store at least the catalog plus one specific database.
 * Multi-databases are essentially a map from "magic values" to database buffers.
 */

#define MULTI_DB_IS_DOUBLE 0x0001
#define MULTI_DB_IS_FLOAT 0x0000

class MultiDatabase {
public:
    /// Create a multidatabase from a serialized multidatabase.
    explicit MultiDatabase(const unsigned char *buffer) : buffer(buffer) { };
    const unsigned char *SubDatabasePointer(int32_t magicValue) const;
private:
    const unsigned char *buffer;
};

class MultiDatabaseEntry {
public:
    MultiDatabaseEntry(int32_t magicValue, std::vector<unsigned char> bytes) // I wonder if making `bytes` a reference would avoid making two copies, or maybe it would be worse by preventing copy-elision
        : magicValue(magicValue), bytes(bytes) { }

    int32_t magicValue;
    uint32_t flags;
    std::vector<unsigned char> bytes;
};

typedef std::vector<MultiDatabaseEntry> MultiDatabaseDescriptor;

void SerializeMultiDatabase(SerializeContext *, const MultiDatabaseDescriptor &dbs, uint32_t flags);

}

#endif
