#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <inttypes.h>
#include <vector>

#include "star-utils.hpp"

namespace lost {

const int32_t kCatalogMagicValue = 0xF9A283BC;

// not an instantiable database on its own -- used in other databases
class KVectorIndex {
public:
    // construct from serialized
    KVectorIndex(const unsigned char *);

    // finds at least all the entries containing the given range. Returns the index (starting from
    // zero) of the first value matching the query
    long QueryLiberal(float minQueryDistance, float maxQueryDistance, long *numReturned) const;
    // TODO: QueryConservative, and QueryTrapezoidal which interpolates linearly between endpoints

    long NumValues() const { return numValues; };
    long NumBins() const { return numBins; };
    float Max() const { return max; };
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

unsigned char *BuildPairDistanceKVectorDatabase(const Catalog &catalog, long *length,
                                    float minDistance, float maxDistance, long numBins);

long SerializeLengthPairDistanceKVector(const Catalog &, float minDistance, float maxDistance, long numBins);
void SerializePairDistanceKVector(const Catalog &, float minDistance, float maxDistance, long numBins, unsigned char *buffer);

// stores angular distance between pairs of stars. Sensitive to uncalibrated camera parameters
class PairDistanceKVectorDatabase {
public:
    PairDistanceKVectorDatabase(const unsigned char *databaseBytes);

    // return at least all the stars between min and max
    const int16_t *FindPairsLiberal(
        float min, float max, long *numReturnedPairs) const;
    // TODO: trapezoidal interpolation

    // for debugging purposes. Return the distances from the given star to each other star it's
    // paired with in the database.
    std::vector<float> StarDistances(int16_t star, const Catalog &) const;

    float MaxDistance() const { return index.Max(); };
    float MinDistance() const { return index.Min(); };

    long NumPairs() const;
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
    int16_t *FindTriplesLiberal(float min, float max, long *numReturnedTriples) const;
    // TODO: trapezoidal interpolation
private:
    KVectorIndex index;
    int16_t *triples;
};

// maximum number of databases in a MultiDatabase
const int kMultiDatabaseMaxDatabases = 64;
const long kMultiDatabaseTocLength = 8*kMultiDatabaseMaxDatabases;

// represents a database that contains multiple databases, which is almost always what will be used
// in the real world, since you'll want to store at least the catalog plus one specific database.
class MultiDatabase {
public:
    MultiDatabase(const unsigned char *buffer) : buffer(buffer) { };
    // return a pointer to the start of the database type indicated by the magic value, if such a
    // sub-database is present in the database. Return null if not found.
    const unsigned char *SubDatabasePointer(int32_t magicValue) const;
private:
    const unsigned char *buffer;
};

class MultiDatabaseBuilder {
public:
    MultiDatabaseBuilder()
        // the () after new ensures it's zero-initialized
        : buffer(new unsigned char[kMultiDatabaseTocLength]()), bulkLength(0) { };
    ~MultiDatabaseBuilder();
    // return pointer to the start of the space allocated for said database. Return null if full.
    unsigned char *AddSubDatabase(int32_t magicValue, long length);
    unsigned char *Buffer() { return buffer; };
    long BufferLength() { return kMultiDatabaseTocLength+bulkLength; };
private:
    unsigned char *buffer;
    // how many bytes are presently allocated for databases (excluding map)
    long bulkLength;
};

}

#endif
