#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <stdlib.h>
#include <inttypes.h>
#include <vector>

#include "star-utils.hpp"

namespace lost {

const int32_t kCatalogMagicValue = 0xF9A283BC;

struct KVectorQuad {
    int16_t* stars;
    float* parameters;
};

/**
 * A data structure enabling constant-time range queries into fixed numerica data.
 * 
 * @note Not an instantiable database on its own -- used in other databases
 */
// TODO: QueryConservative, QueryExact, QueryTrapezoidal?
class KVectorIndex {
public:
    explicit KVectorIndex(const unsigned char *);

    long QueryLiberal(float minQueryDistance, float maxQueryDistance, long *upperIndex) const;

    /// The number of data points in the data referred to by the kvector
    long NumValues() const { return numValues; };
    long NumBins() const { return numBins; };
    /// Upper bound on elements
    float Max() const { return max; };
    // Lower bound on elements
    float Min() const { return min; };
private:
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
std::vector<KVectorQuad> CatalogToQuadDistances(const Catalog &catalog, float minDistance, float maxDistance);
void SerializeKVectorND(const Catalog &catalog, std::vector<KVectorQuad> quads, float minDistance, float maxDistance, long numBins, unsigned char *buffer);
long SerializeLengthQuadStarKVectorND(int numEntries, int bins);
float StarParameterA(float centralToOne, float centralToTwo, float oneToTwo);
float StarParameterB(float centralToOne, float centralToTwo, float centralToThree);
/**
 * A database storing distances between pairs of stars.
 * Supports fast range queries to find all pairs of stars separated by approximately a certain distance.
 * @warning Sensitive to uncalibrated camera parameters
 */
class PairDistanceKVectorDatabase {
public:
    explicit PairDistanceKVectorDatabase(const unsigned char *databaseBytes);

    const int16_t *FindPairsLiberal(float min, float max, const int16_t **end) const;
    std::vector<float> StarDistances(int16_t star, const Catalog &) const;

    /// Upper bound on stored star pair distances
    float MaxDistance() const { return index.Max(); };
    /// Lower bound on stored star pair distances
    float MinDistance() const { return index.Min(); };
    /// Number of star pairs stored in the database
    long NumPairs() const { return index.NumValues(); }

    /// Magic value to use when storing inside a MultiDatabase
    static const int32_t kMagicValue = 0x2536f009;
private:
    KVectorIndex index;
    // TODO: endianness
    const int16_t *pairs;
};

/**
 * A data structure enabling range searches across multidimensional data
 * 
 * @note Hard Coded for a dimension of 4, possible to generalize later
*/
class KVectorND {
public:
    explicit KVectorND(const unsigned char *);
    ~KVectorND();

    std::vector<int16_t *> RangeSearch(const float *maxParameter, const float *minParameter) const;
    int16_t *GetEntry(const long index) const;

    /// Upper bound parameter corresponding to axis i
    float MaxParameter(int i) const { return max[i]; };
    /// Lower bound parameter corresponding to axis i
    float MinParameter(int i) const { return min[i]; };
    /// Exact number of stored pairs
    long NumQuads() const { return numValues; };

    /// Magic value to use when storing inside a MultiDatabase
    static const int32_t kMagicValue = 0xff23ab79;

    /// The number of data points in the data referred to by the kvector
    long NumValues() const { return numValues; };
    // Number of bins on each axis
    long NumBins() const { return numBins; };

    
private:
    // Total Entries
    long numValues;
    // Maximum and Minimum Parameter Values on each Axis
    float* min;
    float* max;
    // Range of bins on each axis
    float* binWidth;
    // Bins on each axis
    long numBins;
    // Indexes for each bin
    const int32_t* bins;

    const int16_t* quads;
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

//     void FindTriplesLiberal(float min, float max, long **begin, long **end) const;
// private:
//     KVectorIndex index;
//     int16_t *triples;
// };

/// maximum number of databases in a MultiDatabase
const int kMultiDatabaseMaxDatabases = 64;
/// The size of the table of contents in a multidatabase (stores subdatabase locations)
const long kMultiDatabaseTocLength = 8*kMultiDatabaseMaxDatabases;

/**
 * A database that contains multiple databases
 * This is almost always the database that is actually passed to star-id algorithms in the real world, since you'll want to store at least the catalog plus one specific database.
 * Multi-databases are essentially a map from "magic values" to database buffers.
 */
class MultiDatabase {
public:
    /// Create a multidatabase from a serialized multidatabase.
    explicit MultiDatabase(const unsigned char *buffer) : buffer(buffer) { };
    const unsigned char *SubDatabasePointer(int32_t magicValue) const;
private:
    const unsigned char *buffer;
};

/// Class for easily creating a MultiDatabase
class MultiDatabaseBuilder {
public:
    MultiDatabaseBuilder()
        : buffer((unsigned char *)calloc(1, kMultiDatabaseTocLength)), bulkLength(0) { };
    ~MultiDatabaseBuilder();

    unsigned char *AddSubDatabase(int32_t magicValue, long length);

    /// When done adding databases, use this to get the buffer you should write to disk.
    unsigned char *Buffer() { return buffer; };
    /// The length of the buffer returned by Buffer
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
