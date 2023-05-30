#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <inttypes.h>
#include <stdlib.h>

#include <utility>
#include <vector>

#include "star-utils.hpp"

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

long SerializeLengthPairDistanceKVector(const Catalog &, float minDistance, float maxDistance,
                                        long numBins);
void SerializePairDistanceKVector(const Catalog &, float minDistance, float maxDistance,
                                  long numBins, unsigned char *buffer);

/**
 * A database storing distances between pairs of stars.
 * Supports fast range queries to find all pairs of stars separated by approximately a certain
 * distance.
 * @warning Sensitive to uncalibrated camera parameters
 */
class PairDistanceKVectorDatabase {
   public:
    // TODO: databaseBytes just = buffer?
    explicit PairDistanceKVectorDatabase(const unsigned char *databaseBytes);

    const int16_t *FindPairsLiberal(float min, float max, const int16_t **end) const;
    const int16_t *FindPairsExact(const Catalog &, float min, float max, const int16_t **end) const;
    std::vector<float> StarDistances(int16_t star, const Catalog &) const;

    /// Upper bound on stored star pair distances
    float MaxDistance() const { return index.Max(); };
    /// Lower bound on stored star pair distances
    float MinDistance() const { return index.Min(); };
    /// Exact number of stored pairs
    long NumPairs() const;

    /// Magic value to use when storing inside a MultiDatabase
    static const int32_t kMagicValue = 0x2536f009;

   private:
    KVectorIndex index;
    // TODO: endianness
    const int16_t *pairs;
};

/*
Pre-processing for Tetra star-id algorithm
Return:
(a) List of stars we can use for Tetra
(b) Subset of (a) that we use to generate Tetra star patterns
*/
std::pair<std::vector<uint16_t>, std::vector<uint16_t>> TetraPreparePattCat(const Catalog &,
                                                                            const float maxFovDeg);

// long SerializeTetraDatabase(const Catalog &, float maxFov, unsigned char *buffer,
//                             const std::vector<uint16_t> &,
//                             const std::vector<uint16_t> &, bool ser);

long SerializeLengthTetraDatabase(const Catalog &, float maxFov,
                            const std::vector<uint16_t> &,
                            const std::vector<uint16_t> &);

void SerializeTetraDatabase(const Catalog &, float maxFov, unsigned char *buffer,
                            const std::vector<uint16_t> &,
                            const std::vector<uint16_t> &);

typedef std::vector<int> Pattern;

/*
Layout:

/////////////////// Header //////////////////////////////
- Max FOV (float)
- Number of patterns in pattern catalog (int)
////////////////////////////////////////////////////////
- All patterns (number of patterns * pattSize=4 * sizeof(uint16_t))
- List of Catalog indices to use for Tetra star ID algo
/////////////////////////////////////////////////////////

*/
class TetraDatabase {
   public:
    explicit TetraDatabase(const unsigned char *buffer);

    // Get max angle (in degrees) allowed between stars in the same pattern
    float MaxAngle() const;

    /// Number of rows in pattern catalog
    // With load factor of 0.5, size = number of patterns * 2
    int PattCatSize() const;

    // Get the 4-tuple pattern at row=index, 0-based
    Pattern GetPattern(int index) const;

    /*
    Returns a list of rows (4-tuples) from the Pattern Catalog
    Start from index and does quadratic probing
    */
    std::vector<Pattern> GetPatternMatches(int index) const;

    uint16_t GetTrueCatInd(int tetraIndex) const;
    // TODO: should probably have a field describing number of indices for future updates to db

    /// Magic value to use when storing inside a MultiDatabase
    static const int32_t kMagicValue = 0xDEADBEEF;
    static const int headerSize = sizeof(float) + sizeof(int32_t);

   private:
    const unsigned char *buffer_;
    float maxAngle_;
    uint32_t catalogSize_;
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
const long kMultiDatabaseTocLength = 8 * kMultiDatabaseMaxDatabases;

/**
 * A database that contains multiple databases
 * This is almost always the database that is actually passed to star-id algorithms in the real
 * world, since you'll want to store at least the catalog plus one specific database.
 * Multi-databases are essentially a map from "magic values" to database buffers.
 */
class MultiDatabase {
   public:
    /// Create a multidatabase from a serialized multidatabase.
    explicit MultiDatabase(const unsigned char *buffer) : buffer(buffer){};
    const unsigned char *SubDatabasePointer(int32_t magicValue) const;

   private:
    const unsigned char *buffer;
};

/// Class for easily creating a MultiDatabase
class MultiDatabaseBuilder {
   public:
    MultiDatabaseBuilder()
        : buffer((unsigned char *)calloc(1, kMultiDatabaseTocLength)), bulkLength(0){};
    ~MultiDatabaseBuilder();

    unsigned char *AddSubDatabase(int32_t magicValue, long length);

    /// When done adding databases, use this to get the buffer you should write to disk.
    unsigned char *Buffer() { return buffer; };
    /// The length of the buffer returned by Buffer
    long BufferLength() { return kMultiDatabaseTocLength + bulkLength; };

   private:
    // Throughout LOST, most dynamic memory is managed with `new` and `delete` to make it easier to
    // use unique pointers. Here, however, we use realloc, so C-style memory management.
    unsigned char *buffer;
    // how many bytes are presently allocated for databases (excluding map)
    long bulkLength;
};

}  // namespace lost

#endif
