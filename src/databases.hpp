#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <stdlib.h>
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
        : buffer((unsigned char *)calloc(1, kMultiDatabaseTocLength)), bulkLength(0) { };
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

//Hash table database storing patterns in buckets to reduce number + time of database calls 
class HashTableDatabase {
public: 
    HashTableDatabase();
private: 
    /* Number of stars in star pattern. */
    /* Must be at least 3.  Recommended number is 4. */
    int num_stars_in_pattern;
    /* Minimum star brightness (in magnitude) for inclusion in catalog. */
    /* Note that lower magnitude values are brighter. */
    float min_magnitude;
    /* Maximum Field of View for catalog in radians. */
    /* Also the maximum angle between any two stars in a tetrahedron. */
    /* Typically equal to the angle subtended by the imager's diagonal. */
    /* Must be less than pi, but should be much less to prevent invalidation */
    /* of small angle approximation and to minimize non-linearity of FOV error. */
    /* .247 radians is about 14.1 degrees or the diagonal of a 10 degree FOV */
    float max_fov;
    /* Maximum star coordinate centroiding error as fraction of maximum FOV. */
    /* .001 is .1% of the max FOV or 1.414 pixels in a 1000x1000 image. */
    // 1 / (1024*sqrt(2)) < .00069054
    float max_centroid_error; 
    /* Maximum error in imager FOV estimate as fraction of true imager FOV. */
    /* max_fov*(1+max_fov_error) must be less than pi, but should be much less. */
    /* .01 for a 10 degree FOV imager covers estimates from 9.9 to 10.1 degrees. */
    float max_fov_error; 
    /* Minimum angle between stars in catalog in radians. */
    /* Optionally used to remove double stars, set to 0 otherwise. */
    /* .001 is about 5.7 pixels distance with 1000 pixels over a 10 degree FOV. */
    float min_star_separation;
    /* Ratio of bin size to error range, where bins are the discretization of */
    /* Pattern data values to allow them to be hashed and where the error range covers */
    /* a range of twice the data value's maximum error.  When a data value's error range */
    /* overlaps two bins, it's replicated into both.  By linearity of expectation, */
    /* the expected number of replicas of a given Pattern is: */
    /* (1+1/bin_size_ratio)^num_dims, where num_dims is the number of data values */
    /* stored in Patterns.  The expected ratio of Patterns with matching bins to */
    /* Patterns with matching values is: (1 + bin_size_ratio)^num_dims. */
    /* The bin_size_ratio represents a tradeoff between catalog size and */
    /* runtime, as more replicas means a larger catalog and more mismatching results */
    /* means more time spent checking for a match.  In most cases, Patterns with */
    /* matching values are rare, so a larger bin_size_ratio is better.  Must be */
    /* greater than 1.0 without changes to Pattern struct, recommended value of 2.0. */
    float bin_size_ratio; 
    /* Ratio specifying how much of the catalog is occupied by Patterns rather than */
    /* empty space.  Values below .5 create unnecessarily large catalog sizes, while */
    /* values above .7 cause longer lookup times due to increased collision frequency. */
    float catalog_density; 
    /* Size of catalog cache in number of patterns.  Generating the catalog in cached */
    /* pieces reduces disk I/O and greatly speeds up catalog generation. */
    long max_cat_cache_size; 
    /* Maximum distance matching pattern can be from original offset in catalog. */
    /* Measured in number of patterns.  Determines overlap between cached catalog */
    /* pieces during catalog generation and the number of extra catalog positions */
    /* needed at the end of the catalog due to the catalog hash table being non-cyclic. */
    int max_probe_depth; 
    /* Number of entries in the Hipparchos catalog. */
    int STARN; 
    /* The current calendar year. */
    int current_year;  

    /* Feature struct */
    /* Data format of what's stored in a given feature of a catalog Pattern. */
    /* A coordinate system is constructed by centering the largest edge along the x axis, */
    /* where edges are straight lines between stars in the pattern. */
    /* The x and y coordinates of every star are divided by the length of the largest edge. */
    /* This places two stars at the (x,y) coordinates (-.5,0) and (.5,0) and the remaining */
    /* stars in one of two places due to a 180 degree rotational ambiguity. */
    /* Each Feature encodes each of the remaining stars' bins, coordinates, and id. */
    struct Feature {
        /* Implicitly encoded doubles are used to store the coordinates. */
        /* This cuts down on disc usage by instead using integers with an */
        /* implicit division factor, as all coordinates are between -1 and 1. */
        int x : 15;
        /* Bins are identifying values for the Pattern's catalog position.  */
        /* They discretize Pattern data values, allowing them to be hashed. */
        /* A bin offset is the offset of this Pattern's bin from the */
        /* lowest bin the Pattern can occupy.  With a bin_size_ratio */
        /* greater than 1.0, bin offsets can be stored in a single bit. */
        unsigned int x_bin_offset : 1;
        int y : 15;
        unsigned int y_bin_offset : 1;
        /* ID of the Feature's star.  Unsigned 15 bit integers support 32768 unique stars. */
        unsigned int star_id : 15;
        /* Extra, unused bit. */
        unsigned int pad : 1;
    };

    /* Pattern struct */
    /* Data format of what's stored in a given catalog position for a given star pattern. */
    /* The rotation that results in the largest x bin (y bin as tie-breaker) */
    /* is chosen to deal with the pattern's 180 degree rotational ambiguity. */
    /* The largest edge length is also stored to allow for sanity-checking of FOV/image scale. */
    class Pattern {
        public:
            /* Features are stored in order of increasing x bin, with increasing y bin */
            /* as a tie breaker.  Two Features cannot contain the same bin pair, as this */
            /* would cause ambiguities while matching image stars to catalog stars. */
            std::vector<Feature> features;
            /* Length of the largest edge in the star pattern.  Used to recompute the FOV as well as */
            /* sanity check matches.  Also stored as an implicitly encoded double between 0 and 1, */
            /* as edge lengths must be between 0 and max_le_length, the maximum largest edge length. */
            uint16_t largest_edge;
            /* As the largest edge length is always greater than zero, it is also used as */
            /* a boolean to indicate whether the hash table slot contains a star pattern. */
            /* Note that this is therefore not explicitly set during catalog creation. */
            #define has_pattern largest_edge
            unsigned int le_bin_offset : 1;
            unsigned int fixed_star_id_1 : 15;
            unsigned int fixed_star_id_2 : 15;
            /* Boolean value indicating whether this catalog position contains the last */
            /* matching pattern in the catalog.  This avoids having to probe to the next slot. */
            unsigned int is_last : 1;

            Pattern();
    };

    /* Star struct */
    /* Data format of what's stored in the stars array for a given star. */
    struct Star {
        /* Unit vector x,y,z coordinates pointing to the star from the center of the celestial sphere. */
        double vec[3];
        /* Star magnitude */
        double mag;
        /* Star id, also known as the Hipparchos number */
        unsigned int star_id;
    };

    struct CatalogEntry {
        double RA; 
        double DEC;
        int XNO; 
        char IS; 
        int mag; 
    };
};

}

#endif
