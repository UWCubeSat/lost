#include "databases.hpp"

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <vector>
#include <algorithm>
#include <iostream>

#include "attitude-utils.hpp"
#include "star-utils.hpp"

namespace lost {

struct KVectorPair {
    int16_t index1;
    int16_t index2;
    float distance;
};

bool CompareKVectorPairs(const KVectorPair &p1, const KVectorPair &p2) {
    return p1.distance < p2.distance;
}

// just the index part of the kvector, doesn't store the sorted list it refers to. This makes it
// flexible whether used to store star distance angles, an angle of a triangle, etc.
/**
 K-vector index layout. No magic value because its never used on its own.
 | size          | name       | description                                                 |
 |---------------+------------+-------------------------------------------------------------|
 | 4             | numEntries |                                                             |
 | sizeof float  | min        | minimum value contained in the database                     |
 | sizeof float  | max        | max value contained in index                                |
 | 4             | numBins    |                                                             |
 | 4*(numBins+1) | bins       | The `i'th bin (starting from zero) stores how many pairs of |
 |               |            | stars have a distance lesst han or equal to:                |
 |               |            | min+i*(max-min)/numBins                                     |
 */

/// The number of bytes that a kvector index will take up whe serialized
long SerializeLengthKVectorIndex(long numBins) {
    return 4+sizeof(float)+sizeof(float)+4+4*(numBins+1);
}

// apparently there's no easy way to accept an iterator argument. Hate java all you want, but at
// least it can do that!
// https://stackoverflow.com/questions/5054087/declare-a-function-accepting-generic-iterator or
// rather, the correct way is to use a template and just use the ++ and * operators willy-nilly,
// which will throw an error if the operators are not implemented.

/**
 * Serialize a KVector index to disk.
 * Use SerializeLengthKVectorIndex to determine how long the buffer should be.
 * @param values The actual entries the kvector should be referring to, sorted in ascending order.
 * @pre values must be sorted in ascending order!
 * @param min,max Guaranteed bounds on the entries of values
 * @todo Consider replacing min and max parameters by just calculating the min and max of values?
 * @param numBins the number of "bins" the KVector should use. A higher number makes query results "tighter" but takes up more disk space. Usually should be set somewhat smaller than (max-min) divided by the "width" of the typical query.
 * @param buffer[out] index is written here.
 */
void SerializeKVectorIndex(const std::vector<float> &values, float min, float max, long numBins, unsigned char *buffer) {
    std::vector<int32_t> kVector(numBins+1); // numBins = length, all elements zero
    float binWidth = (max - min) / numBins;

    // generate the k-vector part
    // Idea: When we find the first star that's across any bin boundary, we want to update all the newly sealed bins
    long lastBin = 0; // first bin the last star belonged to
    for (int32_t i = 0; i < (int32_t)values.size(); i++) {
        if (i > 0) {
            assert(values[i] >= values[i-1]);
        }
        assert(values[i] >= min);
        assert(values[i] <= max);
        long thisBin = (long)ceil((values[i] - min) / binWidth); // first bin we belong to
        assert(thisBin >= 0);
        assert(thisBin <= numBins); // thisBin == numBins is acceptable since kvector length == numBins + 1
        // if thisBin > lastBin, then no more stars can be added to those bins between thisBin and lastBin, so set them.
        for (long bin = lastBin; bin < thisBin; bin++) {
            kVector[bin] = i; // our index is the number of pairs with shorter distance
        }
        lastBin = thisBin;
    }
    for (long bin = lastBin; bin <= numBins; bin++) {
        kVector[bin] = values.size();
    }

    // verify kvector
    int32_t lastBinVal = -1;
    for (const int32_t &bin : kVector) {
        assert(bin >= lastBinVal);
        lastBinVal = bin;
    }

    unsigned char *bufferStart = buffer;
    // metadata fields
    *(int32_t *)buffer = values.size();
    buffer += sizeof(int32_t);
    *(float *)buffer = min;
    buffer += sizeof(float);
    *(float *)buffer = max;
    buffer += sizeof(float);
    *(int32_t *)buffer = numBins;
    buffer += sizeof(int32_t);

    // kvector index field
    // you could probably do this with memcpy instead, but the explicit loop is necessary for endian
    // concerns? TODO endianness
    for (const int32_t &bin : kVector) {
        *(int32_t *)buffer = bin;
        buffer += sizeof(int32_t);
    }

    // TODO: verify length
    // assert(buffer - bufferStart == SerializeLengthKVectorIndex(numBins));
}

/// Construct from serialized buffer.
KVectorIndex::KVectorIndex(const unsigned char *buffer) {
    numValues = *(int32_t *)buffer;
    buffer += sizeof(int32_t);
    min = *(float *)buffer;
    buffer += sizeof(float);
    max = *(float *)buffer;
    buffer += sizeof(float);
    numBins = *(int32_t *)buffer;
    buffer += sizeof(int32_t);

    assert(min >= 0.0f);
    assert(max > min);
    binWidth = (max - min) / numBins;

    bins = (const int32_t *)buffer;
}

/**
 * Finds all the entries in the given range, and possibly a few just outside the range on the ends.
 * @param upperIndex[out] Is set to the index of the last returned value +1.
 * @return the index (starting from zero) of the first value matching the query
 */
long KVectorIndex::QueryLiberal(float minQueryDistance, float maxQueryDistance, long *upperIndex) const {
    assert(maxQueryDistance > minQueryDistance);
    if (maxQueryDistance >= max) {
        maxQueryDistance = max - 0.00001; // TODO: better way to avoid hitting the bottom bin
    }
    if (minQueryDistance <= min) {
        minQueryDistance = min + 0.00001;
    }
    if (minQueryDistance > max || maxQueryDistance < min) {
        *upperIndex = 0;
        return 0;
    }
    long lowerBin = BinFor(minQueryDistance);
    long upperBin = BinFor(maxQueryDistance);
    assert(upperBin >= lowerBin);
    assert(upperBin <= numBins);
    // bins[lowerBin-1]=number of pairs <= r < query distance, so it is the index of the
    // first possible item that might be equal to the query distance
    int lowerIndex = bins[lowerBin-1];
    if (lowerIndex >= numValues) {
        // all pairs have distance less than queried. Return value is irrelevant as long as
        // numReturned=0
        return 0;
    }
    // bins[upperBin]=number of pairs <= r >= query distance
    *upperIndex = bins[upperBin];
    return lowerIndex;
}

/// return the lowest-indexed bin that contains the number of pairs with distance <= dist
long KVectorIndex::BinFor(float query) const {
    long result = (long)ceil((query - min) / binWidth);
    assert(result >= 0);
    assert(result <= numBins);
    return result;
}

/**
 pair K-vector database layout. The kvector appears before the bulk pair data because it contains the
 number of pairs, which is necessary to read the bulk pair data.

     | size (bytes)             | name         | description                                                 |
     |--------------------------+--------------+-------------------------------------------------------------|
     | sizeof kvectorIndex      | kVectorIndex | Serialized KVector index                                    |
     | 2*sizeof(int16)*numPairs | pairs        | Bulk pair data                                              |
 */
std::vector<KVectorPair> CatalogToPairDistances(const Catalog &catalog, float minDistance, float maxDistance) {
    std::vector<KVectorPair> result;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t k = i+1; k < (int16_t)catalog.size(); k++) {

            KVectorPair pair = { i, k, AngleUnit(catalog[i].spatial, catalog[k].spatial) };
            assert(isfinite(pair.distance));
            assert(pair.distance >= 0);
            assert(pair.distance <= M_PI);

            if (pair.distance >= minDistance && pair.distance <= maxDistance) {
                // we'll sort later
                result.push_back(pair);
            }
        }
    }
    return result;
}

long SerializeLengthPairDistanceKVector(long numPairs, long numBins) {
    return SerializeLengthKVectorIndex(numBins) + 2*sizeof(int16_t)*numPairs;
}

/// Number of bytes that a serialized KVectorDatabase will take up
long SerializeLengthPairDistanceKVector(const Catalog &catalog, float minDistance, float maxDistance, long numBins) {
    return SerializeLengthPairDistanceKVector(CatalogToPairDistances(catalog, minDistance, maxDistance).size(), numBins);
}

/**
 * Serialize a pair-distance KVector into buffer.
 * Use SerializeLengthPairDistanceKVector to determine how large the buffer needs to be. See command line documentation for other options.
 */
void SerializePairDistanceKVector(const Catalog &catalog, float minDistance, float maxDistance, long numBins, unsigned char *buffer) {
    std::vector<int32_t> kVector(numBins+1); // numBins = length, all elements zero
    std::vector<KVectorPair> pairs = CatalogToPairDistances(catalog, minDistance, maxDistance);

    // sort pairs in increasing order.
    std::sort(pairs.begin(), pairs.end(), CompareKVectorPairs);
    std::vector<float> distances;

    for (const KVectorPair &pair : pairs) {
        distances.push_back(pair.distance);
    }

    unsigned char *bufferStart = buffer;

    // index field
    SerializeKVectorIndex(distances, minDistance, maxDistance, numBins, buffer);
    buffer += SerializeLengthKVectorIndex(numBins);

    // bulk pairs field
    for (const KVectorPair &pair : pairs) {
        *(int16_t *)buffer = pair.index1;
        buffer += sizeof(int16_t);
        *(int16_t *)buffer = pair.index2;
        buffer += sizeof(int16_t);
    }

    // verify length
    assert(buffer - bufferStart == SerializeLengthPairDistanceKVector(pairs.size(), numBins));
}

// Obtains the size of the buffer with just the metadata and KVector number fields
long SerializeLengthKVectorNDIndex(int totalBins) {
    return sizeof(int32_t) + 4 * sizeof(float) + 4 * sizeof(float) 
            + sizeof(int32_t) + sizeof(int32_t) * (totalBins + 1);
}

// Obtains the size of an entire KVectorND database
long SerializeLengthQuadStarKVectorND(int numEntries, int bins) {
    return SerializeLengthKVectorNDIndex(bins) + 4 * sizeof(int16_t) * numEntries;
}

/// Create the database from a serialized buffer.
PairDistanceKVectorDatabase::PairDistanceKVectorDatabase(const unsigned char *buffer)
    : index(KVectorIndex(buffer)) {

    // TODO: errors? (not even sure what i meant by this comment anymore)
    buffer += SerializeLengthKVectorIndex(index.NumBins());
    pairs = (const int16_t *)buffer;
}

KVectorND::KVectorND(const unsigned char *buffer) {
    numValues = *(int32_t *)buffer; // Entries
    buffer += sizeof(int32_t);

    // Bounds
    min = (float *) malloc(4 * sizeof(float));
    max = (float *) malloc(4 * sizeof(float));
    for(int i = 0; i < 4; i++) {
        min[i] = *(float *)buffer;
        buffer += sizeof(float);
        max[i] = *(float *)buffer;
        buffer += sizeof(float);
    }

    long numBins = *(int32_t *)buffer; // quartic root of the real number of bins
    buffer += sizeof(int32_t);

    binWidth = (float *) malloc(4 * sizeof(float));
    for(int i = 0; i < 4; i++) {
        binWidth[i] = (max[i] - min[i]) / numBins;
        buffer += sizeof(float);
    }

    bins = (int32_t *)buffer;

    buffer += SerializeLengthKVectorNDIndex(numBins * numBins * numBins * numBins);
    quads = (int16_t *)buffer;
}

KVectorND::~KVectorND() {
    free(min);
    free(max);
    free(binWidth);
}

unsigned char* SerializeQuadKVectorIndex(const int32_t kVector[], float* min, float* max, long numBins, unsigned char *buffer, int entries) {

    unsigned char *bufferStart = buffer;
    // metadata fields
    *(int32_t *)buffer = entries; // Entries
    buffer += sizeof(int32_t);

    // Bounds
    for(int i = 0; i < 4; i++) {
        *(float *)buffer = min[i];
        buffer += sizeof(float);
        *(float *)buffer = max[i];
        buffer += sizeof(float);
    }

    *(int32_t *)buffer = numBins; // quartic root of the real number of bins
    buffer += sizeof(int32_t);

    // kvector index field
    for (int i = 0; i < sizeof(kVector); i++) {
        *(int32_t *)buffer = kVector[i];
        buffer += sizeof(int32_t);
    }

    // verify length
    assert(buffer - bufferStart == SerializeLengthKVectorIndex(numBins));

    return buffer;
}

void Insert(const Catalog &catalog, std::vector<int16_t> quad, int16_t centralIndex, int16_t starIndex) {
    int size = quad.size();
    for(int i = 0; i < quad.size(); i++) {
        int dist1 = AngleUnit(catalog[centralIndex].spatial, catalog[i].spatial);
        int dist2 = AngleUnit(catalog[centralIndex].spatial, catalog[starIndex].spatial);
        if(dist2 < dist1) {
            quad.insert(quad.cbegin() + i, starIndex);
            break;
        }
    }
    if(quad.size() - size == 0) {
        quad.push_back(starIndex);
    }
}



/** 
 * Produces a parameter characterising a star pattern of 3 stars given 3 dot products
 * @param centralToOne, centralToTwo, and oneToTwo are the dot products between 3 stars
 * @param centralToOneError, centralToTwoError, and oneToTwo error are the approximate errors of the above parameters
 * @return The modified Liebe parameter charactising this star triplet
 */ 
float StarParameterA(float centralToOne, float centralToTwo, float oneToTwo) {
    return (1 - oneToTwo) / (2 - centralToOne - centralToTwo);
}

/**
 * Produces a parameter characterising a star pattern of 4 stars given 3 dot products
 * @param centralToOne, centralToTwo, and centralToThree are dot products between 4 stars
 * @param centralToOneError, centralToTwoError, and centralToThree error are the approximate errors of the above parameters
 * @note Using this error style does not give the true errors, but rather give a larger error radius than is real
 * @return The modified Liebe parameter characterising this star quadruple
*/
float StarParameterB(float centralToOne, float centralToTwo, float centralToThree) {
    return (centralToTwo - centralToOne) / (centralToOne - centralToThree);
}

float* StarParameterSet(Vec3 central, std::vector<Vec3> stars) {
    float result[4];
    for(int i = 0; i <= 1; i++) {
        for(int j = i + 1; j <= 2; i++) {
                result[i + j - 1] = StarParameterA(central * stars.at(i), 
                        central * stars.at(j), 
                        stars.at(i) * stars.at(j));
        }
    }
    result[3] = StarParameterB(central * stars.at(0), 
                    central * stars.at(1), 
                    central * stars.at(2));
    return result;
}


std::vector<KVectorQuad> CatalogToQuadDistances(const Catalog &catalog, float minDistance, float maxDistance) {
    std::vector<KVectorQuad> result;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        std::vector<int16_t> quad;
        for (int16_t j = 0; j < (int16_t)catalog.size(); j++) {
            if(i != j) {
                Insert(catalog, quad, i, j);
                if(sizeof(quad) > 3) {
                    quad.resize(3);
                }
            }
        }
        assert(sizeof(quad) == 3);
        if(AngleUnit(catalog[i].spatial, catalog[quad[0]].spatial) > minDistance && AngleUnit(catalog[i].spatial, catalog[quad[2]].spatial) < maxDistance) {
            int16_t stars[] = {i, quad.at(0), quad.at(1), quad.at(2)};
            result.push_back({stars, StarParameterSet(catalog[i].spatial, {catalog[quad[0]].spatial, catalog[quad[1]].spatial, catalog[quad[2]].spatial})});
        }
            
    }
    return result;
}


int BinFunctionND(float* parameters, float* offset, float* scale, int* product) {
    int bin = 0;
    for(int i = 0; i < 4; i++) {
        bin += std::floor(scale[i] * (parameters[i] - offset[i])) * product[i];
    }
    return bin;
}

void SerializeKVectorND(const Catalog &catalog, std::vector<KVectorQuad> quads, float minDistance, float maxDistance, long numBins, unsigned char *buffer) {
    assert(sizeof(quads) > 0);

    // TODO: Statistics Printout

    // Preprocessing: Finding boundaries along axes
    float min[] = {quads[0].parameters[0], quads[0].parameters[1], quads[0].parameters[2], quads[0].parameters[3]};
    float max[] = {quads[0].parameters[0], quads[0].parameters[1], quads[0].parameters[2], quads[0].parameters[3]};
    for(int i = 1; i < sizeof(quads); i++) {
        KVectorQuad quad = quads[i];
        for(int j = 0; j < 4; j++) {
            float param = quad.parameters[i];
            if(param < min[i]) {
                min[i] = param;
            }
            if(param > max[i]) {
                max[i] = param;
            }
        }
    }
    
    // Preprocessing: Calculating scaling factors
    float scale[4];
    int product[] = {1, 1, 1, 1, 1};
    for(int i = 0; i < 4; i++) {
        scale[i] = numBins / (max[i] - min[i]);
        product[i + 1] = product[i] * numBins;
    }
    
    // Preprocessing: Finding bins for all elements
    KVectorQuad* pointers[sizeof(quads)];
    int bins[sizeof(quads)];
    for(int i = 0; i < sizeof(quads); i++) {
        pointers[i] = &quads[i];
        bins[i] = BinFunctionND(quads[i].parameters, min, scale, product);
    }

    // KVector Positions
    int32_t kVector[product[4] + 1];
    for(int i = 0; i < product[4]; i++) {
        kVector[bins[i]]++;
    }
    int32_t sum = 0;
    for(int i = 0; i < sizeof(kVector); i++) {
        int32_t temp = kVector[i];
        kVector[i] = sum;
        sum += temp;
    }
    kVector[product[4]] = sum;

    // Sorting all elements
    KVectorQuad* sortedPointers[sizeof(pointers)];
    int positions[product[4]];
    for(int i = 0; i < sizeof(sortedPointers); i++) {
        int bin = bins[i];
        sortedPointers[kVector[bin] + positions[bin]] = pointers[i];
        positions[bin]++;
    }
    
    // Metadata and KVector Positions
    unsigned char *bufferStart = buffer;
    buffer = SerializeQuadKVectorIndex(kVector, min, max, numBins, buffer, sizeof(quads));

    // indexes for each quad entry
    for (int i = 0; i < sizeof(sortedPointers); i++) {
        for(int j = 0; j < 4; j++) {
            *(int16_t *)buffer = sortedPointers[i]->stars[i];
            buffer += sizeof(int16_t);
        }
    }

    // TODO: verify length
    // assert(buffer - bufferStart == SerializeLengthPairDistanceKVector(pairs.size(), numBins));

}

int16_t *KVectorND::GetEntry(const long index) const {
    int16_t stars[] = {0, 0, 0, 0};
    int i = 0;
    for(const int16_t *ptr = quads + 4 * index; ptr < quads + 4 + 4 * index; ptr++) {
        stars[i] = *ptr;
        i++;
    }
    return stars;
}

std::vector<int16_t *> KVectorND::RangeSearch(const float *maxParameter, const float *minParameter) const {
    // Somehow assert or guarentee the size of each array
    long* minBins = new long[4];
    long* maxBins = new long[4];
    std::vector<int16_t *> result;

    // Preprocessing: Calculating scaling factors
    float scale[4];
    int product[] = {1, 1, 1, 1, 1};
    for(int i = 0; i < 4; i++) {
        scale[i] = numBins / (max[i] - min[i]);
        product[i + 1] = product[i] * numBins;
    }

    long b0 = 0;
    for(int i = 0; i < 4; i++) {
        minBins[i] = std::floor(scale[i] * (minParameter[i] - min[i])) * product[i];
        maxBins[i] = std::floor(scale[i] * (maxParameter[i] - min[i])) * product[i];
        b0 += minBins[i];
    }
    const long difference = maxBins[0] - minBins[0];

    for(int k = 0; k <= maxBins[3] - minBins[3]; k++) {
        for(int j = 0; j <= maxBins[2] - minBins[2]; j++) {
            for(int i = 0; i <= maxBins[1] - minBins[1]; i++) {
                long begin = b0 + i * numBins + j * numBins * numBins + k * numBins * numBins * numBins;
                long end = begin + difference;
                for(; begin < end; begin++) {
                    result.push_back(GetEntry(begin));
                }
            }
        }
    }
    delete[] minBins;
    delete[] maxBins;
    return result;
}

/// Return the value in the range [low,high] which is closest to num
float Clamp(float num, float low, float high) {
    return num < low ? low : num > high ? high : num;
}

/**
 * Return at least all the star pairs whose inter-star distance is between min and max
 * @param end[out] Is set to an "off-the-end" pointer, one past the last pair being returned by the query.
 * @return A pointer to the start of the matched pairs. Each pair is stored as simply two 16-bit integers, each of which is a catalog index. (you must increment the pointer twice to get to the next pair).
 */
const int16_t *PairDistanceKVectorDatabase::FindPairsLiberal(
    float minQueryDistance, float maxQueryDistance, const int16_t **end) const {

    long upperIndex;
    long lowerIndex = index.QueryLiberal(minQueryDistance, maxQueryDistance, &upperIndex);
    *end = &pairs[upperIndex * 2];
    return &pairs[lowerIndex * 2];
}

/// Return the distances from the given star to each star it's paired with in the database (for debugging).
std::vector<float> PairDistanceKVectorDatabase::StarDistances(int16_t star, const Catalog &catalog) const {
    std::vector<float> result;
    for (int i = 0; i < NumPairs(); i++) {
        if (pairs[i*2] == star || pairs[i*2+1] == star) {
            result.push_back(AngleUnit(catalog[pairs[i*2]].spatial, catalog[pairs[i*2+1]].spatial));
        }
    }
    return result;
}

/**
   MultiDatabase memory layout:

   | size           | name              | description                                             |
   |----------------+-------------------+---------------------------------------------------------|
   | 8*maxDatabases | table of contents | each 8-byte entry is the 4-byte magic value followed by |
   |                |                   | a 4-byte index into the bulk where that db begins       |
   | Large          | databases         | the database contents                                   |
 */

/**
 * @brief return a pointer to the start of the database type indicated by the magic value, if such
 * a sub-database is present in the database
 * @param magicValue
 * @return Returns a pointer to the start of the database type indicated by the magic value, null if not found
 */
const unsigned char *MultiDatabase::SubDatabasePointer(int32_t magicValue) const {
    long databaseIndex = -1;
    int32_t *toc = (int32_t *)buffer;
    for (int i = 0; i < kMultiDatabaseMaxDatabases; i++) {
        int32_t curMagicValue = *toc;
        toc++;
        if (curMagicValue == magicValue) {
            databaseIndex = *toc;
            break;
        }
        toc++;
    }
    // the database was not found
    if (databaseIndex < 0) {
        return NULL;
    }

    return buffer+kMultiDatabaseTocLength+databaseIndex;
}

/**
 * Add a database to a MultiDatabase
 * @param magicValue A value unique to this type of database which is used to extract it out of the database later.
 * @param length The number of bytes to allocate for this database.
 * @return Pointer to the start of the space allocated for said database. Return null if full (too many databases).
 */
unsigned char *MultiDatabaseBuilder::AddSubDatabase(int32_t magicValue, long length) {
    // find unused spot in toc and take it!
    int32_t *toc = (int32_t *)buffer;
    bool foundSpot = false;
    for (int i = 0; i < kMultiDatabaseMaxDatabases; i++) {
        if (*toc == 0) {
            *toc = magicValue;
            toc++;
            *toc = bulkLength;
            foundSpot = true;
            break;
        }
        // skip the entry
        toc += 2;
    }

    // database is full
    if (!foundSpot) {
        return NULL;
    }

    buffer = (unsigned char *)realloc(buffer, kMultiDatabaseTocLength+bulkLength+length);
    // just past the end of the last database
    unsigned char *result = buffer+kMultiDatabaseTocLength+bulkLength;
    bulkLength += length;
    return result;
}

MultiDatabaseBuilder::~MultiDatabaseBuilder() {
    free(buffer);
}

}

// TODO: after creating the database, print more statistics, such as average number of pairs per
// star, stars per bin, space distribution between array and index.
