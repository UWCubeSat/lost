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

/**
 * Obtains the size of the buffer with just the metadata and KVectorND number fields
 * 
 * @param totalBins The total number of bins stored in this KVectorND
 * 
 * @return The size, in bits, of the KVectorND with just the metadata fields
 * as well as the index numbers that map to the beginning of each bin
 * */
long SerializeLengthKVectorNDIndex(int totalBins) {
    return sizeof(int32_t) + 4 * sizeof(float) + 4 * sizeof(float) 
            + sizeof(int32_t) + sizeof(int32_t) * (totalBins + 1);
}

/**
 * Obtains the size of an entire KVectorND in bits
 * 
 * @param numEntries The number of star patterns (quadruples) stored in a given KVectorND
 * @param bins The number of bins that are stored in a given KVectorND
 * 
 * @return The size, in bits, of the KVector with a given number of star patterns (quadruples)
 * and bins given by the parameters
*/
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

/**
 * Creates a new KVectorND from a given buffer that holds that information
 * 
 * @param buffer The string buffer which holds the KVectorND information
 * 
 * @note Creates a KVectorND based off of the information stored in the buffer
 * 
 * @pre Requires that the buffer points to the beginning of the KVectorND
 * @warning Precondition: Will create an incorrect KVectorND if buffer does not point to the beginning
 * of the KVectorND
*/
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

/**
 * Destroys a KVectorND
 * 
 * @post Destroys this
 * 
*/
KVectorND::~KVectorND() {
    free(min);
    free(max);
    free(binWidth);
}

/**
 * Initializes the metadata fields and the index mapping from bins to entries and stores them into a buffer
 * for a given KVectorND
 * 
 * @param kVector The vector which holds the mappings from bins to star pattern entries
 * @param min The array that holds the minimum value of each axis defined by the KVector. In this case, it
 * is the minimum value obtained from each of the 4 star parameters.
 * @param max The array that holds the maximum value of each axis defined by the KVector. In this case, it
 * is the maximum value obtained from each of the 4 star parameters.
 * @param numBins The number of bins at each axis.
 * @param buffer The buffer to store the metadata and index mappings into
 * @param entries The number of entries found in a given KVectorND
 * 
 * @return The buffer position that points to the next empty slot after the metadata fields and index mappings
 * have been stored
 * 
 * @note buffer is mutated, storing the metadata fields (Number of Entries, Axis Bounds, Number of Bins at Each Bound)
 * and then the index mappings
 * 
 * @post buffer now contains metadata fields and index mappings for the given KVectorND
*/
unsigned char* SerializeQuadKVectorIndex(const std::vector<int32_t> &kVector, float* min, float* max, long numBins, unsigned char *buffer, int entries) {

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
    for (int i = 0; i < kVector.size(); i++) {
        *(int32_t *)buffer = kVector[i];
        buffer += sizeof(int32_t);
    }

    // TODO: verify length. Here is code for the normal KVector
    assert(buffer - bufferStart == SerializeLengthKVectorNDIndex(numBins * numBins * numBins * numBins));

    return buffer;
}

/**
 * Assistant function to CatalogToQuadDistances, which inserts a given star into a vector, such that
 * the vector's entries are ordered by increasing distance away from a central star
 * 
 * @param catalog The catalog of stars that the indexes are based off of
 * @param quad The vector containing star indexes, sorted in increasing order from the central star
 * @param centralIndex The index corresponding to the star in catalog that is the central star
 * @param starIndex THe index corresponding to a star in catalog that will be inserted
 * 
 * @note quad is altered such that it includes starIndex, while preserving the original ordering of quad,
 * which is to say that quad still is sorted from increasing distance from the star of centralIndex
*/
void Insert(const Catalog &catalog, std::vector<int16_t> &quad, int16_t centralIndex, int16_t starIndex) {
    assert(centralIndex > 0 && centralIndex < catalog.size());
    assert(starIndex > 0 && starIndex < catalog.size());
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
 * Produces a star parameter characterising a star pattern of 3 stars given 3 dot products
 * 
 * @param centralToOne,centralToTwo,oneToTwo The dot products between 3 stars
 * 
 * @return The modified Liebe star parameter charactising this star triplet
 * 
 * @pre centralToOne + centralToTwo != 2
 * @warning Precondition: Precondition: centralToOne + centralToTwo != 2
 */ 
float StarParameterA(float centralToOne, float centralToTwo, float oneToTwo) {
    assert(centralToOne + centralToTwo != 2);
    return (1 - oneToTwo) / (2 - centralToOne - centralToTwo);
}

/**
 * Produces a star parameter characterising a star pattern of 4 stars given 3 dot products
 * 
 * @param centralToOne,centralToTwo,centralToThree The dot products between 4 stars
 * 
 * @return The modified Liebe star parameter characterising this star quadruple
 * 
 * @pre centralToOne != centralToThree
 * @warning Precondition: centralToOne != centralToThree
*/
float StarParameterB(float centralToOne, float centralToTwo, float centralToThree) {
    assert(centralToOne != centralToThree);
    return (centralToTwo - centralToOne) / (centralToOne - centralToThree);
}

/**
 * Produces the star parameter set for a given star quadruple
 * 
 * @param central The vector of the central star
 * @param stars The vector of the 3 other stars
 * 
 * @return A vector holding the 4 star parameters corresponding to the given star quadruple.
 * If we call the central star 0, and the others 1, 2, and 3, in order of increasing
 * angular distance from 0, then the star parameters originate from the following patterns,
 * in this exact order:
 * (0, 1, 2), (0, 1, 3), (0, 2, 3), (0, 1, 2, 3)
 * 
 * @pre central != null && stars != null && entries in stars are not null && stars.size() == 3
 * @warning Precondition: central != null && stars != null && entries in stars are not null && stars.size() == 3
*/
std::vector<float> StarParameterSet(Vec3 central, std::vector<Vec3> stars) {
    assert(stars.size() == 3);
    std::vector<float> result;
    for(int i = 0; i <= 1; i++) {
        for(int j = i + 1; j <= 2; j++) {
                result.push_back(StarParameterA(central * stars.at(i), 
                        central * stars.at(j), 
                        stars.at(i) * stars.at(j)));
        }
    }
    result.push_back(StarParameterB(central * stars.at(0), 
                    central * stars.at(1), 
                    central * stars.at(2)));
    assert(result.size() == 4);
    return result;
}

/**
 * Calculates all star patterns (quadruples) for a KVectorND that follow the requirements:
 *      1. Each pattern has 3 stars surrounding a central star, where those 3 stars are the 
 *         closest out of all stars to the central one
 *      2. In each pattern, the closest star to the central one has a distance greater than
 *         the specified minimum distance
 *      3. In each pattern, the farthest star from the central one has a distance less than
 *         the specified maximum distance
 * 
 * @param catalog The catalog of stars to generate the patterns from
 * @param minDistance The specified minimum distance used as specified above
 * @param maxDistance The specified maximum distance used as specified above
 * 
 * @return A vector of KVectorQuad structures, which hold all the star patterns obeying the rules
 * above with their indicies corresponding to catalog and their corresponding star parameter set.
 * 
*/
std::vector<KVectorQuad> CatalogToQuadDistances(const Catalog &catalog, float minDistance, float maxDistance) {
    // Holds both the indexes and the corresponding vectors of stars
    std::vector<KVectorQuad> result;

    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {

        // Holds indexes of stars
        std::vector<int16_t> quad;
        for (int16_t j = 0; j < (int16_t)catalog.size(); j++) {
            if(i != j) {
                Insert(catalog, quad, i, j);
                if(quad.size() > 3) {
                    quad.resize(3);
                }
            }
        }
        assert(quad.size() == 3);
        if(AngleUnit(catalog[i].spatial, catalog[quad.at(0)].spatial) > minDistance && AngleUnit(catalog[i].spatial, catalog[quad.at(2)].spatial) < maxDistance) {
            std::vector<float> parameters = StarParameterSet(catalog[i].spatial, {catalog[quad.at(0)].spatial, catalog[quad.at(1)].spatial, catalog[quad.at(2)].spatial});
            result.push_back({{i, quad.at(0), quad.at(1), quad.at(2)}, parameters});
        }
    }
    // Returns a list of star patterns of 4 stars where each pattern contains the indexes and then the 
    // modified Liebe star parameter sets
    return result;
}

/**
 * Provides the bin value for the a given set of values corresponding to a star quadruple
 * 
 * @param values The vector holding the star parameter set for a given star quadruple
 * @param offset The offset to apply to each star parameter in the star parameter set
 * @param scale The scaling value for each axis (corresponding to each star parameter)
 * @param product The vector holding the progressive number of bins as the axis progress
 * 
 * @return The bin value for a star quadruple with the given values
 * 
 * @pre values.size() == 4
 * @warning Precondition: values.size() == 4
 * 
*/
long BinFunctionND(const std::vector<float> &values, float* offset, float* scale, int* product) {
    assert(values.size() == 4);
    int bin = 0;
    for(int i = 0; i < 4; i++) {
        if(std::floor(scale[i] * (values[i] - offset[i])) == product[1]) {
            bin += ((int)product[1] - 1) * product[i];
        } else {
            bin += std::floor(scale[i] * (values[i] - offset[i])) * product[i];
        }
    }
    assert(bin < product[4]);
    assert(bin > 0);
    return bin;
}

/**
 * Serializes an entire KVectorND into a buffer
 * 
 * @param quads The entries of star quadruples for a given KVectorND
 * @param numBins The number of bins along each axis (corresponding to each star parameter type)
 * @param buffer The buffer to place the data into
 * 
 * @pre quads.size() > 0
 * @warning Precondition: quads.size() > 0
 * 
 * @note Modifies buffer by placing the data relavent to the KVectorND in the following order:
 *      1. Number of Entries
 *      2. Axis Bounds
 *      3. Number of Bins at Each Bound
 *      4. Index Mappings from Bin number to Entries
 *      5. Star Quadruple Entries
 * 
 * @post buffer contains all information pertaining to the related KVectorND
*/
void SerializeKVectorND(std::vector<KVectorQuad> &quads, long numBins, unsigned char *buffer) {
    assert(quads.size() > 0);

    // TODO: Statistics Printout

    // Preprocessing: Finding boundaries along axes (min and max parameters A and B)
    float min[4] = {quads.at(0).parameters.at(0), quads.at(0).parameters.at(1), quads.at(0).parameters.at(2), quads.at(0).parameters.at(3)};
    float max[4] = {quads.at(0).parameters.at(0), quads.at(0).parameters.at(1), quads.at(0).parameters.at(2), quads.at(0).parameters.at(3)};
    for(int i = 1; i < quads.size(); i++) {
        for(int j = 0; j < 4; j++) {
            float param = quads.at(i).parameters.at(j);
            if(param < min[j]) {
                min[j] = param;
            }
            if(param > max[j]) {
                max[j] = param;
            }
        }
    }
    
    // Preprocessing: Calculating scaling factors and product,
    // a measure of the progressive amount of bins as the
    // user progresses through the axes
    float scale[4];
    int product[5] = {1, 1, 1, 1, 1};
    for(int i = 0; i < 4; i++) {
        scale[i] = numBins / (max[i] - min[i]);
        product[i + 1] = product[i] * numBins;
    }
    
    // Preprocessing: Finding bins numbers for all elements
    KVectorQuad* pointers[(int)quads.size()];
    long bins[(int)quads.size()];
    for(int i = 0; i < quads.size(); i++) {
        pointers[i] = &quads[i];
        bins[i] = BinFunctionND(quads[i].parameters, min, scale, product);
    }

    // Producing the mappings from bin numbers to
    // index within array of star quadruples
    // that will appear in buffer
    std::vector<int32_t> kVector((int)product[4] + 1);
    for(int i = 0; i < quads.size(); i++) {
        kVector[bins[i]]++;
    }
    int32_t sum = 0;
    for(int i = 0; i < product[4] + 1; i++) {
        int32_t temp = kVector[i];
        kVector[i] = sum;
        sum += temp;
    }
    kVector[product[4]] = sum;

    
    // Sorting all star quadruples in order of their bin
    // number
    KVectorQuad* sortedPointers[(int)product[4]];
    int positions[(int)product[4]];
    for(int i = 0; i < product[4]; i++) {
        positions[i] = 0;
    }
    for(int i = 0; i < quads.size(); i++) {
        int bin = bins[i];
        sortedPointers[kVector[bin] + positions[bin]] = pointers[i];
        positions[bin]++;
    }
    
    // For checking postcondition
    unsigned char *bufferStart = buffer;

    // Store metadata and index mappings into the buffer
    buffer = SerializeQuadKVectorIndex(kVector, min, max, numBins, buffer, quads.size());

    // Store star quadruple entries in order
    for (int i = 0; i < quads.size(); i++) {
        for(int j = 0; j < 4; j++) {
            *(int16_t *)buffer = sortedPointers[i]->stars[j];
            buffer += sizeof(int16_t);
        }
    }

    // Verify postcondition
    assert(buffer - bufferStart == SerializeLengthQuadStarKVectorND(quads.size(), numBins * numBins * numBins * numBins));

}

/**
 * Provides the mapping between a given bin number and a star pattern
 * entry
 * 
 * @param binNumber The bin number for finding the corresponding star
 * pattern in this
 * 
 * @return The index of the star pattern to look for
*/
int32_t KVectorND::KVectorToStarQuadIndex(const long binNumber) const {
    const int32_t targetIndex = *(bins + (int32_t)binNumber);
    assert(targetIndex > 0 && targetIndex < *((int32_t *)quads - 1));
    return targetIndex;
}

/**
 * Gets the complete star pattern for the star pattern at a given index
 * 
 * @param index The index of the star pattern in this
 * 
 * @return An array of length 4 that contains the indexes of all stars
 * in the corresponding catalog for a given star pattern
 * 
*/
int16_t *KVectorND::GetEntry(const long index) const {
    int16_t stars[] = {0, 0, 0, 0};
    int i = 0;
    for(const int16_t *ptr = quads + 4 * index; ptr < quads + 4 + 4 * index; ptr++) {
        stars[i] = *ptr;
        i++;
    }
    assert(i == 4);
    return stars;
}

/**
 * Performs a range search on this KVectorND and finds all possible
 * star patterns meeting the specified criterion
 * 
 * @param minParameter The array expressing all minimum star parameter values
 * that the results must fit
 * @param maxParameter The array expressing all maximum star parameter values
 * that the results must fit
 * 
 * @return A vector containing all star patterns that fit the maximum and minimum
 * star parameter bounds as specified by the parameters
 * 
 * @pre minParameter and maxParameter have lengths of 4
 * @warning minParameter and maxParameter have lengths of 4
*/
std::vector<int16_t *> KVectorND::RangeSearch(const float *minParameter, const float *maxParameter) const {
    // Somehow assert or guarentee the size of each array - We can't :(
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
        if(std::floor(scale[i] * (minParameter[i] - min[i])) == product[1]) {
            minBins[i] = ((int)product[1] - 1) * product[i];
        } else {
            minBins[i] = std::floor(scale[i] * (minParameter[i] - min[i])) * product[i];
        }
        if(std::floor(scale[i] * (maxParameter[i] - min[i])) == product[1]) {
            maxBins[i] = ((int)product[1] - 1) * product[i];
        } else {
            maxBins[i] = std::floor(scale[i] * (maxParameter[i] - min[i])) * product[i];
        }
        assert(maxBins[i] >= minBins[i]);
        b0 += minBins[i];
    }
    assert(b0 <= product[4]);

    const long difference = maxBins[0] - minBins[0];

    for(int k = 0; k <= maxBins[3] - minBins[3]; k++) {
        for(int j = 0; j <= maxBins[2] - minBins[2]; j++) {
            for(int i = 0; i <= maxBins[1] - minBins[1]; i++) {

                long begin = b0 + i * numBins + j * numBins * numBins + k * numBins * numBins * numBins;
                assert(b0 <= product[4]);

                long end = begin + difference;
                assert(end <= product[4]);
                assert(begin <= end);

                int32_t qBegin = KVectorToStarQuadIndex(begin);
                int32_t qEnd = KVectorToStarQuadIndex(end);

                assert(qBegin <= qEnd);

                for(qBegin;  qBegin <= qEnd; qBegin++) {
                    result.push_back(GetEntry(qBegin));
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
