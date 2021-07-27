#include "databases.hpp"
#include "attitude-utils.hpp"
#include "star-utils.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

namespace lost {

struct KVectorPair {
    int16_t index1;
    int16_t index2;
    float distance;
};

// triple of stars 
// small angle and large angle in triangle
struct KVectorTriple {
    int16_t index1;
    int16_t index2;
    int16_t index3;
    float distance;
};

bool CompareKVectorPairs(const KVectorPair &p1, const KVectorPair &p2) {
    return p1.distance < p2.distance;
}

// sort triangles by the small angle in the KVector
bool CompareKVectorTriples(const KVectorTriple &t1, const KVectorTriple &t2) {
    return t1.distance < t2.distance;
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

long SerializeLengthKVectorIndex(long numBins) {
    return 4+sizeof(float)+sizeof(float)+4+4*(numBins+1);
}

// apparently there's no easy way to accept an iterator argument. Hate java all you want, but at
// least it can do that!
// https://stackoverflow.com/questions/5054087/declare-a-function-accepting-generic-iterator

// values should be presorted in ascending order
void SerializeKVectorIndex(const std::vector<float> &values, float min, float max, long numBins, unsigned char *buffer) {
    std::vector<int32_t> kVector(numBins+1); // numBins = length, all elements zero
    float binWidth = (max - min) / numBins;

    // generate the k-vector part
    // Idea: When we find the first star that's across any bin boundary, we want to update all the newly sealed bins
    long lastBin = 0; // first bin the last star belonged to
    for (int32_t i = 0; i < (int32_t)values.size(); i++) {
        assert(values[i] >= min);
        assert(values[i] <= max);
        long thisBin = (long)ceil((values[i] - min) / binWidth); // first bin we belong to
        assert(thisBin >= 0);
        assert(thisBin <= numBins); // thisBin == numBins is acceptable since kvector length == numBins + 1
        // if thisBin > lastBin, then no more stars can be added to those bins between thisBin and lastBin, so set them.
        for (long bin = lastBin; bin < thisBin; bin++) {
            kVector[bin]=i; // our index is the number of pairs with shorter distance
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

    // verify length
    assert(buffer - bufferStart == SerializeLengthKVectorIndex(numBins));
}

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

long KVectorIndex::QueryLiberal(float minQueryDistance, float maxQueryDistance, long *numReturned) const {
    assert(maxQueryDistance > minQueryDistance);
    if (maxQueryDistance >= max) {
        maxQueryDistance = max - 0.00001; // TODO: better way to avoid hitting the bottom bin
    }
    if (minQueryDistance <= min) {
        minQueryDistance = min + 0.00001;
    }
    if (minQueryDistance > max || maxQueryDistance < min) {
        *numReturned = 0;
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
    int upperIndex = bins[upperBin] - 1;
    *numReturned = upperIndex - lowerIndex + 1;
    assert(*numReturned >= 0);
    return lowerIndex;
}

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

/**
 triple K-vector database layout. The kvector appears before the bulk pair data because it contains the
 number of triples, which is necessary to read the bulk triple data.

     | size (bytes)                 | name         | description                                                 |
     |------------------------------+--------------+-------------------------------------------------------------|
     | sizeof kvectorIndex          | kVectorIndex | Serialized KVector index                                    |
     | 3*sizeof(int16)*numTriples   | triples      | Bulk triple data                                            |
 */
std::vector<KVectorTriple> CatalogToTripleDistances(const Catalog &catalog, float minDistance, float maxDistance) {
    std::vector<std::vector<int16_t>> sufficientlyClose((int16_t)catalog.size());
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t j = i+1; j < (int16_t)catalog.size(); j++) {
            float d = AngleUnit(catalog[i].spatial, catalog[j].spatial);
            if (d >= minDistance && d <= maxDistance) {
                sufficientlyClose[i].push_back(j);
            }
        }
    }
    std::vector<KVectorTriple> result;
    for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
        for (int16_t closeIndexJ = 0; closeIndexJ < (int16_t) sufficientlyClose[i].size(); closeIndexJ++) {
            int16_t j = sufficientlyClose[i][closeIndexJ];
            for (int16_t closeIndexK = 0; closeIndexK < (int16_t) sufficientlyClose[j].size(); closeIndexK++) {
                int16_t k = sufficientlyClose[j][closeIndexK];
                float d = AngleUnit(catalog[i].spatial, catalog[k].spatial);
                if (d < minDistance || d > maxDistance) {
                    continue;
                }
                int mindex;
                KVectorTriple triple = { i, j, k, minInnerAngle(catalog, mindex, i, j, k) };
                assert(isfinite(triple.distance));
                assert(triple.distance >= 0);
                assert(triple.distance <= M_PI);
                result.push_back(triple);
            }
        }
    }
    return result;
}

long SerializeLengthPairDistanceKVector(long numPairs, long numBins) {
    return SerializeLengthKVectorIndex(numBins) + 2*sizeof(int16_t)*numPairs;
}

long SerializeLengthTripleDistanceKVector(long numTriples, long numBins) {
    return SerializeLengthKVectorIndex(numBins) + 3*sizeof(int16_t)*numTriples;
}

long SerializeLengthPairDistanceKVector(const Catalog &catalog, float minDistance, float maxDistance, long numBins) {
    return SerializeLengthPairDistanceKVector(CatalogToPairDistances(catalog, minDistance, maxDistance).size(), numBins);
}

long SerializeLengthTripleDistanceKVector(const Catalog &catalog, float minDistance, float maxDistance, long numBins) {
    return SerializeLengthTripleDistanceKVector(CatalogToTripleDistances(catalog, minDistance, maxDistance).size(), numBins);
}

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

void SerializeTripleDistanceKVector(const Catalog &catalog, float minDistance, float maxDistance, long numBins, unsigned char *buffer) {
    std::vector<int32_t> kVector(numBins+1); // numBins = length, all elements zero
    std::vector<KVectorTriple> triples = CatalogToTripleDistances(catalog, minDistance, maxDistance);

    // sort triples in increasing order.
    std::sort(triples.begin(), triples.end(), CompareKVectorTriples);
    std::vector<float> distances;

    for (const KVectorTriple &triple : triples) {
        distances.push_back(triple.distance);
    }

    unsigned char *bufferStart = buffer;

    // index field
    SerializeKVectorIndex(distances, 0, M_PI, numBins, buffer);
    buffer += SerializeLengthKVectorIndex(numBins);

    // bulk triples field
    for (const KVectorTriple &triple : triples) {
        *(int16_t *)buffer = triple.index1;
        buffer += sizeof(int16_t);
        *(int16_t *)buffer = triple.index2;
        buffer += sizeof(int16_t);
        *(int16_t *)buffer = triple.index3;
        buffer += sizeof(int16_t);
    }

    // verify length
    assert(buffer - bufferStart == SerializeLengthTripleDistanceKVector(triples.size(), numBins));
}

PairDistanceKVectorDatabase::PairDistanceKVectorDatabase(const unsigned char *buffer)
    : index(KVectorIndex(buffer)) {
    
    // TODO: errors? (not even sure what i meant by this comment anymore)
    buffer += SerializeLengthKVectorIndex(index.NumBins());
    pairs = (const int16_t *)buffer;
}

float Clamp(float num, float low, float high) {
    return num < low ? low : num > high ? high : num;
}

TripleDistanceKVectorDatabase::TripleDistanceKVectorDatabase(const unsigned char *buffer)
    : index(KVectorIndex(buffer)) {
    
    // TODO: errors? (not even sure what i meant by this comment anymore); me neither...
    buffer += SerializeLengthKVectorIndex(index.NumBins());
    triples = (const int16_t *)buffer;
}

const int16_t *PairDistanceKVectorDatabase::FindPairsLiberal(
    float minQueryDistance, float maxQueryDistance, long *numReturnedPairs) const {

    long lowerIndex = index.QueryLiberal(minQueryDistance, maxQueryDistance, numReturnedPairs);
    return &pairs[lowerIndex * 2];
}

const int16_t *TripleDistanceKVectorDatabase::FindTriplesLiberal(
    float minQueryDistance, float maxQueryDistance, long *numReturnedTriples) const {

    long lowerIndex = index.QueryLiberal(minQueryDistance, maxQueryDistance, numReturnedTriples);
    return &triples[lowerIndex * 3];
}

long PairDistanceKVectorDatabase::NumPairs() const {
    return index.NumValues();
}

long TripleDistanceKVectorDatabase::NumTriples() const {
    return index.NumValues();
}

std::vector<float> PairDistanceKVectorDatabase::StarDistances(int16_t star, const Catalog &catalog) const {
    std::vector<float> result;
    for (int i = 0; i < NumPairs(); i++) {
        if (pairs[i*2] == star || pairs[i*2+1] == star) {
            result.push_back(AngleUnit(catalog[pairs[i*2]].spatial, catalog[pairs[i*2+1]].spatial));
        }
    }
    return result;
}

/*
std::vector<float> TripleDistanceKVectorDatabase::StarDistances(int16_t star, const Catalog &catalog) const {
    std::vector<float> result;
    for (int i = 0; i < NumTriples(); i++) {
        if (triples[i*2] == star || triples[i*2+1] == star) {
            result.push_back(AngleUnit(catalog[triples[i*2]].spatial, catalog[triples[i*2+1]].spatial));
        }
    }
    return result;
}
*/

/**
   MultiDatabase memory layout:

   | size           | name              | description                                             |
   |----------------+-------------------+---------------------------------------------------------|
   | 8*maxDatabases | table of contents | each 8-byte entry is the 4-byte magic value followed by |
   |                |                   | a 4-byte index into the bulk where that db begins       |
   | Large          | databases         | the database contents                                   |
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
    delete[] buffer;
}

}

// TODO: after creating the database, print more statistics, such as average number of pairs per
// star, stars per bin, space distribution between array and index.
