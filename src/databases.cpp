#include "databases.hpp"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

// TODO: check which ones are useless
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

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
  return 4 + sizeof(float) + sizeof(float) + 4 + 4 * (numBins + 1);
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
 * @param numBins the number of "bins" the KVector should use. A higher number makes query results
 * "tighter" but takes up more disk space. Usually should be set somewhat smaller than (max-min)
 * divided by the "width" of the typical query.
 * @param buffer[out] index is written here.
 */
void SerializeKVectorIndex(const std::vector<float> &values, float min, float max, long numBins,
                           unsigned char *buffer) {
  std::vector<int32_t> kVector(numBins + 1);  // numBins = length, all elements zero
  float binWidth = (max - min) / numBins;

  // generate the k-vector part
  // Idea: When we find the first star that's across any bin boundary, we want to update all the
  // newly sealed bins
  long lastBin = 0;  // first bin the last star belonged to
  for (int32_t i = 0; i < (int32_t)values.size(); i++) {
    if (i > 0) {
      assert(values[i] >= values[i - 1]);
    }
    assert(values[i] >= min);
    assert(values[i] <= max);
    long thisBin = (long)ceil((values[i] - min) / binWidth);  // first bin we belong to
    assert(thisBin >= 0);
    assert(thisBin <=
           numBins);  // thisBin == numBins is acceptable since kvector length == numBins + 1
    // if thisBin > lastBin, then no more stars can be added to those bins between thisBin and
    // lastBin, so set them.
    for (long bin = lastBin; bin < thisBin; bin++) {
      kVector[bin] = i;  // our index is the number of pairs with shorter distance
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
long KVectorIndex::QueryLiberal(float minQueryDistance, float maxQueryDistance,
                                long *upperIndex) const {
  assert(maxQueryDistance > minQueryDistance);
  if (maxQueryDistance >= max) {
    maxQueryDistance = max - 0.00001;  // TODO: better way to avoid hitting the bottom bin
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
  int lowerIndex = bins[lowerBin - 1];
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
 pair K-vector database layout. The kvector appears before the bulk pair data because it contains
 the number of pairs, which is necessary to read the bulk pair data.

     | size (bytes)             | name         | description |
     |--------------------------+--------------+-------------------------------------------------------------|
     | sizeof kvectorIndex      | kVectorIndex | Serialized KVector index | |
 2*sizeof(int16)*numPairs | pairs        | Bulk pair data |
 */
std::vector<KVectorPair> CatalogToPairDistances(const Catalog &catalog, float minDistance,
                                                float maxDistance) {
  std::vector<KVectorPair> result;
  for (int16_t i = 0; i < (int16_t)catalog.size(); i++) {
    for (int16_t k = i + 1; k < (int16_t)catalog.size(); k++) {
      KVectorPair pair = {i, k, AngleUnit(catalog[i].spatial, catalog[k].spatial)};
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
  return SerializeLengthKVectorIndex(numBins) + 2 * sizeof(int16_t) * numPairs;
}

/// Number of bytes that a serialized KVectorDatabase will take up
long SerializeLengthPairDistanceKVector(const Catalog &catalog, float minDistance,
                                        float maxDistance, long numBins) {
  return SerializeLengthPairDistanceKVector(
      CatalogToPairDistances(catalog, minDistance, maxDistance).size(), numBins);
}

/**
 * Serialize a pair-distance KVector into buffer.
 * Use SerializeLengthPairDistanceKVector to determine how large the buffer needs to be. See command
 * line documentation for other options.
 */
void SerializePairDistanceKVector(const Catalog &catalog, float minDistance, float maxDistance,
                                  long numBins, unsigned char *buffer) {
  std::vector<int32_t> kVector(numBins + 1);  // numBins = length, all elements zero
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

/// Create the database from a serialized buffer.
PairDistanceKVectorDatabase::PairDistanceKVectorDatabase(const unsigned char *buffer)
    : index(KVectorIndex(buffer)) {
  // TODO: errors? (not even sure what i meant by this comment anymore)
  buffer += SerializeLengthKVectorIndex(index.NumBins());
  pairs = (const int16_t *)buffer;
}

/// Return the value in the range [low,high] which is closest to num
float Clamp(float num, float low, float high) { return num < low ? low : num > high ? high : num; }

/**
 * Return at least all the star pairs whose inter-star distance is between min and max
 * @param end[out] Is set to an "off-the-end" pointer, one past the last pair being returned by the
 * query.
 * @return A pointer to the start of the matched pairs. Each pair is stored as simply two 16-bit
 * integers, each of which is a catalog index. (you must increment the pointer twice to get to the
 * next pair).
 */
const int16_t *PairDistanceKVectorDatabase::FindPairsLiberal(float minQueryDistance,
                                                             float maxQueryDistance,
                                                             const int16_t **end) const {
  long upperIndex;
  long lowerIndex = index.QueryLiberal(minQueryDistance, maxQueryDistance, &upperIndex);
  *end = &pairs[upperIndex * 2];
  return &pairs[lowerIndex * 2];
}

/// Number of star pairs stored in the database
long PairDistanceKVectorDatabase::NumPairs() const { return index.NumValues(); }

/// Return the distances from the given star to each star it's paired with in the database (for
/// debugging).
std::vector<float> PairDistanceKVectorDatabase::StarDistances(int16_t star,
                                                              const Catalog &catalog) const {
  std::vector<float> result;
  for (int i = 0; i < NumPairs(); i++) {
    if (pairs[i * 2] == star || pairs[i * 2 + 1] == star) {
      result.push_back(AngleUnit(catalog[pairs[i * 2]].spatial, catalog[pairs[i * 2 + 1]].spatial));
    }
  }
  return result;
}

///////////////////// Tetra database //////////////////////

int KeyToIndex(std::vector<int> key, int binFactor, int maxIndex) {
  const long long MAGIC_RAND = 2654435761;
  long index = 0;
  for (int i = 0; i < (int)key.size(); i++) {
    index += key[i] * std::pow(binFactor, i);
  }
  return (index * MAGIC_RAND) % maxIndex;
}

typedef std::array<short, 3> ShortVec3;
typedef std::array<float, 3> FloatVec3;

long SerializeTetraDatabase(const Catalog &catalog, float maxFov, unsigned char *buffer,
                            const std::vector<short> &pattStars, bool ser) {
  maxFov = DegToRad(maxFov);
  const short pattBins = 25;
  const int tempBins = 4;


  // TODO: unorderd_map might be better
  std::map<Vec3, std::vector<short>> tempCoarseSkyMap;

  //  // std::map<ShortVec3, std::set<short>> tempCoarseSkyMap;

  for (const short starID : pattStars) {
    // ShortVec3 hash = {
    //   (short)((catalog[starID].spatial.x + 1) * tempBins),
    //   (short)((catalog[starID].spatial.y + 1) * tempBins),
    //   (short)((catalog[starID].spatial.z + 1) * tempBins)
    // };
    Vec3 v(catalog[starID].spatial);

    Vec3 hash{short((v.x + 1) * tempBins), short((v.y + 1) * tempBins),
              short((v.z + 1) * tempBins)};

    // Vec3 hash{(short)((catalog[starID].spatial.x + 1) * tempBins),
    //                   (short)((catalog[starID].spatial.y + 1) * tempBins),
    //                   (short)((catalog[starID].spatial.z + 1) * tempBins)};
    tempCoarseSkyMap[hash].push_back(starID);
  }

  // Return a list of stars that are nearby
  auto tempGetNearbyStars = [&tempCoarseSkyMap, &catalog](const Vec3 &vec, float radius) {
    std::vector<float> components{vec.x, vec.y, vec.z};

    std::vector<std::vector<int>> hcSpace;
    for (float x : components) {
      std::vector<int> range;
      int lo = int((x + 1 - radius) * tempBins);
      lo = std::max(lo, 0);
      int hi = int((x + 1 + radius) * tempBins);
      hi = std::min(hi + 1, 2 * tempBins);
      range.push_back(lo);
      range.push_back(hi);
      hcSpace.push_back(range);
    }

    // hashcode space has 3 ranges, one for each of [x, y, z]

    std::vector<short> nearbyStarIDs;

    for (int a = hcSpace[0][0]; a < hcSpace[0][1]; a++) {
      for (int b = hcSpace[1][0]; b < hcSpace[1][1]; b++) {
        for (int c = hcSpace[2][0]; c < hcSpace[2][1]; c++) {
          //   ShortVec3 code{short(a), short(b), short(c)};
          Vec3 code{short(a), short(b), short(c)};

          // For each star j in partition with key=code,
          // see if our star and j have angle < radius. If so, they are nearby
          for (short starID : tempCoarseSkyMap[code]) {

            float dotProd = vec * catalog[starID].spatial;
            if(dotProd > std::cos(radius)){
                nearbyStarIDs.push_back(starID);
            }

            // float dotProd = vec[0] * catalog[starID].spatial.x +
            //                 vec[1] * catalog[starID].spatial.y + vec[2] * catalog[starID].spatial.z;
            // if (dotProd > std::cos(radius)) {
            //   nearbyStarIDs.push_back(starID);
            // }
          }
        }
      }
    }

    return nearbyStarIDs;
  };

  const int pattSize = 4;
  typedef std::array<short, pattSize> Pattern;
  std::vector<Pattern> pattList;
  Pattern patt{0, 0, 0, 0};

  for (int ind = 0; ind < (int)pattStars.size(); ind++) {
    short firstStarID = pattStars[ind];
    patt[0] = firstStarID;

    FloatVec3 firstStarVec = {catalog[firstStarID].spatial.x, catalog[firstStarID].spatial.y,
                              catalog[firstStarID].spatial.z};

    ShortVec3 hashCode = {
        short((firstStarVec[0] + 1) * (float)tempBins),
        short((firstStarVec[1] + 1) * (float)tempBins),
        short((firstStarVec[2] + 1) * (float)tempBins),
    };

    tempCoarseSkyMap[hashCode].erase(patt[0]);
    auto nearbyStars = tempGetNearbyStars(firstStarVec, maxFov, ind);
    const int n = (int)nearbyStars.size();

    for (int i = 0; i < n; i++) {
      for (int j = i + 1; j < n; j++) {
        for (int k = j + 1; k < n; k++) {
          patt[1] = nearbyStars[i];
          patt[2] = nearbyStars[j];
          patt[3] = nearbyStars[k];

          bool pattFits = true;
          for (int pair1 = 0; pair1 < pattSize; pair1++) {
            for (int pair2 = pair1 + 1; pair2 < pattSize; pair2++) {
              float dotProd = catalog[patt[pair1]].spatial.x * catalog[patt[pair2]].spatial.x +
                              catalog[patt[pair1]].spatial.y * catalog[patt[pair2]].spatial.y +
                              catalog[patt[pair1]].spatial.z * catalog[patt[pair2]].spatial.z;
              if (dotProd <= std::cos(maxFov)) {
                pattFits = false;
                break;
              }
            }
          }
          // end of for loop
          if (pattFits) {
            pattList.push_back(patt);
          }
        }
      }
    }
  }

  std::cout << "Found " << pattList.size() << " patterns" << std::endl;

  int catalogLength = 2 * (int)pattList.size();
  std::vector<Pattern> pattCatalog(catalogLength);

  //
  for (Pattern patt : pattList) {
    std::vector<float> pattEdgeLengths;
    for (int i = 0; i < pattSize; i++) {
      CatalogStar star1 = catalog[patt[i]];
      for (int j = i + 1; j < pattSize; j++) {
        // calculate distance between vectors
        CatalogStar star2 = catalog[patt[j]];
        float edgeLen = std::sqrt(std::pow(star2.spatial.x - star1.spatial.x, 2) +
                                  std::pow(star2.spatial.y - star1.spatial.y, 2) +
                                  std::pow(star2.spatial.z - star1.spatial.z, 2));
        pattEdgeLengths.push_back(edgeLen);
      }
    }
    std::sort(pattEdgeLengths.begin(), pattEdgeLengths.end());

    float pattLargestEdge = pattEdgeLengths[(int)(pattEdgeLengths.size() - 1)];
    std::vector<float> pattEdgeRatios;
    for (int i = 0; i < (int)pattEdgeLengths.size() - 1;
         i++) {  // size()-1, since we ignore the largest edge
      pattEdgeRatios.push_back(pattEdgeLengths[i] / pattLargestEdge);
    }

    std::vector<int> key;
    for (float edgeRatio : pattEdgeRatios) {
      key.push_back(int(edgeRatio * pattBins));
    }

    int hashIndex = KeyToIndex(key, pattBins, catalogLength);
    long long offset = 0;
    while (true) {
      int index = int(hashIndex + std::pow(offset, 2)) % catalogLength;
      offset++;
      if (pattCatalog[index][0] == 0) {
        pattCatalog[index] = patt;
        break;
      }
    }
  }

  // Done with everything, write to buffer

  if (ser) {
    for (Pattern patt : pattCatalog) {
      for (int i = 0; i < pattSize; i++) {
        *((short *)buffer) = patt[i];
        buffer += sizeof(short);
      }
    }
    return -1;
  } else {
    return sizeof(short) * pattCatalog.size() * pattSize;
  }
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
 * @return Returns a pointer to the start of the database type indicated by the magic value, null if
 * not found
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

  return buffer + kMultiDatabaseTocLength + databaseIndex;
}

/**
 * Add a database to a MultiDatabase
 * @param magicValue A value unique to this type of database which is used to extract it out of the
 * database later.
 * @param length The number of bytes to allocate for this database.
 * @return Pointer to the start of the space allocated for said database. Return null if full (too
 * many databases).
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

  buffer = (unsigned char *)realloc(buffer, kMultiDatabaseTocLength + bulkLength + length);
  // just past the end of the last database
  unsigned char *result = buffer + kMultiDatabaseTocLength + bulkLength;
  bulkLength += length;
  return result;
}

MultiDatabaseBuilder::~MultiDatabaseBuilder() { free(buffer); }

}  // namespace lost

// TODO: after creating the database, print more statistics, such as average number of pairs per
// star, stars per bin, space distribution between array and index.
