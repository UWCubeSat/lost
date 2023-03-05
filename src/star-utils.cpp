#include "star-utils.hpp"

#include <assert.h>
#include <math.h>

#include <iostream> // TODO: remove
#include <vector>

#include <algorithm>

namespace lost {

// brightest star first
bool CatalogStarMagnitudeCompare(const CatalogStar &a, const CatalogStar &b) {
  return a.magnitude < b.magnitude;
}

/**
 * Remove dim stars from a catalog.
 * @param catalog Original catalog. Not mutated.
 * @param maxMagnitude Stars dimmer than this (greater numerical magnitude) are removed.
 * @param maxStars If there are greater than `maxStars` many stars left after filtering by
 * `maxMagnitude`, only the `maxStars` brightest of them are kept.
 * @return Newly constructed catalog containing only the sufficiently bright stars.
 */
Catalog NarrowCatalog(const Catalog &catalog, int maxMagnitude, int maxStars) {
  Catalog result;
  for (int i = 0; i < (int)catalog.size(); i++) {
    // Higher magnitude implies the star is dimmer
    if (catalog[i].magnitude <= maxMagnitude) {
      result.push_back(catalog[i]);
    }
  }
  // TODO: sort regardless?
  std::sort(result.begin(), result.end(), CatalogStarMagnitudeCompare);

  if (maxStars < (int)result.size()) {
    // std::sort(result.begin(), result.end(), CatalogStarMagnitudeCompare);
    result.resize(maxStars);
  }

  return result;
}

int KeyToIndex(std::vector<int> key, int binFactor, long long maxIndex) {
  const long long MAGIC_RAND = 2654435761;
  long index = 0;
  for (int i = 0; i < (int)key.size(); i++) {
    index += key[i] * std::pow(binFactor, i);
  }

  return ((index % maxIndex) * (MAGIC_RAND % maxIndex)) % maxIndex;
}

// TODO: make static or declare in header file or move
// std::pair<Catalog, std::vector<short>> TetraPreparePattCat(const Catalog &catalog,
//                                                            const float maxFovDeg) {
std::pair<std::vector<short>, std::vector<short>> TetraPreparePattCat(const Catalog &catalog,
                                                           const float maxFovDeg) {
  // TODO: these should scale based on FOV
  // Larger FOV should allow more patterns
  // 10, 20 for maxFovDeg=20ish seemed to work
  const int pattStarsPerFOV = 10;
  const int verificationStarsPerFOV = 20;
  // 25, 25
  // const int pattStarsPerFOV = 20;
  // const int verificationStarsPerFOV = 30;
  // To eliminate double stars, specify that star must be > 0.05 degrees apart
  const float starMinSep = 0.05;

  const float maxFOV = DegToRad(maxFovDeg);

  int numEntries = catalog.size();
  int keepForPattCount = 1;
  std::vector<bool> keepForPatterns(numEntries, false);
  std::vector<bool> keepForVerifying(numEntries, false);

  // NOTE: pattern set will always be a subset of verification set
  // We technically don't even need the verification set unless we plan
  // to calculate probability of mismatch - (which we probably should)

  // Definitely keep the first star
  keepForPatterns[0] = true;
  keepForVerifying[0] = true;

  for (int i = 1; i < numEntries; i++) {
    // vec representing new star
    Vec3 vec = catalog[i].spatial;

    bool angsPatternsOK = true;
    int numPattStarsInFov = 0;

    // We should test each new star's angular distance to all stars
    // we've already selected to be kept for pattern construction
    // Stop early if:
    // a) double star: angle < min separation allowed
    // b) Number of stars in region maxFov/2 >= pattStarsPerFOV
    // std::vector<float> angsToPatterns;
    for (int j = 0; j < i; j++) {
      if (keepForPatterns[j]) {
        float dotProd = vec * catalog[j].spatial;
        // angsToPatterns.push_back(dotProd);
        if (dotProd >= std::cos(DegToRad(starMinSep))) {
          angsPatternsOK = false;
          break;
        }
        // If angle between new star i and old star j is less than maxFov/2, OK
        if (dotProd > std::cos(maxFOV / 2)) {
          numPattStarsInFov++;
          if (numPattStarsInFov >= pattStarsPerFOV) {
            angsPatternsOK = false;
            break;
          }
        }
      }
    }

    if (angsPatternsOK) {
      keepForPatterns[i] = true;
      keepForVerifying[i] = true;
      keepForPattCount++;
      continue;
    }

    bool angsVerifyingOK = true;
    int numVerStarsInFov = 0;

    // Same thing here, we should test each new star's angular distance to all
    // stars we've already selected to be kept for verification
    // std::vector<float> angsVerifying;
    for (int j = 0; j < i; j++) {
      if (keepForVerifying[j]) {
        float dotProd = vec * catalog[j].spatial;
        if (dotProd >= std::cos(DegToRad(starMinSep))) {
          angsVerifyingOK = false;
          break;
        }
        if (dotProd > std::cos(maxFOV / 2)) {
          numVerStarsInFov++;
          // Not really a bug, more like mistake on my part
          // verificationStarsPerFOV is still low, so when FOV is big,
          // very few patterns are generated
          if (numVerStarsInFov >= verificationStarsPerFOV) {
            angsVerifyingOK = false;
            break;
          }
        }
      }
    }

    if (angsVerifyingOK) {
      keepForVerifying[i] = true;
    }
  }

  // Catalog finalCat;
  std::vector<short> finalCatIndices;
  std::vector<short> pattStars;

  // finalCat is the final version of the star table
  for (int i = 0; i < (int)keepForVerifying.size(); i++) {
    if (keepForVerifying[i]) {
      // finalCat.push_back(catalog[i]);
      finalCatIndices.push_back(i);
    }
  }

  // Pretty clever way of finding which stars in the FINAL star table
  // should be used for pattern construction later in Tetra's database generation step
  short cumulativeSum = -1;
  for (int i = 0; i < (int)keepForVerifying.size(); i++) {
    if (keepForVerifying[i]) {
      cumulativeSum++;
    }
    if (keepForPatterns[i]) {
      pattStars.push_back(cumulativeSum);
    }
  }
  return std::pair<std::vector<short>, std::vector<short>>{finalCatIndices, pattStars};
  // return std::pair<Catalog, std::vector<short>>{finalCat, pattStars};
}

/// Return a pointer to the star with the given name, or NULL if not found.
Catalog::const_iterator FindNamedStar(const Catalog &catalog, int name) {
    for (auto it = catalog.cbegin(); it != catalog.cend(); ++it) {
        if (it->name == name) {
            return it;
        }
    }
    return catalog.cend();
}

// TODO: ok? Seems kinda stupid, anyways here's a function to get index of a CatalogStar in catalog
// given name
int FindCatalogStarIndex(const Catalog &catalog, int name) {
  for (int i = 0; i < (int)catalog.size(); i++) {
    if (catalog[i].name == name) {
      return i;
    }
  }
  return -1;  // no star found
}

/// @sa SerializeCatalogStar
long SerializeLengthCatalogStar(bool inclMagnitude, bool inclName) {
  long starSize = SerializeLengthVec3();
  if (inclMagnitude) {
    starSize += sizeof(float);
  }
  if (inclName) {
    starSize += sizeof(int16_t);
  }
  return starSize;
}

/**
 * Serialize a CatalogStar into a byte buffer.
 * Use SerializeLengthCatalogStar() to determine how many bytes to allocate in `buffer`
 * @param inclMagnitude Whether to include the magnitude of the star.
 * @param inclName Whether to include the (numerical) name of the star.
 * @param buffer[out] Where the serialized star is stored.
 */
// TODO: make inclusion of name/magnitude true by default?
// Actually why give the option in the first place, algos like Tetra need this to work
void SerializeCatalogStar(const CatalogStar &catalogStar, bool inclMagnitude, bool inclName,
                          unsigned char *buffer) {
  SerializeVec3(catalogStar.spatial, buffer);
  buffer += SerializeLengthVec3();
  if (inclMagnitude) {
    *(float *)buffer = catalogStar.magnitude;
    buffer += sizeof(float);
  }
  if (inclName) {
    // TODO: double check that bools aren't some special bitwise thing in C++
    *(int16_t *)buffer = catalogStar.name;
    buffer += sizeof(int16_t);
  }
}

/**
 * Deserialize a catalog star.
 * @warn The `inclMagnitude` and `inclName` parameters must be the same as passed to
 * SerializeCatalogStar()
 * @sa SerializeCatalogStar
 */
CatalogStar DeserializeCatalogStar(const unsigned char *buffer, bool inclMagnitude, bool inclName) {
  CatalogStar result;
  result.spatial = DeserializeVec3(buffer);
  buffer += SerializeLengthVec3();
  if (inclMagnitude) {
    result.magnitude = *(float *)buffer;
    buffer += sizeof(float);
  } else {
    result.magnitude =
        -424242;  // TODO, what to do about special values, since there's no good ones for ints.
  }
  if (inclName) {
    result.name = *(int16_t *)buffer;
    buffer += sizeof(int16_t);
  } else {
    result.name = -1;
  }
  return result;
}

/// @sa SerializeCatalog
long SerializeLengthCatalog(const Catalog &catalog, bool inclMagnitude, bool inclName) {
  return sizeof(int16_t) + sizeof(int8_t) +
         catalog.size() * SerializeLengthCatalogStar(inclMagnitude, inclName);
}

/**
 * Serialize the catalog to `buffer`.
 * Use SerializeLengthCatalog() to determine how many bytes to allocate in `buffer`
 * @param inclMagnitude,inclName See SerializeCatalogStar()
 */
void SerializeCatalog(const Catalog &catalog, bool inclMagnitude, bool inclName,
                      unsigned char *buffer) {
  unsigned char *bufferStart = buffer;

  // size
  *(int16_t *)buffer = catalog.size();
  buffer += sizeof(int16_t);

  // flags
  int8_t flags = (inclMagnitude) | (inclName << 1);
  *(int8_t *)buffer = flags;
  buffer += sizeof(int8_t);

  long catalogStarLength = SerializeLengthCatalogStar(inclMagnitude, inclName);
  for (const CatalogStar &catalogStar : catalog) {
    SerializeCatalogStar(catalogStar, inclMagnitude, inclName, buffer);
    buffer += catalogStarLength;
  }

  assert(buffer - bufferStart == SerializeLengthCatalog(catalog, inclMagnitude, inclName));
}

// TODO (longer term): don't deserialize the catalog, store it on disk using the
// in-memory format so we can just copy it to memory then cast Sus, 333 hw3
// comes to mind...
// https://courses.cs.washington.edu/courses/cse333/22au/hw/hw3/hw3.html

/**
 * Deserialize a catalog.
 * @param[out] inclMagnitudeReturn,inclNameReturn Will store whether `inclMagnitude` and
 * `inclNameReturn` were set in the corresponding SerializeCatalog() call.
 */
Catalog DeserializeCatalog(const unsigned char *buffer, bool *inclMagnitudeReturn,
                           bool *inclNameReturn) {
  bool inclName, inclMagnitude;
  Catalog result;

  int16_t numStars = *(int16_t *)buffer;
  buffer += sizeof(int16_t);

  int8_t flags = *(int8_t *)buffer;
  inclMagnitude = (flags)&1;
  inclName = (flags >> 1) & 1;
  if (inclMagnitudeReturn != NULL) {
    *inclMagnitudeReturn = inclMagnitude;
  }
  if (inclNameReturn != NULL) {
    *inclNameReturn = inclName;
  }
  buffer += sizeof(int8_t);

  int catalogStarLength = SerializeLengthCatalogStar(inclMagnitude, inclName);
  for (int i = 0; i < numStars; i++) {
    result.push_back(DeserializeCatalogStar(buffer, inclMagnitude, inclName));
    buffer += catalogStarLength;
  }

  return result;
}

}  // namespace lost
