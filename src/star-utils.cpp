#include "star-utils.hpp"

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <iostream>  // TODO: remove
#include <vector>
#include <set>

#include "serialize-helpers.hpp"

namespace lost {

const int TetraConstants::numPattStars = 4;
const int TetraConstants::numPattBins = 50;
const float TetraConstants::pattErrorRange = 0.001;
const float TetraConstants::pattMaxError = 0.001;

// brightest star first
bool CatalogStarMagnitudeCompare(const CatalogStar &a, const CatalogStar &b) {
    return a.magnitude < b.magnitude;
}

Catalog NarrowCatalog(const Catalog &catalog, int maxMagnitude, int maxStars, float minSeparation) {
    Catalog result;
    for (int i = 0; i < (int)catalog.size(); i++) {
        // Higher magnitude implies the star is dimmer
        if (catalog[i].magnitude <= maxMagnitude) {
            result.push_back(catalog[i]);
        }
    }
    // TODO: sort regardless?
    std::sort(result.begin(), result.end(), CatalogStarMagnitudeCompare);

    // TODO: bug? Why remove both i and j
    // remove stars that are too close to each other
    std::set<int> tooCloseIndices;
    // filter out stars that are too close together
    // easy enough to n^2 brute force, the catalog isn't that big
    for (int i = 0; i < (int)result.size(); i++) {
        for (int j = i + 1; j < (int)result.size(); j++) {
            if (AngleUnit(result[i].spatial.Normalize(), result[j].spatial.Normalize()) < minSeparation) {
                tooCloseIndices.insert(i);
                tooCloseIndices.insert(j);
            }
        }
    }

    // // Erase all the stars whose indices are in tooCloseIndices from the result.
    // // Loop backwards so indices don't get messed up as we iterate.
    for (auto it = tooCloseIndices.rbegin(); it != tooCloseIndices.rend(); it++) {
        result.erase(result.begin() + *it);
    }

    // and finally limit to n brightest stars
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

/// Return a pointer to the star with the given name, or NULL if not found.
Catalog::const_iterator FindNamedStar(const Catalog &catalog, int name) {
    for (auto it = catalog.cbegin(); it != catalog.cend(); ++it) {
        if (it->name == name) {
            return it;
        }
    }
    return catalog.cend();
}

/**
 * Serialize a CatalogStar into a byte buffer.
 * Use SerializeLengthCatalogStar() to determine how many bytes to allocate in `buffer`
 * @param inclMagnitude Whether to include the magnitude of the star.
 * @param inclName Whether to include the (numerical) name of the star.
 * @param buffer[out] Where the serialized star is stored.
 */
void SerializeCatalogStar(SerializeContext *ser, const CatalogStar &catalogStar, bool inclMagnitude, bool inclName) {
    SerializeVec3(ser, catalogStar.spatial);
    if (inclMagnitude) {
        SerializePrimitive<float>(ser, catalogStar.magnitude);
    }
    if (inclName) {
        // TODO: double check that bools aren't some special bitwise thing in C++
        SerializePrimitive<int16_t>(ser, catalogStar.name);
    }
}

/**
 * Deserialize a catalog star.
 * @warn The `inclMagnitude` and `inclName` parameters must be the same as passed to SerializeCatalogStar()
 * @sa SerializeCatalogStar
 */
CatalogStar DeserializeCatalogStar(DeserializeContext *des, bool inclMagnitude, bool inclName) {
    CatalogStar result;
    result.spatial = DeserializeVec3(des);
    if (inclMagnitude) {
        result.magnitude = DeserializePrimitive<float>(des);
    } else {
        result.magnitude = -424242; // TODO, what to do about special values, since there's no good ones for ints.
    }
    if (inclName) {
        result.name = DeserializePrimitive<int16_t>(des);
    } else {
        result.name = -1;
    }
    return result;
}

/**
 * Serialize the catalog to `buffer`.
 * Use SerializeLengthCatalog() to determine how many bytes to allocate in `buffer`
 * @param inclMagnitude,inclName See SerializeCatalogStar()
 */
void SerializeCatalog(SerializeContext *ser, const Catalog &catalog, bool inclMagnitude, bool inclName) {
    SerializePrimitive<int16_t>(ser, catalog.size());

    // flags
    int8_t flags = (inclMagnitude) | (inclName << 1);
    SerializePrimitive<int8_t>(ser, flags);

    for (const CatalogStar &catalogStar : catalog) {
        SerializeCatalogStar(ser, catalogStar, inclMagnitude, inclName);
    }
}

/**
 * Deserialize a catalog.
 * @param[out] inclMagnitudeReturn,inclNameReturn Will store whether `inclMagnitude` and `inclNameReturn` were set in the corresponding SerializeCatalog() call.
 */
Catalog DeserializeCatalog(DeserializeContext *des, bool *inclMagnitudeReturn, bool *inclNameReturn) {
    bool inclName, inclMagnitude;
    Catalog result;

    int16_t numStars = DeserializePrimitive<int16_t>(des);

    int8_t flags = DeserializePrimitive<int8_t>(des);
    inclMagnitude = (flags) & 1;
    inclName = (flags>>1) & 1;

    if (inclMagnitudeReturn != NULL) {
        *inclMagnitudeReturn = inclMagnitude;
    }
    if (inclNameReturn != NULL) {
        *inclNameReturn = inclName;
    }

    for (int i = 0; i < numStars; i++) {
        result.push_back(DeserializeCatalogStar(des, inclMagnitude, inclName));
    }

    return result;
}

float MagToBrightness(int mag) {
    return pow(10.0, -mag/250.0);
}

}
