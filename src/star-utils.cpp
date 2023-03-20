#include "star-utils.hpp"

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <iostream>  // TODO: remove
#include <vector>
#include <set>

namespace lost {

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
    // std::set<int> tooCloseIndices;
    // // filter out stars that are too close together
    // // easy enough to n^2 brute force, the catalog isn't that big
    // for (int i = 0; i < (int)result.size(); i++) {
    //     for (int j = i + 1; j < (int)result.size(); j++) {
    //         if (AngleUnit(result[i].spatial.Normalize(), result[j].spatial.Normalize()) < minSeparation) {
    //             tooCloseIndices.insert(i);
    //             tooCloseIndices.insert(j);
    //         }
    //     }
    // }

    // // Erase all the stars whose indices are in tooCloseIndices from the result.
    // // Loop backwards so indices don't get messed up as we iterate.
    // for (auto it = tooCloseIndices.rbegin(); it != tooCloseIndices.rend(); it++) {
    //     result.erase(result.begin() + *it);
    // }

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

/// Get index of a CatalogStar in catalog given name
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
