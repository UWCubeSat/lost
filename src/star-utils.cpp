#include "star-utils.hpp"

#include <math.h>
#include <assert.h>
#include <algorithm>

namespace lost {

float StarDistancePixels(Star one, Star two) {
    float distX = one.x - two.x;
    float distY = one.y - two.y;
    return sqrt(distX*distX + distY*distY);
}

// brightest star first
bool CatalogStarMagnitudeCompare(const CatalogStar &a, const CatalogStar &b) {
    return a.magnitude < b.magnitude;
}

// TODO: this function is about to cause a horribly difficult to debug error, because the catalog
// used in the database may not equal the catalog used by star id or something.

// TODO: the maxStars is worthless, it doesn't get the brightest stars
Catalog NarrowCatalog(const Catalog &catalog, int maxMagnitude, int maxStars) {
    Catalog result;
    for (int i = 0; i < (int)catalog.size(); i++) {
        if (catalog[i].magnitude <= maxMagnitude) {
            result.push_back(catalog[i]);
        }
    }
    if (maxStars < (int)result.size()) {
        std::sort(result.begin(), result.end(), CatalogStarMagnitudeCompare);
        result.resize(maxStars);
    }

    return result;
}

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

void SerializeCatalogStar(const CatalogStar &catalogStar, bool inclMagnitude, bool inclName, unsigned char *buffer) {
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

CatalogStar DeserializeCatalogStar(const unsigned char *buffer, bool inclMagnitude, bool inclName) {
    CatalogStar result;
    result.spatial = DeserializeVec3(buffer);
    buffer += SerializeLengthVec3();
    if (inclMagnitude) {
        result.magnitude = *(float *)buffer;
        buffer += sizeof(float);
    } else {
        result.magnitude = -424242; // TODO, what to do about special values, since there's no good ones for ints.
    }
    if (inclName) {
        result.name = *(int16_t *)buffer;
        buffer += sizeof(int16_t);
    } else {
        result.name = -1;
    }
    return result;
}

long SerializeLengthCatalog(const Catalog &catalog, bool inclMagnitude, bool inclName) {
    return sizeof(int16_t) + sizeof(int8_t) + catalog.size()*SerializeLengthCatalogStar(inclMagnitude, inclName);
}

void SerializeCatalog(const Catalog &catalog, bool inclMagnitude, bool inclName, unsigned char *buffer) {
    unsigned char *bufferStart = buffer;

    // size
    *(int16_t *)buffer = catalog.size();
    buffer += sizeof(int16_t);

    // flags
    int8_t flags = (inclMagnitude) | (inclName<<1);
    *(int8_t *)buffer = flags;
    buffer += sizeof(int8_t);

    long catalogStarLength = SerializeLengthCatalogStar(inclMagnitude, inclName);
    for (const CatalogStar &catalogStar : catalog) {
        SerializeCatalogStar(catalogStar, inclMagnitude, inclName, buffer);
        buffer += catalogStarLength;
    }

    assert(buffer-bufferStart == SerializeLengthCatalog(catalog, inclMagnitude, inclName));
}

// TODO (longer term): don't deserialize the catalog, store it on disk using the in-memory format so
// we can just copy it to memory then cast
Catalog DeserializeCatalog(const unsigned char *buffer, bool *inclMagnitudeReturn, bool *inclNameReturn) {
    bool inclName, inclMagnitude;
    Catalog result;

    int16_t numStars = *(int16_t *)buffer;
    buffer += sizeof(int16_t);

    int8_t flags = *(int8_t *)buffer;
    inclMagnitude = (flags) & 1;
    inclName = (flags>>1) & 1;
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

}
