#include "star-utils.hpp"

#include <math.h>
#include <assert.h>
#include <algorithm>

namespace lost {

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

// compute min mid max values and arguments
// vals and args must be of length 3
void tripleArgs(float &min, float &mid, float &max, int &argmin, int &argmid, int &argmax, float vals[], int args[]) {
    if (vals[0] <= vals[1] && vals[0] <= vals[2]) {
        min = vals[0];
        argmin = args[0];
        if (vals[1] <= vals[2]) {
            mid = vals[1];
            argmid = args[1];
            max = vals[2];
            argmax = args[2];
        } else {
            mid = vals[2];
            argmid = args[2];
            max = vals[1];
            argmax = args[1];
        }
    } else if (vals[1] <= vals[0] && vals[1] <= vals[2]) {
        min = vals[1];
        argmin = args[1];
        if (vals[0] <= vals[2]) {
            mid = vals[0];
            argmid = args[0];
            max = vals[2];
            argmax = args[2];
        } else {
            mid = vals[2];
            argmid = args[2];
            max = vals[0];
            argmax = args[0];
        }
    } else {
        min = vals[2];
        argmin = args[2];
        if (vals[0] <= vals[1]) {
            mid = vals[0];
            argmid = args[0];
            max = vals[1];
            argmax = args[1];
        } else {
            mid = vals[1];
            argmid = args[1];
            max = vals[0];
            argmax = args[0];
        }
    }
}

void innerAngles(const Catalog &catalog, float &min, float &mid, float &max, int &mindex, int &middex, int &maxdex, int index1, int index2, int index3) {
    float a1 = Angle(catalog[index2].spatial-catalog[index1].spatial, catalog[index3].spatial-catalog[index1].spatial);
    float a2 = Angle(catalog[index2].spatial-catalog[index3].spatial, catalog[index2].spatial-catalog[index1].spatial);
    float a3 = Angle(catalog[index3].spatial-catalog[index1].spatial, catalog[index3].spatial-catalog[index2].spatial);
    float vals[] = {a1, a2, a3};
    int indices[] = {index1, index2, index3};
    tripleArgs(min, mid, max, mindex, middex, maxdex, vals, indices);
}

void focalPlaneAngles(const Stars &stars, float &min, float &mid, float &max, int &mindex, int &middex, int &maxdex, int i_index, int j_index, int k_index) {
    Star i = stars[i_index];
    Star j = stars[j_index];
    Star k = stars[k_index];
    float b1 = Angle(j.position-i.position, k.position-i.position);
    float b2 = Angle(i.position-j.position, k.position-j.position);
    float b3 = Angle(j.position-k.position, i.position-k.position);
    float vals[] = {b1, b2, b3};
    int indices[] = {i_index, j_index, k_index};
    tripleArgs(min, mid, max, mindex, middex, maxdex, vals, indices);
}

float minInnerAngle(const Catalog &catalog, int &arg, int index1, int index2, int index3) {
    float a1 = Angle(catalog[index2].spatial-catalog[index1].spatial, catalog[index3].spatial-catalog[index1].spatial);
    float a2 = Angle(catalog[index2].spatial-catalog[index3].spatial, catalog[index2].spatial-catalog[index1].spatial);
    float a3 = Angle(catalog[index3].spatial-catalog[index1].spatial, catalog[index3].spatial-catalog[index2].spatial);
    if (a1 <= a2 && a1 <= a3) {
        arg = index1;
    } else if (a2 <= a1 && a2 <= a3) {
        arg = index2;
    } else {
        arg = index3;
    }
    return std::min(std::min(a1, a2), a3);
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
