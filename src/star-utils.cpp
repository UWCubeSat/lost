#include "star-utils.hpp"

#include <math.h>
#include <assert.h>
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
 * @param maxStars If there are greater than `maxStars` many stars left after filtering by `maxMagnitude`, only the `maxStars` brightest of them are kept.
 * @return Newly constructed catalog containing only the sufficiently bright stars.
 */
Catalog NarrowCatalog(const Catalog &catalog, int maxMagnitude, int maxStars) {
    Catalog result;
    for (int i = 0; i < (int)catalog.size(); i++) {
        // Somewhat unintuitively, higher magnitude => dimmmer and v.v.
        if (catalog[i].magnitude <= maxMagnitude) {
            result.push_back(catalog[i]);
        }
    }
    if (maxStars < (int)result.size()) {
        // TODO: could change to stable sort
        std::sort(result.begin(), result.end(), CatalogStarMagnitudeCompare);
        result.resize(maxStars);
    }

    return result;
}

// TODO: maxFovDeg is possibly misleading, more accurately maxAovDeg
std::pair<Catalog, std::vector<short>> TetraPreparePattCat(const Catalog &catalog, const float maxFovDeg){

    // TODO: pass through constructor or keep constant?
    const int pattStarsPerFOV = 10;
    const int verificationStarsPerFOV = 20;
    // const float starMaxMag = 7.0;
    const float starMinSep = 0.05;
    // const float pattMaxError = 0.005;

    // Should change to maxAOV
    const float maxFOV = DegToRad(maxFovDeg);
    // const short pattSize = 4;
    // const short pattBins = 25;

    int numEntries = catalog.size();
    int keepForPattCount = 1;
    std::vector<bool> keepForPatterns(numEntries);
    std::vector<bool> keepForVerifying(numEntries);

    keepForPatterns[0] = true;
    keepForVerifying[0] = true;

    for(int ind = 1; ind < numEntries; ind++){
        Vec3 vec = catalog[ind].spatial;

        std::vector<float> angsPatterns;
        for(int j = 0; j < numEntries; j++){
            if(keepForPatterns[j]){
                // float dotProd = vec *
                float dotProd = vec * catalog[j].spatial;
                angsPatterns.push_back(dotProd);
            }
        }

        std::vector<float> angsVerifying;
        for(int j = 0; j < numEntries; j++){
            if(keepForVerifying[j]){
                float dotProd = vec * catalog[j].spatial;
                angsVerifying.push_back(dotProd);
            }
        }

        bool angsPatternsOK = true;
        for(float angPatt : angsPatterns){
            // if(angPatt >= std::cos())
            if(angPatt >= std::cos(DegToRad(starMinSep))){
                angsPatternsOK = false;
                break;
            }
        }
        if (angsPatternsOK) {
            int numStarsInFOV = 0;
            for(float angPatt : angsPatterns){
                if(angPatt > std::cos(maxFOV / 2)){
                    numStarsInFOV++;
                }
            }
            if(numStarsInFOV < pattStarsPerFOV){
                keepForPatterns[ind] = true;
                keepForVerifying[ind] = true;
                keepForPattCount++;
            }
        }

        bool angsVerifyingOK = true;
        for(float angVer: angsVerifying){
            if(angVer >= std::cos(DegToRad(starMinSep))){
                angsVerifyingOK = false;
                break;
            }
        }
        if(angsVerifyingOK){
            int numStarsInFOV = 0;
            for(float angVer : angsVerifying){
                if(angVer > std::cos(maxFOV / 2)){
                    numStarsInFOV++;
                }
            }
            if(numStarsInFOV < verificationStarsPerFOV){
                keepForVerifying[ind] = true;
            }
        }
    }

    Catalog finalCat;
    std::vector<short> pattStars;
    short cumulativeSum = -1;

    // TODO: feels like there's a smarter way to do this
    for(int i = 0; i < (int)keepForVerifying.size(); i++){
        if(keepForVerifying[i]){
            cumulativeSum++;
            finalCat.push_back(catalog[i]);
        }
        if(keepForPatterns[i]){
            pattStars.push_back(cumulativeSum);
        }
    }

    std::pair<Catalog, std::vector<short>> res{finalCat, pattStars};
    return res;
}


/// Return a pointer to the star with the given name, or NULL if not found.
const CatalogStar *FindNamedStar(const Catalog &catalog, int name) {
    for (const CatalogStar &catalogStar : catalog) {
        if (catalogStar.name == name) {
            return &catalogStar;
        }
    }
    return NULL;
}

// TODO: ok? Seems kinda stupid, anyways here's a function to get index of a CatalogStar in catalog given name
int FindCatalogStarIndex(const Catalog &catalog, int name){
    for(int i = 0; i < (int)catalog.size(); i++){
        if(catalog[i].name == name){
            return i;
        }
    }
    return -1; // no star found
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

/**
 * Deserialize a catalog star.
 * @warn The `inclMagnitude` and `inclName` parameters must be the same as passed to SerializeCatalogStar()
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

/// @sa SerializeCatalog
long SerializeLengthCatalog(const Catalog &catalog, bool inclMagnitude, bool inclName) {
    return sizeof(int16_t) + sizeof(int8_t) + catalog.size()*SerializeLengthCatalogStar(inclMagnitude, inclName);
}

/**
 * Serialize the catalog to `buffer`.
 * Use SerializeLengthCatalog() to determine how many bytes to allocate in `buffer`
 * @param inclMagnitude,inclName See SerializeCatalogStar()
 */
void SerializeCatalog(const Catalog &catalog, bool inclMagnitude, bool inclName, unsigned char *buffer) {
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

    assert(buffer-bufferStart == SerializeLengthCatalog(catalog, inclMagnitude, inclName));
}

// TODO (longer term): don't deserialize the catalog, store it on disk using the
// in-memory format so we can just copy it to memory then cast Sus, 333 hw3
// comes to mind...
// https://courses.cs.washington.edu/courses/cse333/22au/hw/hw3/hw3.html


/**
 * Deserialize a catalog.
 * @param[out] inclMagnitudeReturn,inclNameReturn Will store whether `inclMagnitude` and `inclNameReturn` were set in the corresponding SerializeCatalog() call.
 */
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
