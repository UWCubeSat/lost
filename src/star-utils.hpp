#ifndef STAR_UTILS_H
#define STAR_UTILS_H

#include <vector>

#include "attitude-utils.hpp"

namespace lost {

/**
 * @brief
 * @details
 */
class CatalogStar {
public:
    /// @brief
    CatalogStar() = default;

    /**
     * @brief
     * @param raj2000
     * @param dej2000
     * @param magnitude
     * @param name
     */
    CatalogStar(float raj2000, float dej2000, int magnitude, int name) :
        spatial(SphericalToSpatial(raj2000, dej2000)), magnitude(magnitude), name(name) {}

    /**
     * @brief
     * @param spatial
     * @param magnitude
     * @param name
     */
    CatalogStar(Vec3 spatial, int magnitude, int name) :
        spatial(spatial), magnitude(magnitude), name(name) {}

    /// @brief
    Vec3 spatial;
    /**
     * @brief
     * @note *10^-2
     */
    int magnitude;

    /// @brief
    int name;
};

/**
 * @brief
 * @details
 */
class Star {
public:
    /**
     * @brief
     * @param x
     * @param y
     * @param radiusX
     * @param radiusY
     * @param magnitude
     */
    Star(float x, float y, float radiusX, float radiusY, int magnitude) :
        position({x, y}), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude) {};

    /**
     * @brief
     * @param x
     * @param y
     * @param radiusX
     */
    Star(float x, float y, float radiusX) : Star(x, y, radiusX, radiusX, 0) {};

    /**
     * @brief
     */
    Star() : Star(0.0, 0.0, 0.0) {};

    /// @brief
    Vec2 position;

    /// @brief
    float radiusX;

    /**
     * @brief
     * @note if omitted, but x is present, assume circular.
     */
    float radiusY;

    /**
     * @brief
     * @details some relative number
     */
    int magnitude;
    // eccentricity?
};

/**
 * @brief
 * @details
 */
class StarIdentifier {
public:
    /**
     * @brief
     * @param starIndex
     * @param catalogIndex
     * @param weight
     */
    StarIdentifier(int starIndex, int catalogIndex, int weight)
        : starIndex(starIndex), catalogIndex(catalogIndex), weight(weight) { };

    /**
     * @brief
     * @param starIndex
     * @param catalogIndex
     */
    StarIdentifier(int starIndex, int catalogIndex)
        : StarIdentifier(starIndex, catalogIndex, 1.0f) { };

    /// @brief
    int starIndex;

    /// @brief
    int catalogIndex;

    /// @brief
    float weight;
};

typedef std::vector<CatalogStar> Catalog;
typedef std::vector<Star> Stars;
typedef std::vector<StarIdentifier> StarIdentifiers;

long SerializeLengthCatalog(const Catalog &, bool inclMagnitude, bool inclName);
void SerializeCatalog(const Catalog &, bool inclMagnitude, bool inclName, unsigned char *buffer);
// sets magnited and name to whether the catalog in the database contained magnitude and name
Catalog DeserializeCatalog(const unsigned char *buffer, bool *inclMagnitudeReturn, bool *inclNameReturn);
const CatalogStar *findNamedStar(const Catalog &, int name);

// TODO: make maxStars work right, need to sort by magnitude before filter! maxMagnitude is 10^-2
// (so 523 = 5.23)
Catalog NarrowCatalog(const Catalog &, int maxMagnitude, int maxStars);

}

#endif
