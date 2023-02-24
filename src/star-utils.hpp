#ifndef STAR_UTILS_H
#define STAR_UTILS_H

#include <vector>

#include "attitude-utils.hpp"

namespace lost {

/// A star from the Bright Star Catalog.
class CatalogStar {
public:
    CatalogStar() = default;

    /**
     * Create a CatalogStar using its celestial/spherical coordinates.
     * @param raj2000,dej2000 The right ascension and declination of the star in the year 2000
     * @param magnitude See ::magnitude
     * @param name See ::name
     */
    CatalogStar(float raj2000, float dej2000, int magnitude, int name) :
        spatial(SphericalToSpatial(raj2000, dej2000)), magnitude(magnitude), name(name) {}

    CatalogStar(Vec3 spatial, int magnitude, int name) :
        spatial(spatial), magnitude(magnitude), name(name) {}

    /**
     * The point on the unit sphere where the star lies.
     * The earth is at the center of the sphere. A star with right ascension and declination both equal to zero will lie at (1,0,0).
     */
    Vec3 spatial;
    /**
     * The magnitude of the star, with the decimal point shifted two places right.
     * I.e., multiply this by 10<sup>-2</sup> to get the magnitude.
     */
    int magnitude;
    /**
     * A unique number which unambiguously identifies the catalog star.
     * The catalog is an array, so in algorithms we usually refer to catalog stars by their index in the catalog. However, if the catalog is filtered (say, to remove stars that are too dim or too close to each other), the indexes change, which makes everything really hard to debug because indexes change based on how the catalog was filtered. That's what names are for: They're fixed and don't change.
     */
    int name;
};

/**
 * A "centroid" detected in an image.
 * This represents a star from an image which has not necessarily been identified.
 */
class Star {
public:
    Star(float x, float y, float radiusX, float radiusY, int magnitude) :
        position({x, y}), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude) {};

    /// Convenience constructor that sets Star.radiusY = radiusX and Star.magnitude = 0
    Star(float x, float y, float radiusX) : Star(x, y, radiusX, radiusX, 0) {};

    /// Create a zeroed-out star. Fields should be set immediately after construction.
    Star() : Star(0.0, 0.0, 0.0) {};

    /// The (x,y) pixel coordinates in the image (top left is 0,0)
    Vec2 position;
    /// Approximate horizontal radius of the bright area in pixels.
    float radiusX;
    /// Approximate vertical radius of the bright area in pixels.
    float radiusY;
    /**
     * A relative measure of magnitude of the star. Larger is brighter.
     * It's impossible to tell the true magnitude of the star from the image, without really good camera calibration. Anyway, this field is not meant to correspond to the usual measurement of magnitude. Instead, it's just some measure of brightness which may be specific to the centroiding algorithm. For example, it might be the total number of bright pixels in the star.
     */
    int magnitude;
    // eccentricity?
};

/**
 * Records that a certain Star (detected in the image) corresponds to a certain CatalogStar.
 * Only stores indexes into an array of Star objects and an array of CatalogStar objects, which are stored elsewhere.
 */
class StarIdentifier {
public:
    StarIdentifier(int starIndex, int catalogIndex, int weight)
        : starIndex(starIndex), catalogIndex(catalogIndex), weight(weight) { };
    /// Sets StarIdentifier.weight = 1
    StarIdentifier(int starIndex, int catalogIndex)
        : StarIdentifier(starIndex, catalogIndex, 1.0f) { };

    // does not check weight
    bool operator==(const StarIdentifier& other) const {
        return starIndex == other.starIndex &&
            catalogIndex == other.catalogIndex;
    }

    /// An index into an array of Star objects.
    int starIndex;
    /// An index into an array of CatalogStar objects.
    int catalogIndex;
    /// A weight indicating the confidence of this idenification. Often just 1.
    float weight;
};

typedef std::vector<CatalogStar> Catalog;
typedef std::vector<Star> Stars;
typedef std::vector<StarIdentifier> StarIdentifiers;

long SerializeLengthCatalog(const Catalog &, bool inclMagnitude, bool inclName);
void SerializeCatalog(const Catalog &, bool inclMagnitude, bool inclName, unsigned char *buffer);
// sets magnited and name to whether the catalog in the database contained magnitude and name
Catalog DeserializeCatalog(const unsigned char *buffer, bool *inclMagnitudeReturn, bool *inclNameReturn);
Catalog::const_iterator FindNamedStar(const Catalog &catalog, int name);

// TODO: make maxStars work right, need to sort by magnitude before filter! maxMagnitude is 10^-2
// (so 523 = 5.23)
Catalog NarrowCatalog(const Catalog &, int maxMagnitude, int maxStars);

}

#endif
