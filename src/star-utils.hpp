#ifndef STAR_UTILS_H
#define STAR_UTILS_H

#include <vector>

#include "attitude-utils.hpp"

namespace lost {

class CatalogStar {
public:
    CatalogStar(float raj2000, float dej2000, int magnitude, bool weird, int name) :
        spatial(SphericalToSpatial(raj2000, dej2000)), magnitude(magnitude), weird(weird), name(name) { }

    Vec3 spatial;
    int  magnitude;         // *10^-2
    bool weird;             // nonzero for binary, etc
    int name;
};

class Star {
public:
    Star(float x, float y, float radiusX, float radiusY, int magnitude) :
        x(x), y(y), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude) { };
    Star(float x, float y, float radiusX) : Star(x, y, radiusX, radiusX, 0) { };
    Star() : Star(0.0, 0.0, 0.0) { };

    // TODO: store a Vec2 instead and then simply point to that Vec2 when passing Star to
    // CameraToSpatial?
    float x; // pixels*10^6
    float y; // pixels*10^6
    float radiusX;
    float radiusY;  // if omitted, but x is present, assume circular.
    int   magnitude; // some relative number
    // eccentricity?
};

class StarIdentifier {
public:
    StarIdentifier(int starIndex, int catalogIndex, int weight)
        : starIndex(starIndex), catalogIndex(catalogIndex), weight(weight) { };
    StarIdentifier(int starIndex, int catalogIndex)
        : StarIdentifier(starIndex, catalogIndex, 1.0f) { };

    int starIndex;
    int catalogIndex;
    float weight;
};

typedef std::vector<CatalogStar> Catalog;
typedef std::vector<Star> Stars;
typedef std::vector<StarIdentifier> StarIdentifiers;

float StarDistancePixels(Star one, Star two);

// TODO: make maxStars work right, need to sort by magnitude before filter! maxMagnitude is 10^-2
// (so 523 = 5.23)
Catalog NarrowCatalog(const Catalog &, int maxMagnitude, int maxStars);

}

#endif
