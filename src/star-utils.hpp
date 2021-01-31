#ifndef STAR_UTILS_H
#define STAR_UTILS_H

#include <vector>

namespace lost {

class CatalogStar {
public:
    CatalogStar(float raj2000, float dej2000, int magnitude, bool weird, int name) :
        raj2000(raj2000), dej2000(dej2000), magnitude(magnitude), weird(weird), name(name) { }
    float raj2000;           // *10^-6, right ascension
    float dej2000;           // *10^-6, declination
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

}

#endif
