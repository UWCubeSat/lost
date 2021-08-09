#ifndef STAR_UTILS_H
#define STAR_UTILS_H

#include <vector>

#include "attitude-utils.hpp"
#include "camera.hpp"

namespace lost {

class CatalogStar {
public:
    CatalogStar() = default;
    CatalogStar(float raj2000, float dej2000, int magnitude, int name) :
        spatial(SphericalToSpatial(raj2000, dej2000)), magnitude(magnitude), name(name) { }

    Vec3 spatial;
    int  magnitude;         // *10^-2
    int  name;
};

class Star {
public:
    Star(float x, float y, float radiusX, float radiusY, int magnitude) :
        position({x, y}), radiusX(radiusX), radiusY(radiusY), magnitude(magnitude) { };
    Star(float x, float y, float radiusX) : Star(x, y, radiusX, radiusX, 0) { };
    Star() : Star(0.0, 0.0, 0.0) { };

    Vec2 position;
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

long SerializeLengthCatalog(const Catalog &, bool inclMagnitude, bool inclName);
void SerializeCatalog(const Catalog &, bool inclMagnitude, bool inclName, unsigned char *buffer);
// sets magnited and name to whether the catalog in the database contained magnitude and name
Catalog DeserializeCatalog(const unsigned char *buffer, bool *inclMagnitudeReturn, bool *inclNameReturn);

float StarDistancePixels(Star one, Star two);

void FocalPlaneAngles(const Stars &stars, const Camera &camera, float &min, float &mid, float &max, int &mindex, int &middex, int &maxdex, int i_index, int j_index, int k_index);

void InnerAngles(const Catalog &catalog, float &min, float &mid, float &max, int &mindex, int &middex, int &maxdex, int index1, int index2, int index3);

float MinInnerAngle(const Catalog &catalog, int &arg, int index1, int index2, int index3);


// TODO: make maxStars work right, need to sort by magnitude before filter! maxMagnitude is 10^-2
// (so 523 = 5.23)
Catalog NarrowCatalog(const Catalog &, int maxMagnitude, int maxStars);

}

#endif
