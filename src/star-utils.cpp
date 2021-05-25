#include "star-utils.hpp"

#include <math.h>

namespace lost {

float StarDistancePixels(Star one, Star two) {
    float distX = one.x - two.x;
    float distY = one.y - two.y;
    return sqrt(distX*distX + distY*distY);
}

// TODO: this function is about to cause a horribly difficult to debug error, because the catalog
// used in the database may not equal the catalog used by star id or something.

// TODO: the maxStars is worthless, it doesn't get the brightest stars
Catalog NarrowCatalog(const Catalog &catalog, int maxMagnitude, int maxStars) {
    Catalog result;
    for (int i = 0; i < (((int)catalog.size() > maxStars) ? maxStars : (int)catalog.size()); i++) {
        if (catalog[i].magnitude <= maxMagnitude) {
            result.push_back(catalog[i]);
        }
    }

    return result;
}

}
