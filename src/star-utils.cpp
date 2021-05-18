#include "star-utils.hpp"

#include <math.h>

namespace lost {

float StarDistancePixels(Star one, Star two) {
    float distX = one.x - two.x;
    float distY = one.y - two.y;
    return sqrt(distX*distX + distY*distY);
}

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
