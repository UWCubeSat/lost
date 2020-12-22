#include "star-utils.hpp"

#include <math.h>

namespace lost {

float StarDistancePixels(Star one, Star two) {
    float distX = one.x - two.x;
    float distY = one.y - two.y;
    return sqrt(distX*distX + distY*distY);
}

}
