#ifndef FIXTURES_H
#define FIXTURES_H

#include "star-utils.hpp"
#include "attitude-utils.hpp"
#include "camera.hpp"

namespace lost {

extern Catalog tripleCatalog;
extern Catalog quadrupleCatalog;
extern Catalog harderQuadrupleCatalog;
extern Catalog integralCatalog;

extern Camera smolCamera;
extern Camera smolCameraOff;
extern Quaternion straightAhead;

extern Stars elevenStars;
extern StarIdentifiers elevenStarIds;

}

#endif
