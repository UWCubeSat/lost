// unlike most files in this folder, there /is/ a fixtures.hpp that must be kept up-to-date with this file.

#include <math.h>

#include "star-utils.hpp"
#include "attitude-utils.hpp"
#include "camera.hpp"

using namespace lost; // NOLINT

// three stars that should be easy to identify based on their distances to each other
Catalog tripleCatalog = {
    CatalogStar(DegToRad(2), DegToRad(-3), 3.0, 42),
    CatalogStar(DegToRad(4), DegToRad(7), 2.0, 43),
    CatalogStar(DegToRad(2), DegToRad(6), 4.0, 44),
};

Catalog quadrupleCatalog = {
    CatalogStar(DegToRad(2), DegToRad(-3), 3.0, 42),
    CatalogStar(DegToRad(4), DegToRad(7), 2.0, 43),
    CatalogStar(DegToRad(2), DegToRad(6), 4.0, 44),
    CatalogStar(DegToRad(-1), DegToRad(-4), 1.0, 45),
};

// distance between 42/43 == distance between 44/45, have to use distances to 43 to differentiate
Catalog harderQuadrupleCatalog = {
    CatalogStar(DegToRad(2), DegToRad(-3), 3.0, 42),
    CatalogStar(DegToRad(4), DegToRad(7), 2.0, 43),
    CatalogStar(DegToRad(2), DegToRad(5), 4.0, 44),
    CatalogStar(DegToRad(2), DegToRad(1), 1.0, 45),
};

Catalog integralCatalog = {
    CatalogStar(0, 0, 3.0, 42), // (0,0,0)
    CatalogStar(M_PI_4, 0, 3.0, 43), // (.707,.707,0)
    CatalogStar(M_PI_2, 0, 3.0, 44), // (0,1,0)
    CatalogStar(3 * M_PI_4, 0, 3.0, 45), // (-.707,.707,0)

    CatalogStar(0, M_PI_4, 3.0, 46), // (.707,0,.707)
    CatalogStar(M_PI_4, M_PI_4, 3.0, 47), // (.5,.5,.707)
    CatalogStar(M_PI_2, M_PI_4, 3.0, 48), // (0,.707,.707)
    CatalogStar(3 * M_PI_4, M_PI_4, 3.0, 49), // (-.5,.5,.707)

    CatalogStar(0, M_PI_2, 3.0, 50), // (0,0,1)

    CatalogStar(0, -M_PI_4, 3.0, 54), // (.707,0,-.707)
    CatalogStar(M_PI_4, -M_PI_4, 3.0, 55), // (.5,.5,-.707)
    CatalogStar(M_PI_2, -M_PI_4, 3.0, 56), // (0,.707,-.707)
    CatalogStar(3 * M_PI_4, -M_PI_4, 3.0, 57), // (-.5,.5,-.707)

    CatalogStar(0, -M_PI_2, 3.0, 58), // (0,0,-1)
};

Camera smolCamera(FovToFocalLength(DegToRad(36.0), 256), 256, 256);
Camera smolCameraOff(FovToFocalLength(DegToRad(33.0), 256), 256, 256);
Quaternion straightAhead(1.0, 0.0, 0.0, 0.0);

Stars elevenStars = {
    Star(0, 0, 1), // 0
    Star(1, 0, 1), // 1
    Star(0, 1, 1), // 2
    Star(1, 1, 1), // 3
    Star(2, 0, 1), // 4
    Star(2, 1, 1), // 5
    Star(2, 2, 1), // 6
    Star(3, 0, 1), // 7
    Star(3, 1, 1), // 8
    Star(3, 2, 1), // 9
    Star(3, 3, 1), // 10
    Star(1, 3, 1), // 11
};

StarIdentifiers elevenStarIds = {
    StarIdentifier(0, 0, 1),
    StarIdentifier(1, 1, 1),
    StarIdentifier(2, 2, 1),
    StarIdentifier(3, 3, 1),
    StarIdentifier(4, 4, 1),
    StarIdentifier(5, 5, 1),
    StarIdentifier(6, 6, 1),
    StarIdentifier(7, 7, 1),
    StarIdentifier(8, 8, 1),
    StarIdentifier(9, 9, 1),
    StarIdentifier(10, 10, 1),
    StarIdentifier(11, 11, 1),
};
