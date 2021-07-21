#include "catch.hpp"

#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "databases.hpp"
#include "star-id.hpp"
#include "star-utils.hpp"
#include "camera.hpp"
#include "io.hpp"

using namespace lost;

// three stars that should be easy to identify based on their distances to each other
// lets see how it performs when the angles are more important
Catalog tripleCatalog = {
    CatalogStar(DegToRad(0), DegToRad(0), 1.0, 42),
    CatalogStar(DegToRad(90), DegToRad(0), 0, 43),
    CatalogStar(DegToRad(90), DegToRad(90), 1.0, 44),
};

TEST_CASE("inner angle test simple", "[gv] [fast]") {
    const Catalog &catalog = GENERATE((const Catalog &)tripleCatalog);
    float min;
    float mid;
    float max;
    int mindex;
    int middex;
    int maxdex;
    innerAngles(catalog, min, mid, max, mindex, middex, maxdex, 42, 43, 44);
    std::cout << "mindex: " << mindex << std::endl;
    std::cout << "middex: " << middex << std::endl;
    std::cout << "maxdex: " << maxdex << std::endl;
    CHECK(min == DegToRad(45));
    CHECK(max == DegToRad(90));
    //REQUIRE(comparison.numTotal == catalog.size());
    //CHECK(comparison.numCorrect == catalog.size());
}

TEST_CASE("focal plane angle test simple", "[gv] [fast]") {
    Stars stars;
    Star s1 = Star(0, 0, 1.0);
    Star s2 = Star(1, 0, 1.0);
    Star s3 = Star(0, 1, 1.0);
    stars.push_back(s1);
    stars.push_back(s2);
    stars.push_back(s3);
    float min;
    float mid;
    float max;
    int mindex;
    int middex;
    int maxdex;
    focalPlaneAngles(stars, min, mid, max, mindex, middex, maxdex, 0, 1, 2);
    CHECK(max == DegToRad(90));
    CHECK(mid == DegToRad(45));
    CHECK(min == DegToRad(45));
    CHECK(mindex > 0);
    CHECK(middex > 0);
    CHECK(maxdex == 0);
}