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
    CatalogStar(DegToRad(2), DegToRad(-3), 1.0, 42),
    CatalogStar(DegToRad(4), DegToRad(7), 1.0, 43),
    CatalogStar(DegToRad(2), DegToRad(6), 1.0, 44),
};

TEST_CASE("angle tests", "[gv] [fast]") {
    const Catalog &catalog = GENERATE((const Catalog &)tripleCatalog);
    float min;
    float mid;
    float max;
    int mindex;
    int middex;
    int maxdex;
    innerAngles(catalog, min, mid, max, mindex, middex, maxdex, 42, 43, 44);
    std::cout << "min: " << min << std::endl;
    std::cout << "mid: " << mid << std::endl;
    std::cout << "max: " << max << std::endl;
    std::cout << "mindex: " << mindex << std::endl;
    std::cout << "middex: " << middex << std::endl;
    std::cout << "maxdex: " << maxdex << std::endl;
    //StarIdComparison comparison = StarIdsCompare(*gpi.ExpectedStarIds(), outputStarIds, 0.0, gpi.ExpectedStars(), NULL);
    //REQUIRE(comparison.numTotal == catalog.size());
    //CHECK(comparison.numCorrect == catalog.size());
}