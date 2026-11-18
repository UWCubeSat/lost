
#include <assert.h>
#include <math.h>

#include <cstdint>
#include <algorithm>
#include <set>
#include <vector>

#include <catch.hpp>

#include "star-utils.hpp"

namespace lost {
  bool CatalogStarMagnitudeCompare(const CatalogStar &a, const CatalogStar &b);
}

using namespace lost;

// EXAMPLE TEST CASE
// TEST_CASE("Star A bigger than Star B") {
//     // two stars created with star a bigger than star b
//     const int w = 32, h = 32;
//     const int cx = 15, cy = 15;
//     unsigned char image[w * h];
//     memset(image, 0, sizeof(image));

//     for (int dy = -1; dy <= 1; dy++) {
//         for (int dx = -1; dx <= 1; dx++) {
//             image[(cy + dy) * w + (cx + dx)] = 255;
//         }
//     }

//     REQUIRE(stars.size() == 2);
//     CHECK(stars[0].position.x == Approx(cx + 0.5).margin(0.01));
//     CHECK(stars[0].position.y == Approx(cy + 0.5).margin(0.01));
//     CHECK(stars[0].magnitude == 9);
// }

TEST_CASE("CatalogStarMagnitudeCompare") {
    // lower magnitude = brighter
    CatalogStar bright(0, 0, 100, 1);
    CatalogStar dim(0, 0, 500, 2);

    CHECK(CatalogStarMagnitudeCompare(bright, dim) == true);
    CHECK(CatalogStarMagnitudeCompare(dim, bright) == false);

    // equal magnitudes: neither comes before the other
    CatalogStar same(1, 1, 100, 3);
    CHECK(CatalogStarMagnitudeCompare(bright, same) == false);
}

TEST_CASE("NarrowCatalog") {
    Catalog catalog;
    catalog.push_back(CatalogStar(0, 0, 1, 1));
    catalog.push_back(CatalogStar(1, 0, 2, 2));
    catalog.push_back(CatalogStar(2, 0, 3, 3));
    catalog.push_back(CatalogStar(3, 0, 4, 4));
    catalog.push_back(CatalogStar(4, 0, 5, 5));

    Catalog result = NarrowCatalog(catalog, 10, 3, 0.5);
    CHECK(result.size() == 3);
    CHECK(result[0].name == 1);
}

TEST_CASE("TetraFindPatternIndex") {
    int patternEdgeRatios[] = {1, 2, 3};
    size_t numRatios = 3;
    int numBins = 10;
    size_t tableSize = 101;
    size_t index = TetraFindPatternIndex(patternEdgeRatios, numRatios, numBins, tableSize);

    CHECK(index < tableSize);
    CHECK(index == 39);
}

TEST_CASE("FindNamedStar") {
    Catalog catalog;
    catalog.push_back(CatalogStar(0, 0, 100, 10));
    catalog.push_back(CatalogStar(1, 0, 200, 20));
 
    CHECK(FindNamedStar(catalog, 20)->name == 20);
    CHECK(FindNamedStar(catalog, 99) == catalog.cend());
}
 
TEST_CASE("Serialize and deserialize catalog") {
    Catalog catalog;
    catalog.push_back(CatalogStar(1, 0.5, 200, 7));
 
    SerializeContext ser;
    SerializeCatalog(&ser, catalog, true, true);
    DeserializeContext des(ser.buffer.data());
    Catalog result = DeserializeCatalog(&des, NULL, NULL);
 
    REQUIRE(result.size() == 1);
    CHECK(result[0].magnitude == 200);
    CHECK(result[0].name == 7);
}
 
TEST_CASE("MagToBrightness") {
    CHECK(MagToBrightness(0) == Approx(1.0));
    CHECK(MagToBrightness(0) / MagToBrightness(500) == Approx(100.0));
}