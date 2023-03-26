// Tests for catalog narrowing

#include <catch.hpp>

#include "io.hpp"
#include "star-utils.hpp"

using namespace lost; // NOLINT

TEST_CASE("Narrow catalog, maxStars and maxMagnitude", "[narrow-catalog]") {
    const Catalog &catalog = CatalogRead();
    int maxMagnitude = std::max_element(catalog.begin(), catalog.end(),
                                        [](const CatalogStar &a, const CatalogStar &b) {
                                            return a.magnitude < b.magnitude;
                                        })->magnitude;
    CHECK(maxMagnitude > 600);

    Catalog narrowed;
    SECTION("Narrow by maxStars") {
        narrowed = NarrowCatalog(catalog, 9999, 5000, 0);
        CHECK(narrowed.size() == 5000);
    }
    SECTION("Narrow by maxMagnitude") {
        narrowed = NarrowCatalog(catalog, 600, 9999, 0);
    }

    // assert that there are no stars with magnitude >600
    int narrowedMaxMagnitude = std::max_element(narrowed.begin(), narrowed.end(),
                                        [](const CatalogStar &a, const CatalogStar &b) {
                                            return a.magnitude < b.magnitude;
                                        })->magnitude;
    CHECK(narrowedMaxMagnitude <= 600);
}

TEST_CASE("Narrow catalog, minDistance", "[narrow-catalog]") {
    const Catalog &catalog = CatalogRead();

    // stars 2061 and 1999 are separated by about 0.0352 radians.

    // should remove most of it:
    Catalog narrowed1 = NarrowCatalog(catalog, 9999, 9999, 0.035);
    Catalog narrowed2 = NarrowCatalog(catalog, 9999, 9999, 0.036);
    CHECK(0 < narrowed1.size());
    CHECK(0 < narrowed2.size());
    CHECK(FindNamedStar(narrowed1, 2061) != narrowed1.end());
    CHECK(FindNamedStar(narrowed1, 1999) != narrowed1.end());
    CHECK(FindNamedStar(narrowed2, 2061) == narrowed2.end());
    CHECK(FindNamedStar(narrowed2, 1999) == narrowed2.end());
}
