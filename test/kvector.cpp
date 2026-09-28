#include <catch.hpp>

#include "databases.hpp"
#include "io.hpp"
#include "attitude-utils.hpp"
#include "serialize-helpers.hpp"

#include "utils.hpp"

using namespace lost; // NOLINT

TEST_CASE("Kvector full database stuff", "[kvector]") {
    const Catalog &catalog = CatalogRead();
    std::vector<unsigned char> dbBytes;
    SerializeContext ser;
    SerializePairDistanceKVector(&ser, catalog, DegToRad(SCALAR(1.0)), DegToRad(SCALAR(2.0)), 100);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    SECTION("basic consistency checks") {
        long lastNumReturnedPairs = 999999;
        for (scalar i = SCALAR(1.1); i < SCALAR(1.99); i += SCALAR(0.1)) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsExact(catalog, i * SCALAR_M_PI/SCALAR(180.0), SCALAR(2.0) * SCALAR_M_PI/SCALAR(180.0), &end);
            scalar shortestDistance = INFINITY;
            for (const int16_t *pair = pairs; pair != end; pair += 2) {
                scalar distance = AngleUnit(catalog[pair[0]].spatial, catalog[pair[1]].spatial);
                if (distance < shortestDistance) {
                    shortestDistance = distance;
                }
                CHECK(i * SCALAR_M_PI/SCALAR(180.0) <= distance);
                CHECK(distance <= SCALAR(2.01) * SCALAR_M_PI/SCALAR(180.0));
            }
            long numReturnedPairs = (end - pairs)/2;
            REQUIRE(0 < numReturnedPairs);
            REQUIRE(numReturnedPairs < lastNumReturnedPairs);
            REQUIRE(shortestDistance < (i + SCALAR(0.01)) * SCALAR_M_PI/SCALAR(180.0));
            lastNumReturnedPairs = numReturnedPairs;
        }
    }

    SECTION("form a partition") {
        long totalReturnedPairs = 0;
        for (scalar i = SCALAR(1.1); i < SCALAR(2.01); i+= SCALAR(0.1)) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsLiberal(DegToRad(i-SCALAR(0.1))+SCALAR(0.00001), DegToRad(i)-SCALAR(0.00001), &end);
            long numReturnedPairs = (end-pairs)/2;
            totalReturnedPairs += numReturnedPairs;
        }
        REQUIRE(totalReturnedPairs == db.NumPairs());
    }
}

TEST_CASE("Tighter tolerance test", "[kvector]") {
    const Catalog &catalog = CatalogRead();
    SerializeContext ser;
    SerializePairDistanceKVector(&ser, catalog, DegToRad(SCALAR(0.5)), DegToRad(SCALAR(5.0)), 1000);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);
    // radius we'll request
    scalar delta = SCALAR(0.0001);
    // radius we expect back: radius we request + width of a bin
    scalar epsilon = delta + DegToRad(SCALAR(5.0) - SCALAR(0.5)) / 1000;
    // in the first test_case, the ends of each request pretty much line up with the ends of the
    // buckets (intentionally), so that we can do the "form a partition" test. Here, however, a
    // request may intersect a bucket, in which case things slightly outside the requested range should
    // be returned.
    SECTION("liberal") {
        bool outsideRangeReturned = false;
        for (scalar i = DegToRad(SCALAR(0.6)); i < DegToRad(SCALAR(4.9)); i += DegToRad(SCALAR(0.1228))) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsLiberal(i - delta, i + delta, &end);
            for (const int16_t *pair = pairs; pair != end; pair += 2) {
                scalar distance = AngleUnit(catalog[pair[0]].spatial, catalog[pair[1]].spatial);
                // only need to check one side, since we're only looking for one exception.
                if (i - delta > distance) {
                    outsideRangeReturned = true;
                }
                CHECK(i - epsilon <= distance);
                CHECK(distance<= i + epsilon);
            }
        }
        CHECK(outsideRangeReturned);
    }
    SECTION("exact") {
        bool outsideRangeReturned = false;
        for (scalar i = DegToRad(SCALAR(0.6)); i < DegToRad(SCALAR(4.9)); i += DegToRad(SCALAR(0.1228))) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsExact(catalog, i - delta, i + delta, &end);
            for (const int16_t *pair = pairs; pair != end; pair += 2) {
                scalar distance = AngleUnit(catalog[pair[0]].spatial, catalog[pair[1]].spatial);
                // only need to check one side, since we're only looking for one exception.
                if (i - delta > distance) {
                    outsideRangeReturned = true;
                }
                CHECK(i - epsilon <= distance);
                CHECK(distance <= i + epsilon);
            }
        }
        CHECK(!outsideRangeReturned);
    }
}

TEST_CASE("3-star database, check exact results", "[kvector] [fast]") {
    Catalog tripleCatalog = {
        CatalogStar(DegToRad(2), DegToRad(-3), SCALAR(3.0), 42),
        CatalogStar(DegToRad(4), DegToRad(7), SCALAR(2.0), 43),
        CatalogStar(DegToRad(2), DegToRad(6), SCALAR(4.0), 44),
    };
    SerializeContext ser;
    SerializePairDistanceKVector(&ser, tripleCatalog, DegToRad(SCALAR(0.5)), DegToRad(SCALAR(20.0)), 1000);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);
    REQUIRE(db.NumPairs() == 3);

    scalar distances[] = {0.038825754, 0.15707963, 0.177976474};
    SECTION("liberal") {
        for (scalar distance : distances) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsLiberal(distance - SCALAR(1e-6), distance + SCALAR(1e-6), &end);
            REQUIRE(end - pairs == 2);
            CHECK(AngleUnit(tripleCatalog[pairs[0]].spatial, tripleCatalog[pairs[1]].spatial) == Approx(distance).epsilon(1e-4));
        }
    }

    // also serves as a regression test for an off-by-one error that used to be present in exact, where it assumed the end index was inclusive instead of "off-the-end"
    SECTION("exact") {
        for (scalar distance : distances) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsExact(tripleCatalog, distance - SCALAR(1e-4), distance + SCALAR(1e-4), &end);
            REQUIRE(end - pairs == 2);
            CHECK(AngleUnit(tripleCatalog[pairs[0]].spatial, tripleCatalog[pairs[1]].spatial) == Approx(distance).epsilon(1e-4));
        }
    }
}
