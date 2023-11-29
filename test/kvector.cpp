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
    SerializePairDistanceKVector(&ser, catalog, DegToRad(DECIMAL(1.0)), DegToRad(DECIMAL(2.0)), 100);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    SECTION("basic consistency checks") {
        long lastNumReturnedPairs = 999999;
        for (decimal i = DECIMAL(1.1); i < DECIMAL(1.99); i += DECIMAL(0.1)) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsExact(catalog, i * DECIMAL_M_PI/DECIMAL(180.0), DECIMAL(2.0) * DECIMAL_M_PI/DECIMAL(180.0), &end);
            decimal shortestDistance = INFINITY;
            for (const int16_t *pair = pairs; pair != end; pair += 2) {
                decimal distance = AngleUnit(catalog[pair[0]].spatial, catalog[pair[1]].spatial);
                if (distance < shortestDistance) {
                    shortestDistance = distance;
                }
                CHECK(i * DECIMAL_M_PI/DECIMAL(180.0) <= distance);
                CHECK(distance <= DECIMAL(2.01) * DECIMAL_M_PI/DECIMAL(180.0));
            }
            long numReturnedPairs = (end - pairs)/2;
            REQUIRE(0 < numReturnedPairs);
            REQUIRE(numReturnedPairs < lastNumReturnedPairs);
            REQUIRE(shortestDistance < (i + DECIMAL(0.01)) * DECIMAL_M_PI/DECIMAL(180.0));
            lastNumReturnedPairs = numReturnedPairs;
        }
    }

    SECTION("form a partition") {
        long totalReturnedPairs = 0;
        for (decimal i = DECIMAL(1.1); i < DECIMAL(2.01); i+= DECIMAL(0.1)) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsLiberal(DegToRad(i-DECIMAL(0.1))+DECIMAL(0.00001), DegToRad(i)-DECIMAL(0.00001), &end);
            long numReturnedPairs = (end-pairs)/2;
            totalReturnedPairs += numReturnedPairs;
        }
        REQUIRE(totalReturnedPairs == db.NumPairs());
    }
}

TEST_CASE("Tighter tolerance test", "[kvector]") {
    const Catalog &catalog = CatalogRead();
    SerializeContext ser;
    SerializePairDistanceKVector(&ser, catalog, DegToRad(DECIMAL(0.5)), DegToRad(DECIMAL(5.0)), 1000);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);
    // radius we'll request
    decimal delta = DECIMAL(0.0001);
    // radius we expect back: radius we request + width of a bin
    decimal epsilon = delta + DegToRad(DECIMAL(5.0) - DECIMAL(0.5)) / 1000;
    // in the first test_case, the ends of each request pretty much line up with the ends of the
    // buckets (intentionally), so that we can do the "form a partition" test. Here, however, a
    // request may intersect a bucket, in which case things slightly outside the requested range should
    // be returned.
    SECTION("liberal") {
        bool outsideRangeReturned = false;
        for (decimal i = DegToRad(DECIMAL(0.6)); i < DegToRad(DECIMAL(4.9)); i += DegToRad(DECIMAL(0.1228))) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsLiberal(i - delta, i + delta, &end);
            for (const int16_t *pair = pairs; pair != end; pair += 2) {
                decimal distance = AngleUnit(catalog[pair[0]].spatial, catalog[pair[1]].spatial);
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
        for (decimal i = DegToRad(DECIMAL(0.6)); i < DegToRad(DECIMAL(4.9)); i += DegToRad(DECIMAL(0.1228))) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsExact(catalog, i - delta, i + delta, &end);
            for (const int16_t *pair = pairs; pair != end; pair += 2) {
                decimal distance = AngleUnit(catalog[pair[0]].spatial, catalog[pair[1]].spatial);
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
        CatalogStar(DegToRad(2), DegToRad(-3), DECIMAL(3.0), 42),
        CatalogStar(DegToRad(4), DegToRad(7), DECIMAL(2.0), 43),
        CatalogStar(DegToRad(2), DegToRad(6), DECIMAL(4.0), 44),
    };
    SerializeContext ser;
    SerializePairDistanceKVector(&ser, tripleCatalog, DegToRad(DECIMAL(0.5)), DegToRad(DECIMAL(20.0)), 1000);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);
    REQUIRE(db.NumPairs() == 3);

    decimal distances[] = {0.038825754, 0.15707963, 0.177976474};
    SECTION("liberal") {
        for (decimal distance : distances) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsLiberal(distance - DECIMAL(1e-6), distance + DECIMAL(1e-6), &end);
            REQUIRE(end - pairs == 2);
            CHECK(AngleUnit(tripleCatalog[pairs[0]].spatial, tripleCatalog[pairs[1]].spatial) == Approx(distance).epsilon(1e-4));
        }
    }

    // also serves as a regression test for an off-by-one error that used to be present in exact, where it assumed the end index was inclusive instead of "off-the-end"
    SECTION("exact") {
        for (decimal distance : distances) {
            const int16_t *end;
            const int16_t *pairs = db.FindPairsExact(tripleCatalog, distance - DECIMAL(1e-4), distance + DECIMAL(1e-4), &end);
            REQUIRE(end - pairs == 2);
            CHECK(AngleUnit(tripleCatalog[pairs[0]].spatial, tripleCatalog[pairs[1]].spatial) == Approx(distance).epsilon(1e-4));
        }
    }
}
