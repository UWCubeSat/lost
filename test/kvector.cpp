#include <catch.hpp>

#include "databases.hpp"
#include "io.hpp"
#include "attitude-utils.hpp"

using namespace lost;

TEST_CASE("Kvector full database stuff", "[kvector]") {
    MultiDatabaseBuilder builder = MultiDatabaseBuilder();
    Catalog &catalog = CatalogRead();
    long length = SerializeLengthPairDistanceKVector(catalog, DegToRad(1.0), DegToRad(2.0), 100);
    unsigned char *buffer = builder.AddSubDatabase(PairDistanceKVectorDatabase::kMagicValue, length);
    SerializePairDistanceKVector(catalog, DegToRad(1.0), DegToRad(2.0), 100, buffer);
    REQUIRE(length < 999999);
    PairDistanceKVectorDatabase db(buffer);

    SECTION("basic consistency checks") {
        long lastNumReturnedPairs = 999999;
        for (float i = 1.1; i < 1.99; i += 0.1) {
            long numReturnedPairs;
            const int16_t *pairs = db.FindPairsLiberal(i * M_PI/180.0, 2.0 * M_PI/180.0, &numReturnedPairs);
            float shortestDistance = INFINITY;
            for (long k = 0; k < numReturnedPairs; k += 2) {
                float distance = AngleUnit(catalog[pairs[k]].spatial, catalog[pairs[k+1]].spatial);
                if (distance < shortestDistance) {
                    shortestDistance = distance;
                }
                CHECK(i * M_PI/180.0 <=distance);
                CHECK(distance<= 2.01 * M_PI/180.0);
            }
            REQUIRE(0 < numReturnedPairs);
            REQUIRE(numReturnedPairs < lastNumReturnedPairs);
            REQUIRE(shortestDistance < (i + 0.01) * M_PI/180.0);
            lastNumReturnedPairs = numReturnedPairs;
        }
    }

    SECTION("form a partition") {
        long totalReturnedPairs = 0;
        for (float i = 1.1; i < 2.01; i+= 0.1) {
            long numReturnedPairs;
            db.FindPairsLiberal(DegToRad(i-0.1)+0.00001, DegToRad(i)-0.00001, &numReturnedPairs);
            totalReturnedPairs += numReturnedPairs;
        }
        REQUIRE(totalReturnedPairs == db.NumPairs());
    }
}

TEST_CASE("Tighter tolerance test", "[kvector]") {
    MultiDatabaseBuilder builder = MultiDatabaseBuilder();
    Catalog &catalog = CatalogRead();
    long length = SerializeLengthPairDistanceKVector(catalog, DegToRad(0.5), DegToRad(5.0), 1000);
    unsigned char *buffer = builder.AddSubDatabase(PairDistanceKVectorDatabase::kMagicValue, length);
    SerializePairDistanceKVector(catalog, DegToRad(0.5), DegToRad(5.0), 1000, buffer);
    REQUIRE(length < 999999);
    PairDistanceKVectorDatabase db(buffer);
    // radius we'll request
    float delta = 0.0001;
    // radius we expect back: radius we request + width of a bin
    float epsilon = delta + DegToRad(5.0 - 0.5) / 1000;
    // in the first test_case, the ends of each request pretty much line up with the ends of the
    // buckets (intentionally), so that we can do the "form a partition" test. Here, however, a
    // request may intersect a bucket, in which case things slightly outside the requested range should
    // be returned.
    bool outsideRangeReturned = false;
    for (float i = DegToRad(0.6); i < DegToRad(4.9); i += DegToRad(0.1228)) {
        long numReturnedPairs;
        const int16_t *pairs = db.FindPairsLiberal(i - delta, i + delta, &numReturnedPairs);
        for (long k = 0; k < numReturnedPairs; k += 2) {
            float distance = AngleUnit(catalog[pairs[k]].spatial, catalog[pairs[k+1]].spatial);
            // only need to check one side, since we're only looking for one exception.
            if (i - delta > distance) {
                outsideRangeReturned = true;
            }
            CHECK(i - epsilon <=distance);
            CHECK(distance<= i + epsilon);
        }
    }
    CHECK(outsideRangeReturned);
}

TEST_CASE("3-star database, check exact results", "[kvector] [fast]") {
    Catalog tripleCatalog = {
        CatalogStar(DegToRad(2), DegToRad(-3), 3.0, 42),
        CatalogStar(DegToRad(4), DegToRad(7), 2.0, 43),
        CatalogStar(DegToRad(2), DegToRad(6), 4.0, 44),
    };

    MultiDatabaseBuilder builder = MultiDatabaseBuilder();
    long length = SerializeLengthTripleDistanceKVector(tripleCatalog, DegToRad(0.5), DegToRad(20.0), 1000);
    unsigned char *buffer = builder.AddSubDatabase(TripleDistanceKVectorDatabase::kMagicValue, length);
    SerializeTripleDistanceKVector(tripleCatalog, DegToRad(0.5), DegToRad(20.0), 1000, buffer);
    TripleDistanceKVectorDatabase db(buffer);
    REQUIRE(db.NumTriples() == 3);

    float distances[] = {0.038823101, 0.157079488, 0.177976221};
    for (float distance : distances) {
        long numReturnedPairs;
        const int16_t *pairs = db.FindTriplesLiberal(distance - 0.001, distance + 0.001, &numReturnedPairs);
        REQUIRE(numReturnedPairs == 1);
        CHECK(AngleUnit(tripleCatalog[pairs[0]].spatial, tripleCatalog[pairs[1]].spatial) == Approx(distance));
    }
}
