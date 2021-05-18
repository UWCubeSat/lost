#include <catch.hpp>

#include "databases.hpp"
#include "io.hpp"
#include "attitude-utils.hpp"

using namespace lost;

TEST_CASE("Kvector full database stuff", "[kvector]") {
    long length;
    Catalog &catalog = CatalogRead();
    unsigned char *dbBytes = BuildKVectorDatabase(catalog, &length, 1.0 * M_PI/180.0, 2.0 * M_PI/180.0, 100);
    KVectorDatabase db(dbBytes);
    REQUIRE(length < 999999);

    SECTION("basic consistency checks") {
        long lastNumReturnedPairs = 999999;
        for (float i = 1.1; i < 1.99; i += 0.1) {
            long numReturnedPairs;
            int16_t *pairs = db.FindPossibleStarPairsApprox(i * M_PI/180.0, 2.0 * M_PI/180.0, &numReturnedPairs);
            float shortestDistance = INFINITY;
            for (long k = 0; k < numReturnedPairs; k += 2) {
                float distance = AngleUnit(catalog[pairs[k]].spatial, catalog[pairs[k+1]].spatial);
                if (distance < shortestDistance) {
                    shortestDistance = distance;
                }
                REQUIRE(i * M_PI/180.0 <=distance);
                REQUIRE(distance<= 2.01 * M_PI/180.0);
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
            db.FindPossibleStarPairsApprox((i-0.1) * M_PI/180.0, i * M_PI/180.0, &numReturnedPairs);
            totalReturnedPairs += numReturnedPairs;
            printf("%f to %f: %ld\n", i-.1, i, numReturnedPairs);
        }
        REQUIRE(totalReturnedPairs == db.NumPairs());
    }
}
