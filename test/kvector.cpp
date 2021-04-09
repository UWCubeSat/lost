#include <catch.hpp>

#include "databases.hpp"
#include "io.hpp"
#include "attitude-utils.hpp"

using namespace lost;

TEST_CASE("Kvector basics: Create, verify that all stars are actually in range.") {
    long length;
    Catalog &catalog = CatalogRead();
    unsigned char *dbBytes = BuildKVectorDatabase(catalog, &length, 1.0 * M_PI/180.0, 2.0 * M_PI/180.0, 100);
    KVectorDatabase db(dbBytes);
    REQUIRE(length < 999999);

    long lastNumReturnedPairs = 999999;
    for (float i = 1.1; i < 1.99; i += 0.1) {
        long numReturnedPairs;
        int16_t *pairs = db.FindPossibleStarPairsApprox(i, 1.0, &numReturnedPairs);
        for (long k = 0; k < numReturnedPairs; k++) {
            float distance = GreatCircleDistance(
                catalog[pairs[k]].raj2000, catalog[pairs[k]].dej2000,
                catalog[pairs[k+1]].raj2000, catalog[pairs[k]].dej2000);
            REQUIRE(i<=distance);
            REQUIRE(distance<=2.01);
        }
        REQUIRE(0 < numReturnedPairs);
        REQUIRE(numReturnedPairs < lastNumReturnedPairs);
        lastNumReturnedPairs = numReturnedPairs;
        delete pairs;
    }
}
