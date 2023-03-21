#include <catch.hpp>

#include <random>
#include <assert.h>

#include "attitude-utils.hpp"
#include "io.hpp"
#include "star-id.hpp"
#include "star-id-private.hpp"

#include "utils.hpp"

using namespace lost; // NOLINT

TEST_CASE("Never don't identify a pyramid", "[pyramid]") {
    float minDistance = DegToRad(0.5);
    float maxDistance = DegToRad(10.0);
    float tolerance = DegToRad(0.1);
    // What fraction of the pyramids must be /uniquely/ identified. The test always requires that at
    // least one identification be made for each pyramid, but sometimes there are multiple.
    float minFractionUniquelyIdentified = 0.75;

    long length;
    const Catalog catalog = NarrowCatalog(CatalogRead(), 9999, 9999, DegToRad(0.5));
    unsigned char *dbBytes = BuildPairDistanceKVectorDatabase(
        catalog, &length, minDistance, maxDistance, 10000);
    PairDistanceKVectorDatabase db(dbBytes);
    std::cerr << "done narrowing and building" << std::endl;

    // now the fun begins
    int numPyramidsToTry = 20;
    // how do we "randomly" pick first stars without repetition? Pretty easy: Pick a prime modulus
    // >catalog.size(), then choose some other number, which you keep multiplying and modulo-ing,
    // you know it will create every element mod p before it repeats.
    int modulus = 9103; // is prime
    int multiplier1 = 737;
    int multiplier2 = 8822;
    int numPyramidsTried = 0;
    int numUniquelyIdentified = 0;
    for (int i = 1; numPyramidsTried < numPyramidsToTry; i++) {
        int startIndex = i * multiplier1 % modulus;
        // just keep sampling until we're inside the catalog
        if (startIndex > (int)catalog.size()) {
            continue;
        }
        numPyramidsTried++;

        // find other stars, creating a set of stars that are all within range.
        // This loop could be sped up substantially by sorting on x-value, but whatever, it's just a test!
        std::vector<int> catalogIndices{startIndex};
        for (int j = 1; catalogIndices.size() < 4; j++) {
            // There should always be three other stars within 20 degrees!
            assert(j <= (int)catalog.size());

            int otherIndex = j * multiplier2 % modulus;
            if (otherIndex > (int)catalog.size()) {
                continue;
            }
            for (int alreadyChosenIndex : catalogIndices) {
                float dist = AngleUnit(catalog[alreadyChosenIndex].spatial, catalog[otherIndex].spatial);
                if (dist <= minDistance+tolerance || dist >= maxDistance-tolerance) {
                    goto nextOtherIndex;
                }
            }
            // made it through the gauntlet!
            catalogIndices.push_back(otherIndex);

        nextOtherIndex:;
        }

        int a, b, c, d;
        IdentifyPyramidResult matchResult
            = IdentifyPyramid(db,
                              catalog,
                              tolerance,
                              catalog[catalogIndices[0]].spatial,
                              catalog[catalogIndices[1]].spatial,
                              catalog[catalogIndices[2]].spatial,
                              catalog[catalogIndices[3]].spatial,
                              &a, &b, &c, &d);
        CHECK((matchResult == IdentifyPyramidResult::MatchedUniquely
               || matchResult == IdentifyPyramidResult::MatchedAmbiguously));
        if (matchResult == IdentifyPyramidResult::MatchedUniquely) {
            numUniquelyIdentified++;
            CHECK(a == catalogIndices[0]);
            CHECK(b == catalogIndices[1]);
            CHECK(c == catalogIndices[2]);
            CHECK(d == catalogIndices[3]);
        }
    }

    CHECK((float)numUniquelyIdentified / numPyramidsToTry >= minFractionUniquelyIdentified);
}

// TODO: one where spectrality is nearly zero, and thus needs to be ignored. Might be tested by above test already, but unsure.
