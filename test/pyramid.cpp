#include <catch.hpp>
#include <vector>

#include "fixtures.hpp"
#include "attitude-utils.hpp"
#include "io.hpp"
#include "star-id.hpp"
#include "star-id-private.hpp"

#include "utils.hpp"

using namespace lost; // NOLINT

TEST_CASE("Never don't identify a pyramid", "[pyramid]") {
    float minDistance = DegToRad(0.5);
    float maxDistance = DegToRad(10.0);
    float tolerance = DegToRad(0.05);
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
            REQUIRE(j <= (int)catalog.size());

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
    delete[] dbBytes;
}

// TODO: one where spectrality is nearly zero, and thus needs to be ignored. Might be tested by above test already, but unsure.

TEST_CASE("Pyramid selection: Basic strategy only on starCloisters", "[pyramid] [fast]") {
    std::vector<Vec3> starCloistersSpatials;
    for (Star star : starCloisters) {
        starCloistersSpatials.push_back(smolCamera.CameraToSpatial(star.position).Normalize());
    }

    PyramidIterator pyIter(starCloistersSpatials, 0.0, 100.0);
    BestPyramidAtStar py1 = pyIter.Next();
    BestPyramidAtStar py2 = pyIter.Next();
    BestPyramidAtStar py3 = pyIter.Next();
    // CHECK(py1.distancesSum == Approx(3).margin(0.3)); // they're not even close, need to way up the margin
    // CHECK(py2.distancesSum == Approx(6).margin(0.3));
    CHECK(py1.distancesSum < py2.distancesSum);
    // iteration stopped:
    CHECK(py3.distancesSum < 0);

    // of course, you could use a loop here...but it's way easier to read unrolled, and I have the unlimited power of Github copilot!
    CHECK(py1.centroidIndices[0] == 0);
    CHECK(py1.centroidIndices[1] == 1);
    CHECK(py1.centroidIndices[2] == 2);
    CHECK(py1.centroidIndices[3] == 3);

    CHECK(py2.centroidIndices[0] == 4);
    CHECK(py2.centroidIndices[1] == 5);
    CHECK(py2.centroidIndices[2] == 6);
    CHECK(py2.centroidIndices[3] == 7);
}

// TODO:
// TEST_CASE("Pyramids that exist and are unique are always matched immediately", "[pyramid] [fast]") {

// }
