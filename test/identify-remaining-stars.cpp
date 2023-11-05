#include <random>

#include <catch.hpp>

#include "databases.hpp"
#include "star-id.hpp"
#include "star-id-private.hpp"

#include "fixtures.hpp"
#include "utils.hpp"

using namespace lost; // NOLINT

TEST_CASE("IRUnidentifiedCentroid simple orthogonal", "[identify-remaining] [fast]") {
    IRUnidentifiedCentroid centroid(elevenStars[3], 0);
    REQUIRE(centroid.bestAngleFrom90 > 9e9);
    centroid.AddIdentifiedStar(elevenStarIds[1], elevenStars);
    // one star is not enough to get angle
    REQUIRE(centroid.bestAngleFrom90 > 9e9);
    centroid.AddIdentifiedStar(elevenStarIds[2], elevenStars);
    // we've set them up to be almost orthogonal
    REQUIRE(centroid.bestAngleFrom90 == Approx(0).margin(1e-6));
    REQUIRE(((centroid.bestStar1 == elevenStarIds[1] && centroid.bestStar2 == elevenStarIds[2]) ||
             (centroid.bestStar1 == elevenStarIds[2] && centroid.bestStar2 == elevenStarIds[1])));

    // adding another, non-orthogonal one shouldn't break it
    centroid.AddIdentifiedStar(elevenStarIds[0], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(0).margin(1e-6));
    REQUIRE(((centroid.bestStar1 == elevenStarIds[1] && centroid.bestStar2 == elevenStarIds[2]) ||
             (centroid.bestStar1 == elevenStarIds[2] && centroid.bestStar2 == elevenStarIds[1])));
}

TEST_CASE("IRUnidentifiedCentroid not orthogonal until they are", "[identify-remaining] [fast]") {
    IRUnidentifiedCentroid centroid(elevenStars[3], 0);
    centroid.AddIdentifiedStar(elevenStarIds[1], elevenStars);
    centroid.AddIdentifiedStar(elevenStarIds[0], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(M_PI_4));
    centroid.AddIdentifiedStar(elevenStarIds[6], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(M_PI_4));
    centroid.AddIdentifiedStar(elevenStarIds[8], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(0).margin(1e-6));
}

TEST_CASE("IRUnidentifiedCentroid obtuse angle", "[identify-remaining] [fast]") {
    IRUnidentifiedCentroid centroid(elevenStars[3], 0);
    centroid.AddIdentifiedStar(elevenStarIds[1], elevenStars);
    centroid.AddIdentifiedStar(elevenStarIds[6], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(M_PI_4));
}

// TODO: Tests for FindAllInRange if we ever make the logic more complicated

std::vector<int16_t> IdentifyThirdStarTest(const Catalog &catalog, int16_t catalogName1, int16_t catalogName2,
                                           float dist1, float dist2, float tolerance) {
    SerializeContext ser;
    SerializePairDistanceKVector(&ser, integralCatalog, 0, M_PI, 1000);
    DeserializeContext des(ser.buffer.data());
    auto cs1 = FindNamedStar(catalog, catalogName1);
    auto cs2 = FindNamedStar(catalog, catalogName2);

    PairDistanceKVectorDatabase db(&des);
    auto result = IdentifyThirdStar(db,
                                    catalog,
                                    cs1 - catalog.cbegin(), cs2 - catalog.cbegin(),
                                    dist1, dist2, tolerance);
    return result;
}

TEST_CASE("IdentifyThirdStar", "[identify-remaining] [fast]") { // TODO: does it /really/ logically belong with identify-remaining? Maybe we should coin a new term for star pattern identification related functions

    std::vector<int16_t> stars = IdentifyThirdStarTest(integralCatalog,
                                                       42, 44, // (1,0,0), (0,1,0)
                                                       M_PI_2, M_PI_2,
                                                       1e-6);
    REQUIRE(stars.size() == 1);
    REQUIRE(integralCatalog[stars[0]].name == 50);

}

TEST_CASE("IdentifyThirdStar with tolerance", "[identify-remaining] [fast]") {

    // try it again where we actually need the tolerance
    std::vector<int16_t> stars = IdentifyThirdStarTest(integralCatalog,
                                                       42, 44, // (1,0,0), (0,1,0)
                                                       M_PI_2 - DegToRad(1.0), M_PI_2 + DegToRad(1.0),
                                                       0.1);
    REQUIRE(stars.size() == 1);
    REQUIRE(integralCatalog[stars[0]].name == 50);
}

TEST_CASE("IdentifyThirdStar reversed spectrality", "[identify-remaining] [fast]") {
    std::vector<int16_t> stars = IdentifyThirdStarTest(integralCatalog,
                                                       44, 42, // (0,1,0), (1,0,0)
                                                       M_PI_2, M_PI_2,
                                                       1e-6);
    REQUIRE(stars.size() == 1);
    REQUIRE(integralCatalog[stars[0]].name == 58);
}

TEST_CASE("IdentifyThirdStar no third star", "[identify-remaining] [fast]") {
    std::vector<int16_t> stars = IdentifyThirdStarTest(integralCatalog,
                                                       42, 44, // (1,0,0), (0,1,0)
                                                       1, M_PI_2,
                                                       1e-6);
    REQUIRE(stars.size() == 0);
}

TEST_CASE("IdentifyThirdStar just out of tolerance", "[identify-remaining] [fast]") {
    std::vector<int16_t> stars2 = IdentifyThirdStarTest(integralCatalog,
                                                        42, 44, // (1,0,0), (0,1,0)
                                                        M_PI_2 - 2e-6, M_PI_2,
                                                        1e-6);
    REQUIRE(stars2.size() == 0);
}

// This test relies on something marked TODO in star-id.cpp, but really isn't that important. Leaving here for posterity.
// TEST_CASE("IdentifyThirdStar nearly colinear, shouldn't check spectrality", "[identify-remaining] [fast]") {
//     std::vector<int16_t> stars2 = IdentifyThirdStarTest(nearlyColinearCatalog,
//                                                         42, 43,
//                                                         DegToRad(2.0), DegToRad(1.0),
//                                                         1e-2);
//     std::vector<int16_t> stars3 = IdentifyThirdStarTest(nearlyColinearCatalog,
//                                                         43, 42,
//                                                         DegToRad(1.0), DegToRad(2.0),
//                                                         1e-2);
//     REQUIRE(stars2.size() == 1);
//     REQUIRE(nearlyColinearCatalog[stars2[0]].name == 44);
//     REQUIRE(stars3.size() == 1);
//     REQUIRE(nearlyColinearCatalog[stars3[0]].name == 44);
// }

#define kIdentifyRemainingNumImages 10

TEST_CASE("IdentifyRemainingStars fuzz", "[identify-remaining] [fuzz]") {

    // generate a whole bunch of centroids, project all to 3d to create a fake catalog, create
    // kvector,then choose a few randomly to pre-identify, and then ensure the rest are all
    // identified. To ensure there's no ambiguity, we always increment the x-value. This results in
    // a weird triangular distribution of stars, but that's better than nothing!
    int numFakeStars = 100;

    std::default_random_engine rng(GENERATE(take(kIdentifyRemainingNumImages, random(0, 1000000))));
    std::uniform_real_distribution<float> yDist(0.0, 256.0);
    std::uniform_int_distribution<int> starIndexDist(0, numFakeStars - 1);
    std::uniform_int_distribution<int> moreStartingStars(0, 1);

    Stars fakeCentroids;
    StarIdentifiers fakeStarIds;
    Catalog fakeCatalog;
    for (int i = 0; i < numFakeStars; i++) {
        // smolCamera has width 256
        float x = i * 256.0 / numFakeStars;
        float y = yDist(rng);
        fakeCentroids.emplace_back(x, y, 1);
        fakeCatalog.emplace_back(smolCamera.CameraToSpatial({x, y}).Normalize(), 1, i);
        fakeStarIds.emplace_back(i, i);
    }

    // can't just do resize because there's no default constructor for StarIdentifier
    StarIdentifiers fakeStarIdsCopy = fakeStarIds;
    std::shuffle(fakeStarIdsCopy.begin(), fakeStarIdsCopy.end(), rng);
    StarIdentifiers someFakeStarIds;
    int fakePatternSize = moreStartingStars(rng) ? 4 : 2;
    for (int i = 0; i < fakePatternSize; i++) {
        someFakeStarIds.push_back(fakeStarIdsCopy[i]);
    }

    SerializeContext ser;
    SerializePairDistanceKVector(&ser, fakeCatalog, 0, M_PI, 1000);
    DeserializeContext des(ser.buffer.data());
    PairDistanceKVectorDatabase db(&des);

    int numIdentified = IdentifyRemainingStarsPairDistance(&someFakeStarIds, fakeCentroids, db, fakeCatalog, smolCamera, 1e-5);

    REQUIRE(numIdentified == numFakeStars - fakePatternSize);
    REQUIRE(AreStarIdentifiersEquivalent(fakeStarIds, someFakeStarIds));
}

// TODO: Test if kvector MaxDistance is substantially less than FOV

// TODO: Test when some stars can't be identified.

// TODO: Test with false stars.
