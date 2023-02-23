#include <catch.hpp>

#include "star-id.hpp"
#include "star-id-private.hpp"

#include "fixtures.hpp"

using namespace lost; // NOLINT

TEST_CASE("IRUnidentifiedCentroid simple orthogonal", "[identify-remaining] [fast]") {
    IRUnidentifiedCentroid centroid(elevenStars[3]);
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
    IRUnidentifiedCentroid centroid(elevenStars[3]);
    centroid.AddIdentifiedStar(elevenStarIds[1], elevenStars);
    centroid.AddIdentifiedStar(elevenStarIds[0], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(M_PI_4));
    centroid.AddIdentifiedStar(elevenStarIds[6], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(M_PI_4));
    centroid.AddIdentifiedStar(elevenStarIds[8], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(0).margin(1e-6));
}

TEST_CASE("IRUnidentifiedCentroid obtuse angle", "[identify-remaining] [fast]") {
    IRUnidentifiedCentroid centroid(elevenStars[3]);
    centroid.AddIdentifiedStar(elevenStarIds[1], elevenStars);
    centroid.AddIdentifiedStar(elevenStarIds[6], elevenStars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(M_PI_4));
}

// TODO: Tests for FindAllInRange if we ever make the logic more complicated

TEST_CASE("IdentifyThirdStar simple case (does require spectrality)", "[identify-remaining] [fast]") { // TODO: does it /really/ logically belong with identify-remaining? Maybe we should coin a new term for star pattern identification related functions
    std::vector<int16_t> stars = IdentifyThirdStar(integralCatalog,
                                                   Vec3(0,0,0), Vec3(0,1,0),
                                                   1.0, sqrt(2.0),
                                                   1e-6);
    REQUIRE(stars.size() == 1);
    REQUIRE(integralCatalog[stars[0]].name == 50);
}

TEST_CASE("IdentifyThirdStar tolerance", "[identify-remaining] [fast]") {
    // try it again where we actually need the tolerance
    std::vector<int16_t> stars2 = IdentifyThirdStar(integralCatalog,
                                                    Vec3(0,0,0), Vec3(0,1,0),
                                                    0.99, sqrt(2.0) + 0.02,
                                                    0.1);
    REQUIRE(stars.size() == 1);
    REQUIRE(integralCatalog[stars[0]].name == 50);
}

TEST_CASE("IdentifyThirdStar reversed spectrality", "[identify-remaining] [fast]") {
    std::vector<int16_t> stars = IdentifyThirdStar(integralCatalog,
                                                    Vec3(0,1,0), Vec3(0,0,0),
                                                    sqrt(2.0), 1.0,
                                                    1e-6);
    REQUIRE(stars.size() == 1);
    REQUIRE(integralCatalog[stars[0]].name == 58);
}

TEST_CASE("IdentifyThirdStar no third star", "[identify-remaining] [fast]") {
    std::vector<int16_t> stars = IdentifyThirdStar(integralCatalog,
                                                    Vec3(0,0,0), Vec3(0,1,0),
                                                    0.5, sqrt(2.0),
                                                    1e-6);
    REQUIRE(stars.size() == 0);
}

TEST_CASE("IdentifyThirdStar just out of tolerance", "[identify-remaining] [fast]") {
    std::vector<int16_t> stars2 = IdentifyThirdStar(integralCatalog,
                                                    Vec3(0,0,0), Vec3(0,1,0),
                                                    0.99, sqrt(2.0),
                                                    1e-6);
    REQUIRE(stars2.size() == 0);
}

