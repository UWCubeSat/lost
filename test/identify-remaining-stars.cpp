#include <catch.hpp>

#include "star-id.cpp"

using namespace lost; // NOLINT

Stars stars = {
    Star(0, 0, 1), // 0
    Star(1, 0, 1), // 1
    Star(0, 1, 1), // 2
    Star(1, 1, 1), // 3
    Star(2, 0, 1), // 4
    Star(2, 1, 1), // 5
    Star(2, 2, 1), // 6
    Star(3, 0, 1), // 7
    Star(3, 1, 1), // 8
    Star(3, 2, 1), // 9
    Star(3, 3, 1), // 10
    Star(1, 3, 1), // 11
};

StarIdentifiers starIdentifiers = {
    StarIdentifier(0, 0, 1),
    StarIdentifier(1, 1, 1),
    StarIdentifier(2, 2, 1),
    StarIdentifier(3, 3, 1),
    StarIdentifier(4, 4, 1),
    StarIdentifier(5, 5, 1),
    StarIdentifier(6, 6, 1),
    StarIdentifier(7, 7, 1),
    StarIdentifier(8, 8, 1),
    StarIdentifier(9, 9, 1),
    StarIdentifier(10, 10, 1),
    StarIdentifier(11, 11, 1),
};

TEST_CASE("IRUnidentifiedCentroid simple orthogonal", "[identify-remaining]") {
    IRUnidentifiedCentroid centroid(stars[3]);
    REQUIRE(centroid.bestAngleFrom90 == -1);
    centroid.AddIdentifiedStar(starIdentifiers[1], stars);
    // one star is not enough to get angle
    REQUIRE(centroid.bestAngleFrom90 == -1);
    centroid.AddIdentifiedStar(starIdentifiers[2], stars);
    // we've set them up to be almost orthogonal
    REQUIRE(centroid.bestAngleFrom90 == Approx(0));
    REQUIRE(centroid.bestStar1 == starIdentifiers[1] && centroid.bestStar2 == starIdentifiers[2] ||
            centroid.bestStar1 == starIdentifiers[2] && centroid.bestStar2 == starIdentifiers[1]);

    // adding another, non-orthogonal one shouldn't break it
    centroid.AddIdentifiedStar(starIdentifiers[0], stars);
    REQUIRE(centroid.bestAngleFrom90 == Approx(0));
    REQUIRE(centroid.bestStar1 == starIdentifiers[1] && centroid.bestStar2 == starIdentifiers[2] ||
            centroid.bestStar1 == starIdentifiers[2] && centroid.bestStar2 == starIdentifiers[1]);
}

