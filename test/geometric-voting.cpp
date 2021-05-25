#include "catch.hpp"

#include <vector>
#include <stdlib.h>
#include <stdio.h>

#include "databases.hpp"
#include "star-id.hpp"
#include "star-utils.hpp"
#include "camera.hpp"
#include "io.hpp"

using namespace lost;

// three stars that should be easy to identify based on their distances to each other
Catalog tripleCatalog = {
    CatalogStar(DegToRad(2), DegToRad(-3), 3.0, false, 42),
    CatalogStar(DegToRad(4), DegToRad(7), 2.0, false, 43),
    CatalogStar(DegToRad(2), DegToRad(6), 4.0, false, 44),
};

Catalog quadrupleCatalog = {
    CatalogStar(DegToRad(2), DegToRad(-3), 3.0, false, 42),
    CatalogStar(DegToRad(4), DegToRad(7), 2.0, false, 43),
    CatalogStar(DegToRad(2), DegToRad(6), 4.0, false, 44),
    CatalogStar(DegToRad(-1), DegToRad(-4), 1.0, false, 45),
};

// distance between 42/43 == distance between 44/45, have to use distances to 43 to differentiate
Catalog harderQuadrupleCatalog = {
    CatalogStar(DegToRad(2), DegToRad(-3), 3.0, false, 42),
    CatalogStar(DegToRad(4), DegToRad(7), 2.0, false, 43),
    CatalogStar(DegToRad(2), DegToRad(5), 4.0, false, 44),
    CatalogStar(DegToRad(2), DegToRad(1), 1.0, false, 45),
};

Camera smolCamera(FovToFocalLength(DegToRad(36.0), 256), 256, 256);
Camera smolCameraOff(FovToFocalLength(DegToRad(33.0), 256), 256, 256);
Quaternion straightAhead(1.0, 0.0, 0.0, 0.0);

TEST_CASE("Tres Commas", "[gv] [fast]") {
    // TODO: Check that both stars receive exactly 2 votes, there's no reason for them to get 1 vote each!
    const Catalog &catalog = GENERATE((const Catalog &)tripleCatalog,
                                      (const Catalog &)quadrupleCatalog,
                                      (const Catalog &)harderQuadrupleCatalog);
    const Camera &camera = GENERATE((const Camera &)smolCamera,
                                    (const Camera &)smolCameraOff);
    GeneratedPipelineInput gpi(catalog, straightAhead, camera, 0, 0.0, 0.0);
    unsigned char *db = BuildKVectorDatabase(catalog, NULL, DegToRad(.5), DegToRad(20), 1000);
    StarIdentifiers outputStarIds = GeometricVotingStarIdAlgorithm(DegToRad(0.01)).Go(db, *gpi.InputStars(), catalog, camera);

    StarIdComparison comparison = StarIdsCompare(*gpi.ExpectedStarIds(), outputStarIds, 0.0, gpi.ExpectedStars(), NULL);
    REQUIRE(comparison.numTotal == catalog.size());
    CHECK(comparison.numCorrect == catalog.size());
}

TEST_CASE("Real catalog tests", "[gv]") {
    const Catalog &catalog = CatalogRead();
    const Camera &camera = GENERATE((const Camera &)smolCamera,
                                    (const Camera &)smolCameraOff);
    float ra = GENERATE(DegToRad(88.0), DegToRad(24.2));
    float de = GENERATE(DegToRad(3.9), DegToRad(33.9));
    float roll = GENERATE(DegToRad(84.5), DegToRad(66.8));
    Quaternion attitude = SphericalToQuaternion(ra, de, roll);
    printf("GV test: %s, ra=%f, de=%f, roll=%f\n",
        &camera == &smolCamera ? "normal camera" : "offset camera",
        ra, de, roll);
    GeneratedPipelineInput gpi(catalog, attitude, camera, 0, 0.0, 0.0);
    unsigned char *db = BuildKVectorDatabase(catalog, NULL, DegToRad(.5), DegToRad(20), 1000);
    StarIdentifiers outputStarIds = GeometricVotingStarIdAlgorithm(DegToRad(0.01)).Go(db, *gpi.InputStars(), catalog, camera);

    StarIdComparison comparison = StarIdsCompare(*gpi.ExpectedStarIds(), outputStarIds, 0.0, gpi.ExpectedStars(), NULL);
    
    REQUIRE(comparison.numTotal == gpi.InputStars()->size());
    CHECK(comparison.numCorrect >= gpi.InputStars()->size() * 7 / 10);
    CHECK(comparison.numIncorrect == 0);
}

// TEST_CASE("Wrong focal length", "[gv] [fast]") {
//     GeneratedPipelineInput gpi(CatalogRead(), straightAhead, smolCamera, 0, 0.0, 0.0);
//     unsigned char *db = BuildKVectorDatabase(CatalogRead(), NULL, DegToRad(.5), DegToRad(20), 1000);
//     GeometricVotingStarIdAlgorithm *gv = new GeometricVotingStarIdAlgorithm(DegToRad(0.001));
//     gv->Go()
// }
