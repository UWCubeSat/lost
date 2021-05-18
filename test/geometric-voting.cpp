#include "catch.hpp"

#include <vector>

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

Camera smolCamera(FovToFocalLength(DegToRad(36.0), 256), 256, 256);
Quaternion straightAhead(1.0, 0.0, 0.0, 0.0);

TEST_CASE("Tres Commas", "[gv] [fast]") {
    // TODO: Check that both stars receive exactly 2 votes, there's no reason for them to get 1 vote each!
    GeneratedPipelineInput gpi(tripleCatalog, straightAhead, smolCamera, 0, 0.0, 0.0);
    unsigned char *db = BuildKVectorDatabase(tripleCatalog, NULL, DegToRad(.5), DegToRad(20), 1000);
    GeometricVotingStarIdAlgorithm *gv = new GeometricVotingStarIdAlgorithm(DegToRad(0.001));
    Pipeline pipeline(NULL, gv, NULL, db);

    PipelineOutput output = pipeline.Go(gpi);
    StarIdComparison comparison = StarIdsCompare(*gpi.ExpectedStarIds(), *output.starIds, 0.0, gpi.ExpectedStars(), NULL);
    REQUIRE(comparison.numTotal == 3);
    REQUIRE(comparison.numCorrect == 3);
}
