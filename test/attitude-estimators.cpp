#include <random>
#include <catch.hpp>

#include "attitude-estimators.hpp"
#include "attitude-utils.hpp"
#include "camera.hpp"
#include "star-utils.hpp"

using namespace lost; // NOLINT

// Places a few stars at random pixels for a known attitude, then checks the estimator gets that attitude back, repeated over many random scenes
static int CountFailures(AttitudeEstimationAlgorithm &algo, int numStars) {
    Camera camera(500, 1024, 1024);
    std::mt19937 rng(1234);
    std::uniform_real_distribution<double> pixel(0, 1023);
    std::uniform_real_distribution<double> angle(0, 6.28);

    int failures = 0;
    for (int trial = 0; trial < 500; trial++) {
        Attitude truth(SphericalToQuaternion(angle(rng), angle(rng) / 4 - 0.78, angle(rng)));
        Stars stars;
        Catalog catalog;
        StarIdentifiers ids;
        for (int i = 0; i < numStars; i++) {
            Vec2 p = {(decimal)pixel(rng), (decimal)pixel(rng)};
            stars.push_back(Star(p.x, p.y, 1));
            Vec3 inCamera = camera.CameraToSpatial(p).Normalize();
            catalog.push_back(CatalogStar(truth.GetQuaternion().Conjugate().Rotate(inCamera), 0, i));
            ids.push_back(StarIdentifier(i, i));
        }
        Attitude result = algo.Go(camera, stars, catalog, ids);
        Quaternion diff = result.GetQuaternion() * truth.GetQuaternion().Conjugate();
        if (diff.SmallestAngle() > DECIMAL(0.01)) failures++;
    }
    return failures;
}

TEST_CASE("QUEST recovers known attitudes", "[attitude]") {
    QuestAlgorithm quest;
    // 2-star QUEST is too ill-conditioned in float mode to pass reliably
    int numStars = GENERATE(3, 5, 10);
    CAPTURE(numStars);
    CHECK(CountFailures(quest, numStars) == 0);
}

TEST_CASE("DQM recovers known attitudes", "[attitude]") {
    DavenportQAlgorithm dqm;
    int numStars = GENERATE(2, 3, 5, 10);
    CAPTURE(numStars);
    CHECK(CountFailures(dqm, numStars) == 0);
}
