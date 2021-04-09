#include <catch.hpp>

#include <math.h>
#include <iostream>

#include "camera.hpp"
#include "attitude-utils.hpp"

using namespace lost;

TEST_CASE("Convert coordinates: pixel -> spatial -> pixel", "[camera]") {
    Camera camera(100, 512, 1024);

    float expectedX = GENERATE(142, 90, 512, 255);
    float expectedY = GENERATE(18, 512, 0, 800);

    Vec3 spatial = camera.CameraToSpatial({expectedX, expectedY});
    Vec2 actualPixels = camera.SpatialToCamera(spatial);

    CHECK((int)round(actualPixels.x) == expectedX);
    CHECK((int)round(actualPixels.y) == expectedY);
}

TEST_CASE("Convert coordinates explicit: pixel -> spatial", "[camera]") {
    Camera camera(100, 512, 1024);

    Vec3 spatial = camera.CameraToSpatial({300, 728}).Normalize();
    CHECK(spatial.x == Approx(0.413).margin(0.001));
    CHECK(spatial.y == Approx(-0.182).margin(0.001));
    CHECK(spatial.z == Approx(-0.892).margin(0.001));
}
