#include <catch.hpp>

#include <math.h>
#include <iostream>

#include "camera.hpp"
#include "attitude-utils.hpp"

using namespace lost;

TEST_CASE("Convert coordinates: pixel -> spatial -> pixel", "[geometry]") {
    Camera camera(100, 512, 1024);

    float expectedX = GENERATE(142, 90, 512, 255);
    float expectedY = GENERATE(18, 512, 0, 800);

    Vec3 spatial = camera.CameraToSpatial({expectedX, expectedY});
    Vec2 actualPixels = camera.SpatialToCamera(spatial);

    CHECK((int)round(actualPixels.x) == expectedX);
    CHECK((int)round(actualPixels.y) == expectedY);
}

TEST_CASE("Centered coordinates: pixel -> spatial", "[geometry]") {
    Camera camera(100, 512, 1024);

    Vec3 spatial = camera.CameraToSpatial({256, 512});
    CHECK(spatial.y == Approx(0.0));
    CHECK(spatial.z == Approx(0.0));
}

TEST_CASE("Convert coordinates explicit: pixel -> spatial", "[geometry]") {
    Camera camera(100, 512, 1024);

    Vec3 spatial = camera.CameraToSpatial({300, 728}).Normalize();
    CHECK(spatial.x == Approx(0.413).margin(0.001));
    CHECK(spatial.y == Approx(-0.182).margin(0.001));
    CHECK(spatial.z == Approx(-0.892).margin(0.001));
}

// the passing of the above two tests implies spatial -> pixel works too

TEST_CASE("Angle from camera", "[geometry]") {
    Camera camera(128, 256, 256);

    SECTION("horizontal line") {
        Vec3 s1 = camera.CameraToSpatial({0, 128});
        Vec3 s2 = camera.CameraToSpatial({256,128});
        CHECK(Angle(s1, s2) == Approx(M_PI / 2.0));
    }

    SECTION("vertical line") {
        Vec3 s1 = camera.CameraToSpatial({128, 256});
        Vec3 s2 = camera.CameraToSpatial({128,0});
        CHECK(Angle(s1, s2) == Approx(M_PI / 2.0));
    }
}

TEST_CASE("Angle from camera, diagonal", "[geometry]") {
    Camera camera(128, 128, 128);

    Vec3 s1 = camera.CameraToSpatial({0, 0});
    Vec3 s2 = camera.CameraToSpatial({128, 128});
    // it's not as simple as this! TODO
    // CHECK(Angle(s1, s2) == Approx(acos(0.25)));
}
