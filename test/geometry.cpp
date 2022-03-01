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

    CHECK(actualPixels.x == Approx(expectedX).margin(0.000001));
    CHECK(actualPixels.y == Approx(expectedY).margin(0.000001));
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

TEST_CASE("spherical -> quaternion -> spherical", "[geometry]") {
    // 0.1 instead of 0, because at 0 it might sometimes return 2PI, which is fine for most
    // circumstances. Also, at 0, the epsilon=0, so there's no tolerance in Approx by default!
    float ra = DegToRad(GENERATE(28.9, 83.2, 14.0, 0.1, 329.8));
    float de = DegToRad(GENERATE(7.82, 9.88, 88.8, 0.1, -72.0, -9.9));
    float roll = DegToRad(GENERATE(9.38, 300.9, 37.8, 199.9));

    Quaternion quat = SphericalToQuaternion(ra, de, roll);

    EulerAngles angles = quat.ToSpherical();
    // TODO: for small angles, the error is quite large as a fraction of the angle.
    CHECK(angles.ra == Approx(ra).margin(0.00001));
    CHECK(angles.de == Approx(de).margin(0.00001));
    CHECK(angles.roll == Approx(roll).margin(0.00001));
}

TEST_CASE("spherical -> spatial -> spherical", "[geometry]") {
    float ra = DegToRad(GENERATE(28.9, 83.2, 14.0, 0.1, 329.8));
    float de = DegToRad(GENERATE(7.82, 9.88, 88.8, 0.1, -72.0, -9.9));

    float raOut, deOut;
    SpatialToSpherical(SphericalToSpatial(ra, de), &raOut, &deOut);

    CHECK(ra == Approx(raOut));
    CHECK(de == Approx(deOut));
}

TEST_CASE("quat -> dcm -> quat", "[geometry]") {
    float ra = GENERATE(take(5, random(0.1, 3.14*2)));
    float de = GENERATE(take(5, random(-3.14, 3.14)));
    float roll = GENERATE(take(5, random(0.1, 3.14*2)));

    Quaternion quat1 = SphericalToQuaternion(ra, de, roll).Canonicalize();
    Mat3 dcm = QuaternionToDCM(quat1);
    Quaternion quat2 = DCMToQuaternion(dcm).Canonicalize();
    CHECK(quat1.real == Approx(quat2.real));
    CHECK(quat1.i == Approx(quat2.i));
    CHECK(quat1.j == Approx(quat2.j));
    CHECK(quat1.k == Approx(quat2.k));
}

// I know cross product seems simple, perhaps even too simple to be worth testing...but I coded it
// wrong initially and missed the bug for months.
TEST_CASE("cross product", "[geometry]") {
    Vec3 x = { 1, 0, 0 };
    Vec3 y = { 0, 1, 0 };
    Vec3 z = x.crossProduct(y);
    CHECK(z.x == 0);
    CHECK(z.y == 0);
    CHECK(z.z == Approx(1.0));

    Vec3 a = { 5, 8, -2 };
    Vec3 b = { -5.9, 92, 3 };
    Vec3 c = a.crossProduct(b);
    CHECK(c.x == 208);
    CHECK(c.y == Approx(-3.2));
    CHECK(c.z == Approx(507.2));
}
