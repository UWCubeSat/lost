#include <algorithm>
#include <cstring>

#include <catch.hpp>

#include "centroiders.hpp"

using namespace lost; // NOLINT

TEST_CASE("CoG: single symmetric blob", "[centroid]") {
    // Single 3x3 blob of equal-brightness pixels about the center
    const int w = 32, h = 32;
    const int cx = 15, cy = 15;
    unsigned char image[w * h];
    memset(image, 0, sizeof(image));

    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            image[(cy + dy) * w + (cx + dx)] = 255;
        }
    }

    CenterOfGravityAlgorithm algo;
    Stars stars = algo.Go(image, w, h);

    REQUIRE(stars.size() == 1);
    CHECK(stars[0].position.x == Approx(cx + 0.5).margin(0.01));
    CHECK(stars[0].position.y == Approx(cy + 0.5).margin(0.01));
    CHECK(stars[0].magnitude == 9);
}

TEST_CASE("CoG: weighted blob shifts toward brighter pixels", "[centroid]") {
    const int w = 32, h = 32;
    const int cx = 15, cy = 15;
    const int left = 150, mid = 200, right = 255;
    unsigned char image[w * h];
    memset(image, 0, sizeof(image));

    image[cy * w + (cx - 1)] = left;
    image[cy * w + cx] = mid;
    image[cy * w + (cx + 1)] = right;

    CenterOfGravityAlgorithm algo;
    Stars stars = algo.Go(image, w, h);

    decimal expectedX =
        DECIMAL((cx - 1) * left + cx * mid + (cx + 1) * right) / DECIMAL(left + mid + right) +
        DECIMAL(0.5);

    REQUIRE(stars.size() == 1);
    CHECK(stars[0].position.x == Approx(expectedX).margin(0.01));
    CHECK(stars[0].position.y == Approx(cy + 0.5).margin(0.01));
    CHECK(stars[0].position.x > cx + DECIMAL(0.5));
}

TEST_CASE("CoG: two separate blobs", "[centroid]") {
    const int w = 32, h = 32;
    const int ax = 5, ay = 5;
    const int bx = 25, by = 25;
    unsigned char image[w * h];
    memset(image, 0, sizeof(image));

    for (int dy = 0; dy <= 1; dy++) {
        for (int dx = 0; dx <= 1; dx++) {
            image[(ay + dy) * w + (ax + dx)] = 255;
            image[(by + dy) * w + (bx + dx)] = 255;
        }
    }

    CenterOfGravityAlgorithm algo;
    Stars stars = algo.Go(image, w, h);

    REQUIRE(stars.size() == 2);

    std::sort(stars.begin(), stars.end(), [](const Star& a, const Star& b) {
        return a.position.x < b.position.x;
    });

    CHECK(stars[0].position.x == Approx(ax + 1).margin(0.01));
    CHECK(stars[0].position.y == Approx(ay + 1).margin(0.01));
    CHECK(stars[1].position.x == Approx(bx + 1).margin(0.01));
    CHECK(stars[1].position.y == Approx(by + 1).margin(0.01));
}

TEST_CASE("CoG: edge blob rejected", "[centroid]") {
    const int w = 32, h = 32;
    // Corner blob should be rejected since it touches the edge of the image
    const int cornerBlobX = 0, cornerBlobY = 0;
    const int centerBlobX = 15, centerBlobY = 15;
    unsigned char image[w * h];
    memset(image, 0, sizeof(image));

    for (int dy = 0; dy <= 1; dy++) {
        for (int dx = 0; dx <= 1; dx++) {
            image[(cornerBlobY + dy) * w + (cornerBlobX + dx)] = 255;
            image[(centerBlobY + dy) * w + (centerBlobX + dx)] = 255;
        }
    }

    CenterOfGravityAlgorithm algo;
    Stars stars = algo.Go(image, w, h);

    REQUIRE(stars.size() == 1);
    CHECK(stars[0].position.x == Approx(centerBlobX + 1).margin(0.01));
    CHECK(stars[0].position.y == Approx(centerBlobY + 1).margin(0.01));
}

TEST_CASE("CoG: all dark image", "[centroid]") {
    const int w = 32, h = 32;
    unsigned char image[w * h];
    memset(image, 0, sizeof(image));

    CenterOfGravityAlgorithm algo;
    Stars stars = algo.Go(image, w, h);

    CHECK(stars.empty());
}
