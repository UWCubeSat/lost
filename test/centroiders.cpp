#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

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

// WINDOWED CENTER OF GRAVITY

namespace {

/// Add a Gaussian star of the given total brightness centered at (x, y) in pixel coordinates.
void DrawStar(std::vector<unsigned char> *image, int width, int height, decimal x, decimal y, decimal peak) {
    const decimal sigma = DECIMAL(1.0);
    for (int py = 0; py < height; py++) {
        for (int px = 0; px < width; px++) {
            decimal dx = px + DECIMAL(0.5) - x, dy = py + DECIMAL(0.5) - y;
            decimal value = (*image)[py * width + px] + peak * DECIMAL_EXP(-(dx * dx + dy * dy) / (2 * sigma * sigma));
            (*image)[py * width + px] = (unsigned char)std::min(DECIMAL(255.0), DECIMAL_ROUND(value));
        }
    }
}

std::vector<unsigned char> BlankImage(int width, int height) {
    return std::vector<unsigned char>(width * height, 10);
}

}  // namespace

TEST_CASE("Windowed CoG finds a star near its prediction to subpixel accuracy", "[centroid]") {
    const int w = 64, h = 64;
    std::vector<unsigned char> image = BlankImage(w, h);
    DrawStar(&image, w, h, DECIMAL(30.3), DECIMAL(25.7), 200);

    WindowedCenterOfGravity centroider(4, 4, 3);
    std::vector<int> matched;
    Stars stars = centroider.Go(image.data(), w, h, {{DECIMAL(31.5), DECIMAL(24.5)}}, &matched);
    REQUIRE(stars.size() == 1);
    CHECK(matched == std::vector<int>{0});
    CHECK(stars[0].position.x == Approx(30.3).margin(0.05));
    CHECK(stars[0].position.y == Approx(25.7).margin(0.05));
}

TEST_CASE("Windowed CoG ignores stars outside the gate", "[centroid]") {
    const int w = 64, h = 64;
    std::vector<unsigned char> image = BlankImage(w, h);
    DrawStar(&image, w, h, DECIMAL(40.0), DECIMAL(32.0), 250);  // bright, 8 px from the prediction
    DrawStar(&image, w, h, DECIMAL(32.5), DECIMAL(32.5), 60);   // faint, at the prediction

    WindowedCenterOfGravity centroider(4, 4, 3);
    std::vector<int> matched;
    Stars stars = centroider.Go(image.data(), w, h, {{DECIMAL(32.0), DECIMAL(32.0)}}, &matched);
    REQUIRE(stars.size() == 1);
    CHECK(stars[0].position.x == Approx(32.5).margin(0.1));

    std::vector<unsigned char> onlyFar = BlankImage(w, h);
    DrawStar(&onlyFar, w, h, DECIMAL(40.0), DECIMAL(32.0), 250);
    CHECK(centroider.Go(onlyFar.data(), w, h, {{DECIMAL(32.0), DECIMAL(32.0)}}, &matched).empty());

    // the edge of this star reaches into the gate, but its centroid is outside it
    std::vector<unsigned char> justOutside = BlankImage(w, h);
    DrawStar(&justOutside, w, h, DECIMAL(37.5), DECIMAL(32.0), 200);
    CHECK(centroider.Go(justOutside.data(), w, h, {{DECIMAL(32.0), DECIMAL(32.0)}}, &matched).empty());
}

TEST_CASE("Windowed CoG drops a star claimed by two predictions", "[centroid]") {
    const int w = 64, h = 64;
    std::vector<unsigned char> image = BlankImage(w, h);
    DrawStar(&image, w, h, DECIMAL(32.0), DECIMAL(32.0), 200);
    DrawStar(&image, w, h, DECIMAL(12.0), DECIMAL(12.0), 200);

    WindowedCenterOfGravity centroider(4, 4, 3);
    std::vector<int> matched;
    Stars stars = centroider.Go(image.data(), w, h,
                                {{DECIMAL(31.0), DECIMAL(32.0)}, {DECIMAL(12.0), DECIMAL(12.0)}, {DECIMAL(33.0), DECIMAL(32.0)}},
                                &matched);
    CHECK(matched == std::vector<int>{1});
}

TEST_CASE("Windowed CoG rejects a star clipped by the window or the image", "[centroid]") {
    const int w = 64, h = 64;
    std::vector<unsigned char> image = BlankImage(w, h);
    // a long streak: any window around it cuts it off
    for (int x = 10; x < 54; x++) image[32 * w + x] = 200;
    WindowedCenterOfGravity centroider(4, 4, 3);
    std::vector<int> matched;
    CHECK(centroider.Go(image.data(), w, h, {{DECIMAL(32.0), DECIMAL(32.5)}}, &matched).empty());

    std::vector<unsigned char> edge = BlankImage(w, h);
    DrawStar(&edge, w, h, DECIMAL(1.0), DECIMAL(32.0), 200);
    CHECK(centroider.Go(edge.data(), w, h, {{DECIMAL(1.0), DECIMAL(32.0)}}, &matched).empty());
}

TEST_CASE("Windowed CoG drops a match that disagrees with the common shift", "[centroid]") {
    const int w = 128, h = 128;
    std::vector<unsigned char> image = BlankImage(w, h);
    std::vector<Vec2> predictions;
    for (int i = 0; i < 4; i++) {
        // every star 1 px right of its prediction, except the last, which is 3.5 px left of it
        decimal shift = i < 3 ? DECIMAL(1.0) : DECIMAL(-3.5);
        DrawStar(&image, w, h, 20 + 25 * i + shift, 64, 200);
        predictions.push_back({DECIMAL(20.0) + 25 * i, DECIMAL(64.0)});
    }
    WindowedCenterOfGravity centroider(5, 4, 3);
    std::vector<int> matched;
    centroider.Go(image.data(), w, h, predictions, &matched);
    CHECK(matched == std::vector<int>{0, 1, 2});
}

TEST_CASE("Windowed CoG handles black images and predictions off the image", "[centroid]") {
    const int w = 32, h = 32;
    std::vector<unsigned char> black(w * h, 0);
    WindowedCenterOfGravity centroider(4, 4, 3);
    std::vector<int> matched;
    CHECK(centroider.Go(black.data(), w, h, {{DECIMAL(16.0), DECIMAL(16.0)}, {DECIMAL(-50.0), DECIMAL(500.0)}}, &matched).empty());

    // predictions that aren't even numbers, and an image smaller than the window
    const decimal nan = std::numeric_limits<decimal>::quiet_NaN();
    const decimal inf = std::numeric_limits<decimal>::infinity();
    std::vector<unsigned char> tiny(2 * 2, 255);
    CHECK(centroider.Go(tiny.data(), 2, 2, {{nan, nan}, {inf, -inf}, {DECIMAL(1e30), 1}, {1, 1}}, &matched).empty());
}

TEST_CASE("Windowed CoG finds nothing with unusable settings, and accepts a null output", "[centroid]") {
    const int w = 64, h = 64;
    std::vector<unsigned char> image = BlankImage(w, h);
    DrawStar(&image, w, h, DECIMAL(32.0), DECIMAL(32.0), 200);
    const std::vector<Vec2> predictions = {{DECIMAL(32.0), DECIMAL(32.0)}};

    std::vector<int> matched = {7};
    CHECK(WindowedCenterOfGravity(-1, 4, 3).Go(image.data(), w, h, predictions, &matched).empty());
    CHECK(matched.empty());
    CHECK(WindowedCenterOfGravity(4, -1, 3).Go(image.data(), w, h, predictions, &matched).empty());
    CHECK(WindowedCenterOfGravity(1 << 30, 4, 3).Go(image.data(), w, h, predictions, &matched).empty());
    CHECK(WindowedCenterOfGravity(4, 1 << 30, 3).Go(image.data(), w, h, predictions, &matched).empty());
    CHECK(WindowedCenterOfGravity(4, 4, -1).Go(image.data(), w, h, predictions, &matched).empty());
    CHECK(WindowedCenterOfGravity(4, 4, std::numeric_limits<decimal>::quiet_NaN())
          .Go(image.data(), w, h, predictions, &matched).empty());
    CHECK(WindowedCenterOfGravity(4, 4, 3).Go(nullptr, w, h, predictions, &matched).empty());
    CHECK(WindowedCenterOfGravity(4, 4, 3).Go(image.data(), 0, 0, predictions, &matched).empty());
    CHECK(WindowedCenterOfGravity(4, 4, 3).Go(image.data(), -5, 7, predictions, &matched).empty());

    // a window larger than the image is clipped to it
    CHECK(WindowedCenterOfGravity(WindowedCenterOfGravity::kMaxRadius, 4, 3).Go(image.data(), w, h, predictions, &matched).size() == 1);

    CHECK(WindowedCenterOfGravity(4, 4, 3).Go(image.data(), w, h, predictions, nullptr).size() == 1);
}

TEST_CASE("Windowed CoG drops a star claimed by two predictions whose thresholds differ", "[centroid]") {
    const int w = 64, h = 64;
    // background brighter on the left, so that the prediction left of the star gets a higher threshold
    std::vector<unsigned char> image(w * h);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) image[y * w + x] = x < 30 ? 30 : 10;
    }
    for (int y = 31; y <= 33; y++) {
        for (int x = 31; x <= 33; x++) image[y * w + x] = 200;
    }
    image[30 * w + 32] = 25;  // part of the star for the right-hand prediction only

    WindowedCenterOfGravity centroider(4, 4, 3);
    std::vector<int> matched;
    // each prediction alone finds the star, with a different set of pixels
    Stars left = centroider.Go(image.data(), w, h, {{DECIMAL(29.5), DECIMAL(32.5)}}, &matched);
    Stars right = centroider.Go(image.data(), w, h, {{DECIMAL(34.5), DECIMAL(32.5)}}, &matched);
    REQUIRE(left.size() == 1);
    REQUIRE(right.size() == 1);
    CHECK(left[0].magnitude != right[0].magnitude);

    // together, the star is ambiguous, while a star with a single prediction is kept
    for (int y = 49; y <= 51; y++) {
        for (int x = 49; x <= 51; x++) image[y * w + x] = 200;
    }
    Stars stars = centroider.Go(image.data(), w, h,
                                {{DECIMAL(29.5), DECIMAL(32.5)}, {DECIMAL(50.5), DECIMAL(50.5)}, {DECIMAL(34.5), DECIMAL(32.5)}},
                                &matched);
    CHECK(matched == std::vector<int>{1});
}
