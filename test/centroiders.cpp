#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

#include <cmath>
#include <random>
#include <utility>
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

/// Add a Gaussian star of the given peak level centered at (x, y) in pixel coordinates.
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

typedef std::vector<unsigned char> Image;

/// Set the pixels x0..x1, y0..y1 (inclusive) to value.
void Fill(Image *image, int width, int x0, int x1, int y0, int y1, int value) {
    for (int y = y0; y <= y1; y++) {
        for (int x = x0; x <= x1; x++) (*image)[y * width + x] = (unsigned char)value;
    }
}

/// Paint the one-pixel border of the rectangle x0..x1, y0..y1 with the given (level, count) runs.
void PaintBorder(Image *image, int width, int x0, int x1, int y0, int y1,
                 const std::vector<std::pair<int, int>> &runs) {
    std::vector<int> border;
    for (int x = x0; x <= x1; x++) border.push_back(y0 * width + x);
    for (int x = x0; x <= x1; x++) border.push_back(y1 * width + x);
    for (int y = y0 + 1; y < y1; y++) border.push_back(y * width + x0);
    for (int y = y0 + 1; y < y1; y++) border.push_back(y * width + x1);
    size_t next = 0;
    for (const std::pair<int, int> &run : runs) {
        for (int i = 0; i < run.second; i++) (*image)[border.at(next++)] = (unsigned char)run.first;
    }
    REQUIRE(next == border.size());
}

/// Add a Gaussian star, without rounding, to a floating point image.
void AddStar(std::vector<double> *image, int width, int height, double x, double y, double peak, double sigma) {
    for (int py = 0; py < height; py++) {
        for (int px = 0; px < width; px++) {
            const double dx = px + 0.5 - x, dy = py + 0.5 - y;
            (*image)[py * width + px] += peak * std::exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
        }
    }
}

// The window of a prediction at (32.5, 32.5) with gateRadius 4 and maxStarRadius 4 is 24..40 in x and y,
// and its border has 64 pixels.
const Vec2 kCenter = {DECIMAL(32.5), DECIMAL(32.5)};

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

TEST_CASE("Windowed CoG reports the exact centroid, size and pixel count of a known blob", "[centroid]") {
    const int w = 64, h = 64;
    Image image(w * h, 100);
    image[40 * w + 30] = 110;
    image[40 * w + 31] = 200;
    image[40 * w + 32] = 104;  // exactly at the threshold, 4 above the background: part of the star
    image[41 * w + 31] = 150;
    image[39 * w + 30] = 103;  // just below the threshold: not part of it
    image[41 * w + 33] = 120;  // touches the star only diagonally: not part of it

    WindowedCenterOfGravity centroider(4, 4, 3);
    std::vector<int> matched;
    Stars stars = centroider.Go(image.data(), w, h, {{DECIMAL(31.5), DECIMAL(40.5)}}, &matched);
    REQUIRE(stars.size() == 1);
    CHECK(matched == std::vector<int>{0});
    // weights are the levels above the background
    const double weightSum = 10 + 100 + 4 + 50;
    CHECK(stars[0].position.x == Approx((30 * 10 + 31 * 100 + 32 * 4 + 31 * 50) / weightSum + 0.5).margin(1e-4));
    CHECK(stars[0].position.y == Approx((40 * 114 + 41 * 50) / weightSum + 0.5).margin(1e-4));
    CHECK(stars[0].radiusX == Approx(1.5));
    CHECK(stars[0].radiusY == Approx(1.0));
    CHECK(stars[0].magnitude == 4);
}

TEST_CASE("Windowed CoG threshold is 4 levels above a noiseless background", "[centroid]") {
    // MAD 0, so the threshold is the background plus 5 * 1.4826 * 0.5 = 3.7, rounded up
    const int w = 64, h = 64;
    WindowedCenterOfGravity centroider(4, 4, 3);
    Image image(w * h, 100);
    image[32 * w + 32] = 103;
    CHECK(centroider.Go(image.data(), w, h, {kCenter}, nullptr).empty());
    image[32 * w + 32] = 104;
    CHECK(centroider.Go(image.data(), w, h, {kCenter}, nullptr).size() == 1);
}

TEST_CASE("Windowed CoG threshold is the border's median plus 5 * 1.4826 * (MAD + 0.5)", "[centroid]") {
    const int w = 64, h = 64;
    WindowedCenterOfGravity centroider(4, 4, 3);

    SECTION("spread border") {
        // median 20; deviations 0 (13 pixels), 3 (26), 6 (25): MAD 3; threshold 20 + 7.413 * 3.5 = 45.9
        Image image(w * h, 20);
        PaintBorder(&image, w, 24, 40, 24, 40, {{14, 13}, {17, 13}, {20, 13}, {23, 13}, {26, 12}});
        Fill(&image, w, 31, 33, 31, 33, 45);
        CHECK(centroider.Go(image.data(), w, h, {kCenter}, nullptr).empty());
        Fill(&image, w, 31, 33, 31, 33, 46);
        Stars stars = centroider.Go(image.data(), w, h, {kCenter}, nullptr);
        REQUIRE(stars.size() == 1);
        CHECK(stars[0].magnitude == 9);
    }

    SECTION("a star crossing the border does not raise the threshold, also on a nearly black background") {
        // median 3; deviations 0 (10 pixels), 3 (50, the 30 at level 0 among them), 97 (4): MAD 3;
        // threshold 3 + 7.413 * 3.5 = 28.9
        Image image(w * h, 3);
        PaintBorder(&image, w, 24, 40, 24, 40, {{0, 30}, {3, 10}, {6, 20}, {100, 4}});
        Fill(&image, w, 31, 33, 31, 33, 28);
        CHECK(centroider.Go(image.data(), w, h, {kCenter}, nullptr).empty());
        Fill(&image, w, 31, 33, 31, 33, 29);
        CHECK(centroider.Go(image.data(), w, h, {kCenter}, nullptr).size() == 1);
    }

    SECTION("a border split evenly between two levels: the median is their mean, the MAD half their distance") {
        // 32 border pixels at 30 (the left column and 15 of the top row), 32 at 10: median 20, MAD 10;
        // threshold 20 + 7.413 * 10.5 = 97.8
        Image image(w * h, 10);
        Fill(&image, w, 24, 24, 24, 40, 30);
        Fill(&image, w, 26, 40, 24, 24, 30);
        Fill(&image, w, 31, 33, 31, 33, 97);
        CHECK(centroider.Go(image.data(), w, h, {kCenter}, nullptr).empty());
        Fill(&image, w, 31, 33, 31, 33, 98);
        CHECK(centroider.Go(image.data(), w, h, {kCenter}, nullptr).size() == 1);
    }

    SECTION("an eight pixel border with a median and MAD between gray levels") {
        // a 3 by 3 window: four border pixels at 10 and four at 21. Median 15.5, MAD 5.5;
        // threshold 15.5 + 7.413 * 6 = 59.98
        WindowedCenterOfGravity small(1, 0, 3);
        Image image(w * h, 10);
        image[31 * w + 31] = image[31 * w + 32] = image[31 * w + 33] = image[32 * w + 31] = 21;
        image[32 * w + 32] = 59;
        CHECK(small.Go(image.data(), w, h, {kCenter}, nullptr).empty());
        image[32 * w + 32] = 60;
        Stars stars = small.Go(image.data(), w, h, {kCenter}, nullptr);
        REQUIRE(stars.size() == 1);
        CHECK(stars[0].position.x == Approx(32.5).margin(1e-4));
        CHECK(stars[0].position.y == Approx(32.5).margin(1e-4));
    }
}

TEST_CASE("Windowed CoG refuses radii just outside 0..kMaxRadius", "[centroid]") {
    const int w = 64, h = 64;
    Image image(w * h, 10);
    Fill(&image, w, 31, 33, 31, 33, 200);
    const std::vector<Vec2> predictions = {kCenter};
    const int tooLarge = WindowedCenterOfGravity::kMaxRadius + 1;
    CHECK(WindowedCenterOfGravity(tooLarge, 4, 3).Go(image.data(), w, h, predictions, nullptr).empty());
    CHECK(WindowedCenterOfGravity(4, tooLarge, 3).Go(image.data(), w, h, predictions, nullptr).empty());
    // negative star radius with a gate large enough that the star would otherwise fit in the window
    CHECK(WindowedCenterOfGravity(8, -1, 3).Go(image.data(), w, h, predictions, nullptr).empty());
    CHECK(WindowedCenterOfGravity(8, 0, 3).Go(image.data(), w, h, predictions, nullptr).size() == 1);
    // a window of one pixel, placed just off the image
    CHECK(WindowedCenterOfGravity(0, 0, 3).Go(image.data(), w, h, {{DECIMAL(5.0), (decimal)h}}, nullptr).empty());
    CHECK(WindowedCenterOfGravity(0, 0, 3).Go(image.data(), w, h, {{(decimal)w, DECIMAL(5.0)}}, nullptr).empty());
    CHECK(WindowedCenterOfGravity(0, 0, 3).Go(image.data(), w, h, {kCenter}, nullptr).empty());
}

TEST_CASE("Windowed CoG skips a prediction off each side of the image even with a star in reach", "[centroid]") {
    const int w = 64, h = 64;
    WindowedCenterOfGravity centroider(4, 4, 3);
    const decimal nan = std::numeric_limits<decimal>::quiet_NaN();
    const decimal inf = std::numeric_limits<decimal>::infinity();
    struct Case { int starX, starY; Vec2 inside, outside; };
    const Case cases[] = {
        {2, 32, {DECIMAL(0.5), DECIMAL(32.5)}, {DECIMAL(-0.5), DECIMAL(32.5)}},
        {61, 32, {DECIMAL(63.5), DECIMAL(32.5)}, {DECIMAL(64.5), DECIMAL(32.5)}},
        {32, 2, {DECIMAL(32.5), DECIMAL(0.5)}, {DECIMAL(32.5), DECIMAL(-0.5)}},
        {32, 61, {DECIMAL(32.5), DECIMAL(63.5)}, {DECIMAL(32.5), DECIMAL(64.5)}},
    };
    for (const Case &c : cases) {
        INFO("star at " << c.starX << ", " << c.starY);
        Image image(w * h, 10);
        Fill(&image, w, c.starX - 1, c.starX + 1, c.starY - 1, c.starY + 1, 200);
        CHECK(centroider.Go(image.data(), w, h, {c.inside}, nullptr).size() == 1);
        CHECK(centroider.Go(image.data(), w, h, {c.outside}, nullptr).empty());
    }
    Image image(w * h, 10);
    Fill(&image, w, 31, 33, 31, 33, 200);
    std::vector<int> matched;
    Stars stars = centroider.Go(image.data(), w, h,
                                {{nan, DECIMAL(32.5)}, {DECIMAL(32.5), nan}, {inf, DECIMAL(32.5)}, {DECIMAL(32.5), -inf}, kCenter},
                                &matched);
    CHECK(matched == std::vector<int>{4});
}

TEST_CASE("Windowed CoG rejects a star reaching each side of the window or image, and keeps one that stops short", "[centroid]") {
    const int w = 64, h = 64;
    WindowedCenterOfGravity centroider(4, 4, 3);
    auto found = [&](int x0, int x1, int y0, int y1, const Vec2 &prediction) {
        Image image(w * h, 10);
        Fill(&image, w, x0, x1, y0, y1, 200);
        return centroider.Go(image.data(), w, h, {prediction}, nullptr).size();
    };
    // window 24..40: bars that end one pixel inside it, then on it
    CHECK(found(31, 39, 32, 32, kCenter) == 1);
    CHECK(found(31, 40, 32, 32, kCenter) == 0);
    CHECK(found(25, 33, 32, 32, kCenter) == 1);
    CHECK(found(24, 33, 32, 32, kCenter) == 0);
    CHECK(found(32, 32, 31, 39, kCenter) == 1);
    CHECK(found(32, 32, 31, 40, kCenter) == 0);
    CHECK(found(32, 32, 25, 33, kCenter) == 1);
    CHECK(found(32, 32, 24, 33, kCenter) == 0);
    // image edges
    CHECK(found(1, 3, 31, 33, {DECIMAL(2.5), DECIMAL(32.5)}) == 1);
    CHECK(found(0, 2, 31, 33, {DECIMAL(2.5), DECIMAL(32.5)}) == 0);
    CHECK(found(60, 62, 31, 33, {DECIMAL(61.5), DECIMAL(32.5)}) == 1);
    CHECK(found(61, 63, 31, 33, {DECIMAL(61.5), DECIMAL(32.5)}) == 0);
    CHECK(found(31, 33, 1, 3, {DECIMAL(32.5), DECIMAL(2.5)}) == 1);
    CHECK(found(31, 33, 0, 2, {DECIMAL(32.5), DECIMAL(2.5)}) == 0);
    CHECK(found(31, 33, 60, 62, {DECIMAL(32.5), DECIMAL(61.5)}) == 1);
    CHECK(found(31, 33, 61, 63, {DECIMAL(32.5), DECIMAL(61.5)}) == 0);
}

TEST_CASE("Windowed CoG seeds anywhere in the gate circle and nowhere else", "[centroid]") {
    const int w = 64, h = 64;
    WindowedCenterOfGravity centroider(4, 4, 3);
    auto foundPixel = [&](int x, int y, const Vec2 &prediction) {
        Image image(w * h, 10);
        image[y * w + x] = 200;
        return centroider.Go(image.data(), w, h, {prediction}, nullptr).size();
    };
    // a lone pixel 3.9 px from the prediction, in the outermost column or row of the gate on each side
    CHECK(foundPixel(28, 32, {DECIMAL(32.4), DECIMAL(32.5)}) == 1);
    CHECK(foundPixel(36, 32, {DECIMAL(32.6), DECIMAL(32.5)}) == 1);
    CHECK(foundPixel(32, 28, {DECIMAL(32.5), DECIMAL(32.4)}) == 1);
    CHECK(foundPixel(32, 36, {DECIMAL(32.5), DECIMAL(32.6)}) == 1);
    // 4.1 px away
    CHECK(foundPixel(28, 32, {DECIMAL(32.6), DECIMAL(32.5)}) == 0);
    CHECK(foundPixel(36, 32, {DECIMAL(32.4), DECIMAL(32.5)}) == 0);
    // exactly gateRadius away counts as inside
    CHECK(foundPixel(28, 32, kCenter) == 1);
    CHECK(foundPixel(32, 36, kCenter) == 1);

    // a brighter pixel in the corner of the gate's bounding square, 4.24 px away, does not steal the seed
    Image image(w * h, 10);
    Fill(&image, w, 31, 33, 31, 33, 100);
    image[35 * w + 35] = 250;
    Stars stars = centroider.Go(image.data(), w, h, {kCenter}, nullptr);
    REQUIRE(stars.size() == 1);
    CHECK(stars[0].magnitude == 9);
    CHECK(stars[0].position.x == Approx(32.5).margin(1e-4));
}

TEST_CASE("Windowed CoG calls two stars ambiguous exactly when their bounding boxes overlap", "[centroid]") {
    const int w = 64, h = 64;
    // gate 2, so that neither prediction can seed on the other star
    WindowedCenterOfGravity centroider(2, 6, 3);
    for (int transposed = 0; transposed < 2; transposed++) {
        for (int gap = 0; gap < 2; gap++) {
            INFO("transposed " << transposed << ", gap " << gap);
            Image image(w * h, 10);
            auto set = [&](int x, int y) { image[transposed ? x * w + y : y * w + x] = 200; };
            // star A: bounding box 30..32 by 30..32
            for (int y = 30; y <= 32; y++) {
                set(30, y);
                set(31, y);
            }
            set(32, 30);
            // star B, not touching A: bounding box 32..34 by 31..36, or one further along with gap
            for (int y = 31; y <= 36; y++) {
                set(33 + gap, y);
                set(34 + gap, y);
            }
            set(32 + gap, 35);
            set(32 + gap, 36);
            Vec2 predictionA = {DECIMAL(30.5), DECIMAL(31.5)}, predictionB = {DECIMAL(34.5) + gap, DECIMAL(34.5)};
            if (transposed) {
                std::swap(predictionA.x, predictionA.y);
                std::swap(predictionB.x, predictionB.y);
            }
            // each is found alone
            CHECK(centroider.Go(image.data(), w, h, {predictionA}, nullptr).size() == 1);
            CHECK(centroider.Go(image.data(), w, h, {predictionB}, nullptr).size() == 1);
            std::vector<int> matched;
            centroider.Go(image.data(), w, h, {predictionA, predictionB}, &matched);
            if (gap == 0) {
                CHECK(matched.empty());  // boxes share column (or row) 32
            } else {
                CHECK(matched == std::vector<int>{0, 1});
            }
        }
    }

    // stars in the same columns but far apart in rows are not ambiguous
    Image image(w * h, 10);
    Fill(&image, w, 31, 33, 11, 13, 200);
    Fill(&image, w, 31, 33, 41, 43, 200);
    std::vector<int> matched;
    centroider.Go(image.data(), w, h, {{DECIMAL(32.5), DECIMAL(12.5)}, {DECIMAL(32.5), DECIMAL(42.5)}}, &matched);
    CHECK(matched == std::vector<int>{0, 1});
}

TEST_CASE("Windowed CoG outlier rule uses the median offset, both axes, and needs three matches", "[centroid]") {
    const int w = 256, h = 64;
    // star i is a 3x3 block centered at (20.5 + 30 i, 32.5); its prediction is that minus the offset
    auto run = [&](const std::vector<Vec2> &offsets, int gateRadius, decimal outlierRadius) {
        Image image(w * h, 10);
        std::vector<Vec2> predictions;
        for (int i = 0; i < (int)offsets.size(); i++) {
            Fill(&image, w, 19 + 30 * i, 21 + 30 * i, 31, 33, 200);
            predictions.push_back({DECIMAL(20.5) + 30 * i - offsets[i].x, DECIMAL(32.5) - offsets[i].y});
        }
        std::vector<int> matched;
        Stars stars = WindowedCenterOfGravity(gateRadius, 4, outlierRadius).Go(image.data(), w, h, predictions, &matched);
        CHECK(stars.size() == matched.size());
        return matched;
    };
    // two matches 5 px apart: nothing to compare against, both kept
    CHECK(run({{DECIMAL(2.5), 0}, {DECIMAL(-2.5), 0}}, 5, 3) == std::vector<int>{0, 1});
    // three: the odd one out goes
    CHECK(run({{0, 0}, {DECIMAL(4.0), 0}, {0, 0}}, 5, 3) == std::vector<int>{0, 2});
    // one far outlier does not drag the reference: the mean offset is 1.75 px from the three good ones
    CHECK(run({{0, 0}, {0, 0}, {0, 0}, {DECIMAL(7.0), 0}}, 8, DECIMAL(1.5)) == std::vector<int>{0, 1, 2});
    // offsets 1.5 px either side of the median are within 2 px of it; an outlier in y only is dropped
    CHECK(run({{DECIMAL(-1.5), 0}, {DECIMAL(1.5), 0}, {0, 0}, {DECIMAL(1.5), 0}, {DECIMAL(-1.5), 0}, {0, DECIMAL(5.0)}}, 6, 2)
          == std::vector<int>{0, 1, 2, 3, 4});
    // the reference is the (upper) median, not whichever offset happens to come first
    CHECK(run({{DECIMAL(-1.5), 0}, {DECIMAL(-1.5), 0}, {DECIMAL(-1.5), 0}, {DECIMAL(1.5), 0}, {DECIMAL(1.5), 0},
               {DECIMAL(1.5), 0}, {DECIMAL(1.5), 0}}, 5, 2) == std::vector<int>{3, 4, 5, 6});
    // exactly outlierRadius from the median is kept, just beyond is not
    CHECK(run({{0, 0}, {0, 0}, {0, 0}, {DECIMAL(3.0), 0}}, 5, 3) == std::vector<int>{0, 1, 2, 3});
    CHECK(run({{0, 0}, {0, 0}, {0, 0}, {DECIMAL(3.5), 0}}, 5, 3) == std::vector<int>{0, 1, 2});
}

TEST_CASE("Windowed CoG on a noisy 8-bit frame with several stars, a common shift, a hot pixel and an edge star", "[centroid]") {
    const int w = 200, h = 150;
    const double background = 20, noiseSigma = 3;
    struct TrueStar { double x, y, peak, sigma; };
    const TrueStar trueStars[] = {
        {31.37, 28.81, 230, 1.0}, {77.92, 40.15, 120, 1.2}, {140.08, 22.56, 60, 1.0}, {171.63, 61.24, 180, 0.8},
        {25.45, 99.70, 90, 1.1}, {88.19, 121.33, 45, 1.0}, {133.71, 88.02, 250, 1.3}, {180.50, 128.50, 150, 1.0},
        {2.20, 60.40, 200, 1.0},  // clipped by the left edge of the image
    };
    const int starCount = sizeof(trueStars) / sizeof(trueStars[0]);
    std::vector<double> ideal(w * h, background);
    for (const TrueStar &star : trueStars) AddStar(&ideal, w, h, star.x, star.y, star.peak, star.sigma);

    std::mt19937 generator(42);
    Image image(w * h);
    for (int i = 0; i < w * h; i++) {
        const double u1 = (static_cast<double>(generator()) + 1.0) / 4294967296.0;
        const double u2 = (static_cast<double>(generator()) + 1.0) / 4294967296.0;
        const double noise = noiseSigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        image[i] = (unsigned char)std::min(255.0, std::max(0.0, std::round(ideal[i] + noise)));
    }
    image[75 * w + 110] = 255;  // hot pixel, far from every prediction

    // the whole frame moved by (1.3, -0.8) px since the predictions were made; they are also a little off each
    std::vector<Vec2> predictions;
    for (int i = 0; i < starCount; i++) {
        predictions.push_back({(decimal)(trueStars[i].x - 1.3 + 0.2 * std::sin(i * 1.7)),
                               (decimal)(trueStars[i].y + 0.8 + 0.2 * std::cos(i * 2.3))});
    }
    predictions.push_back({DECIMAL(60.0), DECIMAL(75.0)});  // nothing there

    WindowedCenterOfGravity centroider(4, 4, 2);
    std::vector<int> matched;
    Stars stars = centroider.Go(image.data(), w, h, predictions, &matched);
    REQUIRE(matched == std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7});
    REQUIRE(stars.size() == matched.size());
    for (int i = 0; i < 8; i++) {
        INFO("star " << i);
        const double dx = static_cast<double>(stars[i].position.x) - trueStars[i].x;
        const double dy = static_cast<double>(stars[i].position.y) - trueStars[i].y;
        const double error = std::sqrt(dx * dx + dy * dy);
        // The faintest star's peak is 15 noise sigmas. Few of its pixels reach the threshold, and
        // the centroid of those few is less accurate.
        CHECK(error < (trueStars[i].peak < 50 ? 0.35 : 0.15));
        CHECK(stars[i].magnitude >= 4);
    }
}

TEST_CASE("Windowed CoG finds nothing in an image too large to index", "[centroid]") {
    Image image(4, 10);
    const int huge = std::numeric_limits<int>::max() / 2;
    WindowedCenterOfGravity centroider(4, 4, 3);
    CHECK(centroider.Go(image.data(), huge, 3, {{DECIMAL(1.5), DECIMAL(1.5)}}, nullptr).empty());
    CHECK(centroider.Go(image.data(), 3, huge, {{DECIMAL(1.5), DECIMAL(1.5)}}, nullptr).empty());

    // One row or column, few enough pixels to index, but too long for window coordinates to be
    // safe near its far end. The image is never read.
    const int longest = std::numeric_limits<int>::max();
    const int tooLong = WindowedCenterOfGravity::kMaxImageSide + 1;
    for (int side : {longest, tooLong}) {
        INFO("side " << side);
        const decimal farEnd = static_cast<decimal>(side) - DECIMAL(2.5);
        CHECK(centroider.Go(image.data(), side, 1, {{farEnd, DECIMAL(0.5)}, {DECIMAL(1.5), DECIMAL(0.5)}}, nullptr).empty());
        CHECK(centroider.Go(image.data(), 1, side, {{DECIMAL(0.5), farEnd}, {DECIMAL(0.5), DECIMAL(1.5)}}, nullptr).empty());
        CHECK(centroider.Go(image.data(), side, 2, {{farEnd, DECIMAL(0.5)}, {DECIMAL(1.5), DECIMAL(0.5)}}, nullptr).empty());
    }
}

TEST_CASE("Windowed CoG rarely takes noise for a star", "[centroid]") {
    // Pure Gaussian noise at several levels, no stars. A threshold of 5 true standard deviations
    // would give about 0.01 false stars here; before the MAD allowance this found well over a hundred.
    const int w = 320, h = 240;
    std::mt19937 generator(2026);
    WindowedCenterOfGravity centroider(4, 4, std::numeric_limits<decimal>::infinity());
    int falseStars = 0, windows = 0;
    for (double noiseSigma : {1.0, 2.0, 3.0, 5.0}) {
        Image image(w * h);
        for (int i = 0; i < w * h; i++) {
            const double u1 = (static_cast<double>(generator()) + 1.0) / 4294967296.0;
            const double u2 = (static_cast<double>(generator()) + 1.0) / 4294967296.0;
            const double noise = noiseSigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            image[i] = (unsigned char)std::min(255.0, std::max(0.0, std::round(20.3 + noise)));
        }
        // one prediction at a time, so that the windows are judged independently
        for (int y = 12; y < h - 12; y += 9) {
            for (int x = 12; x < w - 12; x += 9) {
                falseStars += centroider.Go(image.data(), w, h, {{x + DECIMAL(0.5), y + DECIMAL(0.5)}}, nullptr).size();
                windows++;
            }
        }
    }
    CAPTURE(falseStars, windows);
    CHECK(windows > 3000);
    CHECK(falseStars <= 6);
}
