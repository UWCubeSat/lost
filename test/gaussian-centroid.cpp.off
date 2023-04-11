#include <catch.hpp>
#include <iostream>

#include "centroiders.hpp"

using namespace lost;  // NOLINT

// TEST_CASE("Convert coordinates: pixel -> spatial -> pixel", "[geometry]") {
//     Camera camera(100, 512, 1024);

//     float expectedX = GENERATE(142, 90, 512, 255);
//     float expectedY = GENERATE(18, 512, 0, 800);

//     Vec3 spatial = camera.CameraToSpatial({expectedX, expectedY});
//     Vec2 actualPixels = camera.SpatialToCamera(spatial);

//     CHECK(actualPixels.x == Approx(expectedX).margin(0.000001));
//     CHECK(actualPixels.y == Approx(expectedY).margin(0.000001));
// }

TEST_CASE("Get pixel value from image", "[img]"){
  const int w = 5;
  const int h = 5;
  unsigned char image[w*h] = {0, 10, 30, 40, 10,
                              10, 10, 50, 255, 255,
                              0, 50, 255, 255, 255,
                              40, 255, 255, 255, 255,
                              35, 255, 255, 255, 255};

  CHECK(Get(3, 0, image, w) == 40);
  CHECK(Get(0, 4, image, w) == 35);
  CHECK(Get(4, 4, image, w) == 255);
}

TEST_CASE("Calculate X-Marginal", "[lsgf1]"){
  const int w = 5;
  const int h = 5;
  unsigned char image[w*h] = {0, 10, 30, 40, 10,
                              10, 10, 50, 255, 255,
                              0, 50, 255, 255, 255,
                              40, 255, 255, 255, 255,
                              35, 255, 255, 255, 255};

  CHECK(XMarginal(2, 1, 1, image, w) == 30 + 50 + 255);
  CHECK(XMarginal(2, 2, 2, image, w) == 30 + 50 + 255 * 3);
  CHECK(XMarginal(3, 3, 1, image, w) == 255 * 3);
}

TEST_CASE("Calculate Y-Marginal", "[lsgf1]"){
  const int w = 5;
  const int h = 5;
  unsigned char image[w*h] = {0, 10, 30, 40, 10,
                              10, 10, 50, 255, 255,
                              0, 50, 255, 255, 255,
                              40, 255, 255, 255, 255,
                              35, 255, 255, 255, 255};

  CHECK(YMarginal(2, 1, 2, image, w) == 10 + 10 + 50 + 255 + 255);
  CHECK(YMarginal(1, 4, 1, image, w) == 35 + 255 + 255);
}

TEST_CASE("Calculate initial guess for Gaussian Fit params", "[lsgf]"){
  const int w = 5;
  const int h = 5;
  /////////////////////////// 0  1   2   3    4
  unsigned char image[w*h] = {0, 10, 30, 40, 10,          // 0
                              10, 10, 50, 255, 255,       // 1
                              0, 50, 255, 255, 255,       // 2
                              40, 255, 255, 255, 255,     // 3
                              35, 255, 255, 255, 255};    // 4

  float a;
  float xb, yb;
  double sigma;
  InitialGuess(2, 1, 1, image, w, &a, &xb, &yb, &sigma);
  CHECK(a == 255);
  std::cout << "xb, yb = " <<  xb << ", " << yb << std::endl;
}