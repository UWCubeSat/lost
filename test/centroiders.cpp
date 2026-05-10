#include <vector>

#include <catch.hpp>

#include "centroiders.hpp"

using namespace lost; // NOLINT

TEST_CASE("Center of gravity handles a large bright component", "[centroiders]") {
    const int width = 512;
    const int height = 512;
    std::vector<unsigned char> image(width * height, 255);

    CenterOfGravityAlgorithm algo;
    Stars stars = algo.Go(image.data(), width, height);

    // A fully bright image touches all edges, so this component is rejected as invalid.
    CHECK(stars.empty());
}
