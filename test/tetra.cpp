#include <cmath>
#include <set>
#include <vector>

#include <catch.hpp>

#include "attitude-utils.hpp"
#include "star-utils.hpp"

using namespace lost;  // NOLINT

// -- CombinationIterator tests --

TEST_CASE("CombinationIterator produces correct count", "[tetra]") {
    // C(6, 3) = 20
    CombinationIterator it(3, 6);
    int count = 0;
    while (it.next())
        count++;
    REQUIRE(count == 20);
}

TEST_CASE("CombinationIterator C(n,k) edge cases", "[tetra]") {
    SECTION("k == n: exactly one combination") {
        CombinationIterator it(4, 4);
        int count = 0;
        while (it.next())
            count++;
        REQUIRE(count == 1);
    }

    SECTION("k > n: no combinations") {
        CombinationIterator it(5, 3);
        REQUIRE(it.next() == false);
    }

    SECTION("k == 1: n combinations") {
        CombinationIterator it(1, 5);
        int count = 0;
        while (it.next())
            count++;
        REQUIRE(count == 5);
    }
}

TEST_CASE("CombinationIterator produces co-lexicographic order", "[tetra]") {
    // All combinations of {0,1,2,3} should be exhausted before index 4 appears
    CombinationIterator it(3, 6);
    bool seenFour = false;
    std::set<std::vector<int>> combosBeforeFour;

    while (it.next()) {
        const std::vector<int>& idx = it.indices();
        bool hasFour = (idx[0] == 4 || idx[1] == 4 || idx[2] == 4);
        if (hasFour && !seenFour) {
            seenFour = true;
            // All C(4,3) = 4 combinations of {0,1,2,3} should already have been visited
            REQUIRE(combosBeforeFour.size() == 4);
        }
        if (!seenFour) {
            combosBeforeFour.insert(std::vector<int>(idx.begin(), idx.end()));
        }
    }
    REQUIRE(seenFour);
}

TEST_CASE("CombinationIterator produces no duplicates", "[tetra]") {
    CombinationIterator it(3, 7);
    std::set<std::vector<int>> seen;
    while (it.next()) {
        std::vector<int> combo(it.indices().begin(), it.indices().end());
        REQUIRE(seen.find(combo) == seen.end());
        seen.insert(combo);
    }
    // C(7,3) = 35
    REQUIRE(seen.size() == 35);
}

// -- TetraEdgeRatios tests --

TEST_CASE("TetraEdgeRatios produces correct number of ratios", "[tetra]") {
    // 4 stars -> C(4,2) = 6 edges -> 5 ratios (largest is the divisor)
    std::vector<Vec3> vecs = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 1, 0)};
    auto ratios = TetraEdgeRatios(4, [&](int i) -> const Vec3& { return vecs[i]; });
    REQUIRE(ratios.size() == 5);
}

TEST_CASE("TetraEdgeRatios are sorted ascending", "[tetra]") {
    std::vector<Vec3> vecs = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 1, 1)};
    auto ratios = TetraEdgeRatios(4, [&](int i) -> const Vec3& { return vecs[i]; });
    for (size_t i = 1; i < ratios.size(); i++) {
        REQUIRE(ratios[i] >= ratios[i - 1]);
    }
}

TEST_CASE("TetraEdgeRatios are all in [0, 1)", "[tetra]") {
    std::vector<Vec3> vecs = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 1, 0)};
    auto ratios = TetraEdgeRatios(4, [&](int i) -> const Vec3& { return vecs[i]; });
    for (decimal r : ratios) {
        REQUIRE(r >= DECIMAL(0.0));
        REQUIRE(r < DECIMAL(1.0));
    }
}

TEST_CASE("TetraEdgeRatios are scale-invariant", "[tetra]") {
    std::vector<Vec3> vecs = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 1, 0)};
    std::vector<Vec3> scaled = {Vec3(3, 0, 0), Vec3(0, 3, 0), Vec3(0, 0, 3), Vec3(3, 3, 0)};
    auto ratios1 = TetraEdgeRatios(4, [&](int i) -> const Vec3& { return vecs[i]; });
    auto ratios2 = TetraEdgeRatios(4, [&](int i) -> const Vec3& { return scaled[i]; });
    REQUIRE(ratios1.size() == ratios2.size());
    for (size_t i = 0; i < ratios1.size(); i++) {
        REQUIRE(ratios1[i] == Approx(ratios2[i]).margin(1e-6));
    }
}

// -- AllPairsWithinFov tests --

TEST_CASE("AllPairsWithinFov accepts stars within threshold", "[tetra]") {
    // Unit vectors along axes, all within 90 degrees of each other
    std::vector<Vec3> vecs = {
        Vec3(1, 0, 0).Normalize(), Vec3(1, 1, 0).Normalize(), Vec3(1, 0, 1).Normalize()};
    decimal cosMaxFov = std::cos(DegToRad(DECIMAL(91.0)));
    REQUIRE(AllPairsWithinFov(3, [&](int i) -> const Vec3& { return vecs[i]; }, cosMaxFov));
}

TEST_CASE("AllPairsWithinFov rejects stars outside threshold", "[tetra]") {
    // Two stars ~180 degrees apart
    std::vector<Vec3> vecs = {Vec3(1, 0, 0), Vec3(-1, 0, 0), Vec3(0, 1, 0)};
    decimal cosMaxFov = std::cos(DegToRad(DECIMAL(90.0)));
    REQUIRE_FALSE(AllPairsWithinFov(3, [&](int i) -> const Vec3& { return vecs[i]; }, cosMaxFov));
}

TEST_CASE("AllPairsWithinFov with tight threshold", "[tetra]") {
    // Stars ~45 degrees apart, threshold at 30 degrees
    std::vector<Vec3> vecs = {Vec3(1, 0, 0).Normalize(), Vec3(1, 1, 0).Normalize()};
    decimal cosMaxFov = std::cos(DegToRad(DECIMAL(30.0)));
    REQUIRE_FALSE(AllPairsWithinFov(2, [&](int i) -> const Vec3& { return vecs[i]; }, cosMaxFov));
}

// -- TetraFindPatternIndex tests --

TEST_CASE("TetraFindPatternIndex is deterministic", "[tetra]") {
    int key[] = {1, 2, 3, 4, 5};
    size_t idx1 = TetraFindPatternIndex(key, 5, 50, 1000);
    size_t idx2 = TetraFindPatternIndex(key, 5, 50, 1000);
    REQUIRE(idx1 == idx2);
}

TEST_CASE("TetraFindPatternIndex stays within table bounds", "[tetra]") {
    for (int a = 0; a < 10; a++) {
        int key[] = {a, a + 1, a + 2, a + 3, a + 4};
        size_t idx = TetraFindPatternIndex(key, 5, 50, 997);
        REQUIRE(idx < 997);
    }
}

TEST_CASE("TetraFindPatternIndex different keys produce different indices", "[tetra]") {
    int key1[] = {1, 2, 3, 4, 5};
    int key2[] = {5, 4, 3, 2, 1};
    size_t idx1 = TetraFindPatternIndex(key1, 5, 50, 10007);
    size_t idx2 = TetraFindPatternIndex(key2, 5, 50, 10007);
    // Not guaranteed by hash functions, but for these specific inputs should differ
    REQUIRE(idx1 != idx2);
}
