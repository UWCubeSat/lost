#include <catch.hpp>
#include <random>
#include <cmath>

#include "io.hpp"
#include "star-utils.hpp"

#include "fixtures.hpp"

using namespace lost; // NOLINT

Stars Perturb(Stars unperturbed, float perturbation, std::default_random_engine *rng) {
    Stars result = unperturbed;
    std::uniform_real_distribution<float> dist(-perturbation, perturbation);
    for (Star &star : result) {
        star.position.x += dist(*rng);
        star.position.y += dist(*rng);
    }
    return result;
}

TEST_CASE("Star IDs compare: All identical", "[comparators] [fast]") {
    StarIdComparison result = StarIdsCompare(elevenStarIds, elevenStarIds,
                                             integralCatalog, integralCatalog,
                                             0.001,
                                             elevenStars, elevenStars);
    // despite the name "elevenStarIds", that actually means the indices go up to 11, but they star
    // at zero, so there are 12 total:
    CHECK(result.numCorrect == 12);
    CHECK(result.numIncorrect == 0);
    CHECK(result.numTotal == 12);
}

TEST_CASE("Star IDs compare: Some perturbations", "[comparators] [fast]") {
    std::default_random_engine rng(GENERATE(take(3, random(1,10000))));
    Stars perturbedEleven = Perturb(elevenStars, 0.4, &rng);
    StarIdComparison result1 = StarIdsCompare(elevenStarIds, elevenStarIds,
                                              integralCatalog, integralCatalog,
                                              sqrt(0.41*0.41 + 0.41*0.41),
                                              elevenStars, perturbedEleven);
    CHECK(result1.numCorrect == 12);
    CHECK(result1.numIncorrect == 0);
    CHECK(result1.numTotal == 12);

    // make sure it still works with higher tolerance
    // this effectively checks that it identifies things okay when there are multiple input stars near each expected star.
    StarIdComparison result2 = StarIdsCompare(elevenStarIds, elevenStarIds,
                                              integralCatalog, integralCatalog,
                                              5,
                                              elevenStars, perturbedEleven);
    CHECK(result2.numCorrect == 12);
    CHECK(result2.numIncorrect == 0);
    CHECK(result2.numTotal == 12);

    // Finally, make sure it fails with lower tolerance
    StarIdComparison result3 = StarIdsCompare(elevenStarIds, elevenStarIds,
                                              integralCatalog, integralCatalog,
                                              sqrt(0.3*0.3 + 0.3*0.3),
                                              elevenStars, perturbedEleven);
    CHECK(result3.numCorrect < 12);
    CHECK(result3.numCorrect > 0);

    CHECK(result3.numIncorrect > 0);
    CHECK(result3.numIncorrect < result3.numTotal);

    CHECK(result3.numTotal > 0);
    CHECK(result3.numTotal < 12);
}

TEST_CASE("Star IDs compare: Permute a few things", "[comparators] [fast]") {
    // swap the first two star-ids
    std::vector<StarIdentifier> permutedStarIds = elevenStarIds;
    std::swap(permutedStarIds[0], permutedStarIds[1]);
    StarIdComparison result1 = StarIdsCompare(elevenStarIds, permutedStarIds,
                                              integralCatalog, integralCatalog,
                                              0.001,
                                              elevenStars, elevenStars);
    CHECK(result1.numCorrect == 12);
    CHECK(result1.numIncorrect == 0);
    CHECK(result1.numTotal == 12);

    // also swap the first two from catalog
    Catalog permutedCatalog = integralCatalog;
    std::swap(permutedCatalog[0], permutedCatalog[1]);
    permutedStarIds[0].catalogIndex = 0;
    permutedStarIds[1].catalogIndex = 1;
    StarIdComparison result2 = StarIdsCompare(permutedStarIds, elevenStarIds,
                                              permutedCatalog, integralCatalog,
                                              0.001,
                                              elevenStars, elevenStars);
    CHECK(result2.numCorrect == 12);
    CHECK(result2.numIncorrect == 0);
    CHECK(result2.numTotal == 12);
}

TEST_CASE("Star IDs compare: Fewer input stars.", "[comparators] [fast]") {
    // remove the first star from the input
    Stars fiveStars = elevenStars;
    fiveStars.erase(fiveStars.begin() + 5, fiveStars.end());
    StarIdentifiers fiveStarIds = elevenStarIds;
    fiveStarIds.erase(fiveStarIds.begin() + 5, fiveStarIds.end());
    StarIdComparison result = StarIdsCompare(elevenStarIds, fiveStarIds,
                                             integralCatalog, integralCatalog,
                                             0.001,
                                             elevenStars, fiveStars);
    CHECK(result.numCorrect == 5);
    CHECK(result.numIncorrect == 0);
    CHECK(result.numTotal == 5);

    // check backwards
    result = StarIdsCompare(fiveStarIds, elevenStarIds,
                            integralCatalog, integralCatalog,
                            0.001,
                            fiveStars, elevenStars);
    CHECK(result.numCorrect == 5);
    CHECK(result.numIncorrect == 12-5);
    // While there are 12 expected stars, there are only 5 input stars, so even an ideal algo can't
    // identify more than that.
    CHECK(result.numTotal == 5);
}
