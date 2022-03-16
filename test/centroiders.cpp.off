#include <catch.hpp>

#include "io.hpp"
#include "centroiders.hpp"

using namespace lost;

TEST_CASE("Center of Gravity: Basic generated input comparison", "[centroid]") {
    int numImages = 1;
    int xResolution = 1024;
    int yResolution = 1024;
    float xFovDeg = 40;
    float xFocalLength = xResolution / (float)2.0 / tan(DegToRad(xFovDeg)/2);
    int referenceBrightness = 8000;
    float brightnessDeviation = 1;
    float noiseDeviation = 0;

    Quaternion attitude = SphericalToQuaternion(0, 0, 0);

    PipelineInputList resultList;

    for (int i = 0; i < numImages; i++) {
        GeneratedPipelineInput *curr = new GeneratedPipelineInput(
            CatalogRead(),
            attitude,
            Camera(xFocalLength, xResolution, yResolution),
            referenceBrightness, brightnessDeviation, noiseDeviation);

        resultList.push_back(std::unique_ptr<PipelineInput>(curr));
    }

    std::vector<PipelineOutput> outputs;
    for (const std::unique_ptr<PipelineInput> &input : resultList) {
        PipelineOutput result = {0};
        const Image *inputImage = (*input).InputImage();
        CenterOfGravityAlgorithm *centroidAlgorithm = new CenterOfGravityAlgorithm();
        result.stars = std::unique_ptr<Stars>(new std::vector<Star>(
            centroidAlgorithm->Go(inputImage->image, inputImage->width, inputImage->height)));
        outputs.push_back(std::move(result));
    }

    // PromptPipelineComparison(resultList, outputs);
    // assert centroider found some kind of stars, output shouldn't be empty
    // assert not found more centroids than present
    float centroidThreshold = 0.0f;

    std::vector<StarIdComparison> comparisons;
    for (int i = 0; i < (int)resultList.size(); i++) {
        // Need to figure out how to modularize StarIdsCompare()
        comparisons.push_back(
            StarIdsCompare(*resultList[i]->ExpectedStarIds(), *outputs[i].starIds,
                           CatalogRead(), CatalogRead(),
                           centroidThreshold, resultList[i]->ExpectedStars(), outputs[i].stars.get()));
    }

    float fractionIncorrectSum = 0;
    float fractionCorrectSum = 0;
    for (const StarIdComparison &comparison : comparisons) {
        fractionIncorrectSum += comparison.fractionIncorrect;
        fractionCorrectSum += comparison.fractionCorrect;
    }

    float fractionIncorrectMean = fractionIncorrectSum / comparisons.size();
    float fractionCorrectMean = fractionCorrectSum / comparisons.size();
}
