#include <catch.hpp>

#include "io.hpp"

using namespace lost;

TEST_CASE("Center of Gravity: Basic generated input comparison", "[centroid]") {
    PipelineInputList result;
    cairo_surface_t *cairoSurface = NULL;
    std::string pngPath;

    while (cairoSurface == NULL || cairo_surface_status(cairoSurface) != CAIRO_STATUS_SUCCESS) {
        // Hard code PNG file path
        cairoSurface = cairo_image_surface_create_from_png(Prompt<std::string>("PNG Path").c_str());
    }
    result.push_back(std::unique_ptr<PipelineInput>(new PngPipelineInput(cairoSurface)));

    Pipeline result;
    // Following code is making VSCode upset
    // result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(
    //     new CenterOfGravityAlgorithm()());
 
    // std::vector<PipelineOutput> outputs = pipeline.Go(input);

    // PromptPipelineComparison(result, outputs);
}
