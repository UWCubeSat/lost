#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

#include "databases.hpp"
#include "centroiders.hpp"
#include "io.hpp"

namespace lost {

static void DatabaseBuild() {
    Catalog narrowedCatalog = PromptNarrowedCatalog(CatalogRead());
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseBuilder builder;
    unsigned char *catalogBuffer = builder.AddSubDatabase(kCatalogMagicValue,
                                                          // TODO: allow magnitude and weird
                                                          SerializeLengthCatalog(narrowedCatalog, false, true));
    SerializeCatalog(narrowedCatalog, false, true, catalogBuffer);

    PromptDatabases(builder, narrowedCatalog);

    std::cerr << "Generated database with " << builder.BufferLength() << " bytes" << std::endl;
    PromptedOutputStream pos;
    pos.Stream().write((char *)builder.Buffer(), builder.BufferLength());
}

static void PipelineRun() {
    PipelineInputList input = PromptPipelineInput();
    Pipeline pipeline = PromptPipeline();
    std::vector<PipelineOutput> outputs = pipeline.Go(input);
    PromptPipelineComparison(input, outputs);
}

static void PipelineBenchmark() {
    PipelineInputList input = PromptPipelineInput();
    Pipeline pipeline = PromptPipeline();
    int iterations = Prompt<int>("Times to run the pipeline");
    std::cerr << "Benchmarking..." << std::endl;

    // TODO: we can do better than this :| maybe include mean time, 99% time, or allow a vector of
    // input and determine which one took the longest
    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; i++) {
        pipeline.Go(input);
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    auto totalTime = std::chrono::duration<double, std::milli>(endTime - startTime);
    std::cout << "total_ms " << totalTime.count() << std::endl;
}

static void EstimateCamera() {
    std::cerr << "Enter estimated camera details when prompted." << std::endl;
    PipelineInputList inputs = PromptPngPipelineInput();
    float baseFocalLength = inputs[0]->InputCamera()->FocalLength();
    float deviationIncrement = Prompt<float>("Focal length increment (base: " + std::to_string(baseFocalLength) + ")");
    float deviationMax = Prompt<float>("Maximum focal length deviation to attempt");
    Pipeline pipeline = PromptPipeline();

    while (inputs[0]->InputCamera()->FocalLength() - baseFocalLength <= deviationMax) {
        std::cerr << "Attempt focal length " << inputs[0]->InputCamera()->FocalLength() << std::endl;
        std::vector<PipelineOutput> outputs = pipeline.Go(inputs);
        if (outputs[0].nice) {
            std::cout << "camera_identified true" << std::endl << *inputs[0]->InputCamera();
            return;
        }

        Camera camera(*inputs[0]->InputCamera());
        if (camera.FocalLength() - baseFocalLength > 0) {
            // yes i know this expression can be simplified shut up
            camera.SetFocalLength(camera.FocalLength() - 2*(camera.FocalLength() - baseFocalLength));
        } else {
            camera.SetFocalLength(camera.FocalLength() + 2*(baseFocalLength - camera.FocalLength()) + deviationIncrement);
        }
        ((PngPipelineInput *)(inputs[0].get()))->SetCamera(camera);
    }
    std::cout << "camera_identified false" << std::endl;
}

}

int main(int argc, char **argv) {
    lost::RegisterCliArgs(argc, argv);
    std::cerr << "LOST: Open-source Star Tracker" << std::endl;
    lost::InteractiveChoice<void (*)()> mainChoices;
    mainChoices.Register("pipeline", "Run a pipeline", &lost::PipelineRun);
    mainChoices.Register("build_database", "Build database from catalog", &lost::DatabaseBuild);
    mainChoices.Register("benchmark", "Benchmark a pipeline", &lost::PipelineBenchmark);
    mainChoices.Register("estimate_camera", "Estimate camera parameters", &lost::EstimateCamera);
    (*mainChoices.Prompt("Choose action"))();
}
