#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

#include "database-builders.hpp"
#include "centroiders.hpp"
#include "io.hpp"

namespace lost {

static void DatabaseBuild() {
    DbBuilder dbBuilder = PromptDbBuilder();
    long length;
    void *db = dbBuilder(CatalogRead(), &length);
    std::cerr << "Generated database with " << length << " bytes" << std::endl;
    PromptedOutputStream pos;
    pos.Stream().write((char *)db, length);
    free(db);
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

}

int main(int argc, char **argv) {
    lost::RegisterCliArgs(argc, argv);
    std::cerr << "LOST: Open-source Star Tracker" << std::endl;
    lost::InteractiveChoice<void (*)()> mainChoices;
    mainChoices.Register("pipeline", "Run a pipeline", &lost::PipelineRun);
    mainChoices.Register("benchmark", "Benchmark a pipeline", &lost::PipelineBenchmark);
    mainChoices.Register("build_database", "Build database from catalog", &lost::DatabaseBuild);
    (*mainChoices.Prompt("Choose action"))();
}
