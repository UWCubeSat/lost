#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <cairo/cairo.h>

#include <string>
#include <iostream>
#include <chrono>

#include "catalog-generic.hpp"
#include "database-builders.hpp"
#include "centroiders.hpp"
#include "io.hpp"

namespace lost {

static void CatalogBuild() {
    // std::string tsvPath;
    // Catalog_t *catalog;
    // int i_builder_choice;
    // db_builder_algorithm_t x_algorithm;
    // void *pv_algorithm_config, *pv_db;

    // puts("Location of tsv file (./bright-star-database.tsv): ");
    // std::getline(std::cin, tsvPath);

    // if (tsvPath.empty()) {
    //     tsvPath = std::string("./bright-star-database.tsv");
    // }

    // catalog = BsdParse(tsvPath);
    // printf("Read %ld stars from catalog.\n", catalog->numStars);

    // for (int i = 0; x_db_builder_algorithm(i).s_name; i++) {
    //     printf("%d) %s\n", i, x_db_builder_algorithm(i).s_name);
    // }
    // printf("Choose database builder: ");
    // scanf("%d", &i_builder_choice);
    // getchar();
    // x_algorithm = x_db_builder_algorithm(i_builder_choice);
    // if (x_algorithm.s_name == NULL) {
    //     puts("Not found...");
    //     return;
    // }
    // pv_algorithm_config = (*x_algorithm.pf_config)();
    // pv_db = (*(x_algorithm.pf_algorithm))(catalog, pv_algorithm_config);
    // (*x_algorithm.pf_stats)(pv_db);
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
    auto totalTime = std::chrono::duration<double, std::milli>(startTime - endTime);
    std::cerr << "Complete in " << totalTime.count() << " milliseconds" << std::endl;
}

}

int main(int argc, char **argv) {
    lost::RegisterCliArgs(argc, argv);
    std::cout << "LOST: Open-source Star Tracker" << std::endl;
    lost::InteractiveChoice<void (*)()> mainChoices;
    mainChoices.Register("pipeline", "Run a pipeline", &lost::PipelineRun);
    mainChoices.Register("benchmark", "Benchmark a pipeline", &lost::PipelineBenchmark);
    mainChoices.Register("build_database", "Build database from catalog", &lost::CatalogBuild);
    (*mainChoices.Prompt("Choose action"))();
}
