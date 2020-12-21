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

// prompts the user for a path, reads it, creates rgb cairo surface, 
static cairo_surface_t *PngRead() {
    std::string pngPath;
    cairo_surface_t *cairoSurface = NULL;

#ifndef CAIRO_HAS_PNG_FUNCTIONS
    puts("Your version of Cairo was compiled without PNG support. Bailing out.");
    exit(0);
#endif
    
    while (cairoSurface == NULL ||
           cairo_surface_status(cairoSurface) != CAIRO_STATUS_SUCCESS) {

        pngPath = Prompt<std::string>("Location of PNG file");

        cairoSurface = cairo_image_surface_create_from_png(pngPath.c_str());

        printf("Reading file: %s\n", cairo_status_to_string(cairo_surface_status(cairoSurface)));
    }

    return cairoSurface;
}

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

static void CentroidsFind() {
    std::string     outputPath;
    cairo_surface_t *cairoSurface;
    unsigned char   *image;

    cairoSurface = PngRead();
    image   = SurfaceToGrayscaleImage(cairoSurface);

    auto factory = makeCentroidAlgorithmChoice().Prompt(std::string("Choose centroiding algo"));
    CentroidAlgorithm *centroidAlgorithm = factory();

    std::vector<Star> stars = centroidAlgorithm->Go(
        image,
        cairo_image_surface_get_width(cairoSurface),
        cairo_image_surface_get_height(cairoSurface)
        );

    delete centroidAlgorithm;
    free(image);

    std::cout << stars.size() << " stars detected." << std::endl;
    outputPath = Prompt<std::string>("Plot output to PNG file");
    // TODO: show exact coordinates

    // plotting
    SurfacePlotCentroids(cairoSurface, stars, 1.0, 0.0, 0.0, 0.3);
    cairo_surface_write_to_png(cairoSurface, outputPath.c_str());
    cairo_surface_destroy(cairoSurface);
}

static void PipelineBenchmark() {
    PipelineInput input = PromptPipelineInput();
    Pipeline pipeline = PromptPipeline();
    int iterations = Prompt<int>("Times to run the pipeline");
    std::cerr << "Benchmarking..." << std::endl;

    // TODO: we can do better than this :| maybe include mean time, 99% time, or allow a vector of
    // input and determine which one took the longest
    auto startTime = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; i++) {
        PipelineRun()
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
