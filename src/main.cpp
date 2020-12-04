#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <cairo/cairo.h>

#include <string>
#include <iostream>

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

        printf("Location of PNG file: ");
        std::cin >> pngPath;

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

    auto factory = makeCentroidAlgorithmChoice().Prompt(std::cout, std::cin, std::string("Choose centroiding algo"));
    CentroidAlgorithm *centroidAlgorithm = factory(std::cout, std::cin);

    std::vector<Star> stars = centroidAlgorithm->Go(
        image,
        cairo_image_surface_get_width(cairoSurface),
        cairo_image_surface_get_height(cairoSurface)
        );

    delete centroidAlgorithm;
    free(image);

    std::cout << stars.size() << " stars detected." << std::endl;
    std::cout << "Plot output to PNG file: ";
    std::cin >> outputPath;
    // TODO: show exact coordinates

    // plotting
    SurfacePlotCentroids(cairoSurface, stars, 1.0, 0.0, 0.0, 0.5);
    cairo_surface_write_to_png(cairoSurface, outputPath.c_str());
    cairo_surface_destroy(cairoSurface);
}

}

int main() {
    puts("LOST: Open-source Star Tracker");

    while (1) {
        std::cout << "Would you like to Build Catalog (b), Find Centroids (c), or Quit (q)? ";
        char choice;
        scanf(" %c", &choice);
        getchar();
        switch (choice) {
        case 'b':
            lost::CatalogBuild();
            break;
        case 'c':
            lost::CentroidsFind();
        default:
            puts("Bye!");
            exit(0);
        }
    }
    puts("");
}
