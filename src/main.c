#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <inttypes.h>
#include <cairo/cairo.h>

#include "catalog-generic.h"
#include "database-builders.h"
#include "centroiders.h"
#include "io.h"

// prompts the user for a path, reads it, creates rgb cairo surface, 
static cairo_surface_t *px_read_png() {
    char s_png_path[256];
    cairo_surface_t *px_surface = NULL;
    int width, height;

#ifndef CAIRO_HAS_PNG_FUNCTIONS
    puts("Your version of Cairo was compiled without PNG support. Bailing out.");
    exit(0);
#endif
    
    while (px_surface == NULL ||
           cairo_surface_status(px_surface) != CAIRO_STATUS_SUCCESS) {

        printf("Location of PNG file: ");
        scanf("%255s", s_png_path);
        getchar();

        px_surface = cairo_image_surface_create_from_png(s_png_path);

        printf("Reading file: %s\n", cairo_status_to_string(cairo_surface_status(px_surface)));
    }

    return px_surface;
}

static void v_build_catalog() {
    FILE            *px_tsv;
    char            s_default_tsv_path[] = "bright-star-database.tsv";
    char            s_tsv_path[256];
    catalog_t       *px_catalog;
    int             i_builder_choice;
    db_builder_algorithm_t x_algorithm;
    void            *pv_algorithm_config, *pv_db;

    printf("Location of tsv file (%s): ", s_default_tsv_path);
    fgets(s_tsv_path, 256 - 1, stdin);
    // get rid of trailing newline
    s_tsv_path[strlen(s_tsv_path)-1] = '\0';

    if (s_tsv_path[0] == '\0') {
        strcpy(s_tsv_path, s_default_tsv_path);
    }

    px_catalog = px_bsd_parse(s_tsv_path);
    printf("Read %ld stars from catalog.\n", px_catalog->l_stars_length);

    for (int i = 0; x_db_builder_algorithm(i).s_name; i++) {
        printf("%d) %s\n", i, x_db_builder_algorithm(i).s_name);
    }
    printf("Choose database builder: ");
    scanf("%d", &i_builder_choice);
    getchar();
    x_algorithm = x_db_builder_algorithm(i_builder_choice);
    if (x_algorithm.s_name == NULL) {
        puts("Not found...");
        return;
    }
    pv_algorithm_config = (*x_algorithm.pf_config)();
    pv_db = (*(x_algorithm.pf_algorithm))(px_catalog, pv_algorithm_config);
    (*x_algorithm.pf_stats)(pv_db);
}

static void v_centroids_find() {
    char            s_output_path[256];
    cairo_surface_t *px_surface;
    unsigned char   *pc_image;

    int                  i_centroid_choice, i_centroids_length;
    centroid_algorithm_t x_algorithm;
    star_t           *px_centroids;
    void                 *pv_algorithm_config;

    cairo_t   *px_cairo;

    px_surface = px_read_png();
    pc_image   = pc_surface_grayscale(px_surface);

    for (int i = 0; x_centroid_algorithm(i).s_name; i++) {
        printf("%d) %s\n", i, x_centroid_algorithm(i).s_name);
    }
    printf("Choose centroiding algorithm: ");
    scanf("%d", &i_centroid_choice);
    getchar();
    x_algorithm = x_centroid_algorithm(i_centroid_choice);
    if (x_algorithm.s_name == NULL) {
        puts("Not found...");
        return;
    }
    pv_algorithm_config = (*x_algorithm.pf_config)();

    px_centroids = (*x_algorithm.pf_algorithm)(
        pc_image,
        cairo_image_surface_get_width(px_surface),
        cairo_image_surface_get_height(px_surface),
        &i_centroids_length,
        pv_algorithm_config
        );

    free(pc_image);

    printf("%d stars detected.\n", i_centroids_length);
    // TODO: show exact coordinates
    printf("Plot output to PNG file: ");
    scanf("%s", s_output_path);
    getchar();

    // plotting
    v_surface_plot_centroids(px_surface, px_centroids, i_centroids_length, 1.0, 0.0, 0.0, 0.5);
    cairo_surface_write_to_png(px_surface, s_output_path);
    cairo_surface_destroy(px_surface);
}

int main() {
    puts("LOST: Open-source Star Tracker");

    while (1) {
        printf("Would you like to Build Catalog (b), Find Centroids (c), or Quit (q)? ");
        char choice;
        scanf(" %c", &choice);
        getchar();
        switch (choice) {
        case 'b':
            v_build_catalog();
            break;
        case 'c':
            v_centroids_find();
        default:
            puts("Bye!");
            exit(0);
        }
    }
    puts("");
}
