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

static catalog_t *px_bsd_parse(char *s_path) {
    catalog_t      *px_result;
    catalog_star_t *px_star_current;
    FILE           *px_file;
    long           l_raj2000_h, l_raj2000_l, // high and low parts
                   l_dej2000_h, l_dej2000_l;
    int            i_magnitude_h, i_magnitude_l;
    char           c_weird;

    px_result = calloc(sizeof(catalog_t), 1);

    px_file = fopen(s_path, "r");
    if (px_file == NULL) {
        printf("Error opening file: %s\n", strerror(errno));
        return NULL;
    }

    while (EOF != fscanf(px_file, "%ld.%ld|%ld.%ld|%*d|%c|%d.%d",
                         &l_raj2000_h, &l_raj2000_l,
                         &l_dej2000_h, &l_dej2000_l,
                         &c_weird,
                         &i_magnitude_h, &i_magnitude_l)) {
           
        if (px_result->l_stars_length % 1000 == 0) {
            px_result->px_stars = realloc(px_result->px_stars,
                                          sizeof(catalog_star_t) *
                                          (px_result->l_stars_length/1000 + 1) * 1000);
        }

        px_star_current = &px_result->px_stars[px_result->l_stars_length];
        px_star_current->l_raj2000 = l_raj2000_h * 1000000 + l_raj2000_l;
        px_star_current->l_dej2000 = l_dej2000_h * 1000000 + l_dej2000_l;
        px_star_current->i_magnitude = i_magnitude_h * 100 + i_magnitude_l;
        px_star_current->i_weird = c_weird != ' ';
        px_result->l_stars_length++;
    }

    fclose(px_file);
    return px_result;
}

#define RED_MASK   0x00FF0000
#define GREEN_MASK 0x0000FF00
#define BLUE_MASK  0x000000FF

// convert an RGB cairo surface into one byte per pixel.
static unsigned char *pc_surface_grayscale(cairo_surface_t *px_surface) {
    int i_width, i_height;
    unsigned char *pc_result;
    uint32_t      *px_surface_data, x_pixel;

    if (cairo_image_surface_get_format(px_surface) != CAIRO_FORMAT_ARGB32 &&
        cairo_image_surface_get_format(px_surface) != CAIRO_FORMAT_RGB24) {
        puts("Can't convert weird image formats to grayscale.");
        return NULL;
    }
    
    i_width  = cairo_image_surface_get_width(px_surface);
    i_height = cairo_image_surface_get_height(px_surface);

    pc_result = malloc(i_width * i_height);
    px_surface_data = (uint32_t *)cairo_image_surface_get_data(px_surface);

    for (int i = 0; i < i_height * i_width; i++) {
        x_pixel = px_surface_data[i];
        // use "luminosity" method of grayscaling
        pc_result[i] = round(
            (x_pixel&RED_MASK)   *0.21 +
            (x_pixel&GREEN_MASK) *0.71 +
            (x_pixel&BLUE_MASK)  *0.07
            );
    }

    return pc_result;
}

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
    void                 *pv_algorithm_config, *pv_centroids;

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

    pv_centroids = (*x_algorithm.pf_algorithm)(
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
    px_cairo = cairo_create(px_surface);
    cairo_set_source_rgb(px_cairo, 1.0, 0.0, 0.0);
    for (int i = 0; i < i_centroids_length; i++) {
        switch (x_algorithm.i_centroid_type) {
        case CENTROID_TYPE_BASE:
            ;
            centroid_base_t *px_centroid = ((centroid_base_t *)pv_centroids) + i;
            cairo_rectangle(px_cairo,
                            round(px_centroid->l_x / 1000000.0),
                            round(px_centroid->l_y / 1000000.0),
                            1, 1);
            cairo_fill(px_cairo);
            break;
        }
    }
    cairo_surface_write_to_png(px_surface, s_output_path);
    cairo_surface_destroy(px_surface);
    cairo_destroy(px_cairo);
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
