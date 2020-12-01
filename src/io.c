#include "io.h"

#include <stdio.h>
#include <inttypes.h>
#include <math.h>

catalog_t *px_bsd_parse(char *s_path) {
    catalog_t      *px_result;
    catalog_star_t *px_star_current;
    FILE           *px_file;
    long           l_raj2000_h, l_raj2000_l, // high and low parts
                   l_dej2000_h, l_dej2000_l;
    int            i_magnitude_h, i_magnitude_l;
    char           c_weird;

    px_result = (catalog_t *)calloc(sizeof(catalog_t), 1);

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
            px_result->px_stars = (catalog_star_t *)realloc(px_result->px_stars,
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

unsigned char *pc_surface_grayscale(cairo_surface_t *px_surface) {
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

    pc_result = (unsigned char *)malloc(i_width * i_height);
    px_surface_data = (uint32_t *)cairo_image_surface_get_data(px_surface);

    for (int i = 0; i < i_height * i_width; i++) {
        x_pixel = px_surface_data[i];
        // use "luminosity" method of grayscaling
        pc_result[i] = round(
            (x_pixel>>16 &0xFF) *0.21 +
            (x_pixel>>8  &0xFF) *0.71 +
            (x_pixel     &0xFF) *0.07
            );
    }

    return pc_result;
}

void v_surface_plot_centroids(cairo_surface_t *px_surface,
                              star_t *px_centroids,
                              int i_centroids_length,
                              double red,
                              double green,
                              double blue,
                              double alpha) {
    cairo_t *px_cairo;

    px_cairo = cairo_create(px_surface);
    cairo_set_source_rgba(px_cairo, red, green, blue, alpha);
    cairo_set_line_width(px_cairo, 1.0); // I wonder what <1.0 does?

    while (i_centroids_length --> 0) {
        if (px_centroids->f_radius_x > 0.0f) {
            float f_radius_x = px_centroids->f_radius_x;
            float f_radius_y = px_centroids->f_radius_y > 0.0f ?
                px_centroids->f_radius_y : f_radius_x;

            // Rectangles should be entirely /outside/ the radius of the star, so the star is fully
            // visible.
            cairo_rectangle(px_cairo,
                            floor(px_centroids->f_x - f_radius_x) - 1,
                            floor(px_centroids->f_y - f_radius_y) - 1,
                            ceil(f_radius_x) * 2 + 2,
                            ceil(f_radius_y) * 2 + 2);
            cairo_stroke(px_cairo);
        } else {
            cairo_rectangle(px_cairo,
                            floor(px_centroids->f_x),
                            floor(px_centroids->f_y),
                            1, 1);
            cairo_fill(px_cairo);
        }
        px_centroids++;
    }
    cairo_destroy(px_cairo);
}

