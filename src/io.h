// I/O stuff, such as Cairo and star catalog interactions.

#ifndef IO_H
#define IO_H

#include <cairo/cairo.h>

#include "catalog-generic.h"
#include "centroiders.h"

// parse the Bright Star Catalog tsv file (see Bash)
catalog_t *
px_bsd_parse(char *s_path);
// Convert a cairo surface to array of grayscale bytes
unsigned char *pc_surface_grayscale(cairo_surface_t *px_surface);
// plot dots at the specified centroids
void v_surface_plot_centroids(cairo_surface_t *px_surface,
                              star_t *px_centroids,
                              int i_centroids_length,
                              double red,
                              double green,
                              double blue,
                              double alpha);
// take an astrometry download from the bash script, and parse it into stuff.
void v_astrometry_parse(char            *s_path,
                        cairo_surface_t **ppx_surface,   // image data
                        star_t      **ppx_centroids, // centroids according to astrometry
                        int             *pi_centroids_length); // TODO: fov, actual angle, etc
#endif
