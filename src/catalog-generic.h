#ifndef CATALOG_GENERIC_H
#define CATALOG_GENERIC_H

#include <stdbool.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

typedef struct catalog_star_t {
    long l_raj2000;           // *10^-6, right ascension
    long l_dej2000;           // *10^-6, declination
    int  i_magnitude;         // *10^-2
    int  i_weird;             // nonzero for binary, etc
} catalog_star_t;

typedef struct catalog_t {
    catalog_star_t *px_stars;
    long           l_stars_length;
} catalog_t;

#endif
