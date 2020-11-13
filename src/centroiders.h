#ifndef CENTROID_H
#define CENTROID_H

// Different centroiding algorithms may return different amounts of data. Because there won't ever
// be /that/ many different centroids in a single image, it's acceptable to waste some space on
// un-used fields. Set unused fields to negative
typedef struct star_t {
    float f_x; // pixels*10^6
    float f_y; // pixels*10^6
    long  l_magnitude; // some relative number
    float f_radius_x;
    float f_radius_y;  // if omitted, but x is present, assume circular.
    // eccentricity?
} star_t;

typedef struct centroid_comparison_t {
    int i_stars_extra_num;    // stars in actual but not expected. Ideally 0
    int i_stars_missing_num;  // stars is expected but not actual. Ideally 0
    float f_mean_error;       // average distance from actual to expected star
    // I would add 99th percentile or something, but the really far away stars should really just
    // count in extra_num
} centroid_comparison_t;

typedef star_t *(*centroid_function_t)(unsigned char *pc_image,
                                           int  i_image_width,
                                           int  i_image_height,
                                           int  *pi_result_length,
                                           void *pv_config);
typedef void *(*centroid_config_function_t)(void);

typedef struct centroid_algorithm_t {
    char *s_name;

    centroid_function_t        pf_algorithm;
    centroid_config_function_t pf_config;
} centroid_algorithm_t;

centroid_algorithm_t x_centroid_algorithm(int i);

float f_centroid_distance(star_t x_one, star_t x_two);
// compare two lists of centroids to find differences.
centroid_comparison_t x_compare_centroids(float f_distance_threshold, // stars further apart than
                                                                      // are different stars.
                                          star_t *px_expected,
                                          int i_expected_length,
                                          star_t *px_actual,
                                          int i_actual_length);
void v_print_centroid_comparison(centroid_comparison_t x_comparison);

#endif
