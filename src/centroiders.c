#include "centroiders.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// DUMMY

typedef struct dummy_centroid_config_t {
    int num_stars;
} dummy_centroid_config_t;

star_t *pv_dummy_centroid(unsigned char *pc_image,
                              int i_image_width,
                              int i_image_height,
                              int *pi_result_length,
                              void *pv_config) {
    dummy_centroid_config_t *px_config = (dummy_centroid_config_t *)pv_config;
    star_t *px_result = calloc(sizeof(star_t), px_config->num_stars);

    *pi_result_length = px_config->num_stars;
    
    for (int i = 0; i < px_config->num_stars; i++) {
        px_result[i].f_x = rand() % i_image_width;
        px_result[i].f_y = rand() % i_image_height;
        px_result[i].f_radius_x = 10.0;
    }

    return (void *)px_result;
}

void *pv_dummy_centroid_config() {
    dummy_centroid_config_t *px_result = malloc(sizeof(dummy_centroid_config_t));

    printf("How many stars to generate: ");
    scanf("%d", &px_result->num_stars);

    return (void *)px_result;
}

// PLUMBING

centroid_algorithm_t x_dummy_centroid_algorithm = {
    .s_name        = "Dummy Centroider",
    .pf_config     = *pv_dummy_centroid_config,
    .pf_algorithm  = *pv_dummy_centroid,
};

centroid_algorithm_t x_empty_centroid_algorithm = { 0 };

centroid_algorithm_t x_centroid_algorithm(int i) {
    switch (i) {
    case 0:
    return x_dummy_centroid_algorithm;
    default:
        return x_empty_centroid_algorithm;
    }
}

// OTHER CENTROID RELATED FUNCTIONS

float f_centroid_distance(star_t x_one, star_t x_two) {
    float f_dist_x = x_one.f_x - x_two.f_x;
    float f_dist_y = x_one.f_y - x_two.f_y;
    return sqrt(f_dist_x*f_dist_x + f_dist_y*f_dist_y);
}

// helper for x_compare_centroids
static void v_closest_centroids(float f_threshold, star_t *px_one,
                                int i_one_length, star_t *px_two,
                                int i_two_length, int *result) {
    for (int i = 0; i < i_one_length; i++) {
        float f_distance_closest = INFINITY;
        int   i_index_closest = -1;

        for (int k = 0; k < i_two_length; k++) {
            float f_distance_curr = f_centroid_distance(px_one[i], px_two[k]);
            if (f_distance_curr < f_threshold && f_distance_curr < f_distance_closest) {
                f_distance_closest = f_distance_curr;
                i_index_closest = k;
            }
        }

        result[i] = i_index_closest;
    }
}

centroid_comparison_t x_compare_centroids(float f_threshold,
                                          star_t *px_expected,
                                          int i_expected_length,
                                          star_t *px_actual,
                                          int i_actual_length) {

    centroid_comparison_t result = { 0 };
    // maps from indexes in each list to the closest centroid from other list
    int *pi_expected_to_actual = malloc(sizeof(int) * i_expected_length),
        *pi_actual_to_expected = malloc(sizeof(int) * i_actual_length);

    v_closest_centroids(f_threshold,
                        px_expected, i_expected_length,
                        px_actual, i_actual_length,
                        pi_expected_to_actual);
    
    v_closest_centroids(f_threshold,
                        px_actual, i_actual_length,
                        px_expected, i_expected_length,
                        pi_actual_to_expected);

    // any expected stars whose closest actual star does not refer to them are missing
    for (int i = 0; i < i_expected_length; i++) {
        if (pi_expected_to_actual[i] == -1 ||
            pi_actual_to_expected[pi_expected_to_actual[i]] != i) {
            result.i_stars_missing_num++;
        } else {
            result.f_mean_error += f_centroid_distance(px_expected[i],
                                                       px_actual[pi_expected_to_actual[i]]);
        }
    }
    result.f_mean_error /= (i_expected_length - result.i_stars_missing_num);


    // any actual star whose closest expected star does not refer to them is extra
    for (int i = 0; i < i_actual_length; i++) {
        if (pi_actual_to_expected[i] == -1 ||
            pi_expected_to_actual[pi_actual_to_expected[i]] != i) {
            result.i_stars_extra_num++;
        }
    }

    free(pi_expected_to_actual);
    free(pi_actual_to_expected);
}

void v_print_centroid_comparison(centroid_comparison_t x_comparison) {
    printf("Extra stars: %d\nMissing stars: %d\nMean error: %f\n",
           x_comparison.i_stars_extra_num,
           x_comparison.i_stars_missing_num,
           x_comparison.f_mean_error);
}
