#include "centroiders.h"

#include <stdio.h>
#include <stdlib.h>

// DUMMY

typedef struct dummy_centroid_config_t {
    int num_stars;
} dummy_centroid_config_t;

void *pv_dummy_centroid(unsigned char *pc_image,
                        int i_image_width,
                        int i_image_height,
                        int *pi_result_length,
                        void *pv_config) {
    dummy_centroid_config_t *px_config = (dummy_centroid_config_t *)pv_config;
    centroid_base_t *px_result = malloc(sizeof(centroid_base_t) * px_config->num_stars);

    *pi_result_length = px_config->num_stars;
    
    for (int i = 0; i < px_config->num_stars; i++) {
        px_result[i].l_x = rand() % i_image_width * 1000000;
        px_result[i].l_y = rand() % i_image_height * 1000000;
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
