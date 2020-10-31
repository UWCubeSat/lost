#include "database-builders.h"

#include <stdio.h>
#include <stdlib.h>

// PAIRWISE

typedef struct pairwise_db_config_t {
    long l_fov; // degrees*10^6
} pairwise_db_config_t;

void *px_pairwise_db_build(catalog_t *px_catalog, void *pv_config) {
    pairwise_db_config_t *px_config = (pairwise_db_config_t *)pv_config;

    return NULL; // TODO: implement
}

void *pv_pairwise_db_config() {
    pairwise_db_config_t *result = malloc(sizeof(pairwise_db_config_t));

    printf("FOV (in millionths of a degree): ");
    scanf("%ld", &result->l_fov);

    return result;
}

void pv_pairwise_db_stats() {
    puts("In the future, this will print, for example, the size of the database.");
}

db_builder_algorithm_t x_pairwise_algorithm = {
    .s_name        = "Pairwise Distances (SLA, Pyramid, Geometric Voting)",
    .i_config_size = sizeof(pairwise_db_config_t),
    .pf_config     = *pv_pairwise_db_config,
    .pf_algorithm  = *px_pairwise_db_build,
    .pf_stats      = *pv_pairwise_db_stats,
};

db_builder_algorithm_t x_empty_algorithm = { 0 };

db_builder_algorithm_t x_db_builder_algorithm(int i) {
    switch (i) {
    case 0:
    return x_pairwise_algorithm;
    default:
        return x_empty_algorithm;
    }
}

