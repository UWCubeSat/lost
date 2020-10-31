#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include <stdio.h>
#include <stdlib.h>

#include "catalog-generic.h"

typedef void *(*db_builder_function_t)(catalog_t *px_catalog, void *pv_config);
typedef void *(*db_config_function_t)(void);
typedef void (*db_stats_function_t)(void *pv_db);

typedef struct db_builder_algorithm_t {
    char *s_name;
    int  i_config_size;
    
    db_config_function_t  pf_config;
    db_builder_function_t pf_algorithm;
    db_stats_function_t   pf_stats;
} db_builder_algorithm_t;

db_builder_algorithm_t x_db_builder_algorithm(int i);

#endif
