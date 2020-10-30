#ifndef DATABASE_BUILDER_H
#define DATABASE_BUILDER_H

#include "catalog-generic.h"

#include "database-builder-pairwise-kvector.h"

static char *s_database_builder_name(int i) {
	switch (i) {
	case 0:
		return "Pairwise K-Vector (Pyramid)";
	default:
		return NULL;
	}
}

typedef void (*database_builder_t)(catalog_t *x_catalog);

static database_builder_t pf_database_builder(int i) {
	switch (i) {
	case 0:
		return (database_builder_t)(*px_pairwise_kvector_build);
	default:
		return NULL;
	}
}

#endif
