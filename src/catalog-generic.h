#ifndef CATALOG_GENERIC_H
#define CATALOG_GENERIC_H

#include <stdbool.h>
#include <stdlib.h>

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

// @param s_path Pathname to read from on disk
// @return A malloc'd array of stars in the catalog
// Each line of the file should be like 001.610417|+64.196111|   7| | 5.59
catalog_t *px_bsd_parse(char *s_path) {
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

#endif
