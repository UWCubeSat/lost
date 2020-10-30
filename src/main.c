#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "catalog-generic.h"
#include "database-builders.h"

static void v_build_catalog() {
	FILE            *px_tsv;
	char            s_default_tsv_path[] = "bright-star-database.tsv";
	char            s_tsv_path[256];
	catalog_t       *px_catalog;
	int             builder_choice;
	void            *database;

	printf("Location of tsv file (%s): ", s_default_tsv_path);
	fgets(s_tsv_path, 256 - 1, stdin);
	// get rid of trailing newline
	s_tsv_path[strlen(s_tsv_path)-1] = '\0';

	if (s_tsv_path[0] == '\0') {
		strcpy(s_tsv_path, s_default_tsv_path);
	}

	px_catalog = px_bsd_parse(s_tsv_path);

	printf("Read %ld stars from catalog.\n", px_catalog->l_stars_length);
	for (int i = 0; s_database_builder_name(i); i++) {
		printf("%d) %s\n", i, s_database_builder_name(i));
	}
	printf("Choose database builder: ");
	scanf("%d", &builder_choice);
	getchar();
	if (!pf_database_builder(builder_choice)) {
		puts("Not found...");
		return;
	}
	(*(pf_database_builder(builder_choice)))(px_catalog);
}

int main() {
	puts("LOST: Open-source Star Tracker");

	while (1) {
		printf("Would you like to Build Catalog (b), or Quit (q)? ");
		char choice;
		scanf(" %c", &choice);
		getchar();
		switch (choice) {
		case 'b':
			v_build_catalog();
			break;
		default:
			puts("Bye!");
			exit(0);
		}
	}
	puts("");
}
