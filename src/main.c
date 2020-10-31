#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "catalog-generic.h"
#include "database-builders.h"

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

static void v_build_catalog() {
    FILE            *px_tsv;
    char            s_default_tsv_path[] = "bright-star-database.tsv";
    char            s_tsv_path[256];
    catalog_t       *px_catalog;
    int             i_builder_choice;
    db_builder_algorithm_t x_algorithm;
    void            *pv_algorithm_config, *pv_db;

    printf("Location of tsv file (%s): ", s_default_tsv_path);
    fgets(s_tsv_path, 256 - 1, stdin);
    // get rid of trailing newline
    s_tsv_path[strlen(s_tsv_path)-1] = '\0';

    if (s_tsv_path[0] == '\0') {
        strcpy(s_tsv_path, s_default_tsv_path);
    }

    px_catalog = px_bsd_parse(s_tsv_path);

    printf("Read %ld stars from catalog.\n", px_catalog->l_stars_length);
    for (int i = 0; x_db_builder_algorithm(i).i_config_size; i++) {
        printf("%d) %s\n", i, x_db_builder_algorithm(i).s_name);
    }
    printf("Choose database builder: ");
    scanf("%d", &i_builder_choice);
    getchar();
    x_algorithm = x_db_builder_algorithm(i_builder_choice);
    if (x_algorithm.i_config_size == 0) {
        puts("Not found...");
        return;
    }
    pv_algorithm_config = (*x_algorithm.pf_config)();
    pv_db = (*(x_algorithm.pf_algorithm))(px_catalog, pv_algorithm_config);
    (*x_algorithm.pf_stats)(pv_db);
}

static void v_centroids_find() {
        
}

int main() {
    puts("LOST: Open-source Star Tracker");

    while (1) {
        printf("Would you like to Build Catalog (b), Find Centroids (c), or Quit (q)? ");
        char choice;
        scanf(" %c", &choice);
        getchar();
        switch (choice) {
        case 'b':
            v_build_catalog();
            break;
        case 'c':
            v_centroids_find();
        default:
            puts("Bye!");
            exit(0);
        }
    }
    puts("");
}
