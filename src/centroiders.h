#ifndef CENTROID_H
#define CENTROID_H

// centroiding algorithms might return additional information, which only certain algorithms can take advantage of. Those centroiding algorithms should return 
typedef struct centroid_base_t {
    long l_x; // pixels*10^6
    long l_y; // pixels*10^6
} centroid_base_t;

typedef struct centroid_magnitude_t {
    centroid_base_t x_base;
    long l_magnitude;
} centroid_magnitude_t;

enum {
    CENTROID_TYPE_BASE,
    CENTROID_TYPE_MAGNITUDE,
};

typedef void *(*centroid_function_t)(unsigned char *pc_image,
                                     int  i_image_width,
                                     int  i_image_height,
                                     int  *pi_result_length,
                                     void *pv_config);
typedef void *(*centroid_config_function_t)(void);

typedef struct centroid_algorithm_t {
    char *s_name;
    int  i_centroid_type;

    centroid_function_t        pf_algorithm;
    centroid_config_function_t pf_config;
} centroid_algorithm_t;

centroid_algorithm_t x_centroid_algorithm(int i);

#endif
