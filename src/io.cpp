#include "io.hpp"

#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <vector>
#include <sstream>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>

namespace lost {

char **argv = NULL;
int argc = 0;

void RegisterCliArgs(int newArgc, char **newArgv) {
    argv = newArgv + 1;
    argc = newArgc - 1;
}

bool HasNextCliArg() {
    return argc > 0;
}

std::string NextCliArg() {
    if (HasNextCliArg()) {
        argc--;
        return std::string(*argv++);
    }
    return std::string("You incompetent fool!");
}

std::vector<CatalogStar> BsdParse(std::string tsvPath) {
    std::vector<CatalogStar> result;
    FILE           *file;
    long           raj2000High, raj2000Low, // high and low parts
                   dej2000High, dej2000Low;
    int            magnitudeHigh, magnitudeLow;
    char           weird;

    file = fopen(tsvPath.c_str(), "r");
    if (file == NULL) {
        printf("Error opening file: %s\n", strerror(errno));
        return result; // TODO
    }

    while (EOF != fscanf(file, "%ld.%ld|%ld.%ld|%*d|%c|%d.%d",
                         &raj2000High, &raj2000Low,
                         &dej2000High, &dej2000Low,
                         &weird,
                         &magnitudeHigh, &magnitudeLow)) {
           
        result.push_back(CatalogStar(
                             raj2000High * 1000000 + raj2000Low,
                             dej2000High * 1000000 + dej2000Low,
                             magnitudeHigh * 100 + magnitudeLow,
                             weird != ' '));
    }

    fclose(file);
    return result;
}

unsigned char *SurfaceToGrayscaleImage(cairo_surface_t *cairoSurface) {
    int width, height;
    unsigned char *result;
    uint32_t *cairoImage, pixel;

    if (cairo_image_surface_get_format(cairoSurface) != CAIRO_FORMAT_ARGB32 &&
        cairo_image_surface_get_format(cairoSurface) != CAIRO_FORMAT_RGB24) {
        puts("Can't convert weird image formats to grayscale.");
        return NULL;
    }
    
    width  = cairo_image_surface_get_width(cairoSurface);
    height = cairo_image_surface_get_height(cairoSurface);

    result = (unsigned char *)malloc(width * height);
    cairoImage = (uint32_t *)cairo_image_surface_get_data(cairoSurface);

    for (int i = 0; i < height * width; i++) {
        pixel = cairoImage[i];
        // use "luminosity" method of grayscaling
        result[i] = round(
            (pixel>>16 &0xFF) *0.21 +
            (pixel>>8  &0xFF) *0.71 +
            (pixel     &0xFF) *0.07
            );
    }

    return result;
}

void SurfacePlotCentroids(cairo_surface_t *cairoSurface,
                          std::vector<Star> centroids,
                          double red,
                          double green,
                          double blue,
                          double alpha) {
    cairo_t *cairoCtx;

    cairoCtx = cairo_create(cairoSurface);
    cairo_set_source_rgba(cairoCtx, red, green, blue, alpha);
    cairo_set_line_width(cairoCtx, 1.0); // I wonder what <1.0 does?

    for (const Star &centroid : centroids) {
        if (centroid.radiusX > 0.0f) {
            float radiusX = centroid.radiusX;
            float radiusY = centroid.radiusY > 0.0f ?
                centroid.radiusY : radiusX;

            // Rectangles should be entirely /outside/ the radius of the star, so the star is fully
            // visible.
            cairo_rectangle(cairoCtx,
                            floor(centroid.x - radiusX) - 1,
                            floor(centroid.y - radiusY) - 1,
                            ceil(radiusX) * 2 + 2,
                            ceil(radiusY) * 2 + 2);
            cairo_stroke(cairoCtx);
        } else {
            cairo_rectangle(cairoCtx,
                            floor(centroid.x),
                            floor(centroid.y),
                            1, 1);
            cairo_fill(cairoCtx);
        }
    }
    cairo_destroy(cairoCtx);
}

// ALGORITHM PROMPTERS

CentroidAlgorithm *DummyCentroidAlgorithmPrompt() {
    int numStars = Prompt<int>("How many stars to generate");
    return new DummyCentroidAlgorithm(numStars);
}

CentroidAlgorithm *CogCentroidAlgorithmPrompt() {
    return new CenterOfGravityAlgorithm();
}

InteractiveChoice<CentroidAlgorithmFactory> makeCentroidAlgorithmChoice() {
    InteractiveChoice<CentroidAlgorithmFactory> result;

    result.Register("dummy", "Random Centroid Algorithm", &DummyCentroidAlgorithmPrompt);
    result.Register("cog", "Center of Gravity Centroid Algorithm", &CogCentroidAlgorithmPrompt);

    return result;
}


}

