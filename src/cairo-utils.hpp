#ifndef CAIRO_UTILS_H
#define CAIRO_UTILS_H

#include <cairo/cairo.h>

#include <string>

#include "star-utils.hpp"
#include "attitude-utils.hpp"

namespace lost {

/// An 8-bit grayscale 2d image
class Image {
public:
    /**
     * The raw pixel data in the image.
     * This is an array of pixels, of length width*height. Each pixel is a single byte. A zero byte is pure black, and a 255 byte is pure white. Support for pixel resolution greater than 8 bits may be added in the future.
     */
    unsigned char *image;

    int width;
    int height;
};

// Convert a cairo surface to array of grayscale bytes
unsigned char *SurfaceToGrayscaleImage(cairo_surface_t *cairoSurface);
cairo_surface_t *GrayscaleImageToSurface(const unsigned char *, const int width, const int height);

void SurfacePlot(std::string description,
                 cairo_surface_t *cairoSurface,
                 const Stars &stars,
                 const StarIdentifiers *starIds,
                 const Catalog *catalog,
                 const Attitude *attitude,
                 double red,
                 double green,
                 double blue,
                 double alpha,
                 bool rawStarIndexes = false);

}

#endif
