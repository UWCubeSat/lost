// I/O stuff, such as Cairo and star catalog interactions.

#ifndef IO_H
#define IO_H

#include <vector>
#include <map>
#include <utility>
#include <string>
#include <iostream>

#include <cairo/cairo.h>

#include "catalog-generic.hpp"
#include "centroiders.hpp"

namespace lost {

// can prompt the user between multiple options.
template <typename S>
class InteractiveChoice {
public:
    // prompt the user until they enter a valid option
    S Prompt(std::ostream &, std::istream &, const std::string &) const;
    void Register(std::string, std::string, S);
private:
    std::map<std::string, std::pair<std::string, S>> options;
};

template <typename S>
void InteractiveChoice<S>::Register(std::string shortKey, std::string longKey, S value) {
    options.emplace(shortKey, std::make_pair(longKey, value));
}

template <typename S>
S InteractiveChoice<S>::Prompt(std::ostream &os, std::istream &is, const std::string &prompt) const {
    std::string userChoice;
    while (1) {
        for (auto it = options.begin(); it != options.end(); it++) {
            os << "(" << it->first << ") " << it->second.first << std::endl;
        }
        os << prompt << ": ";
        is >> userChoice;
        auto found = options.find(userChoice);
        if (found == options.end()) {
            os << "Peace was never an option." << std::endl;
        } else {
            return found->second.second; // I've been programming too much Lisp
        }
    }

}

// parse the Bright Star Catalog tsv file (see Bash)
std::vector<CatalogStar> BsdParse(std::string tsvPath);
// Convert a cairo surface to array of grayscale bytes
unsigned char *SurfaceToGrayscaleImage(cairo_surface_t *cairoSurface);
// plot dots at the specified centroids
void SurfacePlotCentroids(cairo_surface_t *cairoSurface,
                          std::vector<Star> centroids,
                          double red,
                          double green,
                          double blue,
                          double alpha);

// take an astrometry download from the bash script, and parse it into stuff.
// void v_astrometry_parse(std::string 
//                         cairo_surface_t **pcairoSurface,   // image data
//                         Star      **ppx_centroids, // centroids according to astrometry
//                         int             *pi_centroids_length); // TODO: fov, actual angle, etc

// type for functions that create a centroid algorithm (by prompting the user usually)
typedef CentroidAlgorithm *(*CentroidAlgorithmFactory)(std::ostream &, std::istream &);

InteractiveChoice<CentroidAlgorithmFactory> makeCentroidAlgorithmChoice();

}

#endif