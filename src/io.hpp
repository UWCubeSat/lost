// I/O stuff, such as Cairo and star catalog interactions.

#ifndef IO_H
#define IO_H

#include <vector>
#include <map>
#include <utility>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>

#include <cairo/cairo.h>

#ifndef CAIRO_HAS_PNG_FUNCTIONS
#error LOST requires Cairo to be compiled with PNG support
#endif

#include "catalog-generic.hpp"
#include "centroiders.hpp"
#include "star-id.hpp"
#include "camera.hpp"
#include "attitude-utils.hpp"
#include "attitude-estimators.hpp"

namespace lost {

void RegisterCliArgs(int, char **);
bool HasNextCliArg();
std::string NextCliArg();

template <typename S>
S Prompt(const std::string &prompt) {
    S result;
    std::cout << prompt << ": ";
    if (HasNextCliArg()) {
        std::string nextArg = NextCliArg();
        std::cout << nextArg << std::endl;
        std::stringstream(nextArg) >> result;
    } else {
        std::cin >> result;
    }
    return result;
}


// can prompt the user between multiple options.
template <typename S>
class InteractiveChoice {
public:
    // prompt the user until they enter a valid option
    S Prompt(const std::string &) const;
    void Register(std::string, std::string, S);
private:
    std::map<std::string, std::pair<std::string, S>> options;
};

template <typename S>
void InteractiveChoice<S>::Register(std::string shortKey, std::string longKey, S value) {
    options.emplace(shortKey, std::make_pair(longKey, value));
}

template <typename S>
S InteractiveChoice<S>::Prompt(const std::string &prompt) const {
    std::string userChoice;
    bool useCli = HasNextCliArg();
    while (1) {
        if (!useCli) {
            for (auto it = options.begin(); it != options.end(); it++) {
                std::cout << "(" << it->first << ") " << it->second.first << std::endl;
            }
        }
        userChoice = lost::Prompt<std::string>(prompt);
        auto found = options.find(userChoice);
        if (found == options.end()) {
            std::cout << "Peace was never an option." << std::endl;
            if (useCli) {
                exit(1);
            }
        } else {
            return found->second.second; // I've been programming too much Lisp
        }
    }

}

// parse the Bright Star Catalog tsv file (see Bash)
std::vector<CatalogStar> BsdParse(std::string tsvPath);
// Convert a cairo surface to array of grayscale bytes
unsigned char *SurfaceToGrayscaleImage(cairo_surface_t *cairoSurface);
cairo_surface_t *GrayscaleImageToSurface(const unsigned char *, const int width, const int height);
// plot dots at the specified centroids
void SurfacePlotCentroids(cairo_surface_t *cairoSurface,
                          std::vector<Stars> centroids,
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
typedef CentroidAlgorithm *(*CentroidAlgorithmFactory)();

InteractiveChoice<CentroidAlgorithmFactory> makeCentroidAlgorithmChoice();

class Image {
public:
    unsigned char *image;
    int width;
    int height;
};

////////////////////
// PIPELINE INPUT //
////////////////////

// represents the input and expected outputs of a pipeline run.
class PipelineInput {
public:
    virtual const Image *InputImage() const { return NULL; };
    // TODO: a more solid type here
    virtual const void *InputDatabase() const { return NULL; };
    virtual const Stars *InputCentroids() const { return NULL; };
    virtual const Stars *InputStars() const { return NULL; };
    // for tracking
    virtual const Attitude *InputAttitude() const { return NULL; };
    virtual const Camera *InputCamera() const { return NULL; };

    virtual const Stars *ExpectedCentroids() const { return InputCentroids(); };
    virtual const Stars *ExpectedStars() const { return InputStars(); };
    virtual const Attitude *ExpectedAttitude() const { return InputAttitude(); };
};

typedef std::vector<std::unique_ptr<PipelineInput>> PipelineInputList;

PipelineInputList PromptPipelineInput();

/////////////////////
// PIPELINE OUTPUT //
/////////////////////

class PipelineOutput {
public:
    std::unique_ptr<Stars> centroids;
    std::unique_ptr<Stars> stars;
    std::unique_ptr<Attitude> attitude;
};

//////////////
// PIPELINE //
//////////////

// a pipeline is a set of algorithms that describes all or part of the star-tracking "pipeline"

class Pipeline {
    friend Pipeline PromptPipeline();
public:
    PipelineOutput Go(const PipelineInput &);
    std::vector<PipelineOutput> Go(const PipelineInputList &);
private:
    std::unique_ptr<CentroidAlgorithm> centroidAlgorithm;
    std::unique_ptr<StarIdAlgorithm> starIdAlgorithm;
    std::unique_ptr<AttitudeEstimationAlgorithm> attitudeEstimationAlgorithm;
};

Pipeline PromptPipeline();

// ask the user what to do with actual and expected outputs
void PromptPipelineComparison(const PipelineInputList &expected,
                              const std::vector<PipelineOutput> &actual);

}

#endif
