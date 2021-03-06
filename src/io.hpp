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

#include "centroiders.hpp"
#include "star-utils.hpp"
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
    std::cerr << prompt << ": ";
    if (HasNextCliArg()) {
        std::string nextArg = NextCliArg();
        std::cerr << nextArg << std::endl;
        std::stringstream(nextArg) >> result;
    } else {
        std::cin >> result;
    }
    return result;
}

template <typename S>
class InteractiveChoiceOption {
public:
    InteractiveChoiceOption(std::string shortName, std::string longName, S value)
        : shortName(shortName), longName(longName), value(value) { };

    std::string shortName;
    std::string longName;
    S value;
};

// can prompt the user between multiple options.
template <typename S>
class InteractiveChoice {
public:
    // prompt the user until they enter a valid option
    S Prompt(const std::string &) const;
    void Register(std::string, std::string, S);
private:
    std::vector<InteractiveChoiceOption<S>> options;
};

template <typename S>
void InteractiveChoice<S>::Register(std::string shortKey, std::string longKey, S value) {
    options.push_back(InteractiveChoiceOption<S>(shortKey, longKey, value));
}

template <typename S>
S InteractiveChoice<S>::Prompt(const std::string &prompt) const {
    std::string userChoice;
    bool useCli = HasNextCliArg();
    while (1) {
        if (!useCli) {
            for (const auto &option : options) {
                std::cerr << "(" << option.shortName << ") " << option.longName << std::endl;
            }
        }
        userChoice = lost::Prompt<std::string>(prompt);

        auto found = options.begin();
        while (found != options.end() && found->shortName != userChoice) {
            found++;
        }

        if (found == options.end()) {
            std::cerr << "Peace was never an option." << std::endl;
            if (useCli) {
                exit(1);
            }
        } else {
            return found->value;
        }
    }

}

class PromptedOutputStream {
public:
    PromptedOutputStream();
    ~PromptedOutputStream();
    std::ostream &Stream() { return *stream; };
private:
    bool isFstream;
    std::ostream *stream;
};

// prompts for an output stream, then calls the given function with it.
void WithOutputStream(void (*)(std::ostream *));

// use the environment variable LOST_BSC_PATH, or read from ./bright-star-catalog.tsv
std::vector<CatalogStar> &CatalogRead();
// Convert a cairo surface to array of grayscale bytes
unsigned char *SurfaceToGrayscaleImage(cairo_surface_t *cairoSurface);
cairo_surface_t *GrayscaleImageToSurface(const unsigned char *, const int width, const int height);

// take an astrometry download from the bash script, and parse it into stuff.
// void v_astrometry_parse(std::string 
//                         cairo_surface_t **pcairoSurface,   // image data
//                         Star      **ppx_centroids, // centroids according to astrometry
//                         int             *pi_centroids_length); // TODO: fov, actual angle, etc

// type for functions that create a centroid algorithm (by prompting the user usually)

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
    virtual const Stars *InputStars() const { return NULL; };
    // whether the input stars have identification information.
    virtual const StarIdentifiers *InputStarIds() const { return NULL; };
    // for tracking
    virtual const Quaternion *InputAttitude() const { return NULL; };
    virtual const Camera *InputCamera() const { return NULL; };

    virtual const Stars *ExpectedStars() const { return InputStars(); };
    virtual const StarIdentifiers *ExpectedStarIds() const { return InputStarIds(); };
    virtual const Quaternion *ExpectedAttitude() const { return InputAttitude(); };

    cairo_surface_t *InputImageSurface() const;
};

typedef std::vector<std::unique_ptr<PipelineInput>> PipelineInputList;

PipelineInputList PromptPipelineInput();

/////////////////////
// PIPELINE OUTPUT //
/////////////////////

struct PipelineOutput {
    std::unique_ptr<Stars> stars;
    std::unique_ptr<StarIdentifiers> starIds;
    std::unique_ptr<Quaternion> attitude;
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
    std::unique_ptr<unsigned char[]> database;
};

Pipeline PromptPipeline();

// ask the user what to do with actual and expected outputs
void PromptPipelineComparison(const PipelineInputList &expected,
                              const std::vector<PipelineOutput> &actual);

////////////////
// DB BUILDER //
////////////////

// unlike the other algorithm prompters, db builders aren't a 
typedef unsigned char *(*DbBuilder)(const Catalog &, long *);
DbBuilder PromptDbBuilder();

}

#endif
