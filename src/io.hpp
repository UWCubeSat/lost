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
#include "databases.hpp"

namespace lost {

const char kNoDefaultArgument = 0;

class PromptedOutputStream {
public:
    PromptedOutputStream(std::string filePath);
    ~PromptedOutputStream();
    std::ostream &Stream() { return *stream; };

private:
    bool isFstream;
    std::ostream *stream;
};

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

/**
 * @brief
 * @details
 */
class Image {
public:
    /// @brief
    unsigned char *image;

    /// @brief
    int width;

    /// @brief
    int height;
};

////////////////////
// PIPELINE INPUT //
////////////////////

/**
 * @brief
 * @details
 */
class PipelineOptions {
public:
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
    type prop = defaultVal;
#include "./pipeline-options.hpp"
#undef LOST_CLI_OPTION
};

/**
 * @brief Represents the input and expected outputs of a pipeline run.
 * @details
 */
class PipelineInput {
public:
    /// @brief
    virtual ~PipelineInput(){};
    /**
     * @brief
     * @return
     */
    virtual const Image *InputImage() const { return NULL; };

    /**
     * @brief
     * @return
     */
    virtual const Catalog &GetCatalog() const = 0;

    /**
     * @brief
     * @return
     */
    virtual const Stars *InputStars() const { return NULL; };

    /**
     * @brief Whether the input stars have identification information.
     * @return
     */
    virtual const StarIdentifiers *InputStarIds() const { return NULL; };

    /**
     * @brief For tracking
     * @return
     */
    virtual const Attitude *InputAttitude() const { return NULL; };

    /**
     * @brief
     * @return
     */
    virtual const Camera *InputCamera() const { return NULL; };

    /**
     * @brief
     * @return
     */
    virtual const Stars *ExpectedStars() const { return InputStars(); };

    /**
     * @brief
     * @return
     */
    virtual const StarIdentifiers *ExpectedStarIds() const { return InputStarIds(); };

    /**
     * @brief
     * @return
     */
    virtual const Attitude *ExpectedAttitude() const { return InputAttitude(); };

    cairo_surface_t *InputImageSurface() const;
};

class GeneratedPipelineInput : public PipelineInput {
public:
    // TODO: correct params
    GeneratedPipelineInput(const Catalog &, Attitude, Camera,
                           float observedReferenceBrightness, float starSpreadStdDev,
                           float sensitivity, float darkCurrent, float readNoiseStdDev,
                           Attitude motionBlurDirection, float exposureTime, float readoutTime,
                           bool shotNoise, int oversampling,
                           int numFalseStars, int falseMinMagnitude, int falseMaxMagnitude,
                           int seed);
                           

    const Image *InputImage() const { return &image; };
    const Stars *InputStars() const { return &stars; };
    const Camera *InputCamera() const { return &camera; };
    const StarIdentifiers *InputStarIds() const { return &starIds; };
    bool InputStarsIdentified() const { return true; };
    const Attitude *InputAttitude() const { return &attitude; };
    const Catalog &GetCatalog() const { return catalog; };

private:
    // we don't use an Image here because we want to 
    std::vector<unsigned char> imageData;
    Image image;
    Stars stars;
    Camera camera;
    Attitude attitude;
    const Catalog &catalog;
    StarIdentifiers starIds;
};

typedef std::vector<std::unique_ptr<PipelineInput>> PipelineInputList;

PipelineInputList GetPipelineInput(const PipelineOptions &values);

class PngPipelineInput : public PipelineInput {
public:
    PngPipelineInput(cairo_surface_t *, Camera, const Catalog &);

    const Image *InputImage() const { return &image; };
    const Camera *InputCamera() const { return &camera; };
    const Catalog &GetCatalog() const { return catalog; };

private:
    Image image;
    Camera camera;
    const Catalog &catalog;
};

/////////////////////
// PIPELINE OUTPUT //
/////////////////////

/**
 * @brief
 * @details
 */
struct PipelineOutput {
    /// @brief
    std::unique_ptr<Stars> stars;

    /// @brief
    std::unique_ptr<StarIdentifiers> starIds;

    /// @brief
    std::unique_ptr<Attitude> attitude;

    /**
     * @brief The catalog that the indices in starIds refer to
     * @todo Don't store it here
     */
    Catalog catalog;
};

struct StarIdComparison {
    int numCorrect;
    int numIncorrect;
    int numTotal;
    float fractionCorrect;
    float fractionIncorrect;
};

std::ostream &operator<<(std::ostream &, const Camera &);

// actualStars is optional, in which case it's assumed that expectedStars was passed to the star-id
StarIdComparison StarIdsCompare(const StarIdentifiers &expected, const StarIdentifiers &actual,
                                const Catalog &expectedCatalog, const Catalog &actualCatalog,
                                float centroidThreshold,
                                const Stars *expectedStars, const Stars *actualStars);

//////////////
// PIPELINE //
//////////////

/**
 * @brief A pipeline is a set of algorithms that describes all or part of the star-tracking "pipeline"
 * @details
 */
class Pipeline {
    friend Pipeline SetPipeline(const PipelineOptions &values);

public:
    /**
     * @brief
     * @note Pointers just so they're nullable
     */
    Pipeline() = default;
    Pipeline(CentroidAlgorithm *, StarIdAlgorithm *, AttitudeEstimationAlgorithm *, unsigned char *);
    PipelineOutput Go(const PipelineInput &);
    std::vector<PipelineOutput> Go(const PipelineInputList &);

private:
    std::unique_ptr<CentroidAlgorithm> centroidAlgorithm;
    int centroidMinMagnitude = 0;
    std::unique_ptr<StarIdAlgorithm> starIdAlgorithm;
    std::unique_ptr<AttitudeEstimationAlgorithm> attitudeEstimationAlgorithm;
    std::unique_ptr<unsigned char[]> database;
};

Pipeline SetPipeline(const PipelineOptions &values);

// TODO: rename. Do something with the output
void PipelineComparison(const PipelineInputList &expected,
                        const std::vector<PipelineOutput> &actual,
                        const PipelineOptions &values);

////////////////
// DB BUILDER //
////////////////

// TODO: rename
Catalog PromptNarrowedCatalog(const Catalog &);

/**
 * @brief
 * @details
 */
class DatabaseOptions {
public:
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
    type prop = defaultVal;
#include "database-options.hpp"
#undef LOST_CLI_OPTION   
};

// unlike the other algorithm prompters, db builders aren't a
// typedef void (*DbBuilder)(MultiDatabaseBuilder &, const Catalog &);
void GenerateDatabases(MultiDatabaseBuilder &, const Catalog &, const DatabaseOptions &values);
// void PromptDatabases(MultiDatabaseBuilder &, const Catalog &);

/////////////////////
// INSPECT CATALOG //
/////////////////////

void InspectCatalog();

}

#endif
