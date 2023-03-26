// I/O stuff, such as Cairo and star catalog interactions.

#ifndef IO_H
#define IO_H

#include <cairo/cairo.h>

#include <random>
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>


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

/// An output stream which might be a file or stdout
class UserSpecifiedOutputStream {
public:
    explicit UserSpecifiedOutputStream(std::string filePath, bool isBinary);
    ~UserSpecifiedOutputStream();

    /// return the inner output stream, suitable for use with <<
    std::ostream &Stream() { return *stream; };

private:
    bool isFstream;
    std::ostream *stream;
};

// use the environment variable LOST_BSC_PATH, or read from ./bright-star-catalog.tsv
const Catalog &CatalogRead();
// Convert a cairo surface to array of grayscale bytes
unsigned char *SurfaceToGrayscaleImage(cairo_surface_t *cairoSurface);
cairo_surface_t *GrayscaleImageToSurface(const unsigned char *, const int width, const int height);

// take an astrometry download from the bash script, and parse it into stuff.
// void v_astrometry_parse(std::string
//                         cairo_surface_t **pcairoSurface,   // image data
//                         Star      **ppx_centroids, // centroids according to astrometry
//                         int             *pi_centroids_length); // TODO: fov, actual angle, etc

// type for functions that create a centroid algorithm (by prompting the user usually)

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

////////////////////
// PIPELINE INPUT //
////////////////////

/// The command line options passed when running a pipeline
class PipelineOptions {
public:
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
    type prop = defaultVal;
#include "./pipeline-options.hpp"
#undef LOST_CLI_OPTION
};

/**
 * Represents the input and expected outputs of a pipeline run.
 * This is all the data about the pipeline we are about to run that can be gathered without actually running the pipeline. The "input" members return the things that will be fed into the pipeline. The "expected" methods return the "correct" outputs, i.e. what a perfect star tracking algorithm would output (this is only meaningful whe the image is generated, otherwise we don't know the correct output!). The "expected" methods are used to evaluate the quality of our algorithms. Some of the methods (both input and expected) may return NULL for certain subclasses.
 * By default, the "expected" methods return the corresponding inputs, which is reasonable behavior unless you are trying to intentionally introduce error into the inputs.
 */
class PipelineInput {
public:
    virtual ~PipelineInput(){};

    virtual const Image *InputImage() const { return NULL; };
    /// The catalog to which catalog indexes returned from other methods refer.
    virtual const Catalog &GetCatalog() const = 0;
    virtual const Stars *InputStars() const { return NULL; };
    /// The centroid indices in the StarIdentifiers returned from InputStarIds should be indices
    /// into InputStars(), not ExpectedStars(), when present, because otherwise it's useless.
    virtual const StarIdentifiers *InputStarIds() const { return NULL; };
    /// Only used in tracking mode, in which case it is an estimate of the current attitude based on the last attitude, IMU info, etc.
    virtual const Attitude *InputAttitude() const { return NULL; };
    virtual const Camera *InputCamera() const { return NULL; };
    /// Convert the InputImage() output into a cairo surface
    cairo_surface_t *InputImageSurface() const;

    virtual const Stars *ExpectedStars() const { return InputStars(); };
    /// Centroid indices in the StarIdentifiers returned from ExpectedStarIds should be indices into
    /// ExpectedStars(), /not/ InputStars(). This is in contrast to InputStarIds. If you need to
    /// compare ExpectedStarIds against the input stars, then you should use some function which
    /// uses simple heuristics to match the input stars and expected stars (eg based on distance).
    /// Cf how the star-ID comparator works for a reference implementation.
    virtual const StarIdentifiers *ExpectedStarIds() const { return InputStarIds(); };
    virtual const Attitude *ExpectedAttitude() const { return InputAttitude(); };
};

/**
 * A pipeline input which is generated (fake image).
 * Uses a pretty decent image generation algorithm to create a fake image as the basis for a pipeline input. Since we know everything about the generated image, all of the input and expected methods are available.
 */
class GeneratedPipelineInput : public PipelineInput {
public:
    GeneratedPipelineInput(const Catalog &, Attitude, Camera, std::default_random_engine *,

                           bool centroidsOnly,
                           float observedReferenceBrightness, float starSpreadStdDev,
                           float sensitivity, float darkCurrent, float readNoiseStdDev,
                           Attitude motionBlurDirection, float exposureTime, float readoutTime,
                           bool shotNoise, int oversampling,
                           int numFalseStars, int falseMinMagnitude, int falseMaxMagnitude,
                           float perturbationStddev);


    const Image *InputImage() const { return &image; };
    const Stars *InputStars() const { return &stars; };
    const Camera *InputCamera() const { return &camera; };
    const StarIdentifiers *InputStarIds() const { return &starIds; };
    bool InputStarsIdentified() const { return true; };
    const Attitude *InputAttitude() const { return &attitude; };
    const Catalog &GetCatalog() const { return catalog; };

private:
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

/// A pipeline input created by reading a PNG from a file on disk.
class PngPipelineInput : public PipelineInput {
public:
    PngPipelineInput(cairo_surface_t *, Camera, const Catalog &);
    ~PngPipelineInput();

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
 * @brief The result of running a pipeline.
 * @details Also stores intermediate outputs, not just the final attitude.
 */
struct PipelineOutput {
    std::unique_ptr<Stars> stars;
    std::unique_ptr<StarIdentifiers> starIds;
    std::unique_ptr<Attitude> attitude;

    /**
     * @brief The catalog that the indices in starIds refer to
     * @todo Don't store it here
     */
    Catalog catalog;
};

/// The result of comparing an actual star identification with the true star idenification, used for testing and benchmarking.
struct StarIdComparison {
    /// The number of centroids which were identified as the correct catalog star.
    int numCorrect;

    /// The number of centroids which were identified, but as the wrong catalog star.
    int numIncorrect;

    /// The total number of true stars in the image (the number the ideal star-id algorithm would identify)
    int numTotal;
};

std::ostream &operator<<(std::ostream &, const Camera &);

//////////////
// PIPELINE //
//////////////

/**
 * @brief A set of algorithms that describes all or part of the star-tracking "pipeline"
 * @details A centroiding algorithm identifies the (x,y) pixel coordinates of each star detected in the raw image. The star id algorithm then determines which centroid corresponds to which catalog star. Finally, the attitude estimation algorithm determines the orientation of the camera based on the centroids and identified stars.
 */
class Pipeline {
    friend Pipeline SetPipeline(const PipelineOptions &values);

public:
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

/// Commannd line options when using the `database` command.
class DatabaseOptions {
public:
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
    type prop = defaultVal;
#include "database-options.hpp"
#undef LOST_CLI_OPTION
};

// unlike the other algorithm prompters, db builders aren't a
// typedef void (*DbBuilder)(MultiDatabaseBuilder &, const Catalog &);
void GenerateDatabases(MultiDatabaseBuilder *, const Catalog &, const DatabaseOptions &values);
// void PromptDatabases(MultiDatabaseBuilder &, const Catalog &);

/////////////////////
// INSPECT CATALOG //
/////////////////////

void InspectCatalog();

}

#endif
