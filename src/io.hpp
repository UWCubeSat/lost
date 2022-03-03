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
    
    class PromptedOutputStream 
    {
    public:
        PromptedOutputStream(std::string filePath);
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

    class Image
    {
    public:
        unsigned char *image;
        int width;
        int height;
    };

    ////////////////////
    // PIPELINE INPUT //
    ////////////////////

    // TODO maybe make a class/enum also for the different options options
    class PipelineOptions
    {
    public:
        std::string png = "";
        float focalLength;
        float pixelSize = -1;
        float fov = 20; // degtorad will be calculated later in the pipeline
        std::string centroidAlgo = "dummy";
        int dummyCentroidNumStars = 5;
        int centroidMagFilter = -1;  // value that should not be used, to tell whether this was selected or not
        std::string idAlgo = "dummy"; 
        float gvTolerance = 0.04;
        float pyTolerance = 0.04;
        int pyFalseStars = 500;
        float pyMismatchProb = 0.001;
        std::string attitudeAlgo = "dqm";
        int generate = 1;
        int horizontalRes = 1024;
        int verticalRes = 1024;
        int referenceBrightness = 8000;
        float brightnessDeviation = 0.7;
        float noiseDeviation = 10;
        float ra = 88; // degtorad will be calculated later in the pipeline
        float dec = 7; // degtorad will be calculated later in the pipeline
        float roll = 0;
        float centroidCompareThreshold;
        float attitudeCompareThreshold;
        std::string plotRawInput = "";
        std::string plotInput = "";
        std::string plotOutput = "";
        std::string printCentroids = "";
        std::string compareCentroids = "";
        std::string compareStars = "";
        std::string printAttitude = "";
        std::string compareAttitude = "";
        std::string database = "";
    };

    // represents the input and expected outputs of a pipeline run.
    class PipelineInput
    {
    public:
        virtual ~PipelineInput(){};
        virtual const Image *InputImage() const { return NULL; };
        virtual const Catalog &GetCatalog() const = 0;
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

    class GeneratedPipelineInput : public PipelineInput
    {
    public:
        // TODO: correct params
        GeneratedPipelineInput(const Catalog &, Quaternion, Camera,
                               int referenceBrightness, float brightnessDeviation,
                               float noiseDeviation);

        const Image *InputImage() const { return &image; };
        const Stars *InputStars() const { return &stars; };
        const Camera *InputCamera() const { return &camera; };
        const StarIdentifiers *InputStarIds() const { return &starIds; };
        bool InputStarsIdentified() const { return true; };
        const Quaternion *InputAttitude() const { return &attitude; };
        const Catalog &GetCatalog() const { return catalog; };

    private:
        // we don't use an Image here because we want to
        std::unique_ptr<unsigned char[]> imageData;
        Image image;
        Stars stars;
        Camera camera;
        Quaternion attitude;
        const Catalog &catalog;
        StarIdentifiers starIds;
    };

    typedef std::vector<std::unique_ptr<PipelineInput>> PipelineInputList;

    PipelineInputList GetPipelineInput(const PipelineOptions &values);

    class PngPipelineInput : public PipelineInput
    {
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

    struct PipelineOutput
    {
        std::unique_ptr<Stars> stars;
        std::unique_ptr<StarIdentifiers> starIds;
        std::unique_ptr<Quaternion> attitude;
        Catalog catalog; // the catalog that the indices in starIds refer to. TODO: don't store it here
    };

    struct StarIdComparison
    {
        int numCorrect;
        int numIncorrect;
        int numTotal;
        float fractionCorrect;
        float fractionIncorrect;
    };

    // actualStars is optional, in which case it's assumed that expectedStars was passed to the star-id
    StarIdComparison StarIdsCompare(const StarIdentifiers &expected, const StarIdentifiers &actual,
                                    const Catalog &expectedCatalog, const Catalog &actualCatalog,
                                    float centroidThreshold,
                                    const Stars *expectedStars, const Stars *actualStars);

    //////////////
    // PIPELINE //
    //////////////

    // a pipeline is a set of algorithms that describes all or part of the star-tracking "pipeline"

    class Pipeline
    {
        friend Pipeline SetPipeline(const PipelineOptions &values);

    public:
        // pointers just so they're nullable
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

    // ask the user what to do with actual and expected outputs
    void PipelineComparison(const PipelineInputList &expected,
                                  const std::vector<PipelineOutput> &actual,
                                  const PipelineOptions &values);

    ////////////////
    // DB BUILDER //
    ////////////////

    Catalog PromptNarrowedCatalog(const Catalog &);

    // TODO: make not public?
    class DatabaseOptions
    {
    public:
        int maxMagnitude = 1000;
        int maxStars = 10000;
        std::string databaseBuilder = "";
        float kvectorMinDistance = 0.5; // DegToRad will be calculated in a later step (GenerateDatabases)
        float kvectorMaxDistance = 15;  // DegToRad will be calculated in a later step (GenerateDatabases)
        long kvectorDistanceBins = 10000;
        std::string path = "stdout";
    };

    // unlike the other algorithm prompters, db builders aren't a
    // typedef void (*DbBuilder)(MultiDatabaseBuilder &, const Catalog &);
    void GenerateDatabases(MultiDatabaseBuilder &, const Catalog &, const DatabaseOptions &values);
    // void PromptDatabases(MultiDatabaseBuilder &, const Catalog &);

}

#endif
