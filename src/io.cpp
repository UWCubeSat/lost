#include "io.hpp"
#include "attitude-utils.hpp"
#include "database-builders.hpp"
#include "star-id.hpp"

#include <cairo/cairo.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>

#include <vector>
#include <sstream>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <memory>
#include <cstring>
#include <random>

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

PromptedOutputStream::PromptedOutputStream() {
    std::string filePath = Prompt<std::string>("Output file (or - for stdout)");
    if (filePath == "-") {
        stream = &std::cout;
        isFstream = false;
    } else {
        std::fstream *fs = new std::fstream();
        fs->open(filePath, std::fstream::out);
        stream = fs;
        isFstream = true;
    }
}

PromptedOutputStream::~PromptedOutputStream() {
    if (isFstream) {
        delete stream;
    }
}

std::vector<CatalogStar> BscParse(std::string tsvPath) {
    std::vector<CatalogStar> result;
    FILE *file;
    double raj2000, dej2000;
    int magnitudeHigh, magnitudeLow, name;
    char weird;

    file = fopen(tsvPath.c_str(), "r");
    if (file == NULL) {
        printf("Error opening file: %s\n", strerror(errno));
        return result; // TODO
    }

    while (EOF != fscanf(file, "%lf|%lf|%d|%c|%d.%d",
                         &raj2000, &dej2000,
                         &name, &weird,
                         &magnitudeHigh, &magnitudeLow)) {
        char strName[8];
        sprintf(strName, "%d", name);
        result.push_back(CatalogStar(DegToRad(raj2000), DegToRad(dej2000),
                                     magnitudeHigh*100 + magnitudeLow,
                                     weird != ' ',
                                     strName));
    }

    fclose(file);
    return result;
}

#ifndef DEFAULT_BSC_PATH
#define DEFAULT_BSC_PATH "bright-star-catalog.tsv"
#endif

std::vector<CatalogStar> CatalogRead() {
    static bool readYet = false;
    static std::vector<CatalogStar> catalog;

    if (!readYet) {
        readYet = true;
        char *tsvPath = getenv("LOST_BSC_PATH");
        catalog = BscParse(tsvPath ? tsvPath : DEFAULT_BSC_PATH);
    }
    return catalog;
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

cairo_surface_t *GrayscaleImageToSurface(const unsigned char *image,
                                         const int width, const int height) {
    // cairo's 8-bit type isn't BW, it's an alpha channel only, which would look a bit weird lol
    cairo_surface_t *result = cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
    uint32_t *resultData = (uint32_t *)cairo_image_surface_get_data(result);
    // hopefully unnecessary
    cairo_surface_flush(result);
    for (long i = 0; i < width * height; i++) {
        // equal r, g, and b components
        resultData[i] = (image[i] << 16) + (image[i] << 8) + (image[i]);
    }
    cairo_surface_mark_dirty(result);
    return result;
}

void SurfacePlotCentroids(cairo_surface_t *cairoSurface,
                          const Stars &centroids,
                          double red,
                          double green,
                          double blue,
                          double alpha) {
    cairo_t *cairoCtx;

    cairoCtx = cairo_create(cairoSurface);
    cairo_set_source_rgba(cairoCtx, red, green, blue, alpha);
    cairo_set_line_width(cairoCtx, 1.0);
    cairo_set_antialias(cairoCtx, CAIRO_ANTIALIAS_NONE);

    for (const Star &centroid : centroids) {
        if (centroid.radiusX > 0.0f) {
            float radiusX = centroid.radiusX;
            float radiusY = centroid.radiusY > 0.0f ?
                centroid.radiusY : radiusX;

            // Rectangles should be entirely /outside/ the radius of the star, so the star is fully
            // visible.
            cairo_rectangle(cairoCtx,
                            centroid.x+.5 - radiusX,
                            centroid.y+.5 - radiusY,
                            radiusX * 2,
                            radiusY * 2);
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

typedef CentroidAlgorithm *(*CentroidAlgorithmFactory)();

CentroidAlgorithm *DummyCentroidAlgorithmPrompt() {
    int numStars = Prompt<int>("How many stars to generate");
    return new DummyCentroidAlgorithm(numStars);
}

CentroidAlgorithm *CogCentroidAlgorithmPrompt() {
    return new CenterOfGravityAlgorithm();
}

typedef StarIdAlgorithm *(*StarIdAlgorithmFactory)();

StarIdAlgorithm *DummyStarIdAlgorithmPrompt() {
    return new DummyStarIdAlgorithm();
}

StarIdAlgorithm *GeometricVotingStarIdAlgorithmPrompt() {
    return new GeometricVotingStarIdAlgorithm();
}

StarIdAlgorithm *PyramidStarIdAlgorithmPrompt() {
    return new PyramidStarIdAlgorithm();
}

unsigned char *PromptKVectorDatabaseBuilder(const Catalog &catalog, long *length) {
    float minDistance = DegToRad(Prompt<float>("Min distance (deg)"));
    float maxDistance = DegToRad(Prompt<float>("Max distance (deg)"));
    long numBins = Prompt<long>("Number of distance bins");
    return BuildKVectorDatabase(catalog, length, minDistance, maxDistance, numBins);
}

DbBuilder PromptDbBuilder() {
    InteractiveChoice<DbBuilder> dbBuilderChoice;
    dbBuilderChoice.Register("kvector", "K-Vector (geometric voting & pyramid)", PromptKVectorDatabaseBuilder);
    return dbBuilderChoice.Prompt("Choose database builder");
}

// PIPELINE INPUT STUFF

cairo_surface_t *PipelineInput::InputImageSurface() const {
    const Image *inputImage = InputImage();
    return GrayscaleImageToSurface(inputImage->image, inputImage->width, inputImage->height);
}

class PngPipelineInput : public PipelineInput {
public:
    PngPipelineInput(cairo_surface_t *);

    const Image *InputImage() const { return &image; };
private:
    Image image;
};

class AstrometryPipelineInput : public PipelineInput {
public:
    AstrometryPipelineInput(const std::string &path);

    const Image *InputImage() const { return &image; };
    const Quaternion *InputAttitude() const { return &attitude; };
private:
    Image image;
    Quaternion attitude;
};

class GeneratedPipelineInput : public PipelineInput {
public:
    // TODO: correct params
    GeneratedPipelineInput(const std::vector<CatalogStar> &, Quaternion, Camera,
                           unsigned char noise_deviation);

    const Image *InputImage() const { return &image; };
    const Stars *InputCentroids() const { return &stars; };
    const Stars *InputStars() const { return &stars; };
    const Quaternion *InputAttitude() const { return &attitude; };
private:
    // we don't use an Image here because we want to 
    std::unique_ptr<unsigned char[]> imageData;
    Image image;
    Stars stars;
    Quaternion attitude;
};

PngPipelineInput::PngPipelineInput(cairo_surface_t *cairoSurface) {
    image.image = SurfaceToGrayscaleImage(cairoSurface);
    image.width = cairo_image_surface_get_width(cairoSurface);
    image.height = cairo_image_surface_get_height(cairoSurface);
}

PipelineInputList PromptPngPipelineInput() {
    // I'm not sure why, but i can't get an initializer list to work here. Probably something to do
    // with copying unique ptrs
    PipelineInputList result;
    cairo_surface_t *cairoSurface = NULL;
    std::string pngPath;

    while (cairoSurface == NULL || cairo_surface_status(cairoSurface) != CAIRO_STATUS_SUCCESS) {
        cairoSurface = cairo_image_surface_create_from_png(Prompt<std::string>("PNG Path").c_str());
    }
    result.push_back(std::unique_ptr<PipelineInput>(new PngPipelineInput(cairoSurface)));
    return result;
}

AstrometryPipelineInput::AstrometryPipelineInput(const std::string &path) {
    // create from path, TODO
}

PipelineInputList PromptAstrometryPipelineInput() {
    // TODO: why does it let us do a reference to an ephemeral return value?
    std::string path = Prompt<std::string>("Astrometry download directory");
    PipelineInputList result;
    result.push_back(std::unique_ptr<PipelineInput>(new AstrometryPipelineInput(path)));
    return result;
}

GeneratedPipelineInput::GeneratedPipelineInput(const Catalog &catalog,
                                               Quaternion attitude,
                                               Camera camera,
                                               unsigned char noise_deviation) {
    this->attitude = attitude;
    image.width = camera.xResolution;
    image.height = camera.yResolution;
    unsigned char *imageRaw = (unsigned char *)calloc(image.width * image.height, 1);
    imageData = std::unique_ptr<unsigned char[]>(imageRaw);
    image.image = imageData.get();

    for (const CatalogStar &catalogStar : catalog) {
        Vec3 spatial = {
            cos(catalogStar.raj2000)*cos(catalogStar.dej2000),
            sin(catalogStar.raj2000)*cos(catalogStar.dej2000),
            sin(catalogStar.dej2000),
        };
        Vec3 rotated = attitude.Rotate(spatial);
        if (rotated.x < 0) {
            continue;
        }
        Vec2 camCoords = camera.ConvertCoordinates(rotated);
        if (camera.InSensor(camCoords)) {
            stars.push_back(Star(camCoords.x, camCoords.y,
                                 (750-catalogStar.magnitude)/200.0,
                                 (750-catalogStar.magnitude)/200.0,
                                 catalogStar.magnitude));
        }
    }

    for (const Star &star : stars) {
        for (int k = star.y - star.radiusX; k < star.y + star.radiusX; k++) {
            for(int j = star.x - star.radiusX; j < star.x + star.radiusX; j++) {
                if (k*image.width + j < image.width * image.height && k >= 0 && j >= 0) {
                    imageData[k*image.width + j] = rand() % 128 + 128;
                }
            }
        }
    }

    // generate gaussina noise w a normal distribution
    std::normal_distribution<float> dist(0.0, noise_deviation);
    std::default_random_engine generator;
    for (int i = 0; i < image.width * image.height; i++) {
        int noise = int(dist(generator));
        int temp = imageData[i] + noise;
        if (temp < 0) {
            imageData[i] = 0;
        } else if (temp > 255) {
            imageData[i] = 255;
        } else {
            imageData[i] = temp;
        }
    }  
}

Quaternion PromptSphericalAttitude() {
    float ra = Prompt<float>("Boresight right ascension");
    float dec = Prompt<float>("Boresight declination");
    float roll = Prompt<float>("Boresight roll");
    return SphericalToQuaternion(DegToRad(ra), DegToRad(dec), DegToRad(roll));
}

PipelineInputList PromptGeneratedPipelineInput() {
    // TODO: prompt for attitude, imagewidth, etc and then construct a GeneratedPipelineInput
    int numImages = Prompt<int>("Number of images to generate");
    int xResolution = Prompt<int>("Horizontal Resolution");
    int yResolution = Prompt<int>("Vertical Resolution");
    float xFovDeg = Prompt<float>("Horizontal FOV (in degrees)");
    // TODO: allow random angle generation?
    Quaternion attitude = PromptSphericalAttitude();

    PipelineInputList result;

    for (int i = 0; i < numImages; i++) {
        GeneratedPipelineInput *curr = new GeneratedPipelineInput(
            CatalogRead(),
            attitude,
            Camera(DegToRad(xFovDeg), xResolution, yResolution),
            3);

        result.push_back(std::unique_ptr<PipelineInput>(curr));
    }

    return result;
}

typedef PipelineInputList (*PipelineInputFactory)();

PipelineInputList PromptPipelineInput() {
    InteractiveChoice<PipelineInputFactory> inputTypeChoice;
    inputTypeChoice.Register("png", "PNG files", PromptPngPipelineInput);
    inputTypeChoice.Register("generate", "Generated image", PromptGeneratedPipelineInput);
    inputTypeChoice.Register("astrometry", "Astrometry.net", PromptAstrometryPipelineInput);

    return (inputTypeChoice.Prompt("Input from"))();
}

Pipeline PromptPipeline() {
    enum class PipelineStage {
        Centroid, Database, StarId, AttitudeEstimation, Done
    };

    Pipeline result;
    
    InteractiveChoice<PipelineStage> stageChoice;
    stageChoice.Register("centroid", "Centroid", PipelineStage::Centroid);
    stageChoice.Register("database", "Database", PipelineStage::Database);
    stageChoice.Register("starid", "Star-ID", PipelineStage::StarId);
    stageChoice.Register("attitude", "Attitude Estimation", PipelineStage::AttitudeEstimation);
    stageChoice.Register("done", "Done setting up pipeline", PipelineStage::Done);

    while (true) {
        PipelineStage nextStage = stageChoice.Prompt("Which pipeline stage to set");
        switch (nextStage) {

        case PipelineStage::Centroid: {
            InteractiveChoice<CentroidAlgorithmFactory> centroidChoice;
            centroidChoice.Register("dummy", "Random Centroid Algorithm", DummyCentroidAlgorithmPrompt);
            centroidChoice.Register("cog", "Center of Gravity Centroid Algorithm", CogCentroidAlgorithmPrompt);

            result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(
                (centroidChoice.Prompt("Choose centroid algo"))());
            break;
        }

        case PipelineStage::Database: {
            std::string path = Prompt<std::string>("Database file");
            std::fstream fs;
            fs.open(path, std::fstream::in | std::fstream::binary);
            fs.seekg(0, fs.end);
            long length = fs.tellg();
            fs.seekg(0, fs.beg);
            std::cerr << "Reading " << length << " bytes of database" << std::endl;
            result.database = std::unique_ptr<unsigned char[]>(new unsigned char[length]);
            fs.read((char *)result.database.get(), length);
            std::cerr << "Done" << std::endl;
        }

        case PipelineStage::StarId: {
            InteractiveChoice<StarIdAlgorithmFactory> starIdChoice;
            starIdChoice.Register("dummy", "Random", DummyStarIdAlgorithmPrompt);
            starIdChoice.Register("gv", "Geometric Voting", GeometricVotingStarIdAlgorithmPrompt);
            starIdChoice.Register("pyramid", "Pyramid codes", PyramidStarIdAlgorithmPrompt);

            result.starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(
                (starIdChoice.Prompt("Choose Star-ID algo"))());
            break;
        }

        case PipelineStage::AttitudeEstimation: {
            std::cerr << "TODO" << std::endl;
            break;
        }

        case PipelineStage::Done: {
            // fuck style guides
            goto PipelineDone;
        }
        }
    }
    PipelineDone:

    return result;
}

PipelineOutput Pipeline::Go(const PipelineInput &input) {
    // Start executing the pipeline at the first stage that has both input and an algorithm. From
    // there, execute each successive stage of the pipeline using the output of the last stage
    // (human centipede) until there are no more stages set.
    PipelineOutput result = { 0 };

    const Image *inputImage = input.InputImage();
    const Stars *inputCentroids = input.InputCentroids();
    const Stars *inputStars = input.InputStars();

    if (centroidAlgorithm && inputImage) {
        // TODO: we should probably modify Go to just take an image argument
        // TODO: don't copy the vector!
        result.centroids = std::unique_ptr<Stars>(new std::vector<Star>(
            centroidAlgorithm->Go(inputImage->image, inputImage->width, inputImage->height)));
        inputCentroids = result.centroids.get();
    }

    if (starIdAlgorithm && database && inputCentroids) {
        // TODO: don't copy the vector!
        result.stars = std::unique_ptr<Stars>(new std::vector<Star>(*result.centroids));
        starIdAlgorithm->Go(database.get(), result.stars.get());
        inputStars = result.stars.get();
    }

    if (attitudeEstimationAlgorithm && inputStars && input.InputCamera()) {
        result.attitude = std::unique_ptr<Quaternion>(
            new Quaternion(attitudeEstimationAlgorithm->Go(*input.InputCamera(), *inputStars)));
    }

    return result;
}

std::vector<PipelineOutput> Pipeline::Go(const PipelineInputList &inputs) {
    std::vector<PipelineOutput> result;
    
    for (const std::unique_ptr<PipelineInput> &input : inputs) {
        result.push_back(Go(*input));
    }

    return result;
}

////////////////
// COMPARISON //
////////////////

class CentroidComparison {
public:
    CentroidComparison() : meanError(0.0f), numExtraStars(0.0), numMissingStars(0.0) { };
    float meanError;       // average distance from actual to expected star
    // both these are floats because we may average multiple centroid comparisons together:
    float numExtraStars;    // stars in actual but not expected. Ideally 0
    float numMissingStars;  // stars is expected but not actual. Ideally 0
    // I would add 99th percentile or something, but the really far away stars should really just
    // count in extra_num
};

// helper for StarCentroidsCompare
static std::vector<int> FindClosestCentroids(float threshold,
                                              const Stars &one,
                                              const Stars &two) {
    std::vector<int> result;

    for (int i = 0; i < (int)one.size(); i++) {
        float closestDistance = INFINITY;
        int   closestIndex = -1;

        for (int k = 0; k < (int)two.size(); k++) {
            float currDistance = StarDistancePixels(one[i], two[k]);
            if (currDistance < threshold && currDistance < closestDistance) {
                closestDistance = currDistance;
                closestIndex = k;
            }
        }
        
        result.push_back(closestIndex);
    }

    return result;
}

CentroidComparison CentroidsCompare(float threshold,
                                    const Stars &expected,
                                    const Stars &actual) {

    CentroidComparison result;
    // maps from indexes in each list to the closest centroid from other list
    std::vector<int> expectedToActual, actualToExpected;

    expectedToActual = FindClosestCentroids(threshold, expected, actual);
    actualToExpected = FindClosestCentroids(threshold, actual, expected);

    // any expected stars whose closest actual star does not refer to them are missing
    for (int i = 0; i < (int)expectedToActual.size(); i++) {
        if (expectedToActual[i] == -1 ||
            actualToExpected[expectedToActual[i]] != i) {
            result.numMissingStars++;
        } else {
            result.meanError += StarDistancePixels(expected[i], actual[expectedToActual[i]]);
        }
    }
    result.meanError /= (expected.size() - result.numMissingStars);

    // any actual star whose closest expected star does not refer to them is extra
    for (int i = 0; i < (int)actual.size(); i++) {
        if (actualToExpected[i] == -1 ||
            expectedToActual[actualToExpected[i]] != i) {
            result.numExtraStars++;
        }
    }

    return result;
}

CentroidComparison CentroidComparisonsCombine(std::vector<CentroidComparison> comparisons) {
    assert(comparisons.size() > 0);

    CentroidComparison result;
    
    for (const CentroidComparison &comparison : comparisons) {
        result.meanError += comparison.meanError;
        result.numExtraStars += comparison.numExtraStars;
        result.numMissingStars += comparison.numMissingStars;
    }

    result.meanError /= comparisons.size();
    result.numExtraStars /= comparisons.size();
    result.numMissingStars /= comparisons.size();
    
    return result;
}

/////////////////////
// PIPELINE OUTPUT //
/////////////////////

typedef void (*PipelineComparator)(std::ostream &os,
                                   const PipelineInputList &,
                                   const std::vector<PipelineOutput> &);

cairo_status_t OstreamPlotter(void *closure, const unsigned char *data, unsigned int length) {
    std::ostream *os = (std::ostream *)closure;
    os->write((const char *)data, length);
    return CAIRO_STATUS_SUCCESS;
}

void PipelineComparatorPlotInput(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &actual) {
    
    cairo_surface_write_to_png_stream(expected[0]->InputImageSurface(), OstreamPlotter, &os);
}

void PipelineComparatorCentroids(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &actual) {
    int size = (int)expected.size();

    float threshold = Prompt<float>("Threshold to count as the same star (pixels, float)");

    std::vector<CentroidComparison> comparisons;
    for (int i = 0; i < size; i++) {
        comparisons.push_back(CentroidsCompare(threshold,
                                               *(expected[i]->ExpectedCentroids()),
                                               *(actual[i].centroids)));
    }

    CentroidComparison result = CentroidComparisonsCombine(comparisons);
    os << "extra_stars " << result.numExtraStars << std::endl
       << "missing_stars " << result.numMissingStars << std::endl
       << "mean_error " << result.meanError << std::endl;
}

void PipelineComparatorPlotCentroids(std::ostream &os,
                                     const PipelineInputList &expected,
                                     const std::vector<PipelineOutput> &actual) {

    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    SurfacePlotCentroids(cairoSurface, *actual[0].centroids, 1.0, 0.0, 0.0, 0.5);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
}

void PipelineComparatorStars(std::ostream &os,
                                     const PipelineInputList &expected,
                                     const std::vector<PipelineOutput> &actual) {
    os << "Unimplemented" << std::endl;
    // TODO: plot all stars, maybe a different color for unidentified stars, and put text next to
    // properly identified stars.
}

void PipelineComparatorPlotStars(std::ostream &os,
                                     const PipelineInputList &expected,
                                     const std::vector<PipelineOutput> &actual) {
    os << "Unimplemented" << std::endl;
}

void PromptPipelineComparison(const PipelineInputList &expected,
                              const std::vector<PipelineOutput> &actual) {
    assert(expected.size() == actual.size() && expected.size() > 0);

    InteractiveChoice<PipelineComparator> comparatorChoice;

    if (expected[0]->InputImage() && expected.size() == 1) {
        comparatorChoice.Register("plot_input", "Plot raw BW input image to PNG",
                                  PipelineComparatorPlotInput);
    }

    // Centroids
    if (actual[0].centroids != NULL) {
        if (actual.size() == 1) {
            comparatorChoice.Register("plot_centroids", "Plot centroids to PNG",
                                      PipelineComparatorPlotCentroids);
        }
        if (expected[0]->ExpectedCentroids()) {
            comparatorChoice.Register("compare_centroids", "Compare lists of centroids",
                                      PipelineComparatorCentroids);
        }
    }

    // Stars
    if (actual[0].stars != NULL) {
        if (actual.size() == 1) {
            comparatorChoice.Register("plot_stars", "Plot identified stars to PNG",
                                      PipelineComparatorPlotStars);
        }
        if (expected[0]->ExpectedStars()) {
            comparatorChoice.Register("compare_stars", "Compare lists of identified stars",
                                      PipelineComparatorStars);
        }
    }

    // TODO: Attitude

    comparatorChoice.Register("done", "No more comparisons", NULL);

    while (true) {
        PipelineComparator comparator = comparatorChoice.Prompt("What to do with output");
        if (comparator == NULL) {
            break;
        }

        PromptedOutputStream pos;
        comparator(pos.Stream(), expected, actual);
    }
}

}
