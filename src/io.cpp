#include "io.hpp"

#include "attitude-estimators.hpp"
#include "attitude-utils.hpp"
#include "databases.hpp"
#include "star-id.hpp"
#include "star-utils.hpp"

#include <cairo/cairo.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <cstring>
#include <random>
#include <algorithm>

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

PromptedOutputStream::PromptedOutputStream(std::string filePath) {
    //std::string filePath = Prompt<std::string>("Output file (or - for stdout)");
    if (filePath == "stdout") {
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
        if (weird == ' ') {
            result.push_back(CatalogStar(DegToRad(raj2000),
                                         DegToRad(dej2000),
                                         magnitudeHigh*100 + (magnitudeHigh < 0 ? -magnitudeLow : magnitudeLow),
                                         name));
        }
    }

    fclose(file);
    return result;
}

#ifndef DEFAULT_BSC_PATH
#define DEFAULT_BSC_PATH "bright-star-catalog.tsv"
#endif

std::vector<CatalogStar> &CatalogRead() {
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

void SurfacePlot(cairo_surface_t *cairoSurface,
                 const Stars &stars,
                 const StarIdentifiers *starIds,
                 const Catalog *catalog,
                 const Quaternion *attitude,
                 double red,
                 double green,
                 double blue,
                 double alpha) {
    cairo_t *cairoCtx;
    std::string metadata = "";

    cairoCtx = cairo_create(cairoSurface);
    cairo_set_source_rgba(cairoCtx, red, green, blue, alpha);
    cairo_set_line_width(cairoCtx, 1.0);
    cairo_set_antialias(cairoCtx, CAIRO_ANTIALIAS_NONE);
    cairo_font_options_t *cairoFontOptions = cairo_font_options_create();
    cairo_font_options_set_antialias(cairoFontOptions, CAIRO_ANTIALIAS_NONE);
    cairo_set_font_options(cairoCtx, cairoFontOptions);
    cairo_text_extents_t cairoTextExtents;
    cairo_text_extents(cairoCtx, "1234567890", &cairoTextExtents);
    double textHeight = cairoTextExtents.height;

    for (const Star &centroid : stars) {
        // plot the box around the star
        if (centroid.radiusX > 0.0f) {
            float radiusX = centroid.radiusX;
            float radiusY = centroid.radiusY > 0.0f ?
                centroid.radiusY : radiusX;

            // Rectangles should be entirely /outside/ the radius of the star, so the star is
            // fully visible.
            cairo_rectangle(cairoCtx,
                            centroid.position.x - radiusX,
                            centroid.position.y - radiusY,
                            radiusX * 2,
                            radiusY * 2);
            cairo_stroke(cairoCtx);
        } else {
            cairo_rectangle(cairoCtx,
                            floor(centroid.position.x),
                            floor(centroid.position.y),
                            1, 1);
            cairo_fill(cairoCtx);
        }
    }

    metadata += std::to_string(stars.size()) + " centroids   ";

    if (starIds != NULL) {
        assert(catalog != NULL);

        for (const StarIdentifier &starId : *starIds) {
            const Star &centroid = stars[starId.starIndex];
            cairo_move_to(cairoCtx,

                          centroid.radiusX > 0.0f
                          ? centroid.position.x + centroid.radiusX + 3
                          : centroid.position.x + 8,

                          centroid.radiusY > 0.0f
                          ? centroid.position.y - centroid.radiusY + textHeight
                          : centroid.position.y + 10);

            cairo_show_text(cairoCtx, std::to_string((*catalog)[starId.catalogIndex].name).c_str());
        }
        metadata += std::to_string(starIds->size()) + " identified   ";
    }

    if (attitude != NULL) {
        float ra, de, roll;
        attitude->ToSpherical(&ra, &de, &roll);
        metadata +=
            "RA: " + std::to_string(RadToDeg(ra)) + "  " +
            "DE: " + std::to_string(RadToDeg(de)) + "  " +
            "Roll: " + std::to_string(RadToDeg(roll)) + "   ";
    }

    // plot metadata
    cairo_move_to(cairoCtx, 3, 3 + textHeight);
    cairo_show_text(cairoCtx, metadata.c_str());

    cairo_font_options_destroy(cairoFontOptions);
    cairo_destroy(cairoCtx);
}

// ALGORITHM PROMPTERS
typedef CentroidAlgorithm *(*CentroidAlgorithmFactory)();

typedef StarIdAlgorithm *(*StarIdAlgorithmFactory)();

typedef AttitudeEstimationAlgorithm *(*AttitudeEstimationAlgorithmFactory)();

void PromptKVectorDatabaseBuilder(MultiDatabaseBuilder &builder, const Catalog &catalog, float minDistance, float maxDistance, long numBins) {
    // TODO: calculating the length of the vector duplicates a lot of the work, slowing down
    // database generation
    long length = SerializeLengthPairDistanceKVector(catalog, minDistance, maxDistance, numBins);
    unsigned char *buffer = builder.AddSubDatabase(PairDistanceKVectorDatabase::kMagicValue, length);
    if (buffer == NULL) {
        std::cerr << "No room for another database." << std::endl;
    }
    SerializePairDistanceKVector(catalog, minDistance, maxDistance, numBins, buffer);

    // TODO: also parse it and print out some stats before returning
    
}


void GenerateDatabases(MultiDatabaseBuilder &builder, const Catalog &catalog, const DatabaseOptions &values) {

    if (values.databaseBuilder == "kvector") {
        float minDistance = DegToRad(values.kvectorMinDistance);
        float maxDistance = DegToRad(values.kvectorMaxDistance);
        long numBins = values.kvectorDistanceBins;
        PromptKVectorDatabaseBuilder(builder, catalog, minDistance, maxDistance, numBins);
    } else {
        std::cerr << "No database builder selected -- no database generated." << std::endl;
        exit(1);
    }

}

// PIPELINE INPUT STUFF

cairo_surface_t *PipelineInput::InputImageSurface() const {
    const Image *inputImage = InputImage();
    return GrayscaleImageToSurface(inputImage->image, inputImage->width, inputImage->height);
}

class AstrometryPipelineInput : public PipelineInput {
public:
    AstrometryPipelineInput(const std::string &path);

    const Image *InputImage() const { return &image; };
    const Quaternion *InputAttitude() const { return &attitude; };
private:
    Image image;
    Quaternion attitude;
};

PngPipelineInput::PngPipelineInput(cairo_surface_t *cairoSurface, Camera camera, const Catalog &catalog)
    : camera(camera), catalog(catalog) {

    image.image = SurfaceToGrayscaleImage(cairoSurface);
    image.width = cairo_image_surface_get_width(cairoSurface);
    image.height = cairo_image_surface_get_height(cairoSurface);
}

PipelineInputList GetPngPipelineInput(const PipelineOptions &values) {
    // I'm not sure why, but i can't get an initializer list to work here. Probably something to do
    // with copying unique ptrs
    PipelineInputList result;
    cairo_surface_t *cairoSurface = NULL;
    std::string pngPath = values.png;
    
    while (cairoSurface == NULL || cairo_surface_status(cairoSurface) != CAIRO_STATUS_SUCCESS) {
         cairoSurface = cairo_image_surface_create_from_png(pngPath.c_str());
    }
    int xResolution = cairo_image_surface_get_width(cairoSurface);
    int yResolution = cairo_image_surface_get_height(cairoSurface);

    if (!values.focalLength || !values.pixelSize) {
        std::cerr << "Error: No focal length/pixel size given." << std::endl;
        exit(1);
    } 

    float focalLength = values.focalLength;
    float pixelSize = values.pixelSize;
    float focalLengthPixels = focalLength * 1000 / pixelSize;
    Camera cam = Camera(focalLengthPixels, xResolution, yResolution);

    result.push_back(std::unique_ptr<PipelineInput>(new PngPipelineInput(cairoSurface, cam, CatalogRead())));
    return result;
}

AstrometryPipelineInput::AstrometryPipelineInput(const std::string &path) {
    // create from path, TODO
}

// does NOT protect against multiple evaluation of arguments
#define IncrementPixelXY(x, y, amt) imageData[(y)*image.width+(x)] = \
        std::max(0, std::min(255,imageData[(y)*image.width+(x)]+(amt)))
#define IncrementPixelI(i, amt) imageData[i] = std::max(0, std::min(255,imageData[i]+(amt)))

GeneratedPipelineInput::GeneratedPipelineInput(const Catalog &catalog,
                                               Quaternion attitude,
                                               Camera camera,
                                               int referenceBrightness,
                                               float brightnessDeviation,
                                               float noiseDeviation)
    : camera(camera), attitude(attitude), catalog(catalog) {

    image.width = camera.GetXResolution();
    image.height = camera.GetYResolution();
    unsigned char *imageRaw = (unsigned char *)calloc(image.width * image.height, 1);
    imageData = std::unique_ptr<unsigned char[]>(imageRaw);
    image.image = imageData.get();

    for (int i = 0; i < (int)catalog.size(); i++) {
        const CatalogStar &catalogStar = catalog[i];
        Vec3 rotated = attitude.Rotate(catalog[i].spatial);
        if (rotated.x < 0) {
            continue;
        }
        Vec2 camCoords = camera.SpatialToCamera(rotated);

        float radiusX = ceil(brightnessDeviation*2);
        if (camera.InSensor(camCoords)) {
            stars.push_back(Star(camCoords.x, camCoords.y, radiusX, radiusX, catalogStar.magnitude));
            starIds.push_back(StarIdentifier(stars.size() - 1, i));
        }
    }

    for (const Star &star : stars) {
        // "brightness" = number of photons received, for eg
        int totalBrightness = referenceBrightness * pow(100.0f, -star.magnitude/500.0f);

        // the star.x and star.y refer to the pixel whose top left corner the star should appear at
        // (and fractional amounts are relative to the corner). When we color a pixel, we ideally
        // would integrate the intensity of the star over that pixel, but we can make do by sampling
        // the intensity of the star at the /center/ of the pixel, ie, star.x+.5 and star.y+.5
        for(int j = star.position.x - star.radiusX; j >= 0 && j < star.position.x + star.radiusX && j < image.width; j++) {
            for (int k = star.position.y - star.radiusX; k >= 0 && k < star.position.y + star.radiusX && k < image.height; k++) {
                float distanceSquared = pow(k+.5-star.position.y, 2) + pow(j+.5-star.position.x, 2);
                int pixelBrightness = totalBrightness / pow(brightnessDeviation, 2)
                    * exp(-distanceSquared/pow(brightnessDeviation, 2));
                IncrementPixelXY(j, k, pixelBrightness);
            }
        }
    }

    std::normal_distribution<float> readNoiseDist(0.0, noiseDeviation);
    std::default_random_engine generator;
    for (int i = 0; i < image.width * image.height; i++) {
        // dark current

        // shot noise
        std::poisson_distribution<char> shotNoiseDist(imageData[i]);
        // TODO: prompt for sensitivity so we can accurately apply shot noise.

        // read noise
        IncrementPixelI(i, (int)readNoiseDist(generator));
    }  
}

PipelineInputList GetGeneratedPipelineInput(const PipelineOptions &values) {
    // TODO: prompt for attitude, imagewidth, etc and then construct a GeneratedPipelineInput
    int numImages = values.generate;
    int xResolution = values.horizontalRes;
    int yResolution = values.verticalRes;
    float xFovDeg = values.fov;
    float xFocalLength = FovToFocalLength(DegToRad(xFovDeg), xResolution);
    int referenceBrightness = values.referenceBrightness;
    float brightnessDeviation = values.brightnessDeviation;
    float noiseDeviation = values.noiseDeviation;
    float ra = values.ra;
    float dec = values.dec;
    float roll = values.roll;

    // TODO: allow random angle generation?
    Quaternion attitude = SphericalToQuaternion(DegToRad(ra), DegToRad(dec), DegToRad(roll));

    PipelineInputList result;

    for (int i = 0; i < numImages; i++) {
        GeneratedPipelineInput *curr = new GeneratedPipelineInput(
            CatalogRead(),
            attitude,
            Camera(xFocalLength, xResolution, yResolution),
            referenceBrightness, brightnessDeviation, noiseDeviation);

        result.push_back(std::unique_ptr<PipelineInput>(curr));
    }

    return result;
}

typedef PipelineInputList (*PipelineInputFactory)();

PipelineInputList GetPipelineInput(const PipelineOptions &values) {

    if (values.png != "") {
        return GetPngPipelineInput(values);
    } else {
        return GetGeneratedPipelineInput(values);
    }
}

Pipeline::Pipeline(CentroidAlgorithm *centroidAlgorithm,
                   StarIdAlgorithm *starIdAlgorithm,
                   AttitudeEstimationAlgorithm *attitudeEstimationAlgorithm,
                   unsigned char *database)
    : Pipeline() {
    if (centroidAlgorithm) {
        this->centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(centroidAlgorithm);
    }
    if (starIdAlgorithm) {
        this->starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(starIdAlgorithm);
    }
    if (attitudeEstimationAlgorithm) {
        this->attitudeEstimationAlgorithm = std::unique_ptr<AttitudeEstimationAlgorithm>(attitudeEstimationAlgorithm);
    }
    if (database) {
        this->database = std::unique_ptr<unsigned char[]>(database);
    }
}



Pipeline SetPipeline(const PipelineOptions &values) {
    enum class PipelineStage {
        Centroid, CentroidMagnitudeFilter, Database, StarId, AttitudeEstimation, Done
    };

    Pipeline result;
    
    // TODO: more flexible or sth
    // TODO: don't allow setting star-id until database is set, and perhaps limit the star-id
    // choices to those compatible with the database?

    // centroid algorithm stage
    if (values.centroidAlgo == "dummy") {
        result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(new DummyCentroidAlgorithm(values.dummyCentroidNumStars));
    } else if (values.centroidAlgo == "cog") {
        result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(new CenterOfGravityAlgorithm());
    } else if (values.centroidAlgo == "iwcog") {
        result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(new IterativeWeightedCenterOfGravityAlgorithm());
    }

    // centroid magnitude filter stage
    if (values.centroidMagFilter != -1) result.centroidMinMagnitude = values.centroidMagFilter;

    // database stage
    if (values.png != "") {
        std::fstream fs;
        fs.open(values.png, std::fstream::in | std::fstream::binary);
        fs.seekg(0, fs.end);
        long length = fs.tellg();
        fs.seekg(0, fs.beg);
        std::cerr << "Reading " << length << " bytes of database" << std::endl;
        result.database = std::unique_ptr<unsigned char[]>(new unsigned char[length]);
        fs.read((char *)result.database.get(), length);
        std::cerr << "Done" << std::endl;
    } 

    if (values.idAlgo == "dummy") {
        result.starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(new DummyStarIdAlgorithm());
    } else if (values.idAlgo == "gv") {
        result.starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(new GeometricVotingStarIdAlgorithm(DegToRad(values.gvTolerance)));
    } else if (values.idAlgo == "py") {
        result.starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(new PyramidStarIdAlgorithm(DegToRad(values.pyTolerance), values.pyFalseStars, values.pyMismatchProb, 1000));
    }

    if (values.attitudeAlgo == "dqm") {
        result.attitudeEstimationAlgorithm = std::unique_ptr<AttitudeEstimationAlgorithm>(new DavenportQAlgorithm());
    }
    return result;
}

PipelineOutput Pipeline::Go(const PipelineInput &input) {
    // Start executing the pipeline at the first stage that has both input and an algorithm. From
    // there, execute each successive stage of the pipeline using the output of the last stage
    // (human centipede) until there are no more stages set.
    PipelineOutput result = { 0 };

    const Image *inputImage = input.InputImage();
    const Stars *inputStars = input.InputStars();
    const StarIdentifiers *inputStarIds = input.InputStarIds();

    // if database is provided, that's where we get catalog from.
    if (database) {
        MultiDatabase multiDatabase(database.get());
        const unsigned char *catalogBuffer = multiDatabase.SubDatabasePointer(kCatalogMagicValue);
        if (catalogBuffer != NULL) {
            result.catalog = DeserializeCatalog(multiDatabase.SubDatabasePointer(kCatalogMagicValue), NULL, NULL);
        } else {
            std::cerr << "Warning: That database does not include a catalog." << std::endl;
            result.catalog = input.GetCatalog();
        }
    } else {
        result.catalog = input.GetCatalog();
    }

    if (centroidAlgorithm && inputImage) {
        // TODO: we should probably modify Go to just take an image argument
        Stars unfilteredStars = centroidAlgorithm->Go(inputImage->image, inputImage->width, inputImage->height);
        Stars *filteredStars = new std::vector<Star>();
        for (const Star &star : unfilteredStars) {
            if (star.magnitude >= centroidMinMagnitude) {
                filteredStars->push_back(star);
            }
        }
        result.stars = std::unique_ptr<Stars>(filteredStars);
        inputStars = filteredStars;

        // any starid set up to this point needs to be discarded, because it's based on input
        // centroids instead of our new centroids.
        inputStarIds = NULL;
        result.starIds = NULL;
    }

    if (starIdAlgorithm && database && inputStars && input.InputCamera()) {
        // TODO: don't copy the vector!
        result.starIds = std::unique_ptr<StarIdentifiers>(new std::vector<StarIdentifier>(
            starIdAlgorithm->Go(database.get(), *inputStars, result.catalog, *input.InputCamera())));
        inputStarIds = result.starIds.get();
    }

    if (attitudeEstimationAlgorithm && inputStarIds && input.InputCamera()) {
        assert(inputStars); // ensure that starIds doesn't exist without stars
        result.attitude = std::unique_ptr<Quaternion>(
            new Quaternion(attitudeEstimationAlgorithm->Go(*input.InputCamera(), *inputStars, result.catalog, *inputStarIds)));
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
            float currDistance = Distance(one[i].position, two[k].position);
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
            result.meanError += Distance(expected[i].position, actual[expectedToActual[i]].position);
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

bool StarIdCompare(const StarIdentifier &a, const StarIdentifier &b) {
    return a.starIndex < b.starIndex;
}

StarIdComparison StarIdsCompare(const StarIdentifiers &expected, const StarIdentifiers &actual,
                                // use these to map indices to names for the respective lists of StarIdentifiers
                                const Catalog &expectedCatalog, const Catalog &actualCatalog,
                                // stars are ignored if threshold<0
                                float centroidThreshold,
                                const Stars *expectedStars, const Stars *actualStars) {

    StarIdComparison result = {
        0, // correct
        0, // incorrect
        (int)expected.size() // total
    };
    StarIdentifiers expectedSorted = expected;
    StarIdentifiers actualSorted;
    if (expectedStars && actualStars) {
        // map from actual centroid index -> expected centroid index
        std::vector<int> closestCentroids =
            FindClosestCentroids(centroidThreshold, *actualStars, *expectedStars);
        // ensure that even if two actual centroids are close to the same expected centroid, we
        // don't double count.
        std::vector<bool> expectedCentroidsUsed(expected.size());

        for (const StarIdentifier &starId : actual) {
            int closestIndex = closestCentroids[starId.starIndex];

            if (closestIndex == -1) {
                result.numIncorrect++;
            } else if (!expectedCentroidsUsed[closestIndex]) {
                expectedCentroidsUsed[closestIndex] = true;
                actualSorted.push_back(
                    StarIdentifier(closestIndex, starId.catalogIndex, starId.weight));
            }
        }
    } else {
        actualSorted = actual;
    }

    sort(expectedSorted.begin(), expectedSorted.end(), StarIdCompare);
    sort(actualSorted.begin(), actualSorted.end(), StarIdCompare);

    auto currActual = actualSorted.cbegin();
    for (const StarIdentifier &currExpected : expectedSorted) {
        if (currActual != actualSorted.cend() &&
            currActual->starIndex == currExpected.starIndex) {

            if (actualCatalog[currActual->catalogIndex].name == expectedCatalog[currExpected.catalogIndex].name) {
                result.numCorrect++;
            } else {
                result.numIncorrect++;
            }

            currActual++;
        }
    }

    result.fractionCorrect = (float)result.numCorrect / result.numTotal;
    result.fractionIncorrect = (float)result.numIncorrect / result.numTotal;

    return result;
}

/////////////////////
// PIPELINE OUTPUT //
/////////////////////

typedef void (*PipelineComparator)(std::ostream &os,
                                   const PipelineInputList &,
                                   const std::vector<PipelineOutput> &,
                                   const PipelineOptions &);

cairo_status_t OstreamPlotter(void *closure, const unsigned char *data, unsigned int length) {
    std::ostream *os = (std::ostream *)closure;
    os->write((const char *)data, length);
    return CAIRO_STATUS_SUCCESS;
}

void PipelineComparatorPlotRawInput(std::ostream &os,
                                    const PipelineInputList &expected,
                                    const std::vector<PipelineOutput> &actual,
                                    const PipelineOptions &values) {
    
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

void PipelineComparatorPlotInput(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &actual,
                                 const PipelineOptions &values) {
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    assert(expected[0]->InputStars() != NULL);
    SurfacePlot(cairoSurface,
                *expected[0]->InputStars(),
                expected[0]->InputStarIds(),
                &expected[0]->GetCatalog(),
                expected[0]->InputAttitude(),
                // green
                0.0, 1.0, 0.0, 0.6);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

void PipelineComparatorCentroids(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &actual,
                                 const PipelineOptions &values) {
    int size = (int)expected.size();

    float threshold = values.threshold;

    std::vector<CentroidComparison> comparisons;
    for (int i = 0; i < size; i++) {
        comparisons.push_back(CentroidsCompare(threshold,
                                               *(expected[i]->ExpectedStars()),
                                               *(actual[i].stars)));
    }

    CentroidComparison result = CentroidComparisonsCombine(comparisons);
    os << "extra_stars " << result.numExtraStars << std::endl
       << "missing_stars " << result.numMissingStars << std::endl
       << "mean_error " << result.meanError << std::endl;
}

void PipelineComparatorPrintCentroids(std::ostream &os,
                                      const PipelineInputList &expected,
                                      const std::vector<PipelineOutput> &actual,
                                      const PipelineOptions &values) {
    assert(actual.size() == 1);
    assert(actual[0].stars);

    os << "num_centroids " << actual[0].stars->size() << std::endl;
    for (int i = 0; i < (int)actual[0].stars->size(); i++) {
        const Star &star = actual[0].stars->at(i);
        os << "centroid_" << i << "_x " << star.position.x << std::endl;
        os << "centroid_" << i << "_y " << star.position.y << std::endl;
        // TODO: print other stats too?
    }
}

void PipelineComparatorPlotOutput(std::ostream &os,
                                     const PipelineInputList &expected,
                                     const std::vector<PipelineOutput> &actual,
                                     const PipelineOptions &values) {
    // don't need to worry about mutating the surface; InputImageSurface returns a fresh one
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    SurfacePlot(cairoSurface,
                actual[0].stars ? *actual[0].stars : *expected[0]->ExpectedStars(),
                actual[0].starIds.get(),
                &actual[0].catalog,
                actual[0].attitude.get(),
                // red
                1.0, 0.0, 0.0, 0.5);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

void PipelineComparatorStars(std::ostream &os,
                             const PipelineInputList &expected,
                             const std::vector<PipelineOutput> &actual,
                             const PipelineOptions &values) {
    float centroidThreshold = actual[0].stars
        ? values.threshold
        : 0.0f;

    std::vector<StarIdComparison> comparisons;
    for (int i = 0; i < (int)expected.size(); i++) {
        comparisons.push_back(
            StarIdsCompare(*expected[i]->ExpectedStarIds(), *actual[i].starIds,
                           expected[i]->GetCatalog(), actual[i].catalog,
                           centroidThreshold, expected[i]->ExpectedStars(), actual[i].stars.get()));
    }

    if (comparisons.size() == 1) {
        os << "starid_num_correct " << comparisons[0].numCorrect << std::endl;
        os << "starid_num_incorrect " << comparisons[0].numIncorrect << std::endl;
        os << "starid_num_total " << comparisons[0].numTotal << std::endl;
    }
 
    float fractionIncorrectSum = 0;
    float fractionCorrectSum = 0;
    for (const StarIdComparison &comparison : comparisons) {
        fractionIncorrectSum += comparison.fractionIncorrect;
        fractionCorrectSum += comparison.fractionCorrect;
    }

    float fractionIncorrectMean = fractionIncorrectSum / comparisons.size();
    float fractionCorrectMean = fractionCorrectSum / comparisons.size();

    os << "starid_fraction_correct " << fractionCorrectMean << std::endl;
    os << "starid_fraction_incorrect " << fractionIncorrectMean << std::endl;
}

void PipelineComparatorPrintAttitude(std::ostream &os,
                                     const PipelineInputList &expected,
                                     const std::vector<PipelineOutput> &actual,
                                     const PipelineOptions &values) {
    assert(actual.size() == 1);
    assert(actual[0].attitude);

    os << "attitude_real " << actual[0].attitude->real << std::endl;
    os << "attitude_i " << actual[0].attitude->i << std::endl;
    os << "attitude_j " << actual[0].attitude->j << std::endl;
    os << "attitude_k " << actual[0].attitude->k << std::endl;
    float ra, de, roll;
    actual[0].attitude->ToSpherical(&ra, &de, &roll);
    os << "attitude_ra " << ra << std::endl;
    os << "attitude_de " << de << std::endl;
    os << "attitude_roll " << roll << std::endl;
}

void PipelineComparatorAttitude(std::ostream &os,
                                const PipelineInputList &expected,
                                const std::vector<PipelineOutput> &actual,
                                const PipelineOptions &values) {

    // TODO: use Wahba loss function (maybe average per star) instead of just angle. Also break
    // apart roll error from boresight error. This is just quick and dirty for testing
    float angleThreshold = DegToRad(values.threshold);

    float attitudeErrorSum = 0.0f;
    int numCorrect = 0;

    for (int i = 0; i < (int)expected.size(); i++) {
        float attitudeError = (*expected[i]->ExpectedAttitude() * actual[0].attitude->Conjugate()).Angle();
        assert(attitudeError >= 0);
        attitudeErrorSum += attitudeError;
        if (attitudeError <= angleThreshold) {
            numCorrect++;
        }
    }

    float attitudeErrorMean = attitudeErrorSum / expected.size();
    float fractionCorrect = (float)numCorrect / expected.size();

    os << "attitude_error_mean " << attitudeErrorMean << std::endl;
    os << "attitude_num_correct " << numCorrect << std::endl;
    os << "attitude_fraction_correct " << fractionCorrect << std::endl;
}

void PipelineComparison(const PipelineInputList &expected,
                              const std::vector<PipelineOutput> &actual, 
                              const PipelineOptions &values) {
    assert(expected.size() == actual.size() && expected.size() > 0);

    PipelineComparator comparator;
    std::string path;

    if (values.plotRawInput != "") {
        assert(expected[0]->InputImage() && expected.size() == 1);
        comparator = PipelineComparatorPlotRawInput;
        path = values.plotRawInput;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } 
    if (values.plotInput != "") {
        assert(expected[0]->InputImage() && expected.size() == 1 && expected[0]->InputStars());
        comparator = PipelineComparatorPlotInput;
        path = values.plotInput;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } 
    if (values.plotOutput != "") {
        assert(actual.size() == 1 && (actual[0].stars || actual[0].starIds));
        comparator = PipelineComparatorPlotOutput;
        path = values.plotOutput;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } 
    if (values.printCentroids != "") {
        assert(actual[0].stars && actual.size() == 1);
        comparator = PipelineComparatorPrintCentroids;
        path = values.printCentroids;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } 
    if (values.compareCentroids != "") {
        assert(actual[0].stars && expected[0]->ExpectedStars() && values.threshold);
        comparator = PipelineComparatorCentroids;
        path = values.compareCentroids;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } 
    if (values.compareStars != "") {
        assert(expected[0]->ExpectedStars() && actual[0].starIds && values.threshold);
        comparator = PipelineComparatorStars;
        path = values.compareStars;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } 
    if (values.printAttitude != "") {
        assert(actual[0].attitude && actual.size() == 1);
        comparator = PipelineComparatorPrintAttitude;
        path = values.printAttitude;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } 
    if (values.compareAttitude != "") {
        assert(actual[0].attitude && expected[0]->ExpectedAttitude() && values.threshold);
        comparator = PipelineComparatorAttitude;
        path = values.compareAttitude;
        PromptedOutputStream pos(path);
        comparator(pos.Stream(), expected, actual, values);
    } else {
        return;
    }
}

}
