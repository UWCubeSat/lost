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
#include <limits.h>

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
bool isDebug;

void RegisterCliArgs(int newArgc, char **newArgv) {
    argv = newArgv + 1;
    argc = newArgc - 1;
    // any value of the environment variable LOST_DEBUG enters debug mode
    isDebug = getenv("LOST_DEBUG") != NULL;
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

std::string PromptLine(const std::string &prompt) {
    std::cerr << prompt << ": ";
    if (HasNextCliArg()) {
        return NextCliArg();
    }
    std::string result;
    std::getline(std::cin, result);
    return result;
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

    result = new unsigned char[width*height];
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
                 const Attitude *attitude,
                 double red,
                 double green,
                 double blue,
                 double alpha,
                 // if true, don't use catalog name
                 bool rawStarIndexes = false) {
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

            int plotName = rawStarIndexes ? starId.catalogIndex : (*catalog)[starId.catalogIndex].name;
            cairo_show_text(cairoCtx, std::to_string(plotName).c_str());
        }
        metadata += std::to_string(starIds->size()) + " identified   ";
    }

    if (attitude != NULL) {
        EulerAngles spherical = attitude->ToSpherical();
        metadata +=
            "RA: " + std::to_string(RadToDeg(spherical.ra)) + "  " +
            "DE: " + std::to_string(RadToDeg(spherical.de)) + "  " +
            "Roll: " + std::to_string(RadToDeg(spherical.roll)) + "   ";
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
typedef Santa *(*SantaFactory)();

CentroidAlgorithm *DummyCentroidAlgorithmPrompt() {
    int numStars = Prompt<int>("How many stars to generate");
    return new DummyCentroidAlgorithm(numStars);
}

CentroidAlgorithm *CoGCentroidAlgorithmPrompt() {
    return new CenterOfGravityAlgorithm();
}

CentroidAlgorithm *IWCoGCentroidAlgorithmPrompt() {
    return new IterativeWeightedCenterOfGravityAlgorithm();
}

StarIdAlgorithm *DummyStarIdAlgorithmPrompt() {
    return new DummyStarIdAlgorithm();
}

StarIdAlgorithm *GeometricVotingStarIdAlgorithmPrompt() {
    float tolerance = Prompt<float>("Angular tolerance? (degrees)");
    return new GeometricVotingStarIdAlgorithm(DegToRad(tolerance));
}

StarIdAlgorithm *PyramidStarIdAlgorithmPrompt() {
    float tolerance = Prompt<float>("Angular tolerance? (degrees)");
    int numFalseStars = Prompt<int>("Estimated # of false stars (whole sphere)");
    float maxMismatchProbability = Prompt<float>("Maximum mismatch probability");
    return new PyramidStarIdAlgorithm(DegToRad(tolerance), numFalseStars, maxMismatchProbability, 1000);
}

AttitudeEstimationAlgorithm *DavenportQAlgorithmPrompt() {
    return new DavenportQAlgorithm();
}

AttitudeEstimationAlgorithm *TriadAlgorithmPrompt() {
    return new TriadAlgorithm();
}

class SantaStarId : public Santa {
public:
    SantaStarId(int minStars) : minStars(minStars) { };
    bool Go(const PipelineOutput &output) const override {
        if (output.starIds.get() == NULL) {
            std::cerr << "Warning: Star-id santa in use but star-ids is null." << std::endl;
            return false;
        }

        return (int)output.starIds->size() >= minStars;
    };
private:
    int minStars;
};

Santa *SantaStarIdPrompt() {
    int minStars = Prompt<int>("How many stars to count as nice?");
    return new SantaStarId(minStars);
}

Catalog PromptNarrowedCatalog(const Catalog &catalog) {
    float maxSomething = Prompt<float>("Max magnitude or # of stars");
    int maxMagnitude = 1000;
    int maxStars = 10000;
    assert(maxSomething > 0);
    if (maxSomething < 10) {
        maxMagnitude = 100 * (int)floor(maxSomething);
    } else {
        maxStars = (int)floor(maxSomething);
    }
    return NarrowCatalog(catalog, maxMagnitude, maxStars);
}

void PromptKVectorDatabaseBuilder(MultiDatabaseBuilder &builder, const Catalog &catalog) {
    float minDistance = DegToRad(Prompt<float>("Min distance (deg)"));
    float maxDistance = DegToRad(Prompt<float>("Max distance (deg)"));
    long numBins = Prompt<long>("Number of distance bins");

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

void PromptDatabases(MultiDatabaseBuilder &builder, const Catalog &catalog) {
    InteractiveChoice<DbBuilder> dbBuilderChoice;
    dbBuilderChoice.Register("kvector", "K-Vector (geometric voting & pyramid)", PromptKVectorDatabaseBuilder);
    dbBuilderChoice.Register("done", "Exit", NULL);
    while (true) {
        DbBuilder choice = dbBuilderChoice.Prompt("Choose database builder");
        if (choice == NULL) {
            break;
        }
        (*choice)(builder, catalog);
    }
}

std::ostream &operator<<(std::ostream &os, const Camera &camera) {
    os << "camera_focal_length " << camera.FocalLength() << std::endl
       << "camera_fov " << camera.Fov() << std::endl
       << "camera_resolution_x " << camera.XResolution() << std::endl
       << "camera_resolution_y " << camera.YResolution() << std::endl;
    // TODO: principal point
    return os;
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
    const Attitude *InputAttitude() const { return &attitude; };
private:
    Image image;
    Attitude attitude;
};

Camera PromptCamera(int xResolution, int yResolution) {
    float pixelSize = Prompt<float>("Pixel size (µm) for focal length or zero for FOV");
    float focalLengthOrFov = Prompt<float>("Focal Length (mm) or horizontal FOV (degrees)");
    float focalLengthPixels = pixelSize == 0.0f
        ? FovToFocalLength(DegToRad(focalLengthOrFov), xResolution)
        : focalLengthOrFov * 1000 / pixelSize;

    return Camera(focalLengthPixels, xResolution, yResolution);
}

PngPipelineInput::PngPipelineInput(cairo_surface_t *cairoSurface, Camera camera, const Catalog &catalog)
    : camera(camera), catalog(catalog) {

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
        std::cout << "Read status: " << cairo_status_to_string(cairo_surface_status(cairoSurface)) << std::endl;
    }
    int xResolution = cairo_image_surface_get_width(cairoSurface);
    int yResolution = cairo_image_surface_get_height(cairoSurface);
    result.push_back(std::unique_ptr<PipelineInput>(
                         new PngPipelineInput(cairoSurface, PromptCamera(xResolution, yResolution), CatalogRead())));
    return result;
}

AstrometryPipelineInput::AstrometryPipelineInput(const std::string &path) {
    // create from path, TODO
}

// PipelineInputList PromptAstrometryPipelineInput() {
//     // TODO: why does it let us do a reference to an ephemeral return value?
//     std::string path = Prompt<std::string>("Astrometry download directory");
//     PipelineInputList result;
//     result.push_back(std::unique_ptr<PipelineInput>(new AstrometryPipelineInput(path)));
//     return result;
// }

class GeneratedStar : public Star {
public:
    GeneratedStar(Star star, float peakBrightness, Vec2 motionBlurDelta)
        : Star(star), peakBrightness(peakBrightness), delta(motionBlurDelta) { };
    
    float peakBrightness;
    // vector points to where the star will be visibly after one time unit.
    Vec2 delta;
};

// In the equations for pixel brightness both with motion blur enabled and disabled, we don't need
// any constant scaling factor outside the integral because when d0=0, the brightness at the center
// will be zero without any scaling. The scaling factor you usually see on a Normal distribution is
// so that the Normal distribution integrates to one over the real line, making it a probability
// distribution. But we want the /peak/ to be one, not the integral. motion blur enabled

// pixel is the pixel we want the brightness of, p0 is "center" star position, delta is how star
// moves after one time unit. Call with final t and initial t value to get the integral. Brightness may not be the same as star's brightness because of oversampling. See
// https://wiki.huskysat.org/wiki/index.php/Motion_Blur_and_Rolling_Shutter#Motion_Blur_Math
// for details.
static float motionBlurredPixelBrightness(const Vec2 &pixel, const GeneratedStar &generatedStar,
                                          float t, float stddev, float brightness) {
    const Vec2 &p0 = generatedStar.position;
    const Vec2 &delta = generatedStar.delta;
    const Vec2 d0 = p0 - pixel;
    return brightness
        * stddev*sqrt(M_PI) / (sqrt(2)*delta.Magnitude())
        * exp(pow(d0.x*delta.x + d0.y*delta.y, 2) / (2*stddev*stddev*delta.MagnitudeSq())
              - d0.MagnitudeSq() / (2*stddev*stddev))
        * erf((t*delta.MagnitudeSq() + d0.x*delta.x + d0.y*delta.y) / (stddev*sqrt(2)*delta.Magnitude()));
}

static float staticPixelBrightness(const Vec2 &pixel, const GeneratedStar &generatedStar,
                                   float stddev, float brightness) {
    const Vec2 d0 = generatedStar.position - pixel;
    return brightness * exp(-d0.MagnitudeSq() / (2*stddev*stddev));
}

const int kMaxBrightness = 255;

GeneratedPipelineInput::GeneratedPipelineInput(const Catalog &catalog,
                                               Attitude attitude,
                                               Camera camera,
                                               float observedReferenceBrightness,
                                               float starSpreadStdDev,
                                               float sensitivity,
                                               float darkCurrent,
                                               float readNoiseStdDev,
                                               Attitude motionBlurDirection, // applied on top of the attitude
                                               float exposureTime, // zero for no motion blur
                                               float readoutTime, // zero for no rolling shutter
                                               bool shotNoise,
                                               int oversampling
    ) : camera(camera), attitude(attitude), catalog(catalog) {

    // in photons
    float referenceBrightness = observedReferenceBrightness / sensitivity;

    image.width = camera.XResolution();
    image.height = camera.YResolution();
    // number of true photons each pixel receives.

    assert(oversampling >= 1);
    int oversamplingPerAxis = ceil(sqrt(oversampling));
    if (oversamplingPerAxis*oversamplingPerAxis != oversampling) {
        std::cerr << "WARNING: oversampling was not a perfect square. Rounding up to "
                  << oversamplingPerAxis*oversamplingPerAxis << "." << std::endl;
    }
    std::vector<float> photonsBuffer(image.width*image.height, 0);

    bool motionBlurEnabled = exposureTime > 0;
    // TODO: ensure works when motion blur disabled
    if (!motionBlurEnabled) {
        exposureTime = 1.0;
        motionBlurDirection = Attitude(Quaternion(0, 1, 0, 0));
    }
    Quaternion motionBlurDirectionQ = motionBlurDirection.GetQuaternion();
    // attitude at the middle of exposure time
    Quaternion currentAttitude = attitude.GetQuaternion();
    // attitude 1 time unit after middle of exposure
    Quaternion futureAttitude = motionBlurDirectionQ*currentAttitude;
    std::vector<GeneratedStar> generatedStars;

    for (int i = 0; i < (int)catalog.size(); i++) {
        const CatalogStar &catalogStar = catalog[i];
        Vec3 rotated = attitude.Rotate(catalog[i].spatial);
        if (rotated.x < 0) {
            continue;
        }
        Vec2 camCoords = camera.SpatialToCamera(rotated);

        if (camera.InSensor(camCoords)) {
            Vec3 futureSpatial = futureAttitude.Rotate(catalog[i].spatial);
            Vec2 delta = camera.SpatialToCamera(futureSpatial) - camCoords;
            if (!motionBlurEnabled) {
                delta = {0, 0}; // avoid floating point funny business
            }
            // radiant intensity, in photons per time unit per pixel, at the center of the star.
            float peakBrightness = referenceBrightness * pow(10.0, catalogStar.magnitude/-250.0);
            float interestingThreshold = 0.01; // we don't need to check pixels that are expected to
                                               // receive this many photons or fewer.
            // inverse of the function defining the Gaussian distribution: Find out how far from the
            // mean we'll have to go until the number of photons is less than interestingThreshold
            float radius = ceil(sqrt(-log(interestingThreshold/peakBrightness/exposureTime)*2*starSpreadStdDev*starSpreadStdDev));
            Star star = Star(camCoords.x, camCoords.y,
                             radius, radius,
                             catalogStar.magnitude);
            stars.push_back(star);
            generatedStars.push_back(GeneratedStar(star, peakBrightness, delta));
            starIds.push_back(StarIdentifier(stars.size() - 1, i));
        }
    }

    for (const GeneratedStar &star : generatedStars) {
        // delta will be exactly (0,0) when motion blur disabled
        Vec2 earliestPosition = star.position - star.delta*(exposureTime/2.0 + readoutTime/2.0);
        Vec2 latestPosition = star.position + star.delta*(exposureTime/2.0 + readoutTime/2.0);
        int xMin = std::max(0, (int)std::min(earliestPosition.x - star.radiusX, latestPosition.x - star.radiusX));
        int xMax = std::min(image.width-1, (int)std::max(earliestPosition.x + star.radiusX, latestPosition.x + star.radiusX));
        int yMin = std::max(0, (int)std::min(earliestPosition.y - star.radiusX, latestPosition.y - star.radiusX));
        int yMax = std::min(image.height-1, (int)std::max(earliestPosition.y + star.radiusX, latestPosition.y + star.radiusX));

        // peak brightness is measured in photons per time unit per pixel, so if oversampling, we
        // need to convert units to photons per time unit per sample
        float peakBrightnessPerSample = star.peakBrightness / (oversamplingPerAxis*oversamplingPerAxis);

        // the star.x and star.y refer to the pixel whose top left corner the star should appear at
        // (and fractional amounts are relative to the corner). When we color a pixel, we ideally
        // would integrate the intensity of the star over that pixel, but we can make do by sampling
        // the intensity of the star at the /center/ of the pixel, ie, star.x+.5 and star.y+.5
        for(int xPixel = xMin; xPixel <= xMax; xPixel++) {
            for (int yPixel = yMin; yPixel <= yMax; yPixel++) {
                // offset of beginning & end of readout compared to beginning & end of readout for
                // center row
                float readoutOffset = readoutTime * (yPixel - image.height/2.0) / image.height;
                float tStart = -exposureTime/2.0 + readoutOffset;
                float tEnd = exposureTime/2.0 + readoutOffset;

                // loop through all samples in the current pixel
                for (int xSample = 0; xSample < oversamplingPerAxis; xSample++) {
                    for (int ySample = 0; ySample < oversamplingPerAxis; ySample++) {
                        float x = xPixel + (xSample+0.5)/oversamplingPerAxis;
                        float y = yPixel + (ySample+0.5)/oversamplingPerAxis;

                        float curPhotons;
                        if (motionBlurEnabled) {
                            curPhotons =
                                motionBlurredPixelBrightness({x, y}, star, tEnd, starSpreadStdDev, peakBrightnessPerSample)
                                - motionBlurredPixelBrightness({x, y}, star, tStart, starSpreadStdDev, peakBrightnessPerSample);
                        } else {
                            curPhotons = staticPixelBrightness({x, y}, star, starSpreadStdDev, peakBrightnessPerSample);
                        }

                        assert(0.0 <= curPhotons);

                        photonsBuffer[xPixel + yPixel*image.width] += curPhotons;
                    }
                }
            }
        }
    }

    std::default_random_engine rng;
    std::normal_distribution<float> readNoiseDist(0.0, readNoiseStdDev);

    // convert from photon counts to observed pixel brightnesses, applying noise and such.
    imageData = std::vector<unsigned char>(image.width*image.height);
    image.image = imageData.data();
    for (int i = 0; i < image.width * image.height; i++) {
        float curBrightness = 0;

        // dark current (Constant)
        curBrightness += darkCurrent;

        // read noise (Gaussian)
        curBrightness += readNoiseDist(rng);

        // shot noise (Poisson), and quantize
        long quantizedPhotons;
        if (shotNoise) {
            // with GNU libstdc++, it keeps sampling from the distribution until it's within the min-max
            // range. This is problematic if the mean is far above the max long value, because then it
            // might have to sample many many times (and furthermore, the results won't be useful
            // anyway)
            float photons = photonsBuffer[i];
            if (photons > (float)LONG_MAX - 3.0*sqrt(LONG_MAX)) {
                std::cout << "ERROR: One of the pixels had too many photons. Generated image would not be physically accurate, exiting." << std::endl;
                exit(1);
            }
            std::poisson_distribution<long> shotNoiseDist(photonsBuffer[i]);
            quantizedPhotons = shotNoiseDist(rng);
        } else {
            quantizedPhotons = round(photonsBuffer[i]);
        }
        curBrightness += quantizedPhotons * sensitivity;
        
        // std::clamp not introduced until C++17, so we avoid it.
        float clampedBrightness = std::max(std::min(curBrightness, (float)1.0), (float)0.0);
        imageData[i] = floor(clampedBrightness * kMaxBrightness); // TODO: off-by-one, 256?
    }
}

Quaternion PromptSphericalAttitude(const std::string &name) {
    float ra = Prompt<float>(name + " right ascension");
    float dec = Prompt<float>(name + " declination");
    float roll = Prompt<float>(name + " roll");
    return SphericalToQuaternion(DegToRad(ra), DegToRad(dec), DegToRad(roll));
}

PipelineInputList PromptGeneratedPipelineInput() {
    // TODO: prompt for attitude, imagewidth, etc and then construct a GeneratedPipelineInput
    int numImages = Prompt<int>("Number of images to generate");
    int xResolution = Prompt<int>("Horizontal Resolution");
    int yResolution = Prompt<int>("Vertical Resolution");
    Camera camera = PromptCamera(xResolution, yResolution);
    float observedReferenceBrightness = Prompt<int>("Observed reference star brightness");
    float starSpreadStdDev = Prompt<float>("Star spread stddev");
    float sensitivity = Prompt<float>("Camera sensitivity");
    float darkCurrent = Prompt<float>("Dark current");
    float readNoiseStdDev = Prompt<float>("Noise stddev");
    float exposureTime = Prompt<float>("Exposure time");
    float readoutTime = Prompt<float>("Readout time");
    bool shotNoise = Prompt<bool>("Enable shot noise");
    int oversampling = Prompt<int>("Oversampling");

    // TODO: allow random angle generation?
    Quaternion attitude = PromptSphericalAttitude("Boresight");
    Quaternion motionBlurDirection;
    if (exposureTime > 0) {
        motionBlurDirection = PromptSphericalAttitude("Motion blur direction");
    }

    PipelineInputList result;

    for (int i = 0; i < numImages; i++) {
        GeneratedPipelineInput *curr = new GeneratedPipelineInput(
            CatalogRead(),
            attitude,
            camera,
            observedReferenceBrightness, starSpreadStdDev,
            sensitivity, darkCurrent, readNoiseStdDev, motionBlurDirection,
            exposureTime, readoutTime, shotNoise, oversampling);

        result.push_back(std::unique_ptr<PipelineInput>(curr));
    }

    return result;
}

typedef PipelineInputList (*PipelineInputFactory)();

PipelineInputList PromptPipelineInput() {
    InteractiveChoice<PipelineInputFactory> inputTypeChoice;
    inputTypeChoice.Register("png", "PNG files", PromptPngPipelineInput);
    inputTypeChoice.Register("generate", "Generated image", PromptGeneratedPipelineInput);
    // inputTypeChoice.Register("astrometry", "Astrometry.net", PromptAstrometryPipelineInput);

    return (inputTypeChoice.Prompt("Input from"))();
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

Pipeline PromptPipeline() {
    enum class PipelineStage {
        Centroid, CentroidMagnitudeFilter, Database, StarId, AttitudeEstimation, Santa, Done
    };

    Pipeline result;
    
    InteractiveChoice<PipelineStage> stageChoice;
    stageChoice.Register("centroid", "Centroid", PipelineStage::Centroid);
    // TODO: more flexible or sth
    stageChoice.Register("centroid_magnitude_filter", "Centroid Magnitude Filter", PipelineStage::CentroidMagnitudeFilter);
    // TODO: don't allow setting star-id until database is set, and perhaps limit the star-id
    // choices to those compatible with the database?
    stageChoice.Register("database", "Database", PipelineStage::Database);
    stageChoice.Register("starid", "Star-ID", PipelineStage::StarId);
    stageChoice.Register("attitude", "Attitude Estimation", PipelineStage::AttitudeEstimation);
    stageChoice.Register("santa", "Santa (knows if the output is naughty or nice)", PipelineStage::Santa);
    stageChoice.Register("done", "Done setting up pipeline", PipelineStage::Done);

    while (true) {
        PipelineStage nextStage = stageChoice.Prompt("Which pipeline stage to set");
        switch (nextStage) {

        case PipelineStage::Centroid: {
            InteractiveChoice<CentroidAlgorithmFactory> centroidChoice;
            centroidChoice.Register("cog", "Center of Gravity Centroid Algorithm",
                                    CoGCentroidAlgorithmPrompt);
            centroidChoice.Register("iwcog", "Iterative Weighted Center of Gravity Algorithm",
                                    IWCoGCentroidAlgorithmPrompt);
            if (isDebug) {
                centroidChoice.Register("dummy", "Random Centroid Algorithm",
                                        DummyCentroidAlgorithmPrompt);
            }

            result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(
                (centroidChoice.Prompt("Choose centroid algo"))());
            break;
        }

        case PipelineStage::CentroidMagnitudeFilter: {
            result.centroidMinMagnitude = Prompt<int>("Minimum magnitude");
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
            break;
        }

        case PipelineStage::StarId: {
            InteractiveChoice<StarIdAlgorithmFactory> starIdChoice;
            if (isDebug) {
                starIdChoice.Register("dummy", "Random", DummyStarIdAlgorithmPrompt);
            }
            starIdChoice.Register("gv", "Geometric Voting", GeometricVotingStarIdAlgorithmPrompt);
            starIdChoice.Register("pyramid", "Pyramid Scheme", PyramidStarIdAlgorithmPrompt);

            result.starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(
                (starIdChoice.Prompt("Choose Star-ID algo"))());
            break;
        }

        case PipelineStage::AttitudeEstimation: {
            InteractiveChoice<AttitudeEstimationAlgorithmFactory> attitudeEstimationAlgorithmChoice;
            attitudeEstimationAlgorithmChoice.Register("dqm", "Davenport Q Method", DavenportQAlgorithmPrompt);
            attitudeEstimationAlgorithmChoice.Register("triad", "Triad", TriadAlgorithmPrompt);
            result.attitudeEstimationAlgorithm = std::unique_ptr<AttitudeEstimationAlgorithm>(
                (attitudeEstimationAlgorithmChoice.Prompt("Choose Attitude algo"))());
            break;
        }

        case PipelineStage::Santa: {
            InteractiveChoice<SantaFactory> santaChoice;
            santaChoice.Register("star_id_count", "Number of stars identified", SantaStarIdPrompt);
            result.santa = std::unique_ptr<Santa>(
                (santaChoice.Prompt("Choose Santa"))());
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
    const Stars *inputStars = input.InputStars();
    const StarIdentifiers *inputStarIds = input.InputStarIds();

    // if database is provided, that's where we get catalog from.
    if (database) {
        MultiDatabase multiDatabase(database.get());
        const unsigned char *catalogBuffer = multiDatabase.SubDatabasePointer(kCatalogMagicValue);
        if (catalogBuffer != NULL) {
            result.catalog = DeserializeCatalog(multiDatabase.SubDatabasePointer(kCatalogMagicValue), NULL, NULL);
        } else {
            std::cerr << "Warning: That database does not include a catalog. Proceeding with the full catalog." << std::endl;
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
    }

    if (starIdAlgorithm && database && inputStars && input.InputCamera()) {
        // TODO: don't copy the vector!
        result.starIds = std::unique_ptr<StarIdentifiers>(new std::vector<StarIdentifier>(
            starIdAlgorithm->Go(database.get(), *inputStars, result.catalog, *input.InputCamera())));
        inputStarIds = result.starIds.get();
    }

    if (attitudeEstimationAlgorithm && inputStarIds && input.InputCamera()) {
        assert(inputStars); // ensure that starIds doesn't exist without stars
        result.attitude = std::unique_ptr<Attitude>(
            new Attitude(attitudeEstimationAlgorithm->Go(*input.InputCamera(), *inputStars, result.catalog, *inputStarIds)));
    }

    if (santa) {
        result.nice = santa->Go(result);
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
#if LOST_DEBUG > 2
                std::cout << "Expected: " << expectedCatalog[currExpected.catalogIndex].name
                          << " Actual: " << actualCatalog[currActual->catalogIndex].name << std::endl;
#endif
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
                                   const std::vector<PipelineOutput> &);

cairo_status_t OstreamPlotter(void *closure, const unsigned char *data, unsigned int length) {
    std::ostream *os = (std::ostream *)closure;
    os->write((const char *)data, length);
    return CAIRO_STATUS_SUCCESS;
}

void PipelineComparatorPlotRawInput(std::ostream &os,
                                    const PipelineInputList &expected,
                                    const std::vector<PipelineOutput> &actual) {
    
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

void PipelineComparatorPlotInput(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &actual) {
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
                                 const std::vector<PipelineOutput> &actual) {
    int size = (int)expected.size();

    float threshold = Prompt<float>("Threshold to count as the same star (pixels, float)");

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
                                      const std::vector<PipelineOutput> &actual) {
    assert(actual.size() == 1);
    assert(actual[0].stars);

    os << "num_centroids " << actual[0].stars->size() << std::endl;
    for (int i = 0; i < (int)actual[0].stars->size(); i++) {
        const Star &star = actual[0].stars->at(i);
        os << "centroid_" << i << "_x " << star.position.x << std::endl;
        os << "centroid_" << i << "_y " << star.position.y << std::endl;
        if (actual[0].starIds) {
            for (const StarIdentifier &starId : *actual[0].starIds) {
                if (starId.starIndex == i) {
                    os << "centroid_" << i << "_id " << actual[0].catalog[starId.catalogIndex].name << std::endl;
                }
            }
        }
        if (expected[0]->ExpectedStarIds()) {
            for (const StarIdentifier &starId : *expected[0]->ExpectedStarIds()) {
                if (starId.starIndex == i) {
                    os << "centroid_" << i << "_expected_id " << actual[0].catalog[starId.catalogIndex].name << std::endl;
                }
            }
        }
        // TODO: print other stats too?
    }
}

void PipelineComparatorPlotIndexes(std::ostream &os,
                                   const PipelineInputList &expected,
                                   const std::vector<PipelineOutput> &actual) {
    StarIdentifiers identifiers;
    for (int i = 0; i < (int)actual[0].stars->size(); i++) {
        identifiers.push_back(StarIdentifier(i, i));
    }
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    SurfacePlot(cairoSurface,
                actual[0].stars ? *actual[0].stars : *expected[0]->ExpectedStars(),
                &identifiers,
                &actual[0].catalog,
                NULL,
                // orange
                1.0, 0.5, 0.0, 0.5,
                // don't resolve names
                true);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

void PipelineComparatorPlotOutput(std::ostream &os,
                                  const PipelineInputList &expected,
                                  const std::vector<PipelineOutput> &actual) {
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
                             const std::vector<PipelineOutput> &actual) {
    float centroidThreshold = actual[0].stars
        ? Prompt<float>("Threshold to count centroids as the same star (pixels)")
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
                                     const std::vector<PipelineOutput> &actual) {
    assert(actual.size() == 1);
    assert(actual[0].attitude);

    // os << "attitude_real " << actual[0].attitude->real << std::endl;
    // os << "attitude_i " << actual[0].attitude->i << std::endl;
    // os << "attitude_j " << actual[0].attitude->j << std::endl;
    // os << "attitude_k " << actual[0].attitude->k << std::endl;
    EulerAngles spherical = actual[0].attitude->ToSpherical();
    os << "attitude_ra " << RadToDeg(spherical.ra) << std::endl;
    os << "attitude_de " << RadToDeg(spherical.de) << std::endl;
    os << "attitude_roll " << RadToDeg(spherical.roll) << std::endl;
}

void PipelineComparatorAttitude(std::ostream &os,
                                const PipelineInputList &expected,
                                const std::vector<PipelineOutput> &actual) {

    // TODO: use Wahba loss function (maybe average per star) instead of just angle. Also break
    // apart roll error from boresight error. This is just quick and dirty for testing
    float angleThreshold = DegToRad(
        Prompt<float>("Threshold to count two attitudes as the same (deg)"));

    float attitudeErrorSum = 0.0f;
    int numCorrect = 0;

    for (int i = 0; i < (int)expected.size(); i++) {
        Quaternion expectedQuaternion = expected[i]->ExpectedAttitude()->GetQuaternion();
        Quaternion actualQuaternion = actual[i].attitude->GetQuaternion();
        float attitudeError = (expectedQuaternion * actualQuaternion.Conjugate()).Angle();
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

void PipelineComparatorSanta(std::ostream &os,
                             const PipelineInputList &expected,
                             const std::vector<PipelineOutput> &actual) {
    if (actual.size() == 1) {
        os << "nice " << (actual[0].nice ? "true" : "false") << std::endl;
    } else {
        int numNice = 0;
        for (const PipelineOutput &output : actual) {
            if (output.nice) {
                numNice++;
            }
        }
        os << "num_nice " << numNice << std::endl;
    }
}

void PipelineComparatorPrintPairDistance(std::ostream &os,
                                         const PipelineInputList &expected,
                                         const std::vector<PipelineOutput> &actual) {
    int index1 = Prompt<int>("Index of first star");
    int index2 = Prompt<int>("Index of second star");

    const Camera &camera = *expected[0]->InputCamera();
    const Stars &stars = *actual[0].stars;

    assert(index1 >= 0 && index2 >= 0 && index1 < (int)stars.size() && index2 < (int)stars.size());
    os << "pair_distance " << Angle(camera.CameraToSpatial(stars[index1].position),
                                    camera.CameraToSpatial(stars[index2].position))
       << std::endl;
}

void PipelineComparatorPrintPyramidDistances(std::ostream &os,
                                             const PipelineInputList &expected,
                                             const std::vector<PipelineOutput> &actual) {
    int index1 = Prompt<int>("Catalog name/index of first star");
    int index2 = Prompt<int>("Catalog name/index of second star");
    int index3 = Prompt<int>("Catalog name/index of third star");
    int index4 = Prompt<int>("Catalog name/index of four star");

    const Camera &camera = *expected[0]->InputCamera();
    const Stars &stars = *actual[0].stars;

    Vec3 spatial1 = camera.CameraToSpatial(stars[index1].position);
    Vec3 spatial2 = camera.CameraToSpatial(stars[index2].position);
    Vec3 spatial3 = camera.CameraToSpatial(stars[index3].position);
    Vec3 spatial4 = camera.CameraToSpatial(stars[index4].position);

    std::cout << "pair_distance_12 " << Angle(spatial1, spatial2) << std::endl;
    std::cout << "pair_distance_13 " << Angle(spatial1, spatial3) << std::endl;
    std::cout << "pair_distance_14 " << Angle(spatial1, spatial4) << std::endl;
    std::cout << "pair_distance_23 " << Angle(spatial2, spatial3) << std::endl;
    std::cout << "pair_distance_24 " << Angle(spatial2, spatial4) << std::endl;
    std::cout << "pair_distance_34 " << Angle(spatial3, spatial4) << std::endl;
}

void PipelineComparatorPrintTripleAngle(std::ostream &os,
                                        const PipelineInputList &expected,
                                        const std::vector<PipelineOutput> &actual) {
    int index1 = Prompt<int>("Index of first star");
    int index2 = Prompt<int>("Index of second star");
    int index3 = Prompt<int>("Index of third star");

    const Camera &camera = *expected[0]->InputCamera();
    const Stars &stars = *actual[0].stars;

    assert(index1 >= 0 && index1 < (int)stars.size());
    assert(index2 >= 0 && index2 < (int)stars.size());
    assert(index3 >= 0 && index3 < (int)stars.size());

    // TODO, when merging with nondimensional branch
}

void PromptPipelineComparison(const PipelineInputList &expected,
                              const std::vector<PipelineOutput> &actual) {
    assert(expected.size() == actual.size() && expected.size() > 0);

    InteractiveChoice<PipelineComparator> comparatorChoice;

    if (isDebug) {
        if (actual[0].stars && expected[0]->InputCamera() != NULL) {
            comparatorChoice.Register("print_pair_distance", "Angular distance between two stars",
                                      PipelineComparatorPrintPairDistance);

            comparatorChoice.Register("print_pyramid_distances", "Distances between all pairs in a pyramid",
                                      PipelineComparatorPrintPyramidDistances);

            comparatorChoice.Register("print_triple_angle", "Inner angle of a star triangle",
                                      PipelineComparatorPrintTripleAngle);
        }
    }

    if (expected[0]->InputImage() && expected.size() == 1) {
        comparatorChoice.Register("plot_raw_input", "Plot raw BW input image to PNG",
                                  PipelineComparatorPlotRawInput);

        if (expected[0]->InputStars()) {
            comparatorChoice.Register("plot_input", "Plot annotated input image to PNG",
                                      PipelineComparatorPlotInput);
        }
    }

    if (isDebug && actual.size() == 1 && actual[0].stars) {
        comparatorChoice.Register("plot_indexes", "Plot centroid indexes",
                                  PipelineComparatorPlotIndexes);
    }

    if (actual.size() == 1 && (actual[0].stars || actual[0].starIds)) {
        comparatorChoice.Register("plot_output", "Plot output to PNG",
                                  PipelineComparatorPlotOutput);
    }

    // Centroids
    if (actual[0].stars) {
        if (actual.size() == 1) {
            comparatorChoice.Register("print_centroids", "Print list of centroids",
                                      PipelineComparatorPrintCentroids);
        }

        if (expected[0]->ExpectedStars()) {
            comparatorChoice.Register("compare_centroids", "Compare lists of centroids",
                                      PipelineComparatorCentroids);
        }
    }

    // Star-IDs
    if (expected[0]->ExpectedStars() && actual[0].starIds) {
                comparatorChoice.Register("compare_stars", "Compare lists of identified stars",
                                          PipelineComparatorStars);
    }

    if (actual[0].attitude) {
        if (actual.size() == 1) {
            comparatorChoice.Register("print_attitude", "Print the determined ra, de, and roll",
                                      PipelineComparatorPrintAttitude);
        }

        if (expected[0]->ExpectedAttitude()) {
            comparatorChoice.Register("compare_attitude", "Compare expected to actual attitude",
                                      PipelineComparatorAttitude);
        }
    }

    comparatorChoice.Register("print_nice", "Print if output is naughty or nice",
                              PipelineComparatorSanta);

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

typedef void (*CatalogInspector)(const Catalog &);

static std::vector<const CatalogStar *> PromptCatalogStars(const Catalog &catalog, int howMany) {
    std::vector<const CatalogStar *> result;
    for (int i = 0; i < howMany; i++) {
        int name = Prompt<int>("Catalog name of " + std::to_string(i) + "-th star");
        const CatalogStar *star = findNamedStar(catalog, name);
        if (star == NULL) {
            std::cerr << "Star not found!" << std::endl;
            exit(1);
        }
        result.push_back(star);
    }
    return result;
}

void InspectPairDistance(const Catalog &catalog) {
    auto stars = PromptCatalogStars(catalog, 2);

    // TODO: not cout, prompt for an ostream in inspect and pass argument
    std::cout << Angle(stars[0]->spatial, stars[1]->spatial) << std::endl;
}

void InspectPyramidDistances(const Catalog &catalog) {
    auto stars = PromptCatalogStars(catalog, 4);

    std::cout << "pair_distance_01 " << Angle(stars[0]->spatial, stars[1]->spatial) << std::endl;
    std::cout << "pair_distance_02 " << Angle(stars[0]->spatial, stars[2]->spatial) << std::endl;
    std::cout << "pair_distance_03 " << Angle(stars[0]->spatial, stars[3]->spatial) << std::endl;
    std::cout << "pair_distance_12 " << Angle(stars[1]->spatial, stars[2]->spatial) << std::endl;
    std::cout << "pair_distance_13 " << Angle(stars[1]->spatial, stars[3]->spatial) << std::endl;
    std::cout << "pair_distance_23 " << Angle(stars[2]->spatial, stars[3]->spatial) << std::endl;
}

void InspectTripleAngle(const Catalog &catalog) {
    auto stars = PromptCatalogStars(catalog, 3);

    // TODO
}

void InspectFindStar(const Catalog &catalog) {
    std::string raStr = PromptLine("Right Ascension");

    float raRadians;

    int raHours, raMinutes;
    float raSeconds;
    int raFormatTime = sscanf(raStr.c_str(), "%dh %dm %fs", &raHours, &raMinutes, &raSeconds);
    
    float raDeg;
    int raFormatDeg = sscanf(raStr.c_str(), "%f", &raDeg);

    if (raFormatTime == 3) {
        raRadians = (raHours * 2*M_PI/24) + (raMinutes * 2*M_PI/24/60) + (raSeconds * 2*M_PI/24/60/60);
    } else if (raFormatDeg == 1) {
        raRadians = DegToRad(raFormatDeg);
    } else {
        std::cerr << "Invalid right ascension format. Do \"09h 38m 29.8754s\" or a number of degrees." << std::endl;
        exit(1);
    }

    std::string deStr = PromptLine("Declination");

    float deRadians;

    int deDegPart, deMinPart;
    float deSecPart;
    char dummy[8];
    int deFormatParts = sscanf(deStr.c_str(), "%d%s %d%s %f%s", &deDegPart, dummy, &deMinPart, dummy, &deSecPart, dummy);

    float deDeg;
    int deFormatDeg = sscanf(deStr.c_str(), "%f", &deDeg);

    if (deFormatParts == 6) {
        deRadians = DegToRad(deDegPart + (float)deMinPart/60 + (float)deSecPart/60/60);
    } else if (deFormatDeg == 1) {
        deRadians = DegToRad(deFormatDeg);
    } else {
        std::cerr << "Invalid declination format." << std::endl;
        exit(1);
    }

    // find the star

    float tolerance = 0.001;
    Vec3 userSpatial = SphericalToSpatial(raRadians, deRadians);
    int i = 0;
    for (const CatalogStar &curStar : catalog) {
        if ((curStar.spatial - userSpatial).Magnitude() < tolerance) {
            std::cout << "found_star_" << i << " "  << curStar.name << std::endl;
            std::cout << "fonud_star_magnitude_" << i << " " << curStar.magnitude << std::endl;
            i++;
        }
    }
    if (i == 0) {
        std::cerr << "No stars found" << std::endl;
    }
}

void InspectPrintStar(const Catalog &catalog) {
    auto stars = PromptCatalogStars(catalog, 1);
    float ra, de;
    SpatialToSpherical(stars[0]->spatial, &ra, &de);

    std::cout << "star_ra " << RadToDeg(ra) << std::endl;
    std::cout << "star_de " << RadToDeg(de) << std::endl;
}

void InspectCatalog() {
    InteractiveChoice<CatalogInspector> inspectorChoice;
    inspectorChoice.Register("pair_distance", "pair distance angle", InspectPairDistance);
    inspectorChoice.Register("pyramid_distances", "all pair distances in pyramid", InspectPyramidDistances);
    inspectorChoice.Register("triple_angle", "inner angle of a triangle", InspectTripleAngle);
    inspectorChoice.Register("find_star", "find a star name based on ra/de", InspectFindStar);
    inspectorChoice.Register("print_star", "print coordinates of a star", InspectPrintStar);
    (*inspectorChoice.Prompt("Inspect the catalog"))(CatalogRead());
}

} // namespace lost
