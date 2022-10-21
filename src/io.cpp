#include "io.hpp"

#include <cairo/cairo.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <cstring>
#include <random>
#include <algorithm>

#include "attitude-estimators.hpp"
#include "attitude-utils.hpp"
#include "databases.hpp"
#include "star-id.hpp"
#include "star-utils.hpp"

namespace lost {

/// Create a PromptedOutputStream which will output to the given file.
PromptedOutputStream::PromptedOutputStream(std::string filePath) {
    if (isatty(fileno(stdout)) && (filePath == "stdout" || filePath == "-")) {
        std::cerr << "WARNING: output contains binary contents. Not printed to terminal." << std::endl;
        filePath = "/dev/null";
    }

    if (filePath == "stdout" || filePath == "-") {
        stream = &std::cout;
    } else {
        std::fstream *fs = new std::fstream();
        fs->open(filePath, std::fstream::out);
        stream = fs;
        
    }
}

/// Parse the bright star catalog from the TSV file on disk.
std::vector<CatalogStar> BscParse(std::string tsvPath) {
    std::vector<CatalogStar> result;
    FILE *file;
    double raj2000, dej2000;
    int magnitudeHigh, magnitudeLow, name;
    char weird; // dead variable

    file = fopen(tsvPath.c_str(), "r");
    if (file == NULL) {
        printf("Error opening file: %s\n", strerror(errno));
        exit(1);
        return result;
    }

    while (EOF != fscanf(file, "%lf|%lf|%d|%c|%d.%d",
                         &raj2000, &dej2000,
                         &name, &weird,
                         &magnitudeHigh, &magnitudeLow)) {
        result.push_back(CatalogStar(DegToRad(raj2000),
                                     DegToRad(dej2000),
                                     magnitudeHigh*100 + abs(magnitudeLow),
                                     name));
    }

    fclose(file);
    return result;
}

#ifndef DEFAULT_BSC_PATH
#define DEFAULT_BSC_PATH "bright-star-catalog.tsv"
#endif

/// Read and parse the full catalog from disk. If called multiple times, will re-use the first result.
std::vector<CatalogStar> &CatalogRead() {
    static std::vector<CatalogStar> catalog;

    // An apporximation of whether its been initialized, but its accurate in a practical sense
    if (catalog.size() == 0) {
        char *tsvPath = getenv("LOST_BSC_PATH");
        catalog = BscParse(tsvPath ? tsvPath : DEFAULT_BSC_PATH);
    }
    return catalog;
}

/// Convert a colored Cairo image surface into a row-major array of grayscale pixels
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
            (pixel     &0xFF) *0.07);
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

/**
 * Plot information about an image onto the image using Cairo.
 * Puts a box around each centroid, writes the attitude in the top left, etc.
 * @param red,green,blue The color to use when annotating the image. 0 represents none of that color, 1 represents that color in full.
 * @param alpha The transparency of annotations. 0 is completely transparent, 1 is completely opaque.
 * @param rawStarIndexes If true, print the catalog index. This is in contrast with the default behavior, which is to print the name of the catalog star.
 */
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

typedef StarIdAlgorithm *(*StarIdAlgorithmFactory)();

typedef AttitudeEstimationAlgorithm *(*AttitudeEstimationAlgorithmFactory)();

/// Add a pair-distance KVector database to the given builder.
void BuildPairDistanceKVectorDatabase(MultiDatabaseBuilder *builder, const Catalog &catalog, float minDistance, float maxDistance, long numBins) {
    
    std::vector<lost::KVectorPair> pairs = CatalogToPairDistances(catalog, minDistance, maxDistance);
    long length = SerializeLengthPairDistanceKVector(pairs.size(), numBins);
    unsigned char *buffer = builder->AddSubDatabase(PairDistanceKVectorDatabase::kMagicValue, length);
    if (buffer == NULL) {
        std::cerr << "No room for another database." << std::endl;
    }
    SerializePairDistanceKVector(pairs, minDistance, maxDistance, numBins, buffer);
}

/// Generate and add databases to the given multidatabase builder according to the command line options in `values`
void GenerateDatabases(MultiDatabaseBuilder *builder, const Catalog &catalog, const DatabaseOptions &values) {

    if (values.kvector) {
        float minDistance = DegToRad(values.kvectorMinDistance);
        float maxDistance = DegToRad(values.kvectorMaxDistance);
        long numBins = values.kvectorNumDistanceBins;
        BuildPairDistanceKVectorDatabase(builder, catalog, minDistance, maxDistance, numBins);
    } else {
        std::cerr << "No database builder selected -- no database generated." << std::endl;
        exit(1);
    }

}

/// Print information about the camera in machine and human-readable form.
std::ostream &operator<<(std::ostream &os, const Camera &camera) {
    os << "camera_focal_length " << camera.FocalLength() << std::endl
       << "camera_fov " << camera.Fov() << std::endl
       << "camera_resolution_x " << camera.XResolution() << std::endl
       << "camera_resolution_y " << camera.YResolution() << std::endl;
    // TODO: principal point
    return os;
}

// PIPELINE INPUT STUFF

/**
 * Calculate the focal length, in pixels, based on the given command line options.
 * This function exists because there are two ways to specify how "zoomed-in" the camera is. One way is using just FOV, which is useful when generating false images. Another is a combination of pixel size and focal length, which is useful for physical cameras.
 */
float FocalLengthFromOptions(const PipelineOptions &values) {
    if (values.pixelSize == -1) {
        return FovToFocalLength(DegToRad(values.fov), values.generateXRes);
    } else {
        return values.focalLength * 1000 / values.pixelSize;
    }
}

/// Convert the result of InputImage() into a cairo surface.
cairo_surface_t *PipelineInput::InputImageSurface() const {
    const Image *inputImage = InputImage();
    return GrayscaleImageToSurface(inputImage->image, inputImage->width, inputImage->height);
}

// /**
//  * Pipeline input for an image that's already been identified using Astrometry.net.
//  * @todo Not implemented yet.
//  */
// class AstrometryPipelineInput : public PipelineInput {
// public:
//     explicit AstrometryPipelineInput(const std::string &path);

//     const Image *InputImage() const { return &image; };
//     const Attitude *InputAttitude() const { return &attitude; };
// private:
//     Image image;
//     Attitude attitude;
// };

/**
 * A pipeline input coming from an image with no extra metadata. Only InputImage will be available.
 * @param cairoSurface A cairo surface from the image file.
 * @todo should rename, not specific to PNG.
 */
PngPipelineInput::PngPipelineInput(cairo_surface_t *cairoSurface, Camera *camera, const Catalog &catalog)
    : camera(*camera), catalog(catalog) {

    image.image = SurfaceToGrayscaleImage(cairoSurface);
    image.width = cairo_image_surface_get_width(cairoSurface);
    image.height = cairo_image_surface_get_height(cairoSurface);
}

/// Create a PngPipelineInput using command line options.
PipelineInputList GetPngPipelineInput(const PipelineOptions &values) {
    // I'm not sure why, but i can't get an initializer list to work here. Probably something to do
    // with copying unique ptrs
    PipelineInputList result;
    cairo_surface_t *cairoSurface = NULL;
    std::string pngPath = values.png;

    cairoSurface = cairo_image_surface_create_from_png(pngPath.c_str());
    std::cerr << "PNG Read status: " << cairo_status_to_string(cairo_surface_status(cairoSurface)) << std::endl;
    if (cairoSurface == NULL || cairo_surface_status(cairoSurface) != CAIRO_STATUS_SUCCESS) {
        exit(1);
    }

    int xResolution = cairo_image_surface_get_width(cairoSurface);
    int yResolution = cairo_image_surface_get_height(cairoSurface);
    float focalLengthPixels = FocalLengthFromOptions(values);
    Camera cam = Camera(focalLengthPixels, xResolution, yResolution);

    result.push_back(std::unique_ptr<PipelineInput>(new PngPipelineInput(cairoSurface, &cam, CatalogRead())));
    return result;
}

// AstrometryPipelineInput::AstrometryPipelineInput(const std::string &path) {
//     // create from path, TODO
// }

// PipelineInputList PromptAstrometryPipelineInput() {
//     // TODO: why does it let us do a reference to an ephemeral return value?
//     std::string path = Prompt<std::string>("Astrometry download directory");
//     PipelineInputList result;
//     result.push_back(std::unique_ptr<PipelineInput>(new AstrometryPipelineInput(path)));
//     return result;
// }

/// A star used in simulated image generation. Contains extra data about how to simulate the star.
class GeneratedStar : public Star {
public:
    GeneratedStar(Star star, float peakBrightness, Vec2 motionBlurDelta)
        : Star(star), peakBrightness(peakBrightness), delta(motionBlurDelta) { };

    /// the brightness density at the center of the star. 0.0 is black, 1.0 is white.
    float peakBrightness;

    /// (only meaningful with motion blur) Where the star will appear one time unit in the future.
    Vec2 delta;
};

// In the equations for pixel brightness both with motion blur enabled and disabled, we don't need
// any constant scaling factor outside the integral because when d0=0, the brightness at the center
// will be zero without any scaling. The scaling factor you usually see on a Normal distribution is
// so that the Normal distribution integrates to one over the real line, making it a probability
// distribution. But we want the /peak/ to be one, not the integral. motion blur enabled

/**
 * Calculates the indefinite integral of brightness density at a point due to a single star.
 * When oversampling is disabled, this is called only at pixel centers. When oversampling is enabled, it's called at multiple points in each pixel and then averaged. If multiple stars are near each other, the brightnesses can just be added then clamped.
 * See https://wiki.huskysat.org/wiki/index.php/Motion_Blur_and_Rolling_Shutter#Motion_Blur_Math to learn how these equations were derived.
 * @param pixel The point to calculate brightness density at. If only calculating per-pixel, should be the center of the pixel. ie, Vec2(10.5,5.5) would be appropriate for the pixel (10,5)
 * @param generatedStar the star to calculate brightness based on.
 * @param t Since this function computes the indefinite integral, this is the value it's evaluated at. To calculate the definite integral, which is what you want, call this function twice with different values for the `t` parameter then find the difference.
 * @param stddev The standard deviation of spread of the star. Higher values make stars more spread out. See command line documentation.
 * @return Indefinite integral of brightness density.
 */
static float motionBlurredPixelBrightness(const Vec2 &pixel, const GeneratedStar &generatedStar,
                                          float t, float stddev) {
    const Vec2 &p0 = generatedStar.position;
    const Vec2 &delta = generatedStar.delta;
    const Vec2 d0 = p0 - pixel;
    return generatedStar.peakBrightness
        * stddev*sqrt(M_PI) / (sqrt(2)*delta.Magnitude())
        * exp(pow(d0.x*delta.x + d0.y*delta.y, 2) / (2*stddev*stddev*delta.MagnitudeSq())
              - d0.MagnitudeSq() / (2*stddev*stddev))
        * erf((t*delta.MagnitudeSq() + d0.x*delta.x + d0.y*delta.y) / (stddev*sqrt(2)*delta.Magnitude()));
}

/// Like motinoBlurredPixelBrightness, but for when motion blur is disabled.
static float staticPixelBrightness(const Vec2 &pixel, const GeneratedStar &generatedStar,
                                   float stddev) {
    const Vec2 d0 = generatedStar.position - pixel;
    return generatedStar.peakBrightness * exp(-d0.MagnitudeSq() / (2*stddev*stddev));
}

const int kMaxBrightness = 255;

/**
 * Create a generated pipeline input.
 * The parameters correspond directly to command line options. See the command line documentation for more details. This constructor performs the actual image generation.
 */
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
                                               int oversampling,
                                               int numFalseStars,
                                               int falseStarMinMagnitude,
                                               int falseStarMaxMagnitude,
                                               int seed)
    : camera(camera), attitude(attitude), catalog(catalog) {

    assert(falseStarMaxMagnitude <= falseStarMinMagnitude);
    std::default_random_engine rng(seed);

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

    bool motionBlurEnabled = exposureTime > 0 && abs(motionBlurDirection.GetQuaternion().Angle()) > 0.001;
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

    Catalog catalogWithFalse = catalog;

    std::uniform_real_distribution<float> uniformDistribution(0.0, 1.0);
    std::uniform_int_distribution<int> magnitudeDistribution(falseStarMaxMagnitude, falseStarMinMagnitude);
    for (int i = 0; i < numFalseStars; i++) {
        float ra = uniformDistribution(rng) * 2*M_PI;
        // to be uniform around sphere. Borel-Kolmogorov paradox is calling
        float de = asin(uniformDistribution(rng)*2 - 1);
        float magnitude = magnitudeDistribution(rng);

        catalogWithFalse.push_back(CatalogStar(ra, de, magnitude, -1));
    }

    for (int i = 0; i < (int)catalogWithFalse.size(); i++) {
        bool isTrueStar = i < (int)catalog.size();

        const CatalogStar &catalogStar = catalogWithFalse[i];
        Vec3 rotated = attitude.Rotate(catalogWithFalse[i].spatial);
        if (rotated.x <= 0) {
            continue;
        }
        Vec2 camCoords = camera.SpatialToCamera(rotated);

        if (camera.InSensor(camCoords)) {
            Vec3 futureSpatial = futureAttitude.Rotate(catalogWithFalse[i].spatial);
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
            generatedStars.push_back(GeneratedStar(star, peakBrightness, delta));

            // don't add false stars to centroids or star ids
            if (isTrueStar) {
                stars.push_back(star);
                starIds.push_back(StarIdentifier(stars.size() - 1, i));
            }
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
        float oversamplingBrightnessFactor = oversamplingPerAxis*oversamplingPerAxis;

        // the star.x and star.y refer to the pixel whose top left corner the star should appear at
        // (and fractional amounts are relative to the corner). When we color a pixel, we ideally
        // would integrate the intensity of the star over that pixel, but we can make do by sampling
        // the intensity of the star at the /center/ of the pixel, ie, star.x+.5 and star.y+.5
        for (int xPixel = xMin; xPixel <= xMax; xPixel++) {
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
                                (motionBlurredPixelBrightness({x, y}, star, tEnd, starSpreadStdDev)
                                 - motionBlurredPixelBrightness({x, y}, star, tStart, starSpreadStdDev))
                                / oversamplingBrightnessFactor;
                        } else {
                            curPhotons = staticPixelBrightness({x, y}, star, starSpreadStdDev)
                                / oversamplingBrightnessFactor;
                        }

                        assert(0.0 <= curPhotons);

                        photonsBuffer[xPixel + yPixel*image.width] += curPhotons;
                    }
                }
            }
        }
    }

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

/// Create a GeneratedPipelineInput based on the command line options in `values`
PipelineInputList GetGeneratedPipelineInput(const PipelineOptions &values) {
    // TODO: prompt for attitude, imagewidth, etc and then construct a GeneratedPipelineInput

    // TODO: allow random angle generation?
    Attitude attitude = Attitude(SphericalToQuaternion(DegToRad(values.generateRa),
                                                         DegToRad(values.generateDe),
                                                         DegToRad(values.generateRoll)));
    Attitude motionBlurDirection = Attitude(SphericalToQuaternion(DegToRad(values.generateBlurRa),
                                                                  DegToRad(values.generateBlurDe),
                                                                  DegToRad(values.generateBlurRoll)));

    PipelineInputList result;

    float focalLength = FocalLengthFromOptions(values);

    for (int i = 0; i < values.generate; i++) {
        GeneratedPipelineInput *curr = new GeneratedPipelineInput(
            CatalogRead(),
            attitude,
            Camera(focalLength, values.generateXRes, values.generateYRes),
            values.generateRefBrightness,
            values.generateSpreadStdDev,
            values.generateSensitivity,
            values.generateDarkCurrent,
            values.generateReadNoiseStdDev,
            motionBlurDirection,
            values.generateExposure,
            values.generateReadoutTime,
            values.generateShotNoise,
            values.generateOversampling,
            values.generateNumFalseStars,
            values.generateFalseMinMag,
            values.generateFalseMaxMag,
            values.generateSeed);

        result.push_back(std::unique_ptr<PipelineInput>(curr));
    }

    return result;
}

typedef PipelineInputList (*PipelineInputFactory)();

/// Come up with a list of pipeline inputs based on command line options.
PipelineInputList GetPipelineInput(const PipelineOptions &values) {

    if (values.png != "") {
        return GetPngPipelineInput(values);
    } else {
        return GetGeneratedPipelineInput(values);
    }
}

/**
 * Construct a pipeline using the given algorithms, some of which may be null.
 * @param database A pointer to the raw bytes of the database the star ID algorithm expects. If the database is NULL or not the type of database the star ID algorithm expects (almost always a multi-database), you'll get an error trying to identify stars later.
 */
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


/// Create a pipeline from command line options.
Pipeline SetPipeline(const PipelineOptions &values) {

    Pipeline result;

    // TODO: more flexible or sth
    // TODO: don't allow setting star-id until database is set, and perhaps limit the star-id
    // choices to those compatible with the database?

    // centroid algorithm stage
    if (values.centroidAlgo == "dummy") {
        result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(new DummyCentroidAlgorithm(values.centroidDummyNumStars));
    } else if (values.centroidAlgo == "cog") {
        result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(new CenterOfGravityAlgorithm());
    } else if (values.centroidAlgo == "iwcog") {
        result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(new IterativeWeightedCenterOfGravityAlgorithm());
    } else if (values.centroidAlgo != "") {
        std::cout << "Illegal centroid algorithm." << std::endl;
        exit(1);
    }

    // centroid magnitude filter stage
    if (values.centroidMagFilter != -1) result.centroidMinMagnitude = values.centroidMagFilter;

    // database stage
    if (values.databasePath != "") {
        std::fstream fs;
        fs.open(values.databasePath, std::fstream::in | std::fstream::binary);
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
        result.starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(new GeometricVotingStarIdAlgorithm(DegToRad(values.angularTolerance)));
    } else if (values.idAlgo == "py") {
        result.starIdAlgorithm = std::unique_ptr<StarIdAlgorithm>(new PyramidStarIdAlgorithm(DegToRad(values.angularTolerance), values.estimatedNumFalseStars, values.maxMismatchProb, 1000));
    } else if (values.idAlgo != "") {
        std::cout << "Illegal id algorithm." << std::endl;
        exit(1);
    }

    // attitude estimation stage
    if (values.attitudeAlgo == "dqm") {
        result.attitudeEstimationAlgorithm = std::unique_ptr<AttitudeEstimationAlgorithm>(new DavenportQAlgorithm());
    } else if (values.attitudeAlgo == "triad") {
        result.attitudeEstimationAlgorithm = std::unique_ptr<AttitudeEstimationAlgorithm>(new TriadAlgorithm());
    } else if (values.attitudeAlgo == "quest") {
        result.attitudeEstimationAlgorithm = std::unique_ptr<AttitudeEstimationAlgorithm>(new QuestAlgorithm());
    } else if (values.attitudeAlgo != "") {
        std::cout << "Illegal attitude algorithm." << std::endl;
        exit(1);
    }

    return result;
}

/**
 * Run all stages of a pipeline. This is the "main" method for pipelines.
 * In space (or when using an image file as input), the PipelineInput will contain only an InputImage. In this case, `Go` runs each star tracking algorithm in turn, passing the result of each step into the next one.
 * When running on a generated image (or any pipeline input where methods other than InputImage are available), or using a Pipeline where some algorithms are not set, the behavior is more nuanced. Each algorithm will be run on the return value of the corresponding input method from the PipelineInput object, unless an earlier algorithm in the Pipeline returned a result, in which case that intermediate value is used instead of the value from the PipelineInput.
 */
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
        result.attitude = std::unique_ptr<Attitude>(
            new Attitude(attitudeEstimationAlgorithm->Go(*input.InputCamera(), *inputStars, result.catalog, *inputStarIds)));
    }

    return result;
}

/// Convenience function to run the main `Pipeline::Go` function on each input
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

/**
 * The result of comparing actual and expected centroids
 * Used for debugging and benchmarking
 */
class CentroidComparison {
public:
    CentroidComparison() : meanError(0.0f), numExtraStars(0.0), numMissingStars(0.0) {};
    /**
     * Average distance from actual to expected centroids (in pixels)
     * Only correct centroids are considered in this average.
     */
    float meanError;

    /**
     * Stars in actual but not expected. Ideally 0
     * This is a float because we may average multiple centroid comparisons together.
     */
    float numExtraStars;

    /**
     * Stars in expected but not actual. Ideally 0
     * This is a float because we may average multiple centroid comparisons together.
     */
    float numMissingStars;
};

// helper for CentroidsCompare
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

/**
 * Compare expected and actual centroids.
 * Useful for debugging and benchmarking.
 * @param threshold The maximum number of pixels apart two centroids can be to be considered the same.
 */
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

/**
 * Compare expected and actual star identifications.
 * Useful for debugging and benchmarking.
 */
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
                                   const std::vector<PipelineOutput> &,
                                   const PipelineOptions &);

/// Plotter suitable for `cairo_surface_write_to_png_stream` which simply writes to an std::ostream
cairo_status_t OstreamPlotter(void *closure, const unsigned char *data, unsigned int length) {
    std::ostream *os = (std::ostream *)closure;
    os->write((const char *)data, length);
    return CAIRO_STATUS_SUCCESS;
}

/// Plots the input image with no annotation to `os`
void PipelineComparatorPlotRawInput(std::ostream &os,
                                    const PipelineInputList &expected,
                                    const std::vector<PipelineOutput> &,
                                    const PipelineOptions &) {

    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Plot the annotated input image to `os`
// TODO: should probably use Expected methods, not Input methods, because future PipelineInputs could add noise to the result of the Input methods.
void PipelineComparatorPlotInput(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &,
                                 const PipelineOptions &) {
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

/// Compare the actual and expected centroids, printing key stats to `os`
void PipelineComparatorCentroids(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &actual,
                                 const PipelineOptions &values) {
    int size = (int)expected.size();

    float threshold = values.centroidCompareThreshold;

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

/// Print a list of centroids to `os`
void PipelineComparatorPrintCentroids(std::ostream &os,
                                      const PipelineInputList &expected,
                                      const std::vector<PipelineOutput> &actual,
                                      const PipelineOptions &) {
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

// TODO: add a CLI option to use this!
/// Plot an annotated image where centroids are annotated with their centroid index. For debugging.
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

/// Plot the image annotated with output data computed by the star tracking algorithms.
void PipelineComparatorPlotOutput(std::ostream &os,
                                  const PipelineInputList &expected,
                                  const std::vector<PipelineOutput> &actual,
                                  const PipelineOptions &) {
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

/// Compare the expected and actual star identifiers.
void PipelineComparatorStarIds(std::ostream &os,
                             const PipelineInputList &expected,
                             const std::vector<PipelineOutput> &actual,
                             const PipelineOptions &values) {
    float centroidThreshold = actual[0].stars
        ? values.centroidCompareThreshold
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

/// Print the identifed attitude to `os` in Euler angle format.
void PipelineComparatorPrintAttitude(std::ostream &os,
                                     const PipelineInputList &,
                                     const std::vector<PipelineOutput> &actual,
                                     const PipelineOptions &) {
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

/// Compare the actual and expected attitudes.
void PipelineComparatorAttitude(std::ostream &os,
                                const PipelineInputList &expected,
                                const std::vector<PipelineOutput> &actual,
                                const PipelineOptions &values) {

    // TODO: use Wahba loss function (maybe average per star) instead of just angle. Also break
    // apart roll error from boresight error. This is just quick and dirty for testing
    float angleThreshold = DegToRad(values.attitudeCompareThreshold);

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

// TODO: add these debug comparators back in!
// void PipelineComparatorPrintPairDistance(std::ostream &os,
//                                          const PipelineInputList &expected,
//                                          const std::vector<PipelineOutput> &actual) {
//     int index1 = Prompt<int>("Index of first star");
//     int index2 = Prompt<int>("Index of second star");

//     const Camera &camera = *expected[0]->InputCamera();
//     const Stars &stars = *actual[0].stars;

//     assert(index1 >= 0 && index2 >= 0 && index1 < (int)stars.size() && index2 < (int)stars.size());
//     os << "pair_distance " << Angle(camera.CameraToSpatial(stars[index1].position),
//                                     camera.CameraToSpatial(stars[index2].position))
//        << std::endl;
// }

// void PipelineComparatorPrintPyramidDistances(std::ostream &os,
//                                              const PipelineInputList &expected,
//                                              const std::vector<PipelineOutput> &actual) {
//     int index1 = Prompt<int>("Catalog name/index of first star");
//     int index2 = Prompt<int>("Catalog name/index of second star");
//     int index3 = Prompt<int>("Catalog name/index of third star");
//     int index4 = Prompt<int>("Catalog name/index of four star");

//     const Camera &camera = *expected[0]->InputCamera();
//     const Stars &stars = *actual[0].stars;

//     Vec3 spatial1 = camera.CameraToSpatial(stars[index1].position);
//     Vec3 spatial2 = camera.CameraToSpatial(stars[index2].position);
//     Vec3 spatial3 = camera.CameraToSpatial(stars[index3].position);
//     Vec3 spatial4 = camera.CameraToSpatial(stars[index4].position);

//     std::cout << "pair_distance_12 " << Angle(spatial1, spatial2) << std::endl;
//     std::cout << "pair_distance_13 " << Angle(spatial1, spatial3) << std::endl;
//     std::cout << "pair_distance_14 " << Angle(spatial1, spatial4) << std::endl;
//     std::cout << "pair_distance_23 " << Angle(spatial2, spatial3) << std::endl;
//     std::cout << "pair_distance_24 " << Angle(spatial2, spatial4) << std::endl;
//     std::cout << "pair_distance_34 " << Angle(spatial3, spatial4) << std::endl;
// }

// void PipelineComparatorPrintTripleAngle(std::ostream &os,
//                                         const PipelineInputList &expected,
//                                         const std::vector<PipelineOutput> &actual) {
//     int index1 = Prompt<int>("Index of first star");
//     int index2 = Prompt<int>("Index of second star");
//     int index3 = Prompt<int>("Index of third star");

//     const Camera &camera = *expected[0]->InputCamera();
//     const Stars &stars = *actual[0].stars;

//     assert(index1 >= 0 && index1 < (int)stars.size());
//     assert(index2 >= 0 && index2 < (int)stars.size());
//     assert(index3 >= 0 && index3 < (int)stars.size());

//     // TODO, when merging with nondimensional branch
// }

/**
 * Print or otherwise analyze the results of (perhaps multiple) runs of a star tracking pipeline.
 * Uses the command line options in `values` to determine which analyses to run. Examples include plotting an annotated output image to a png file, comparing the actual and expected centroids, etc
 */
void PipelineComparison(const PipelineInputList &expected,
                        const std::vector<PipelineOutput> &actual,
                        const PipelineOptions &values) {
    assert(expected.size() == actual.size() && expected.size() > 0);

    // TODO: Remove the asserts and print out more reasonable error messages.

#define LOST_PIPELINE_COMPARE(comparator, path) do {            \
        PromptedOutputStream pos(path);                         \
        comparator(pos.Stream(), expected, actual, values);     \
    } while (0)

    if (values.plotRawInput != "") {
        assert(expected[0]->InputImage() && expected.size() == 1);
        LOST_PIPELINE_COMPARE(PipelineComparatorPlotRawInput, values.plotRawInput);
    }

    if (values.plotInput != "") {
        assert(expected[0]->InputImage() && expected.size() == 1 && expected[0]->InputStars());
        LOST_PIPELINE_COMPARE(PipelineComparatorPlotInput, values.plotInput);
    }
    if (values.plotOutput != "") {
        assert(actual.size() == 1 && (actual[0].stars || actual[0].starIds));
        LOST_PIPELINE_COMPARE(PipelineComparatorPlotOutput, values.plotOutput);
    }
    if (values.printCentroids != "") {
        assert(actual[0].stars && actual.size() == 1);
        LOST_PIPELINE_COMPARE(PipelineComparatorPrintCentroids, values.printCentroids);
    }
    if (values.compareCentroids != "") {
        assert(actual[0].stars && expected[0]->ExpectedStars() && values.centroidCompareThreshold);
        LOST_PIPELINE_COMPARE(PipelineComparatorCentroids, values.compareCentroids);
    }
    if (values.compareStarIds != "") {
        assert(expected[0]->ExpectedStarIds() && actual[0].starIds);
        LOST_PIPELINE_COMPARE(PipelineComparatorStarIds, values.compareStarIds);
    }
    if (values.printAttitude != "") {
        assert(actual[0].attitude && actual.size() == 1);
        LOST_PIPELINE_COMPARE(PipelineComparatorPrintAttitude, values.printAttitude);
    }
    if (values.compareAttitudes != "") {
        assert(actual[0].attitude && expected[0]->ExpectedAttitude() && values.attitudeCompareThreshold);
        LOST_PIPELINE_COMPARE(PipelineComparatorAttitude, values.compareAttitudes);
    }

#undef LOST_PIPELINE_COMPARE
}

// TODO: Add CLI options for all the inspectors!

// typedef void (*CatalogInspector)(const Catalog &);

// static std::vector<const CatalogStar *> PromptCatalogStars(const Catalog &catalog, int howMany) {
//     std::vector<const CatalogStar *> result;
//     for (int i = 0; i < howMany; i++) {
//         int name = Prompt<int>("Catalog name of " + std::to_string(i) + "-th star");
//         const CatalogStar *star = findNamedStar(catalog, name);
//         if (star == NULL) {
//             std::cerr << "Star not found!" << std::endl;
//             exit(1);
//         }
//         result.push_back(star);
//     }
//     return result;
// }

// void InspectPairDistance(const Catalog &catalog) {
//     auto stars = PromptCatalogStars(catalog, 2);

//     // TODO: not cout, prompt for an ostream in inspect and pass argument
//     std::cout << Angle(stars[0]->spatial, stars[1]->spatial) << std::endl;
// }

// void InspectPyramidDistances(const Catalog &catalog) {
//     auto stars = PromptCatalogStars(catalog, 4);

//     std::cout << "pair_distance_01 " << Angle(stars[0]->spatial, stars[1]->spatial) << std::endl;
//     std::cout << "pair_distance_02 " << Angle(stars[0]->spatial, stars[2]->spatial) << std::endl;
//     std::cout << "pair_distance_03 " << Angle(stars[0]->spatial, stars[3]->spatial) << std::endl;
//     std::cout << "pair_distance_12 " << Angle(stars[1]->spatial, stars[2]->spatial) << std::endl;
//     std::cout << "pair_distance_13 " << Angle(stars[1]->spatial, stars[3]->spatial) << std::endl;
//     std::cout << "pair_distance_23 " << Angle(stars[2]->spatial, stars[3]->spatial) << std::endl;
// }

// void InspectTripleAngle(const Catalog &catalog) {
//     auto stars = PromptCatalogStars(catalog, 3);

//     // TODO
// }

// void InspectFindStar(const Catalog &catalog) {
//     std::string raStr = PromptLine("Right Ascension");

//     float raRadians;

//     int raHours, raMinutes;
//     float raSeconds;
//     int raFormatTime = sscanf(raStr.c_str(), "%dh %dm %fs", &raHours, &raMinutes, &raSeconds);

//     float raDeg;
//     int raFormatDeg = sscanf(raStr.c_str(), "%f", &raDeg);

//     if (raFormatTime == 3) {
//         raRadians = (raHours * 2*M_PI/24) + (raMinutes * 2*M_PI/24/60) + (raSeconds * 2*M_PI/24/60/60);
//     } else if (raFormatDeg == 1) {
//         raRadians = DegToRad(raFormatDeg);
//     } else {
//         std::cerr << "Invalid right ascension format. Do \"09h 38m 29.8754s\" or a number of degrees." << std::endl;
//         exit(1);
//     }

//     std::string deStr = PromptLine("Declination");

//     float deRadians;

//     int deDegPart, deMinPart;
//     float deSecPart;
//     char dummy[8];
//     int deFormatParts = sscanf(deStr.c_str(), "%d%s %d%s %f%s", &deDegPart, dummy, &deMinPart, dummy, &deSecPart, dummy);

//     float deDeg;
//     int deFormatDeg = sscanf(deStr.c_str(), "%f", &deDeg);

//     if (deFormatParts == 6) {
//         deRadians = DegToRad(deDegPart + (float)deMinPart/60 + (float)deSecPart/60/60);
//     } else if (deFormatDeg == 1) {
//         deRadians = DegToRad(deFormatDeg);
//     } else {
//         std::cerr << "Invalid declination format." << std::endl;
//         exit(1);
//     }

//     // find the star

//     float tolerance = 0.001;
//     Vec3 userSpatial = SphericalToSpatial(raRadians, deRadians);
//     int i = 0;
//     for (const CatalogStar &curStar : catalog) {
//         if ((curStar.spatial - userSpatial).Magnitude() < tolerance) {
//             std::cout << "found_star_" << i << " "  << curStar.name << std::endl;
//             std::cout << "fonud_star_magnitude_" << i << " " << curStar.magnitude << std::endl;
//             i++;
//         }
//     }
//     if (i == 0) {
//         std::cerr << "No stars found" << std::endl;
//     }
// }

// void InspectPrintStar(const Catalog &catalog) {
//     auto stars = PromptCatalogStars(catalog, 1);
//     float ra, de;
//     SpatialToSpherical(stars[0]->spatial, &ra, &de);

//     std::cout << "star_ra " << RadToDeg(ra) << std::endl;
//     std::cout << "star_de " << RadToDeg(de) << std::endl;
// }

// void InspectCatalog() {
//     InteractiveChoice<CatalogInspector> inspectorChoice;
//     inspectorChoice.Register("pair_distance", "pair distance angle", InspectPairDistance);
//     inspectorChoice.Register("pyramid_distances", "all pair distances in pyramid", InspectPyramidDistances);
//     inspectorChoice.Register("triple_angle", "inner angle of a triangle", InspectTripleAngle);
//     inspectorChoice.Register("find_star", "find a star name based on ra/de", InspectFindStar);
//     inspectorChoice.Register("print_star", "print coordinates of a star", InspectPrintStar);
//     (*inspectorChoice.Prompt("Inspect the catalog"))(CatalogRead());
// }

} // namespace lost
