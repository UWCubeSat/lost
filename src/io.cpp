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
#include <map>
#include <chrono>

#include "attitude-estimators.hpp"
#include "attitude-utils.hpp"
#include "databases.hpp"
#include "star-id.hpp"
#include "star-utils.hpp"

namespace lost {

/// Create a PromptedOutputStream which will output to the given file.
UserSpecifiedOutputStream::UserSpecifiedOutputStream(std::string filePath, bool isBinary) {
    if (isBinary && isatty(fileno(stdout)) && (filePath == "stdout" || filePath == "-")) {
        std::cerr << "WARNING: output contains binary contents. Not printed to terminal." << std::endl;
        filePath = "/dev/null";
    }

    if (filePath == "stdout" || filePath == "-") {
        stream = &std::cout;
        isFstream = false;
    } else {
        std::fstream *fs = new std::fstream();
        fs->open(filePath, std::fstream::out);
        stream = fs;
        isFstream = true;
    }
}

UserSpecifiedOutputStream::~UserSpecifiedOutputStream() {
    if (isFstream) {
        delete stream;
    }
}

/// Parse the bright star catalog from the TSV file on disk.
std::vector<CatalogStar> BscParse(std::string tsvPath) {
    std::vector<CatalogStar> result;
    FILE *file;
    double raj2000, dej2000;
    int magnitudeHigh, magnitudeLow, name;
    char weird;

    file = fopen(tsvPath.c_str(), "r");
    if (file == NULL) {
        printf("Error opening file: %s\n", strerror(errno));
        exit(1); // TODO: do we want any other error handling?
        return result;
    }

    while (EOF != fscanf(file, "%lf|%lf|%d|%c|%d.%d",
                         &raj2000, &dej2000,
                         &name, &weird,
                         &magnitudeHigh, &magnitudeLow)) {
        result.push_back(CatalogStar(DegToRad(raj2000),
                                     DegToRad(dej2000),
                                     magnitudeHigh*100 + (magnitudeHigh < 0 ? -magnitudeLow : magnitudeLow),
                                     name));
    }

    fclose(file);
    assert(result.size() > 9000); // basic sanity check
    return result;
}

#ifndef DEFAULT_BSC_PATH
#define DEFAULT_BSC_PATH "bright-star-catalog.tsv"
#endif

/// Read and parse the full catalog from disk. If called multiple times, will re-use the first result.
const Catalog &CatalogRead() {
    static bool readYet = false;
    static std::vector<CatalogStar> catalog;

    if (!readYet) {
        readYet = true;
        char *tsvPath = getenv("LOST_BSC_PATH");
        catalog = BscParse(tsvPath ? tsvPath : DEFAULT_BSC_PATH);
        // perform essential narrowing
        // remove all stars with exactly the same position as another, keeping the one with brighter magnitude
        std::sort(catalog.begin(), catalog.end(), [](const CatalogStar &a, const CatalogStar &b) {
            return a.spatial.x < b.spatial.x;
        });
        for (int i = catalog.size(); i > 0; i--) {
            if ((catalog[i].spatial - catalog[i-1].spatial).Magnitude() < 5e-5) { // 70 stars removed at this threshold.
                if (catalog[i].magnitude > catalog[i-1].magnitude) {
                    catalog.erase(catalog.begin() + i);
                } else {
                    catalog.erase(catalog.begin() + i - 1);
                }
            }
        }
    }
    return catalog;
}

/// Convert a colored Cairo image surface into a row-major array of grayscale pixels.
/// Result is allocated with new[]
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
void SurfacePlot(std::string description,
                 cairo_surface_t *cairoSurface,
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
    std::string metadata = description + " ";

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
    // TODO: calculating the length of the vector duplicates a lot of the work, slowing down
    // database generation
    long length = SerializeLengthPairDistanceKVector(catalog, minDistance, maxDistance, numBins);
    unsigned char *buffer = builder->AddSubDatabase(PairDistanceKVectorDatabase::kMagicValue, length);
    if (buffer == NULL) {
        std::cerr << "No room for another database." << std::endl;
    }
    SerializePairDistanceKVector(catalog, minDistance, maxDistance, numBins, buffer);

    // TODO: also parse it and print out some stats before returning
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
    if ((values.pixelSize != -1) ^ (values.focalLength != 0)) {
        std::cerr << "ERROR: Exactly one of --pixel-size or --focal-length were set." << std::endl;
        exit(1);
    }

    // no surefire way to see if the fov was set on command line, so we just check if it was changed from default.
    if (values.pixelSize != -1 && values.fov != 20) {
        std::cerr << "ERROR: Both focal length and FOV were provided. We only need one of the two methods (pixel size + focal length, or fov) to determine fov, please only provide one." << std::endl;
        exit(1);
    }

    if (values.pixelSize == -1) {
        return FovToFocalLength(DegToRad(values.fov), values.generateXRes);
    } else {
        return values.focalLength * 1000 / values.pixelSize;
    }
}

/// Convert the result of InputImage() into a cairo surface.
/// Allocates a new surface, whih must be destroyed with cairo_surface_destroy
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
 * No references to the surface are kept in the class and it may be freed after construction.
 * @param cairoSurface A cairo surface from the image file.
 * @todo should rename, not specific to PNG.
 */
PngPipelineInput::PngPipelineInput(cairo_surface_t *cairoSurface, Camera camera, const Catalog &catalog)
    : camera(camera), catalog(catalog) {

    image.image = SurfaceToGrayscaleImage(cairoSurface);
    image.width = cairo_image_surface_get_width(cairoSurface);
    image.height = cairo_image_surface_get_height(cairoSurface);
}

PngPipelineInput::~PngPipelineInput() {
    delete[] image.image;
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

    result.push_back(std::unique_ptr<PipelineInput>(new PngPipelineInput(cairoSurface, cam, CatalogRead())));
    cairo_surface_destroy(cairoSurface);
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
static float MotionBlurredPixelBrightness(const Vec2 &pixel, const GeneratedStar &generatedStar,
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

/// Like motionBlurredPixelBrightness, but for when motion blur is disabled.
static float StaticPixelBrightness(const Vec2 &pixel, const GeneratedStar &generatedStar,
                                   float stddev) {
    const Vec2 d0 = generatedStar.position - pixel;
    return generatedStar.peakBrightness * exp(-d0.MagnitudeSq() / (2*stddev*stddev));
}

/**
 * Compute how likely a star is to be imaged, given a "cutoff" magnitude that the camera can see half of.
 *
 * The theory is that there's a threshold of total light energy that must be received to image a
 * star. The main random factors are shot noise and read noise, but to simplify things so we don't
 * need to think about photons, we only focus on read noise, which we assume has a standard
 * deviation 1/5th of the cutoff brightness. We compute the probability that, taking read noise into
 * account, the observed energy would be less than the cutoff energy.
 */
static float CentroidImagingProbability(float mag, float cutoffMag) {
    float brightness = MagToBrightness(mag);
    float cutoffBrightness = MagToBrightness(cutoffMag);
    float stddev = cutoffBrightness/5.0;
    // CDF of Normal distribution with given mean and stddev
    return 1 - (0.5 * (1 + erf((cutoffBrightness-brightness)/(stddev*sqrt(2.0)))));
}

const int kMaxBrightness = 255;

/**
 * Create a generated pipeline input.
 * The parameters correspond directly to command line options. See the command line documentation for more details. This constructor performs the actual image generation.
 */
GeneratedPipelineInput::GeneratedPipelineInput(const Catalog &catalog,
                                               Attitude attitude,
                                               Camera camera,
                                               std::default_random_engine *rng,

                                               bool centroidsOnly,
                                               float zeroMagTotalPhotons,
                                               float starSpreadStdDev,
                                               float saturationPhotons,
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
                                               int cutoffMag,
                                               float perturbationStddev)
    : camera(camera), attitude(attitude), catalog(catalog) {

    assert(falseStarMaxMagnitude <= falseStarMinMagnitude);
    assert(perturbationStddev >= 0.0);

    image.width = camera.XResolution();
    image.height = camera.YResolution();
    // number of true photons each pixel receives.

    assert(oversampling >= 1);
    int oversamplingPerAxis = ceil(sqrt(oversampling));
    if (oversamplingPerAxis*oversamplingPerAxis != oversampling) {
        std::cerr << "WARNING: oversampling was not a perfect square. Rounding up to "
                  << oversamplingPerAxis*oversamplingPerAxis << "." << std::endl;
    }
    assert(exposureTime > 0);
    bool motionBlurEnabled = abs(motionBlurDirection.GetQuaternion().Angle()) > 0.001;
    Quaternion motionBlurDirectionQ = motionBlurDirection.GetQuaternion();
    // attitude at the middle of exposure time
    Quaternion currentAttitude = attitude.GetQuaternion();
    // attitude 1 time unit after middle of exposure
    Quaternion futureAttitude = motionBlurDirectionQ*currentAttitude;
    std::vector<GeneratedStar> generatedStars;

    // a star with 1 photon has peak density 1/(2pi sigma^2), because 2d gaussian formula. Then just
    // multiply up proportionally!
    float zeroMagPeakPhotonDensity = zeroMagTotalPhotons / (2*M_PI * starSpreadStdDev*starSpreadStdDev);

    // TODO: Is it 100% correct to just copy the standard deviation in both dimensions?
    std::normal_distribution<float> perturbation1DDistribution(0.0, perturbationStddev);

    Catalog catalogWithFalse = catalog;

    std::uniform_real_distribution<float> uniformDistribution(0.0, 1.0);
    std::uniform_int_distribution<int> magnitudeDistribution(falseStarMaxMagnitude, falseStarMinMagnitude);
    for (int i = 0; i < numFalseStars; i++) {
        float ra = uniformDistribution(*rng) * 2*M_PI;
        // to be uniform around sphere. Borel-Kolmogorov paradox is calling
        float de = asin(uniformDistribution(*rng)*2 - 1);
        float magnitude = magnitudeDistribution(*rng);

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
            float peakBrightnessPerTime = zeroMagPeakPhotonDensity * MagToBrightness(catalogStar.magnitude);
            float peakBrightness = peakBrightnessPerTime * exposureTime;
            float interestingThreshold = 0.05; // we don't need to check pixels that are expected to
                                               // receive this many photons or fewer.
            // inverse of the function defining the Gaussian distribution: Find out how far from the
            // mean we'll have to go until the number of photons is less than interestingThreshold
            float radius = ceil(sqrt(-log(interestingThreshold/peakBrightness)*2*M_PI*starSpreadStdDev*starSpreadStdDev));
            Star star = Star(camCoords.x, camCoords.y,
                             radius, radius,
                             // important to invert magnitude here, so that centroid magnitude becomes larger for brighter stars.
                             // It's possible to make it so that the magnitude is always positive too, but allowing weirder magnitudes helps keep star-id algos honest about their assumptions on magnitude.
                             // we don't use its magnitude anywhere else in generation; peakBrightness was already calculated.
                             -catalogStar.magnitude);
            generatedStars.push_back(GeneratedStar(star, peakBrightness, delta));

            // Now add the star to the input and expected lists.
            // We do actually want to add false stars as well, because:
            // a) A centroider isn't any worse because it picks up a false star that looks exactly like a normal star, so why should we exclude them from compare-centroids?
            // b) We want to feed false stars into star-ids to evaluate their false-star resilience without running centroiding.

            // Add all stars to expected, cause that's how we roll
            expectedStars.push_back(star);
            // and provide expected identifications for all of them
            if (isTrueStar) {
                expectedStarIds.push_back(StarIdentifier(expectedStars.size()-1, i));
            }

            // for input, though, add perturbation and stuff.
            Star inputStar = star;
            if (perturbationStddev > 0.0) {
                // clamp to within 2 standard deviations for some reason:
                inputStar.position.x += std::max(std::min(perturbation1DDistribution(*rng), 2*perturbationStddev), -2*perturbationStddev);
                inputStar.position.y += std::max(std::min(perturbation1DDistribution(*rng), 2*perturbationStddev), -2*perturbationStddev);
            }
            // If it got perturbed outside of the sensor, don't add it.
            if (camera.InSensor(inputStar.position)
                // and also don't add it if it's too dim.
                && (cutoffMag >= 10000
                    || std::bernoulli_distribution(CentroidImagingProbability(catalogStar.magnitude, cutoffMag))(*rng))) {
                inputStars.push_back(inputStar);
                if (isTrueStar) {
                    inputStarIds.push_back(StarIdentifier(inputStars.size()-1, i));
                }
            }
        }
    }

    if (centroidsOnly) {
        // we're outta here
        return;
    }

    std::vector<float> photonsBuffer(image.width*image.height, 0);

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
                                (MotionBlurredPixelBrightness({x, y}, star, tEnd, starSpreadStdDev)
                                 - MotionBlurredPixelBrightness({x, y}, star, tStart, starSpreadStdDev))
                                / oversamplingBrightnessFactor;
                        } else {
                            curPhotons = StaticPixelBrightness({x, y}, star, starSpreadStdDev)
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
        curBrightness += readNoiseDist(*rng);

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
            quantizedPhotons = shotNoiseDist(*rng);
        } else {
            quantizedPhotons = round(photonsBuffer[i]);
        }
        curBrightness += quantizedPhotons / saturationPhotons;

        // std::clamp not introduced until C++17, so we avoid it.
        float clampedBrightness = std::max(std::min(curBrightness, (float)1.0), (float)0.0);
        imageData[i] = floor(clampedBrightness * kMaxBrightness); // TODO: off-by-one, 256?
    }
}

/**
 * Generates a random attitude for coverage testing.
 * Takes a random engine as a parameter.
 */
static Attitude RandomAttitude(std::default_random_engine* pReng) {
    std::uniform_real_distribution<float> randomAngleDistribution(0, 1);

    // normally the ranges of the Ra and Dec are:
    // Dec: [-90 deg, 90 deg] --> [-pi/2 rad, pi/2 rad], where negative means south
    // of the celestial equator and positive means north
    // Ra: [0 deg, 360 deg] --> [0 rad, 2pi rad ]
    // Roll: [0 rad, 2 pi rad]

    float randomRa = 2 *  M_PI * randomAngleDistribution(*pReng);
    float randomDec = (M_PI / 2) - acos(1 - 2 * randomAngleDistribution(*pReng)); //acos returns a float in range [0, pi]
    float randomRoll = 2 *  M_PI * randomAngleDistribution(*pReng);

    Attitude randAttitude = Attitude(SphericalToQuaternion(randomRa, randomDec, randomRoll));

    return randAttitude;
}

/// Create a GeneratedPipelineInput based on the command line options in `values`
PipelineInputList GetGeneratedPipelineInput(const PipelineOptions &values) {
    // TODO: prompt for attitude, imagewidth, etc and then construct a GeneratedPipelineInput

    int seed;

    // time based seed if option specified
    if (values.timeSeed) {
        seed = time(0);
    } else {
        seed = values.generateSeed;
    }


    std::default_random_engine rng(seed);

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


        Attitude inputAttitude;
        if (values.generateRandomAttitudes) {
            inputAttitude = RandomAttitude(&rng);
        } else {
            inputAttitude = attitude;
        }

        GeneratedPipelineInput *curr = new GeneratedPipelineInput(
                CatalogRead(),
                inputAttitude,
                Camera(focalLength, values.generateXRes, values.generateYRes),
                &rng,

                values.generateCentroidsOnly,
                values.generateZeroMagPhotons,
                values.generateSpreadStdDev,
                values.generateSaturationPhotons,
                values.generateDarkCurrent,
                values.generateReadNoiseStdDev,
                motionBlurDirection,
                values.generateExposure,
                values.generateReadoutTime,
                values.generateShotNoise,
                values.generateOversampling,
                values.generateNumFalseStars,
                (values.generateFalseMinMag * 100),
                (values.generateFalseMaxMag * 100),
                (values.generateCutoffMag * 100),
                values.generatePerturbationStddev);

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
    enum class PipelineStage {
        Centroid, CentroidMagnitudeFilter, Database, StarId, AttitudeEstimation, Done
    };

    Pipeline result;

    // TODO: more flexible or sth
    // TODO: don't allow setting star-id until database is set, and perhaps limit the star-id
    // choices to those compatible with the database?

    // centroid algorithm stage
    if (values.centroidAlgo == "dummy") {
        result.centroidAlgorithm = std::unique_ptr<CentroidAlgorithm>(
            new DummyCentroidAlgorithm(values.centroidDummyNumStars));
    } else if (values.centroidAlgo == "cog") {
        result.centroidAlgorithm =
            std::unique_ptr<CentroidAlgorithm>(new CenterOfGravityAlgorithm());
    } else if (values.centroidAlgo == "iwcog") {
        result.centroidAlgorithm =
            std::unique_ptr<CentroidAlgorithm>(new IterativeWeightedCenterOfGravityAlgorithm());
    } else if (values.centroidAlgo == "lsgf1d") {
        result.centroidAlgorithm =
            std::unique_ptr<CentroidAlgorithm>(new LeastSquaresGaussianFit1D(values.centroidFitRadius));
    } else if (values.centroidAlgo == "lsgf2d") {
        result.centroidAlgorithm =
            std::unique_ptr<CentroidAlgorithm>(new LeastSquaresGaussianFit2D(values.centroidFitRadius));
    } else if (values.centroidAlgo == "ggrid"){
        result.centroidAlgorithm =
            std::unique_ptr<CentroidAlgorithm>(new GaussianGrid());
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
    PipelineOutput result;

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
            std::cerr << "WARNING: That database does not include a catalog. Proceeding with the full catalog." << std::endl;
            result.catalog = input.GetCatalog();
        }
    } else {
        result.catalog = input.GetCatalog();
    }

    if (centroidAlgorithm && inputImage) {

        // run centroiding, keeping track of the time it takes
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

        // TODO: we should probably modify Go to just take an image argument
        Stars unfilteredStars = centroidAlgorithm->Go(inputImage->image, inputImage->width, inputImage->height);

        std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
        result.centroidingTimeNs = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

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
    } else if (centroidAlgorithm) {
        std::cerr << "ERROR: Centroid algorithm specified, but no input image to run it on." << std::endl;
        exit(1);
    }

    if (starIdAlgorithm && database && inputStars && input.InputCamera()) {
        // TODO: don't copy the vector!
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

        result.starIds = std::unique_ptr<StarIdentifiers>(new std::vector<StarIdentifier>(
            starIdAlgorithm->Go(database.get(), *inputStars, result.catalog, *input.InputCamera())));

        std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
        result.starIdTimeNs = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        inputStarIds = result.starIds.get();
    } else if (starIdAlgorithm) {
        std::cerr << "ERROR: Star ID algorithm specified but cannot run because database, centroids, or camera are missing." << std::endl;
        exit(1);
    }

    if (attitudeEstimationAlgorithm && inputStarIds && input.InputCamera()) {
        assert(inputStars); // ensure that starIds doesn't exist without stars
        std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

        result.attitude = std::unique_ptr<Attitude>(
            new Attitude(attitudeEstimationAlgorithm->Go(*input.InputCamera(), *inputStars, result.catalog, *inputStarIds)));

        std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();
        result.attitudeEstimationTimeNs = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    } else if (attitudeEstimationAlgorithm) {
        std::cerr << "ERROR: Attitude estimation algorithm set, but either star IDs or camera are missing. One reason this can happen: Setting a centroid algorithm and attitude algorithm, but no star-id algorithm -- that can't work because the input star-ids won't properly correspond to the output centroids!" << std::endl;
        exit(1);
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
    CentroidComparison() : meanError(0.0), numCorrectCentroids(0), numExtraCentroids(0) { }
    /**
     * Average distance from actual to expected centroids (in pixels)
     * Only correct centroids are considered in this average.
     */
    float meanError;

    /**
     * Number of actual stars within the centroiding threshold of an expected star.
     */
    float numCorrectCentroids;

    /**
     * Stars in actual but not expected. Ideally 0
     * This is a float because we may average multiple centroid comparisons together.
     */
    float numExtraCentroids;

    // We no longer have a num missing because often generated stars have too low of a signal-to-noise ratio and the centroid algo won't pick them up.
};

/// Create a mapping, where keys are indices into `one` and values are indices of all centroids in
/// `two` whose distance to the current `one` star is <=threshold. In a (ordered) multimap, the
/// insertion order is preserved for elements with the same key, and indeed we'll sort the elements
/// corresponding to each key by distance from the corresponding `one` star.
static std::multimap<int, int> FindClosestCentroids(float threshold,
                                                    const Stars &one,
                                                    const Stars &two) {
    std::multimap<int, int> result;

    for (int i = 0; i < (int)one.size(); i++) {
        std::vector<std::pair<float, int>> closest;
        for (int k = 0; k < (int)two.size(); k++) {
            float currDistance = (one[i].position - two[k].position).Magnitude();
            if (currDistance <= threshold) {
                closest.emplace_back(currDistance, k);
            }
        }
        std::sort(closest.begin(), closest.end());
        for (const std::pair<float, int> &pair : closest) {
            result.emplace(i, pair.second);
        }
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

    // TODO: Somehow penalize when multiple centroids correspond to the same expected star (i.e.,
    // one star turned into multiple centroids). That should probably be considered an extra
    // centroid, but rn it isn't.

    CentroidComparison result;
    // maps from indexes in each list to the closest centroid from other list
    std::multimap<int, int> actualToExpected = FindClosestCentroids(threshold, actual, expected);

    for (int i = 0; i < (int)actualToExpected.size(); i++) {
        auto closest = actualToExpected.find(i);
        if (closest == actualToExpected.end()) {
            result.numExtraCentroids++;
        } else {
            result.meanError += (actual[i].position - expected[closest->second].position).Magnitude();
            result.numCorrectCentroids++;
        }
    }
    result.meanError /= result.numCorrectCentroids;

    return result;
}

CentroidComparison CentroidComparisonsCombine(std::vector<CentroidComparison> comparisons) {
    assert(comparisons.size() > 0);

    CentroidComparison result;

    for (const CentroidComparison &comparison : comparisons) {
        result.meanError += comparison.meanError;
        result.numExtraCentroids += comparison.numExtraCentroids;
    }

    result.meanError /= comparisons.size();
    result.numExtraCentroids /= comparisons.size();

    return result;
}

// (documentation in hpp)
StarIdComparison StarIdsCompare(const StarIdentifiers &expected, const StarIdentifiers &actual,
                                // use these to map indices to names for the respective lists of StarIdentifiers
                                const Catalog &expectedCatalog, const Catalog &actualCatalog,
                                float centroidThreshold,
                                const Stars &expectedStars, const Stars &inputStars) {

    StarIdComparison result = {
        0, // correct
        0, // incorrect
        0, // total
    };

    // EXPECTED STAR IDS

    // map from expected star indices to expected catalog indices (basically flattening the expected star-ids)
    std::vector<int> expectedCatalogIndices(expectedStars.size(), -1);
    for (const StarIdentifier &starId : expected) {
        assert(0 <= starId.starIndex && starId.starIndex <= (int)expectedStars.size());
        assert(0 <= starId.catalogIndex && starId.catalogIndex <= (int)expectedCatalog.size());
        expectedCatalogIndices[starId.starIndex] = starId.catalogIndex;
    }

    // FIND NEAREST CENTROIDS

    std::multimap<int, int> inputToExpectedCentroids = FindClosestCentroids(centroidThreshold, inputStars, expectedStars);
    // std::multimap<int, int> expectedToInputCentroids = FindClosestCentroids(centroidThreshold, expectedStars, inputStars);

    // COMPUTE TOTAL
    // Count the number of expected stars with at least one input star near them
    for (int i = 0; i < (int)inputStars.size(); i++) {
        // make sure there's at least one expected star near this input star which has an identification
        auto closestRange = inputToExpectedCentroids.equal_range(i);
        bool found = false;
        for (auto it = closestRange.first; it != closestRange.second; it++) {
            if (expectedCatalogIndices[it->second] != -1) {
                found = true;
                break;
            }
        }
        if (found) {
            result.numTotal++;
        }
    }

    // COMPUTE CORRECT AND INCORRECT

    std::vector<bool> identifiedInputCentroids(inputStars.size(), false);
    for (const StarIdentifier &starId : actual) {
        // as later, there shouldn't be duplicate starIndex. This indicates a bug in the star-id algorithm, not comparison code.
        assert(!identifiedInputCentroids[starId.starIndex]);
        identifiedInputCentroids[starId.starIndex] = true;
        assert(0 <= starId.starIndex && starId.starIndex <= (int)inputStars.size());
        assert(0 <= starId.catalogIndex && starId.catalogIndex <= (int)actualCatalog.size());

        // Check that there's at least one expected centroid in range which agrees with your identification.
        auto expectedCentroidsInRange = inputToExpectedCentroids.equal_range(starId.starIndex);
        bool found = false;
        for (auto it = expectedCentroidsInRange.first; it != expectedCentroidsInRange.second; it++) {
            int expectedCatalogIndex = expectedCatalogIndices[it->second];
            if (expectedCatalogIndex != -1
                && expectedCatalog[expectedCatalogIndex].name == actualCatalog[starId.catalogIndex].name) {

                result.numCorrect++;
                found = true;
                break;
            }
        }

        // Either there's no expected centroid in range, or none of them agree with the identification.
        if (!found) {
            result.numIncorrect++;
        }
    }

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
static cairo_status_t OstreamPlotter(void *closure, const unsigned char *data, unsigned int length) {
    std::ostream *os = (std::ostream *)closure;
    os->write((const char *)data, length);
    return CAIRO_STATUS_SUCCESS;
}

/// Plots the input image with no annotation to `os`
static void PipelineComparatorPlotRawInput(std::ostream &os,
                                    const PipelineInputList &expected,
                                    const std::vector<PipelineOutput> &,
                                    const PipelineOptions &) {

    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Plot the annotated input image to `os`
// TODO: should probably use Expected methods, not Input methods, because future PipelineInputs could add noise to the result of the Input methods.
static void PipelineComparatorPlotInput(std::ostream &os,
                                 const PipelineInputList &expected,
                                 const std::vector<PipelineOutput> &,
                                 const PipelineOptions &) {
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    assert(expected[0]->InputStars() != NULL);
    SurfacePlot("pipeline input",
                cairoSurface,
                *expected[0]->InputStars(),
                expected[0]->InputStarIds(),
                &expected[0]->GetCatalog(),
                expected[0]->InputAttitude(),
                // green
                0.0, 1.0, 0.0, 0.6);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

static void PipelineComparatorPlotExpected(std::ostream &os,
                                    const PipelineInputList &expected,
                                    const std::vector<PipelineOutput> &,
                                    const PipelineOptions &) {
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    assert(expected[0]->ExpectedStars() != NULL);
    SurfacePlot("expected output",
                cairoSurface,
                *expected[0]->ExpectedStars(),
                expected[0]->ExpectedStarIds(),
                &expected[0]->GetCatalog(),
                expected[0]->ExpectedAttitude(),
                // blu
                0.2, 0.5, 1.0, 0.7);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Compare the actual and expected centroids, printing key stats to `os`
static void PipelineComparatorCentroids(std::ostream &os,
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
    os << "centroids_num_correct " << result.numCorrectCentroids << std::endl
       << "centroids_num_extra " << result.numExtraCentroids << std::endl
       << "centroids_mean_error " << result.meanError << std::endl;
}

static void PrintCentroids(const std::string &prefix,
                           std::ostream &os,
                           const Catalog &catalog,
                           const Stars &stars,
                           // may be null:
                           const StarIdentifiers *starIds) {
    os << "num_" << prefix << "_centroids " << stars.size() << std::endl;
    for (int i = 0; i < (int)stars.size(); i++) {
        os << prefix << "_centroid_" << i << "_x " << stars[i].position.x << std::endl;
        os << prefix << "_centroid_" << i << "_y " << stars[i].position.y << std::endl;
        if (starIds) {
            for (const StarIdentifier &starId : *starIds) {
                if (starId.starIndex == i) {
                    os << prefix << "_centroid_" << i << "_id " << catalog[starId.catalogIndex].name << std::endl;
                }
            }
        }
    }
}

/// Print a list of centroids to `os`
static void PipelineComparatorPrintExpectedCentroids(std::ostream &os,
                                                     const PipelineInputList &expected,
                                                     const std::vector<PipelineOutput> &, // actual
                                                     const PipelineOptions &) {
    assert(expected.size() == 1);
    assert(expected[0]->ExpectedStars());

    PrintCentroids("expected",
                   os,
                   expected[0]->GetCatalog(),
                   *expected[0]->ExpectedStars(),
                   expected[0]->ExpectedStarIds());
}

static void PipelineComparatorPrintActualCentroids(std::ostream &os,
                                                   const PipelineInputList &, // expected
                                                   const std::vector<PipelineOutput> &actual,
                                                   const PipelineOptions &) {
    assert(actual.size() == 1);
    assert(actual[0].stars);

    PrintCentroids("actual",
                   os,
                   actual[0].catalog,
                   *actual[0].stars,
                   actual[0].starIds.get());
}

/// Plot an annotated image where centroids are annotated with their centroid index. For debugging.
/// Use whatever stars were input into the star-id algo (so either actual centroids, or inputstars)
void PipelineComparatorPlotCentroidIndices(std::ostream &os,
                                           const PipelineInputList &expected,
                                           const std::vector<PipelineOutput> &actual,
                                           const PipelineOptions &) {
    const Stars &stars = actual[0].stars ? *actual[0].stars : *expected[0]->InputStars();
    StarIdentifiers identifiers;
    for (int i = 0; i < (int)stars.size(); i++) {
        identifiers.push_back(StarIdentifier(i, i));
    }
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    SurfacePlot("centroid indices (input)",
                cairoSurface,
                stars,
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
static void PipelineComparatorPlotOutput(std::ostream &os,
                                         const PipelineInputList &expected,
                                         const std::vector<PipelineOutput> &actual,
                                         const PipelineOptions &) {
    // don't need to worry about mutating the surface; InputImageSurface returns a fresh one
    cairo_surface_t *cairoSurface = expected[0]->InputImageSurface();
    SurfacePlot("pipeline output",
                cairoSurface,
                actual[0].stars ? *actual[0].stars : *expected[0]->InputStars(),
                actual[0].starIds.get(),
                &actual[0].catalog,
                actual[0].attitude.get(),
                // red
                1.0, 0.0, 0.0, 0.5);
    cairo_surface_write_to_png_stream(cairoSurface, OstreamPlotter, &os);
    cairo_surface_destroy(cairoSurface);
}

/// Compare the expected and actual star identifiers.
static void PipelineComparatorStarIds(std::ostream &os,
                                      const PipelineInputList &expected,
                                      const std::vector<PipelineOutput> &actual,
                                      const PipelineOptions &values) {
    int numImagesCorrect = 0;
    int numImagesIncorrect = 0;
    int numImagesTotal = expected.size();
    for (int i = 0; i < numImagesTotal; i++) {
        // since the actual star IDs exist, it must have gotten input from somewhere!
        // TODO: overhaul: It seems that in these comparators there should be a more fundamental way to figure out the input that was actually sent to a stage.
        // I.e., instead of having expected and actual arguments, have some sort of PipelineRunSummary object, where the InputStars method looks at actual, then input.
        assert(actual[i].stars.get() || expected[i]->InputStars());

        const Stars &inputStars = actual[i].stars.get()
            ? *actual[i].stars.get()
            : *expected[i]->InputStars();
        StarIdComparison comparison =
            StarIdsCompare(*expected[i]->ExpectedStarIds(), *actual[i].starIds,
                           expected[i]->GetCatalog(), actual[i].catalog,
                           values.centroidCompareThreshold, *expected[i]->ExpectedStars(), inputStars);

        if (numImagesTotal == 1) {
            os << "starid_num_correct " << comparison.numCorrect << std::endl;
            os << "starid_num_incorrect " << comparison.numIncorrect << std::endl;
            os << "starid_num_total " << comparison.numTotal << std::endl;
        }

        if (comparison.numCorrect > 0 && comparison.numIncorrect == 0) {
            numImagesCorrect++;
        }
        if (comparison.numIncorrect > 0) {
            numImagesIncorrect++;
        }
    }

    // A "correct" image is one where at least two stars are correctly id'd and none are incorrectly id'd
    os << "starid_num_images_correct " << numImagesCorrect << std::endl;
    os << "starid_num_images_incorrect " << numImagesIncorrect << std::endl;
}

/// Print the identifed attitude to `os` in Euler angle format.
static void PipelineComparatorPrintAttitude(std::ostream &os,
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
static void PipelineComparatorAttitude(std::ostream &os,
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

static void PrintTimeStats(std::ostream &os, const std::string &prefix, const std::vector<long> &times) {
    assert(times.size() > 0);

    // print average, min, max, and 95% max
    long sum = 0;
    long min = LONG_MAX;
    long max = 0;
    for (int i = 0; i < (int)times.size(); i++) {
        assert(times[i] > 0);
        sum += times[i];
        min = std::min(min, times[i]);
        max = std::max(max, times[i]);
    }
    long average = sum / times.size();
    std::vector<long> sortedTimes = times;
    std::sort(sortedTimes.begin(), sortedTimes.end());
    // what really is the 95th percentile? Being conservative, we want to pick a value that at least
    // 95% of the times are less than. This means: (1) finding the number of times, (2) Finding
    // Math.ceil(0.95 * numTimes), and (3) subtracting 1 to get the index.
    int ninetyFiveIndex = (int)std::ceil(0.95 * times.size()) - 1;
    assert(ninetyFiveIndex >= 0);
    long ninetyFifthPercentile = sortedTimes[ninetyFiveIndex];

    os << prefix << "_average_ns " << average << std::endl;
    os << prefix << "_min_ns " << min << std::endl;
    os << prefix << "_max_ns " << max << std::endl;
    os << prefix << "_95%_ns " << ninetyFifthPercentile << std::endl;
}

/// For each stage of the pipeline, print statistics about how long it took to run.
static void PipelineComparatorPrintSpeed(std::ostream &os,
                                    const PipelineInputList &,
                                    const std::vector<PipelineOutput> &actual,
                                    const PipelineOptions &) {
    std::vector<long> centroidingTimes;
    std::vector<long> starIdTimes;
    std::vector<long> attitudeTimes;
    std::vector<long> totalTimes;
    for (int i = 0; i < (int)actual.size(); i++) {
        long totalTime = 0;
        if (actual[i].centroidingTimeNs > 0) {
            centroidingTimes.push_back(actual[i].centroidingTimeNs);
            totalTime += actual[i].centroidingTimeNs;
        }
        if (actual[i].starIdTimeNs > 0) {
            starIdTimes.push_back(actual[i].starIdTimeNs);
            totalTime += actual[i].starIdTimeNs;
        }
        if (actual[i].attitudeEstimationTimeNs > 0) {
            attitudeTimes.push_back(actual[i].attitudeEstimationTimeNs);
            totalTime += actual[i].attitudeEstimationTimeNs;
        }
        totalTimes.push_back(totalTime);
    }
    if (centroidingTimes.size() > 0) {
        PrintTimeStats(os, "centroiding", centroidingTimes);
    }
    if (starIdTimes.size() > 0) {
        PrintTimeStats(os, "starid", starIdTimes);
    }
    if (attitudeTimes.size() > 0) {
        PrintTimeStats(os, "attitude", attitudeTimes);
    }
    PrintTimeStats(os, "total", totalTimes);
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
    if (actual.size() == 0) {
        std::cerr << "ERROR: No \"comparator\"/output action was specified. I.e., the star identification is all done, but you didn't specify how to return the results to you! Try --plot-output <filepath>, perhaps." << std::endl;
        exit(1);
    }

    assert(expected.size() == actual.size() && expected.size() > 0);

    // TODO: Remove the asserts and print out more reasonable error messages.

#define LOST_PIPELINE_COMPARE(precondition, errmsg, comparator, path, isBinary) do { \
        if (precondition) {                                             \
            UserSpecifiedOutputStream pos(path, isBinary);              \
            comparator(pos.Stream(), expected, actual, values);         \
        } else {                                                        \
            std::cerr << "ERROR: Comparator not applicable: " << errmsg << std::endl; \
            exit(1);                                                    \
        }                                                               \
    } while (0)

    if (values.plotRawInput != "") {
        LOST_PIPELINE_COMPARE(expected[0]->InputImage() && expected.size() == 1,
                              "--plot-raw-input requires exactly 1 input image, but " + std::to_string(expected.size()) + " many were provided.",
                              PipelineComparatorPlotRawInput, values.plotRawInput, true);
    }

    if (values.plotInput != "") {
        LOST_PIPELINE_COMPARE(expected[0]->InputImage() && expected.size() == 1 && expected[0]->InputStars(),
                              "--plot-input requires exactly 1 input image, and for centroids to be available on that input image. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPlotInput, values.plotInput, true);
    }
    if (values.plotExpected != "") {
        LOST_PIPELINE_COMPARE(expected[0]->InputImage() && expected.size() == 1 && expected[0]->ExpectedStars(),
                              "--plot-expected-input requires exactly 1 input image, and for expected centroids to be available on that input image. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPlotExpected, values.plotExpected, true);
    }
    if (values.plotOutput != "") {
        LOST_PIPELINE_COMPARE(actual.size() == 1 && (actual[0].stars || actual[0].starIds),
                              "--plot-output requires exactly 1 output image, and for either centroids or star IDs to be available on that output image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPlotOutput, values.plotOutput, true);
    }
    if (values.printExpectedCentroids != "") {
        LOST_PIPELINE_COMPARE(expected.size() == 1 && expected[0]->ExpectedStars(),
                              "--print-expected-centroids requires exactly 1 input image, and for expected centroids to be available on that input image. " + std::to_string(expected.size()) + " many input images were provided.",
                              PipelineComparatorPrintExpectedCentroids, values.printExpectedCentroids, false);
    }
    if (values.printActualCentroids != "") {
        LOST_PIPELINE_COMPARE(actual.size() == 1 && actual[0].stars,
                              "--print-actual-centroids requires exactly 1 output image, and for centroids to be available on that output image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPrintActualCentroids, values.printActualCentroids, false);
    }
    if (values.plotCentroidIndices != "") {
        LOST_PIPELINE_COMPARE(expected.size() == 1 && expected[0]->InputImage(),
                              "--plot-centroid-indices requires exactly 1 input with image. " + std::to_string(expected.size()) + " many inputs were provided.",
                              PipelineComparatorPlotCentroidIndices, values.plotCentroidIndices, true);
    }
    if (values.compareCentroids != "") {
        LOST_PIPELINE_COMPARE(actual[0].stars && expected[0]->ExpectedStars() && values.centroidCompareThreshold,
                              "--compare-centroids requires at least 1 output image, and for expected centroids to be available on the input image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorCentroids, values.compareCentroids, false);
    }
    if (values.compareStarIds != "") {
        LOST_PIPELINE_COMPARE(expected[0]->ExpectedStarIds() && actual[0].starIds && expected[0]->ExpectedStars(),
                              "--compare-star-ids requires at least 1 output image, and for expected star IDs and centroids to be available on the input image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorStarIds, values.compareStarIds, false);
    }
    if (values.printAttitude != "") {
        LOST_PIPELINE_COMPARE(actual[0].attitude && actual.size() == 1,
                              "--print-attitude requires exactly 1 output image, and for attitude to be available on that output image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPrintAttitude, values.printAttitude, false);
    }
    if (values.compareAttitudes != "") {
        LOST_PIPELINE_COMPARE(actual[0].attitude && expected[0]->ExpectedAttitude() && values.attitudeCompareThreshold,
                              "--compare-attitudes requires at least 1 output image, and for expected attitude to be available on the input image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorAttitude, values.compareAttitudes, false);
    }
    if (values.printSpeed != "") {
        LOST_PIPELINE_COMPARE(actual.size() > 0,
                              // I don't think this should ever actually happen??
                              "--print-speed requires at least 1 output image. " + std::to_string(actual.size()) + " many output images were provided.",
                              PipelineComparatorPrintSpeed, values.printSpeed, false);
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
