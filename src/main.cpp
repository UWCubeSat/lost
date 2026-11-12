/**
 * LOST starting point
 *
 * Reads in CLI arguments/flags and starts the appropriate pipelines
 */

#include <bitset>
#include <iostream>
#include <string>
#include <vector>

#include "CLI11.hpp"
#include "databases.hpp"
#include "decimal.hpp"
#include "io.hpp"

namespace lost {

static void DatabaseBuild(const DatabaseOptions& values) {
    Catalog narrowedCatalog = NarrowCatalog(
        CatalogRead(), (int)(values.minMag * 100), values.maxStars, DegToRad(values.minSeparation));
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseDescriptor dbEntries = GenerateDatabases(narrowedCatalog, values);
    SerializeContext ser = serFromDbValues(values);

    uint32_t dbFlags = 0;
    dbFlags |= typeid(decimal) == typeid(float) ? MULTI_DB_FLOAT_FLAG : 0;

    SerializeMultiDatabase(&ser, dbEntries, dbFlags);

    std::cerr << "Generated database with " << ser.buffer.size() << " bytes" << std::endl;
    std::cerr << "Database flagged with " << std::bitset<8 * sizeof(dbFlags)>(dbFlags) << std::endl;

    UserSpecifiedOutputStream pos = UserSpecifiedOutputStream(values.outputPath, true);
    pos.Stream().write((char*)ser.buffer.data(), ser.buffer.size());
}

static void PipelineRun(const PipelineOptions& values) {
    PipelineInputList input = GetPipelineInput(values);
    Pipeline pipeline = SetPipeline(values);
    std::vector<PipelineOutput> outputs = pipeline.Go(input);
    PipelineComparison(input, outputs, values);
}

////////////////////////
// CLI OPTION HELPERS //
////////////////////////

/**
 * Add an output path option whose value may be omitted, e.g. both `--compare-star-ids` and
 * `--compare-star-ids out.txt` are valid. No argument supplied means stdout ("-").
 */
static CLI::Option* AddOutputPathOption(CLI::App* cmd,
                                        const std::string& name,
                                        std::string& variable,
                                        const std::string& description) {
    return cmd->add_option(name, variable, description)
        ->expected(0, 1)
        ->default_str("-")
        ->option_text("[PATH] (no value: stdout)");
}

//////////////////////
// DATABASE OPTIONS //
//////////////////////

static void AddCatalogNarrowingOptions(CLI::App* cmd, DatabaseOptions& opts) {
    const std::string group = "CATALOG NARROWING OPTIONS";
    cmd->add_option("--min-mag",
                    opts.minMag,
                    "Narrows catalog to stars with magnitude less than minimum-magnitude (brighter "
                    "stars have lower magnitudes).")
        ->group(group);
    cmd->add_option("--max-stars",
                    opts.maxStars,
                    "Limits catalog to the num-stars brightest stars in the catalog.")
        ->check(CLI::PositiveNumber)
        ->group(group);
    cmd->add_option("--min-separation",
                    opts.minSeparation,
                    "Removes stars from the catalog that are within minimum-separation degrees of "
                    "another catalog star. The default option is probably good enough for you.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
}

static void AddKVectorOptions(CLI::App* cmd, DatabaseOptions& opts) {
    const std::string group = "PAIR-DISTANCE KVECTOR DATABASE OPTIONS";
    // `--kvector=false` is rejected: passing --kvector always means true
    CLI::Option* kvector = cmd->add_flag("--kvector", opts.kvector, "Generate a KVector database")
                               ->disable_flag_override()
                               ->group(group);
    cmd->add_option("--kvector-min-distance",
                    opts.kvectorMinDistance,
                    "Sets the minimum distance bin in the kvector building method to min "
                    "(degrees). Defaults to 0.5 if option is not selected.")
        ->check(CLI::NonNegativeNumber)
        ->needs(kvector)
        ->group(group);
    cmd->add_option("--kvector-max-distance",
                    opts.kvectorMaxDistance,
                    "Sets the maximum distance max (degrees). Defaults to 15 if option is not "
                    "selected. Should be set to no larger than your camera's FOV, otherwise you're "
                    "just wasting time and space.")
        ->check(CLI::PositiveNumber)
        ->needs(kvector)
        ->group(group);
    cmd->add_option("--kvector-distance-bins",
                    opts.kvectorNumDistanceBins,
                    "Sets the number of distance bins in the kvector building method to num-bins.  "
                    "Defaults to 10000 if option is not selected, which is pretty reasonable for "
                    "most cases.")
        ->check(CLI::PositiveNumber)
        ->needs(kvector)
        ->group(group);
}

static void AddTetraDatabaseOptions(CLI::App* cmd, DatabaseOptions& opts) {
    const std::string group = "TETRA DATABASE OPTIONS";
    // `--tetra=false` is rejected: passing --tetra always means true
    CLI::Option* tetra = cmd->add_flag("--tetra", opts.tetra, "Generate a Tetra database.")
                             ->disable_flag_override()
                             ->group(group);
    cmd->add_option("--tetra-max-fov",
                    opts.tetraMaxFov,
                    "Sets the maximum field of view for Tetra patterns to fov degrees. All stars "
                    "in a pattern must be within this angle of each other. Should be at least as "
                    "large as the camera's FOV. Defaults to 12.")
        ->check(CLI::PositiveNumber)
        ->needs(tetra)
        ->group(group);
    cmd->add_option("--tetra-edge-ratio-bins",
                    opts.tetraEdgeRatioBins,
                    "Number of bins used to discretize edge ratios when hashing patterns. Higher "
                    "values give finer discretization and fewer hash collisions. Defaults to 50.")
        ->check(CLI::PositiveNumber)
        ->needs(tetra)
        ->group(group);
    cmd->add_option("--tetra-pattern-stars-per-fov",
                    opts.tetraPatternStarsPerFov,
                    "Maximum number of pattern stars allowed within half the FOV of any other "
                    "given pattern star during catalog preprocessing. Be careful adjusting this "
                    "value. Increasing it can cause large growth in database size while lowering "
                    "it risks gaps in sky coverage. Defaults to 10.")
        ->check(CLI::PositiveNumber)
        ->needs(tetra)
        ->group(group);
}

static void AddOtherDatabaseOptions(CLI::App* cmd, DatabaseOptions& opts) {
    const std::string group = "OTHER DATABASE OPTIONS";
    cmd->add_flag("--swap-integer-endianness",
                  opts.swapIntegerEndianness,
                  "If true, generate databases with all integer values having opposite endianness "
                  "than the "
                  "generating machine. It will not be possible to use the generated databases on "
                  "the system "
                  "they were generated on.")
        ->group(group);
    cmd->add_flag(
           "--swap-decimal-endianness",
           opts.swapDecimalEndianness,
           "If true, generate databases with all floating point values having opposite endianness "
           "than the generating machine. It will not be possible to use the generated databases on "
           "the system they were generated on.")
        ->group(group);
    cmd->add_option(
           "--output", opts.outputPath, "The file to output the database to. Defaults to stdout.")
        ->group(group);
}

static void AddDatabaseOptions(CLI::App* cmd, DatabaseOptions& opts) {
    AddCatalogNarrowingOptions(cmd, opts);
    AddKVectorOptions(cmd, opts);
    AddTetraDatabaseOptions(cmd, opts);
    AddOtherDatabaseOptions(cmd, opts);
}

//////////////////////
// PIPELINE OPTIONS //
//////////////////////

static CLI::Option* AddCameraOptions(CLI::App* cmd, PipelineOptions& opts) {
    const std::string group = "CAMERA OPTIONS";
    CLI::Option* png =
        cmd->add_option("--png", opts.png, "Identify the png image at the given filepath.")
            ->group(group);
    cmd->add_option("--focal-length",
                    opts.focalLength,
                    "The focal length of the camera that took the picture (in mm).")
        ->check(CLI::PositiveNumber)
        ->group(group);
    cmd->add_option("--pixel-size",
                    opts.pixelSize,
                    "The pixel size of the camera that took the picture (in micrometers).")
        ->check(CLI::PositiveNumber)
        ->group(group);
    cmd->add_option("--fov",
                    opts.fov,
                    "The field of view of the camera that took the picture (in degrees). Defaults "
                    "to 20 deg if option is not selected.")
        ->check(CLI::PositiveNumber)
        ->group(group);
    return png;
}

static void AddCentroidOptions(CLI::App* cmd, PipelineOptions& opts) {
    const std::string group = "CENTROID OPTIONS";
    cmd->add_option(
           "--centroid-algo",
           opts.centroidAlgo,
           "Runs the algo centroiding algorithm. Recognized options are: dummy (random centroid "
           "algorithm), cog (center of gravity), and iwcog (iterative weighted center of gravity). "
           " "
           "Defaults to dummy if option is not selected.")
        ->check(CLI::IsMember({"dummy", "cog", "iwcog"}))
        ->group(group);
    cmd->add_option("--centroid-dummy-stars",
                    opts.centroidDummyNumStars,
                    "Runs the dummy centroiding algorithm (random centroid algorithm) with "
                    "num-stars stars centroided. Defaults to 5 if option is not selected.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--centroid-mag-filter",
                    opts.centroidMagFilter,
                    "Will not consider centroids with magnitude below min-mag.")
        ->group(group);
    cmd->add_option(
           "--centroid-filter-brightest",
           opts.centroidFilterBrightest,
           "Remove all but the brightest num-stars many stars from the list of centroids before "
           "sending to star-id. Often a better choice than --centroid-mag-filter, because you can "
           "ensure that you keep enough stars to do star-id. If both this option and "
           "--centroid-mag-filter are provided, then all stars satisfying both criteria are kept "
           "(intersection).")
        ->group(group);
}

static void AddStarIdOptions(CLI::App* cmd, PipelineOptions& opts) {
    const std::string group = "STAR IDENTIFICATION OPTIONS";
    cmd->add_option("--database",
                    opts.databasePath,
                    "Chooses filename as the database to use during star identification.")
        ->group(group);
    cmd->add_option(
           "--star-id-algo",
           opts.idAlgo,
           "Runs the algo star identification algorithm. Current options are \"dummy\", \"gv\", "
           "\"tetra\", and \"pyramid\". Defaults to \"dummy\" if option is not selected.")
        ->check(CLI::IsMember({"dummy", "gv", "py", "tetra"}))
        ->group(group);
    cmd->add_option("--angular-tolerance",
                    opts.angularTolerance,
                    "Sets the estimated angular centroiding error tolerance, used in some star id "
                    "algorithms, to tolerance degrees. Defaults to 0.04 degrees.")
        ->check(CLI::PositiveNumber)
        ->group(group);
    cmd->add_option("--false-stars-estimate,--false-stars",
                    opts.estimatedNumFalseStars,
                    "num is the estimated number of false stars in the whole sphere for the "
                    "pyramid scheme star identification algorithm. Defaults to 500 if option is "
                    "not selected.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--max-mismatch-probability,--max-mismatch-prob",
                    opts.maxMismatchProb,
                    "probability is the maximum allowable probability of an incorrect star "
                    "identification, for star id algorithms which support it. Defaults to 0.001.")
        ->check(CLI::Range(0.0, 1.0))
        ->group(group);
    cmd->add_option("--tetra-pairing",
                    opts.tetraPairingMode,
                    "Sets the method used to pair observed centroids with catalog stars after a "
                    "Tetra pattern match. Options are \"triad\" (try all 24 permutations, pick the "
                    "one with lowest TRIAD attitude residual) and \"centroid-radius\" (sort by "
                    "distance from pattern centroid, faster but less robust for symmetric "
                    "patterns). Defaults to \"triad\".")
        ->check(CLI::IsMember({"triad", "centroid-radius"}))
        ->group(group);
    cmd->add_option("--tetra-pattern-search-range",
                    opts.tetraPatternSearchRange,
                    "How far from the observed edge ratios to search when looking up patterns in "
                    "the database. A larger value is more tolerant of centroiding error but slower "
                    "and more likely to produce false matches. Defaults to 0.001.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--tetra-pattern-match-tolerance",
                    opts.tetraPatternMatchTolerance,
                    "Maximum allowable difference between observed and catalog edge ratios when "
                    "verifying a candidate match from the hash table. Defaults to 0.001.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
}

static void AddAttitudeOptions(CLI::App* cmd, PipelineOptions& opts) {
    cmd->add_option(
           "--attitude-algo",
           opts.attitudeAlgo,
           "Runs the algo algorithm for the attitude stage of the pipeline. Current options are "
           "\"dqm\" (Davenport Q), \"triad\", and \"quest\". Defaults to dqm.")
        ->check(CLI::IsMember({"dqm", "triad", "quest"}))
        ->group("ATTITUDE DETERMINATION OPTIONS");
}

static void AddImageGenerationOptions(CLI::App* cmd, PipelineOptions& opts, CLI::Option* png) {
    const std::string group = "IMAGE GENERATION OPTIONS";
    cmd->add_option("--generate",
                    opts.generate,
                    "Generates num-images many images. Mutually exclusive with --png.")
        ->check(CLI::NonNegativeNumber)
        ->excludes(png)
        ->group(group);
    cmd->add_option("--generate-x-resolution",
                    opts.generateXRes,
                    "Sets the horizontal resolution of the generated image(s) to pixels. Defaults "
                    "to 1024.")
        ->check(CLI::PositiveNumber)
        ->group(group);
    cmd->add_option("--generate-y-resolution",
                    opts.generateYRes,
                    "Sets the vertical resolution of the generated image(s) to pixels. Defaults to "
                    "1024.")
        ->check(CLI::PositiveNumber)
        ->group(group);
    cmd->add_option(
           "--generate-zero-mag-photons",
           opts.generateZeroMagPhotons,
           "A star with magnitude 0 will cause num-photoelectrons many photoelectrons to "
           "be received by the sensor per time unit. A default value of 20,000 is chosen. "
           "See Liebe's, \"tutorial on star tracker accuracy\" for theoretical information "
           "on how to calculate this.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--generate-saturation-photons",
                    opts.generateSaturationPhotons,
                    "When a pixel receives at least this many photoelectrons, it will appear "
                    "completely white (at least before other noise is applied). Note that, because "
                    "of noise, a pixel may still appear completely white if it receives less than "
                    "saturation-photoelectrons many photoelectrons.")
        ->check(CLI::PositiveNumber)
        ->group(group);
    cmd->add_option("--generate-spread-stddev",
                    opts.generateSpreadStdDev,
                    "Sets the standard deviation of the Gaussian distribution we use to "
                    "approximate the distribution of photons received from each star. Defaults to "
                    "1.0.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--generate-shot-noise",
                    opts.generateShotNoise,
                    "Enables or disables shot noise simulation in generated images. Defaults to "
                    "true.")
        ->group(group);
    cmd->add_option("--generate-dark-current",
                    opts.generateDarkCurrent,
                    "Set observed brightness of dark current in the image, from 0 (no dark noise) "
                    "to 1 (whole image pure white). Defaults to 0.1. cf --generate-sensitivity to "
                    "control shot noise intensity.")
        ->check(CLI::Range(0.0, 1.0))
        ->group(group);
    cmd->add_option("--generate-read-noise-stddev",
                    opts.generateReadNoiseStdDev,
                    "Sets the standard deviation of Gaussian noise in the generated image(s) to "
                    "stddev. Noise is measured in observed brightness, where 1 is the difference "
                    "between pure white and pure black. Defaults to 0.05.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_flag("--generate-centroids-only",
                  opts.generateCentroidsOnly,
                  "If set (argument is not needed, but can be a bool), will not actually generate "
                  "any pixel "
                  "data -- only centroid positions. Much faster if only evaluating a star-id "
                  "algorithm (eg "
                  "sky coverage testing).")
        ->group(group);

    const std::string attitudeDesc =
        "Set the right ascension, declination, and roll of the generated image.";
    cmd->add_option("--generate-ra", opts.generateRa, attitudeDesc)->group(group);
    cmd->add_option("--generate-de", opts.generateDe, attitudeDesc)->group(group);
    cmd->add_option("--generate-roll", opts.generateRoll, attitudeDesc)->group(group);
    cmd->add_flag("--generate-random-attitudes",
                  opts.generateRandomAttitudes,
                  "Generate a random attitude for each generated image. Ignores specified "
                  "ra, de, and roll if set.")
        ->group(group);

    cmd->add_option("--generate-exposure",
                    opts.generateExposure,
                    "Controls the exposure time of generated images. Default is 0.2. "
                    "Affects both motion blur and number of photons received per star.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);

    const std::string blurDesc =
        "Set the motion blur direction. By default, there is no motion blur (i.e., all of these "
        "are 0.0)";
    cmd->add_option("--generate-blur-ra", opts.generateBlurRa, blurDesc)->group(group);
    cmd->add_option("--generate-blur-de", opts.generateBlurDe, blurDesc)->group(group);
    cmd->add_option("--generate-blur-roll", opts.generateBlurRoll, blurDesc)->group(group);

    cmd->add_option(
           "--generate-readout-time",
           opts.generateReadoutTime,
           "Rolling shutter readout time for generated images. This is the time it takes from "
           "starting to capture the first row until finishing capturing the last row. Defaults to "
           "0 "
           "(rolling shutter simulation disabled).")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--generate-false-stars",
                    opts.generateNumFalseStars,
                    "The number of false stars generated. All are placed uniformly randomly "
                    "around the celestial sphere. Defaults to 0 (no false stars generated)")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--generate-false-min-mag",
                    opts.generateFalseMinMag,
                    "The minimum (weakest) magnitude that generated false stars should have. "
                    "Defaults to 8.0.")
        ->group(group);
    cmd->add_option("--generate-false-max-mag",
                    opts.generateFalseMaxMag,
                    "The maximum (strongest) magnitude that generated false stars should have. "
                    "Defaults to 1.0.")
        ->group(group);
    cmd->add_option("--generate-perturb-centroids",
                    opts.generatePerturbationStddev,
                    "Introduce artificial centroiding error. If provided, all the input and "
                    "expected centroids "
                    "will be shifted randomly according to a 2D Gaussian distribution with the "
                    "given standard "
                    "deviation (defaults to 0.2 pixel stddev). For evaluating star-id performance "
                    "vs centroid "
                    "error.")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option(
           "--generate-cutoff-mag",
           opts.generateCutoffMag,
           "When evaluating star-id algos on generated images, we don't want every single catalog "
           "star in the FOV to be fed into the star-id algorithm. This option lets you control "
           "that! Every catalog star has some probability of being turned into an input centroid, "
           "and this specifies the magnitude that will have a 50% probability of being included in "
           "the input centroid list. An overly complicated mathematical function determines how "
           "that probability changes as you increase and decrease the magnitude from there, but it "
           "falls to 0 or 1 pretty quickly at either side.\n"
           "Set to a value of 100 or greater to disable cutoff and include all centroids.")
        ->group(group);
    cmd->add_option("--generate-oversampling",
                    opts.generateOversampling,
                    "Perform num-samples many samples per pixel when generating images. "
                    "Effectively performs antialiasing. num-samples must be a perfect square. "
                    "Default 4.")
        ->group(group);
    cmd->add_option("--generate-seed",
                    opts.generateSeed,
                    "An integer seed to pass to the random number generator used during image "
                    "generation.")
        ->group(group);
    cmd->add_flag("--generate-time-based-seed",
                  opts.timeSeed,
                  "If true (the default), uses the current timestamp as a random seed. "
                  "--generate-seed is ignored if time based seed is set.")
        ->group(group);
}

static void AddComparatorOptions(CLI::App* cmd, PipelineOptions& opts) {
    const std::string group = "COMPARATOR OPTIONS";
    auto addOutput = [&](const std::string& name, std::string& variable, const std::string& desc) {
        AddOutputPathOption(cmd, name, variable, desc)->group(group);
    };

    addOutput("--plot-raw-input",
              opts.plotRawInput,
              "Argument is optional. Plots raw BW input image to PNG. Prints to path. Defaults to "
              "stdout.");
    addOutput("--plot-input",
              opts.plotInput,
              "Argument is optional. Plots annotated input image to PNG. Prints to path. Defaults "
              "to stdout.");
    addOutput("--plot-expected",
              opts.plotExpected,
              "Argument is optional. Plots annotated \"expected\" image to PNG. The difference "
              "between \"expected\" and \"input\" is that \"input\" is what's passed to the actual "
              "algorithms, while \"expected\" is compared against the actual output of the "
              "algorithms for evaluation purposes.");
    addOutput("--plot-output",
              opts.plotOutput,
              "Argument is optional. Plot output to PNG. Prints to path. Defaults to stdout.");
    addOutput("--plot-centroid-indices",
              opts.plotCentroidIndices,
              "Mainly for debugging. Argument is optional. Plot output to PNG, but instead of "
              "labeling the stars with the catalog star they were identified as, labels them with "
              "their centroid index.\n"
              "If a centroid algo was provided, prints those indices. Otherwise, prints the "
              "\"input centrtoid\" indices, which is whatever would be passed to the star-id algo "
              "(so for generated images, after perturbations and magnitude cutoff is applied).");
    addOutput("--print-actual-centroids",
              opts.printActualCentroids,
              "Argument is optional. Print list of centroids in the pipeline output to path. If "
              "star identification was also performed, print which catalog star each centroid was "
              "identified as. If input is from a generated image, also print the expected id for "
              "each centroid (note that there can sometimes be multiple acceptable star-ids for a "
              "star, if it's within the centroid comparison threshold of multiple catalog stars. "
              "Usually setting --min-separation sufficiently high prevents this). Defaults to "
              "stdout.");
    addOutput("--print-input-centroids",
              opts.printInputCentroids,
              "Argument is optional. Print list of centroids in the pipeline input to path. This "
              "only works for generated images. Prints whatever centroids would normally be sent "
              "to the star-id algorithm if provided in the pipeline. Also prints star-ids, but "
              "don't rely on them, because if the input centroids were perturbed the generated "
              "star-id may not actually be the closest catalog star anymore! Defaults to stdout.");
    addOutput("--print-expected-centroids",
              opts.printExpectedCentroids,
              "Argument is optional. Print list of \"expected\" centroids from the pipeline input "
              "to path. The \"expected\" centroids are all centroids used to generate the real "
              "image, including false stars and stars below the cutoff magnitude. Also prints out "
              "expected star identifications for those stars. This only works for generated "
              "images. Generally you should use --print-actual-centroids instead, and use "
              "--compare-centroids if you need to compare actual and expected centroids.");
    addOutput("--print-attitude",
              opts.printAttitude,
              "Argument is optional. Print the determined right ascension, declination, and roll "
              "to path. TODO: Option to print as quaternion. Defaults to stdout.");
    // Not documented in the old man page.
    addOutput("--print-expected-attitude",
              opts.printExpectedAttitude,
              "Argument is optional. Print the expected right ascension, declination, and roll to "
              "path. Only works for generated images. Defaults to stdout.");
    addOutput("--print-speed",
              opts.printSpeed,
              "Print the average, min, max, and upper 95th percentile of how long each stage of "
              "the pipeline took to run.");
    addOutput("--compare-centroids",
              opts.compareCentroids,
              "Argument is option. Compare expected to actual centroids. Prints to path. Defaults "
              "to stdout.");
    addOutput("--compare-star-ids",
              opts.compareStarIds,
              "Argument is optional. Compare expected to actual star identifiers. Prints to path. "
              "Defaults to stdout. In the output, starid_num_correct, starid_num_incorrect, and "
              "starid_num_total count the number of individual stars identified in an image. The "
              "starid_num_images_correct and starid_num_images_incorrect outputs are more "
              "aggregated: starid_num_images_correct records the number of images where at least "
              "two stars were correctly id'd and no stars were incorrectly id'd, while "
              "starid_num_images_incorrect records the number of images with at least one "
              "incorrectly id'd star (used when generating more than 1 image at a time, eg for "
              "evaluation).");
    addOutput("--compare-attitudes",
              opts.compareAttitudes,
              "Argument is optional. Compare expected to actual attitude. Prints to path. Defaults "
              "to stdout.");

    cmd->add_option("--centroid-compare-threshold",
                    opts.centroidCompareThreshold,
                    "Sets the threshold (in pixels) to consider two centroids equal (cf "
                    "--compare-centroids).")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
    cmd->add_option("--attitude-compare-threshold",
                    opts.attitudeCompareThreshold,
                    "Sets the threshold (in degrees) to consider two attitudes equal (cf "
                    "--compare-attitudes).")
        ->check(CLI::NonNegativeNumber)
        ->group(group);
}

static void AddPipelineOptions(CLI::App* cmd, PipelineOptions& opts) {
    CLI::Option* png = AddCameraOptions(cmd, opts);
    AddCentroidOptions(cmd, opts);
    AddStarIdOptions(cmd, opts);
    AddAttitudeOptions(cmd, opts);
    AddImageGenerationOptions(cmd, opts, png);
    AddComparatorOptions(cmd, opts);
}

static int LostMain(int argc, char** argv) {
    CLI::App app{"LOST: Open-source Star Tracker"};
    argv = app.ensure_utf8(argv);
    app.require_subcommand(0, 1);
    app.option_defaults()->always_capture_default();

    DatabaseOptions dbOptions;
    CLI::App* dbCmd = app.add_subcommand("database", "Build a star database");
    AddDatabaseOptions(dbCmd, dbOptions);

    PipelineOptions pipelineOptions;
    CLI::App* pipelineCmd = app.add_subcommand("pipeline", "Run a star-tracking pipeline");
    AddPipelineOptions(pipelineCmd, pipelineOptions);

    CLI11_PARSE(app, argc, argv);

    if (app.got_subcommand(pipelineCmd)) {
        PipelineRun(pipelineOptions);
    } else if (app.got_subcommand(dbCmd)) {
        DatabaseBuild(dbOptions);
    } else {
        std::cout << app.help() << std::endl;
    }

    return 0;
}

}  // namespace lost

int main(int argc, char** argv) {
    return lost::LostMain(argc, argv);
}
