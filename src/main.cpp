/**
 * LOST starting point
 *
 * Reads in CLI arguments/flags and starts the appropriate pipelines
 */

#include <bitset>
#include <iostream>

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

static int LostMain(int argc, char** argv) {
    CLI::App app{"LOST: Open-source Star Tracker"};
    argv = app.ensure_utf8(argv);
    app.require_subcommand(0, 1);

    // --- database subcommand ---
    DatabaseOptions dbOptions;
    auto* dbCmd = app.add_subcommand("database", "Build a star database");

    dbCmd->add_option("--min-mag", dbOptions.minMag)->capture_default_str();
    dbCmd->add_option("--max-stars", dbOptions.maxStars)->capture_default_str();
    dbCmd->add_option("--min-separation", dbOptions.minSeparation)->capture_default_str();
    dbCmd->add_option("--kvector", dbOptions.kvector)->expected(0, 1)->default_str("true");
    dbCmd->add_option("--kvector-min-distance", dbOptions.kvectorMinDistance)
        ->capture_default_str();
    dbCmd->add_option("--kvector-max-distance", dbOptions.kvectorMaxDistance)
        ->capture_default_str();
    dbCmd->add_option("--kvector-distance-bins", dbOptions.kvectorNumDistanceBins)
        ->capture_default_str();
    dbCmd->add_option("--swap-integer-endianness", dbOptions.swapIntegerEndianness)
        ->expected(0, 1)
        ->default_str("true");
    dbCmd->add_option("--swap-decimal-endianness", dbOptions.swapDecimalEndianness)
        ->expected(0, 1)
        ->default_str("true");
    dbCmd->add_option("--output", dbOptions.outputPath)->capture_default_str();

    // --- pipeline subcommand ---
    PipelineOptions pipelineOptions;
    auto* pipelineCmd = app.add_subcommand("pipeline", "Run a star-tracking pipeline");
    pipelineCmd->validate_optional_arguments();

    // Camera
    pipelineCmd->add_option("--png", pipelineOptions.png);
    pipelineCmd->add_option("--focal-length", pipelineOptions.focalLength)->capture_default_str();
    pipelineCmd->add_option("--pixel-size", pipelineOptions.pixelSize)->capture_default_str();
    pipelineCmd->add_option("--fov", pipelineOptions.fov)->capture_default_str();

    // Pipeline stages
    pipelineCmd->add_option("--centroid-algo", pipelineOptions.centroidAlgo)
        ->expected(0, 1)
        ->default_str("cog");
    pipelineCmd->add_option("--centroid-dummy-stars", pipelineOptions.centroidDummyNumStars)
        ->capture_default_str();
    pipelineCmd->add_option("--centroid-mag-filter", pipelineOptions.centroidMagFilter)
        ->expected(0, 1)
        ->default_str("5");
    pipelineCmd->add_option("--centroid-filter-brightest", pipelineOptions.centroidFilterBrightest)
        ->expected(0, 1)
        ->default_str("10");
    pipelineCmd->add_option("--database", pipelineOptions.databasePath);
    pipelineCmd->add_option("--star-id-algo", pipelineOptions.idAlgo)
        ->expected(0, 1)
        ->default_str("pyramid");
    pipelineCmd->add_option("--angular-tolerance", pipelineOptions.angularTolerance)
        ->capture_default_str();
    pipelineCmd
        ->add_option("--false-stars-estimate,--false-stars", pipelineOptions.estimatedNumFalseStars)
        ->capture_default_str();
    pipelineCmd
        ->add_option("--max-mismatch-probability,--max-mismatch-prob",
                     pipelineOptions.maxMismatchProb)
        ->capture_default_str();
    pipelineCmd->add_option("--attitude-algo", pipelineOptions.attitudeAlgo)
        ->expected(0, 1)
        ->default_str("dqm");

    // Output comparison
    pipelineCmd
        ->add_option("--centroid-compare-threshold", pipelineOptions.centroidCompareThreshold)
        ->capture_default_str();
    pipelineCmd
        ->add_option("--attitude-compare-threshold", pipelineOptions.attitudeCompareThreshold)
        ->capture_default_str();
    pipelineCmd->add_option("--plot-raw-input", pipelineOptions.plotRawInput)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--plot-input", pipelineOptions.plotInput)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--plot-expected", pipelineOptions.plotExpected)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--plot-centroid-indices", pipelineOptions.plotCentroidIndices)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--plot-output", pipelineOptions.plotOutput)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--print-expected-centroids", pipelineOptions.printExpectedCentroids)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--print-input-centroids", pipelineOptions.printInputCentroids)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--print-actual-centroids", pipelineOptions.printActualCentroids)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--print-attitude", pipelineOptions.printAttitude)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--print-expected-attitude", pipelineOptions.printExpectedAttitude)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--print-speed", pipelineOptions.printSpeed)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--compare-centroids", pipelineOptions.compareCentroids)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--compare-star-ids", pipelineOptions.compareStarIds)
        ->expected(0, 1)
        ->default_str("-");
    pipelineCmd->add_option("--compare-attitudes", pipelineOptions.compareAttitudes)
        ->expected(0, 1)
        ->default_str("-");

    // Image generation
    pipelineCmd->add_option("--generate", pipelineOptions.generate)
        ->expected(0, 1)
        ->default_str("1");
    pipelineCmd->add_option("--generate-x-resolution", pipelineOptions.generateXRes)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-y-resolution", pipelineOptions.generateYRes)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-centroids-only", pipelineOptions.generateCentroidsOnly)
        ->expected(0, 1)
        ->default_str("true");
    pipelineCmd->add_option("--generate-zero-mag-photons", pipelineOptions.generateZeroMagPhotons)
        ->capture_default_str();
    pipelineCmd
        ->add_option("--generate-saturation-photons", pipelineOptions.generateSaturationPhotons)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-spread-stddev", pipelineOptions.generateSpreadStdDev)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-shot-noise", pipelineOptions.generateShotNoise)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-dark-current", pipelineOptions.generateDarkCurrent)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-read-noise-stddev", pipelineOptions.generateReadNoiseStdDev)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-ra", pipelineOptions.generateRa)->capture_default_str();
    pipelineCmd->add_option("--generate-de", pipelineOptions.generateDe)->capture_default_str();
    pipelineCmd->add_option("--generate-roll", pipelineOptions.generateRoll)->capture_default_str();
    pipelineCmd->add_option("--generate-random-attitudes", pipelineOptions.generateRandomAttitudes)
        ->expected(0, 1)
        ->default_str("true");
    pipelineCmd->add_option("--generate-blur-ra", pipelineOptions.generateBlurRa)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-blur-de", pipelineOptions.generateBlurDe)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-blur-roll", pipelineOptions.generateBlurRoll)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-exposure", pipelineOptions.generateExposure)
        ->expected(0, 1)
        ->default_str("0.1");
    pipelineCmd->add_option("--generate-readout-time", pipelineOptions.generateReadoutTime)
        ->expected(0, 1)
        ->default_str("0.01");
    pipelineCmd->add_option("--generate-oversampling", pipelineOptions.generateOversampling)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-false-stars", pipelineOptions.generateNumFalseStars)
        ->expected(0, 1)
        ->default_str("50");
    pipelineCmd->add_option("--generate-false-min-mag", pipelineOptions.generateFalseMinMag)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-false-max-mag", pipelineOptions.generateFalseMaxMag)
        ->capture_default_str();
    pipelineCmd
        ->add_option("--generate-perturb-centroids", pipelineOptions.generatePerturbationStddev)
        ->expected(0, 1)
        ->default_str("0.2");
    pipelineCmd->add_option("--generate-cutoff-mag", pipelineOptions.generateCutoffMag)
        ->capture_default_str();
    pipelineCmd->add_option("--generate-seed", pipelineOptions.generateSeed)->capture_default_str();
    pipelineCmd->add_option("--generate-time-based-seed", pipelineOptions.timeSeed)
        ->expected(0, 1)
        ->default_str("true");

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
