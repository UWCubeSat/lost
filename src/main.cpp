/**
 * LOST starting point
 *
 * Reads in CLI arguments/flags and starts the appropriate pipelines
 */

#include <assert.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstring>
#include <map>

#include "databases.hpp"
#include "centroiders.hpp"
#include "decimal.hpp"
#include "io.hpp"
#include "man-database.h"
#include "man-pipeline.h"

namespace lost {

/// Create a database and write it to a file based on the command line options in \p values
static void DatabaseBuild(const DatabaseOptions &values) {
    Catalog narrowedCatalog = NarrowCatalog(CatalogRead(), (int) (values.minMag * 100), values.maxStars, DegToRad(values.minSeparation));
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseDescriptor dbEntries = GenerateDatabases(narrowedCatalog, values);
    SerializeContext ser = serFromDbValues(values);
    // Inject flags into the Serialized Database.
    uint32_t dbFlags = typeid(decimal) == typeid(double) ? MULTI_DB_IS_DOUBLE : MULTI_DB_IS_FLOAT;
    SerializeMultiDatabase(&ser, dbEntries, dbFlags);

    std::cerr << "Generated database with " << ser.buffer.size() << " bytes" << std::endl;
    std::cerr << "Database flagged with " << dbFlags << std::endl;

    UserSpecifiedOutputStream pos = UserSpecifiedOutputStream(values.outputPath, true);
    pos.Stream().write((char *) ser.buffer.data(), ser.buffer.size());

}

/// Run a star-tracking pipeline (possibly including generating inputs and analyzing outputs) based on command line options in \p values.
static void PipelineRun(const PipelineOptions &values) {
    PipelineInputList input = GetPipelineInput(values);
    Pipeline pipeline = SetPipeline(values);
    std::vector<PipelineOutput> outputs = pipeline.Go(input);
    PipelineComparison(input, outputs, values);
}

// DO NOT DELETE
// static void PipelineBenchmark() {
//     PipelineInputList input = PromptPipelineInput();
//     Pipeline pipeline = PromptPipeline();
//     int iterations = Prompt<int>("Times to run the pipeline");
//     std::cerr << "Benchmarking..." << std::endl;

//     // TODO: we can do better than this :| maybe include mean time, 99% time, or allow a vector of
//     // input and determine which one took the longest
//     auto startTime = std::chrono::high_resolution_clock::now();
//     for (int i = 0; i < iterations; i++) {
//         pipeline.Go(input);
//     }
//     auto endTime = std::chrono::high_resolution_clock::now();
//     auto totalTime = std::chrono::duration<double, std::milli>(endTime - startTime);
//     std::cout << "total_ms " << totalTime.count() << std::endl;
// }

// static void EstimateCamera() {
//     std::cerr << "Enter estimated camera details when prompted." << std::endl;
//     PipelineInputList inputs = PromptPngPipelineInput();
//     float baseFocalLength = inputs[0]->InputCamera()->FocalLength();
//     float deviationIncrement = Prompt<float>("Focal length increment (base: " + std::to_string(baseFocalLength) + ")");
//     float deviationMax = Prompt<float>("Maximum focal length deviation to attempt");
//     Pipeline pipeline = PromptPipeline();

//     while (inputs[0]->InputCamera()->FocalLength() - baseFocalLength <= deviationMax) {
//         std::cerr << "Attempt focal length " << inputs[0]->InputCamera()->FocalLength() << std::endl;
//         std::vector<PipelineOutput> outputs = pipeline.Go(inputs);
//         if (outputs[0].nice) {
//             std::cout << "camera_identified true" << std::endl << *inputs[0]->InputCamera();
//             return;
//         }

//         Camera camera(*inputs[0]->InputCamera());
//         if (camera.FocalLength() - baseFocalLength > 0) {
//             // yes i know this expression can be simplified shut up
//             camera.SetFocalLength(camera.FocalLength() - 2*(camera.FocalLength() - baseFocalLength));
//         } else {
//             camera.SetFocalLength(camera.FocalLength() + 2*(baseFocalLength - camera.FocalLength()) + deviationIncrement);
//         }
//         ((PngPipelineInput *)(inputs[0].get()))->SetCamera(camera);
//     }
//     std::cout << "camera_identified false" << std::endl;
// }

/// Convert string to boolean
bool atobool(const char *cstr) {
    std::string str(cstr);
    if (str == "1" || str == "true") {
        return true;
    }
    if (str == "0" || str == "false") {
        return false;
    }
    assert(false);
}

/**
 * Handle optional CLI arguments
 * https://stackoverflow.com/a/69177115
 */
#define LOST_OPTIONAL_OPTARG()                                   \
    ((optarg == NULL && optind < argc && argv[optind][0] != '-') \
     ? (bool) (optarg = argv[optind++])                          \
     : (optarg != NULL))

// This is separate from `main` just because it's in the `lost` namespace
static int LostMain(int argc, char **argv) {

    if (argc == 1) {
        std::cout << "Usage: ./lost database or ./lost pipeline" << std::endl
                  << "Use --help flag on those commands for further help" << std::endl;
        return 0;
    }

    std::string command(argv[1]);
    optind = 2;

    if (command == "database") {

        enum class DatabaseCliOption {
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) prop,
#include "database-options.hpp"
#undef LOST_CLI_OPTION
            help
        };

        static struct option long_options[] = {
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
            {name,                                                      \
             defaultArg == 0 ? required_argument : optional_argument, \
             0,                                                         \
             (int)DatabaseCliOption::prop},
#include "database-options.hpp" // NOLINT
#undef LOST_CLI_OPTION
                {"help", no_argument, 0, (int) DatabaseCliOption::help},
                {0}
        };

        DatabaseOptions databaseOptions;
        int index;
        int option;

        while ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
            switch (option) {
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
                case (int)DatabaseCliOption::prop :                     \
                    if (defaultArg == 0) {     \
                        databaseOptions.prop = converter;       \
                    } else {                                    \
                        if (LOST_OPTIONAL_OPTARG()) {           \
                            databaseOptions.prop = converter;   \
                        } else {                                \
                            databaseOptions.prop = defaultArg;  \
                        }                                       \
                    }                                           \
            break;
#include "database-options.hpp" // NOLINT
#undef LOST_CLI_OPTION
                case (int) DatabaseCliOption::help :std::cout << documentation_database_txt << std::endl;
                    return 0;
                    break;
                default :std::cout << "Illegal flag" << std::endl;
                    exit(1);
            }
        }

        lost::DatabaseBuild(databaseOptions);

    } else if (command == "pipeline") {

        enum class PipelineCliOption {
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) prop,
#include "pipeline-options.hpp"
#undef LOST_CLI_OPTION
            help
        };

        static struct option long_options[] = {
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
            {name,                                                      \
             defaultArg == 0 ? required_argument : optional_argument, \
             0,                                                         \
             (int)PipelineCliOption::prop},
#include "pipeline-options.hpp" // NOLINT
#undef LOST_CLI_OPTION

                // DATABASES
                {"help", no_argument, 0, (int) PipelineCliOption::help},
                {0, 0, 0, 0}
        };

        lost::PipelineOptions pipelineOptions;
        int index;
        int option;

        while ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
            switch (option) {
#define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
                case (int)PipelineCliOption::prop :                         \
                    if (defaultArg == 0) {    \
                        pipelineOptions.prop = converter;       \
                    } else {                                    \
                        if (LOST_OPTIONAL_OPTARG()) {           \
                            pipelineOptions.prop = converter;   \
                        } else {                                \
                            pipelineOptions.prop = defaultArg;  \
                        }                                       \
                    }                                           \
            break;
#include "pipeline-options.hpp" // NOLINT
#undef LOST_CLI_OPTION
                case (int) PipelineCliOption::help :std::cout << documentation_pipeline_txt << std::endl;
                    return 0;
                    break;
                default :std::cout << "Illegal flag" << std::endl;
                    exit(1);
            }
        }

        lost::PipelineRun(pipelineOptions);

    } else {
        std::cout << "Usage: ./lost database or ./lost pipeline" << std::endl
                  << "Use --help flag on those commands for further help" << std::endl;
    }
    return 0;
}

}

int main(int argc, char **argv) {
    return lost::LostMain(argc, argv);
}
