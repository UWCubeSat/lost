/**
 * LOST entry point.
 *
 * This file parses CLI arguments and dispatches to the appropriate subsystem:
 *   - "database": Build a star catalog database and serialize it to a file.
 *   - "pipeline": Run a star-tracking pipeline (image generation, centroiding,
 *                 star identification, attitude estimation, and output comparison).
 */

#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>

#include <bitset>
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

/// Build a star database from the catalog and write it to the output path.
static void DatabaseBuild(const DatabaseOptions &values) {
    Catalog narrowedCatalog = NarrowCatalog(CatalogRead(), (int) (values.minMag * 100), values.maxStars, DegToRad(values.minSeparation));
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseDescriptor dbEntries = GenerateDatabases(narrowedCatalog, values);
    SerializeContext ser = serFromDbValues(values);

    // Build the flags word. Currently the only flag indicates whether the
    // database was built with single-precision (float) or double-precision decimals.
    uint32_t dbFlags = 0;
    dbFlags |= typeid(decimal) == typeid(float) ? MULTI_DB_FLOAT_FLAG : 0;

    SerializeMultiDatabase(&ser, dbEntries, dbFlags);

    std::cerr << "Generated database with " << ser.buffer.size() << " bytes" << std::endl;
    std::cerr << "Database flagged with " << std::bitset<8*sizeof(dbFlags)>(dbFlags) << std::endl;

    UserSpecifiedOutputStream pos = UserSpecifiedOutputStream(values.outputPath, true);
    pos.Stream().write((char *) ser.buffer.data(), ser.buffer.size());
}

/// Run a star-tracking pipeline and compare outputs against expected values.
static void PipelineRun(const PipelineOptions &values) {
    PipelineInputList input = GetPipelineInput(values);
    Pipeline pipeline = SetPipeline(values);
    std::vector<PipelineOutput> outputs = pipeline.Go(input);
    PipelineComparison(input, outputs, values);
}

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

#define LOST_OPTIONAL_OPTARG()                                   \
    ((optarg == NULL && optind < argc && argv[optind][0] != '-') \
     ? (bool) (optarg = argv[optind++])                          \
     : (optarg != NULL))

// ---------------------------------------------------------------------------
// CLI Parsing Helpers
// ---------------------------------------------------------------------------

static void PrintUsage() {
    std::cout << "Usage: ./lost database or ./lost pipeline" << std::endl
              << "Use --help flag on those commands for further help" << std::endl;
}


// Parse "database" subcommand options from the command line.
static int ParseDatabaseOptions(int argc, char **argv, DatabaseOptions &databaseOptions) {
    // Generate an enum with one value per CLI option
    // Expands to:
    // enum class DatabaseCliOption {
    //     minMag,
    //     maxStars,
    //     ...
    //     help
    // };
    enum class DatabaseCliOption {
        #define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) prop,
        #include "database-options.hpp"
        #undef LOST_CLI_OPTION
        help
    };

    // Generate options array
    // Expands to:
    // static struct option long_options[] = {
    //     {"min-mag", required_argument, 0, (int)DatabaseCliOption::minMag},
    //     {"max-stars", required_argument, 0, (int)DatabaseCliOption::maxStars},
    //     ...
    //     {0}
    // };
    static struct option long_options[] = {
        #define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
                {name,                                                             \
                defaultArg == 0 ? required_argument : optional_argument,            \
                0,                                                                  \
                (int)DatabaseCliOption::prop},
        #include "database-options.hpp" // NOLINT
        #undef LOST_CLI_OPTION
        {"help", no_argument, 0, (int) DatabaseCliOption::help},
        {0}
    };

    // Parse options
    // Expands to:
    // switch (option) {
    //     case (int)DatabaseCliOption::minMag:
    //         databaseOptions.minMag = STR_TO_DECIMAL(optarg);
    //         break;
    //     case (int)DatabaseCliOption::maxStars:
    //         databaseOptions.maxStars = atoi(optarg);
    //         break;
    //     ...
    // }
    int index;
    int option;
    while ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
        switch (option) {
            #define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
                case (int)DatabaseCliOption::prop:                                       \
                    if (defaultArg == 0) {                                               \
                        databaseOptions.prop = converter;                                \
                    } else {                                                             \
                        if (LOST_OPTIONAL_OPTARG()) {                                    \
                            databaseOptions.prop = converter;                            \
                        } else {                                                         \
                            databaseOptions.prop = defaultArg;                           \
                        }                                                                \
                    }                                                                    \
                    break;
            #include "database-options.hpp" // NOLINT
            #undef LOST_CLI_OPTION

            case (int) DatabaseCliOption::help:
                std::cout << documentation_database_txt << std::endl;
                return -1;
            default:
                std::cout << "Illegal flag" << std::endl;
                return 1;
        }
    }

    return 0;
}

// Parse "pipeline" subcommand options from the command line.
static int ParsePipelineOptions(int argc, char **argv, PipelineOptions &pipelineOptions) {
    // Generate an enum with one value per CLI option
    // Expands to:
    // enum class PipelineCliOption {
    //     png,
    //     focalLength,
    //     ...
    //     help
    // };
    enum class PipelineCliOption {
        #define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) prop,
        #include "pipeline-options.hpp"
        #undef LOST_CLI_OPTION
        help
    };

    // Generate options array
    // Expands to:
    // static struct option long_options[] = {
    //     {"png", required_argument, 0, (int)PipelineCliOption::png},
    //     {"focal-length", required_argument, 0, (int)PipelineCliOption::focalLength},
    //     ...
    //     {0}
    // };
    static struct option long_options[] = {
        #define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
                {name,                                                             \
                defaultArg == 0 ? required_argument : optional_argument,            \
                0,                                                                  \
                (int)PipelineCliOption::prop},
        #include "pipeline-options.hpp" // NOLINT
        #undef LOST_CLI_OPTION
        {"help", no_argument, 0, (int) PipelineCliOption::help},
        {0, 0, 0, 0}
    };

    // Parse options
    // Expands to:
    // switch (option) {
    //     case (int)PipelineCliOption::png:
    //         pipelineOptions.png = optarg;
    //         break;
    //     case (int)PipelineCliOption::focalLength:
    //         pipelineOptions.focalLength = atof(optarg);
    //         break;
    //     ...
    // }
    int index;
    int option;
    while ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
        switch (option) {
            #define LOST_CLI_OPTION(name, type, prop, defaultVal, converter, defaultArg) \
                case (int)PipelineCliOption::prop:                                       \
                    if (defaultArg == 0) {                                               \
                        pipelineOptions.prop = converter;                                \
                    } else {                                                             \
                        if (LOST_OPTIONAL_OPTARG()) {                                    \
                            pipelineOptions.prop = converter;                            \
                        } else {                                                         \
                            pipelineOptions.prop = defaultArg;                           \
                        }                                                                \
                    }                                                                    \
                    break;
            #include "pipeline-options.hpp" // NOLINT
            #undef LOST_CLI_OPTION

            case (int) PipelineCliOption::help:
                std::cout << documentation_pipeline_txt << std::endl;
                return -1;
            default:
                std::cout << "Illegal flag" << std::endl;
                return 1;
        }
    }

    return 0;
}

// ---------------------------------------------------------------------------
// Entry Point
// ---------------------------------------------------------------------------

/// Top-level dispatcher. Parses the subcommand and delegates to the appropriate handler.
static int LostMain(int argc, char **argv) {
    if (argc == 1) {
        PrintUsage();
        return 0;
    }

    std::string command(argv[1]);
    optind = 2; // skip program name and subcommand for getopt_long

    if (command == "database") {
        DatabaseOptions databaseOptions;
        int result = ParseDatabaseOptions(argc, argv, databaseOptions);
        if (result == -1) return 0;  // --help was printed
        if (result != 0)  return result;
        DatabaseBuild(databaseOptions);

    } else if (command == "pipeline") {
        PipelineOptions pipelineOptions;
        int result = ParsePipelineOptions(argc, argv, pipelineOptions);
        if (result == -1) return 0;  // --help was printed
        if (result != 0)  return result;
        PipelineRun(pipelineOptions);

    } else {
        PrintUsage();
    }

    return 0;
}

}

int main(int argc, char **argv) {
    return lost::LostMain(argc, argv);
}
