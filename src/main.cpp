#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <getopt.h>
#include <cstring>
#include <map>

#include "databases.hpp"
#include "centroiders.hpp"
#include "io.hpp"

namespace lost {

static void DatabaseBuild(std::map<std::string,std::string> values) {
    Catalog narrowedCatalog = NarrowCatalog(CatalogRead(),stoi(values["maxMagnitude"]),stoi(values["maxStars"]));
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseBuilder builder;
    // TODO: allow magnitude and weird
    unsigned char *catalogBuffer = builder.AddSubDatabase(kCatalogMagicValue,SerializeLengthCatalog(narrowedCatalog, false, true));
    SerializeCatalog(narrowedCatalog, false, true, catalogBuffer);

    GenerateDatabases(builder, narrowedCatalog, values);

    std::cerr << "Generated database with " << builder.BufferLength() << " bytes" << std::endl;
    PromptedOutputStream pos = PromptedOutputStream(values["path"]);
    pos.Stream().write((char *)builder.Buffer(), builder.BufferLength());
}


// static void DatabaseBuild(int maxMagnitude, int maxStars) {
//     Catalog narrowedCatalog = NarrowCatalog(CatalogRead(),maxMagnitude,maxStars);
//     std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

//     MultiDatabaseBuilder builder;
//     unsigned char *catalogBuffer = builder.AddSubDatabase(kCatalogMagicValue,
//                                                           // TODO: allow magnitude and weird
//                                                           SerializeLengthCatalog(narrowedCatalog, false, true));
//     SerializeCatalog(narrowedCatalog, false, true, catalogBuffer);

//     PromptDatabases(builder, narrowedCatalog);

//     std::cerr << "Generated database with " << builder.BufferLength() << " bytes" << std::endl;
//     PromptedOutputStream pos;
//     pos.Stream().write((char *)builder.Buffer(), builder.BufferLength());
// }

// static void PipelineRun() {
//     PipelineInputList input = PromptPipelineInput();
//     Pipeline pipeline = PromptPipeline();
//     std::vector<PipelineOutput> outputs = pipeline.Go(input);
//     PromptPipelineComparison(input, outputs);
// }

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

}

int main(int argc, char **argv) {

    if (strcmp(argv[1], "database") == 0) {

        static struct option long_options[] =
        {
            {"mag",     required_argument, 0, 'm'},
            {"stars",   required_argument, 0, 's'},
            {"kvector", no_argument,       0, 'k'},
            {"kvector-min-distance", required_argument, 0, 'a'},
            {"kvector-max-distance",  required_argument, 0, 'z'},
            {"kvector-distance-bins",  required_argument, 0, 'b'},
            {"help",  no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        // default values/flags
        std::map<std::string,std::string> parsedValues = {
            {"maxMagnitude","1000"},
            {"maxStars","10000"},
            {"databaseBuilder",""},
            {"kvector-min-distance",std::to_string(lost::DegToRad(0.5))},
            {"kvector-max-distance",std::to_string(lost::DegToRad(15))},
            {"kvector-distance-bins","10000"},
            {"path",""}
        };

        int index;
        int option;
        while (optind < argc) {
            if ((option = getopt_long(argc, argv, "m:s:ka:z:b:h", long_options, &index)) != -1) {
                switch (option) {
                    case 'm' :
                        std::cout << "You raised the maginude to " << optarg << std::endl;
                        parsedValues["maxMagnitude"] = optarg;
                        break;
                    case 's' :
                        std::cout << "You lowered the stars to " << optarg << std::endl;
                        parsedValues["maxStars"] = optarg;
                        break;
                    case 'k' :
                        std::cout << "You picked the kvector version! " << std::endl;
                        parsedValues["databaseBuilder"] = "kvector";
                        break;
                    case 'a' :
                        std::cout << "You set the min distance to  " << optarg << std::endl;
                        parsedValues["kvector-min-distance"] = optarg;
                        break;
                    case 'z' :
                        std::cout << "You set the max distance to  " << optarg << std::endl;
                        parsedValues["kvector-max-distance"] = optarg;
                        break;
                    case 'b' :
                        std::cout << "You set the number of bins to  " << optarg << std::endl;
                        parsedValues["kvector-distance-bins"] = optarg;
                        break;
                    case 'h' :
                        system("man documentation/database.man");
                        break;
                    default :
                        std::cout << "Illegal flag" << std::endl;
                        exit(1);
                }
            } else {
                parsedValues["path"] = optarg;
                std::cout << "You set the path to  " << optarg << std::endl;

                optind++;
            }
        }

        lost::DatabaseBuild(parsedValues);
        
    }
    return 0;
}

    // lost::RegisterCliArgs(argc, argv);
    // std::cerr << "LOST: Open-source Star Tracker" << std::endl;
    // lost::InteractiveChoice<void (*)()> mainChoices;
    // mainChoices.Register("pipeline", "Run a pipeline", &lost::PipelineRun);
    // mainChoices.Register("benchmark", "Benchmark a pipeline", &lost::PipelineBenchmark);
    // mainChoices.Register("build_database", "Build database from catalog", &lost::DatabaseBuild);
    // (*mainChoices.Prompt("Choose action"))();