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

static void DatabaseBuild(DatabaseOptions values) {
    Catalog narrowedCatalog = NarrowCatalog(CatalogRead(),values.maxMagnitude,values.maxStars);
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseBuilder builder;
    // TODO: allow magnitude and weird
    unsigned char *catalogBuffer = builder.AddSubDatabase(kCatalogMagicValue,SerializeLengthCatalog(narrowedCatalog, false, true));
    SerializeCatalog(narrowedCatalog, false, true, catalogBuffer);

    GenerateDatabases(builder, narrowedCatalog, values);

    std::cerr << "Generated database with " << builder.BufferLength() << " bytes" << std::endl;
    PromptedOutputStream pos = PromptedOutputStream(values.path);
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

static void PipelineRun(std::map<std::string,std::string> values) {
    PipelineInputList input = GetPipelineInput(values);
    Pipeline pipeline = PromptPipeline();
    std::vector<PipelineOutput> outputs = pipeline.Go(input);
    PromptPipelineComparison(input, outputs);
}

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
            {"path", required_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        lost::DatabaseOptions databaseOptions;
        int index;
        int option;

        while (optind < argc) {
            if ((option = getopt_long(argc, argv, "m:s:ka:z:b:hp:", long_options, &index)) != -1) {
                switch (option) {
                    case 'm' :
                        std::cout << "You raised the maginude to " << optarg << std::endl;
                        databaseOptions.maxMagnitude = atoi(optarg);
                        break;
                    case 's' :
                        std::cout << "You lowered the stars to " << optarg << std::endl;
                        databaseOptions.maxStars = atoi(optarg);
                        break;
                    case 'k' :
                        std::cout << "You picked the kvector version! " << std::endl;
                        databaseOptions.databaseBuilder = "kvector";
                        break;
                    case 'a' :
                        databaseOptions.kvectorMinDistance = atof(optarg);
                        std::cout << "You set the min distance to " << optarg << std::endl;
                        break;
                    case 'z' :
                        std::cout << "You set the max distance to " << optarg << std::endl;
                        databaseOptions.kvectorMaxDistance = atof(optarg);
                        break;
                    case 'b' :
                        std::cout << "You set the number of bins to " << optarg << std::endl;
                        databaseOptions.kvectorDistanceBins = atol(optarg);
                        break;
                    case 'p' :
                        std::cout << "You set the path to " << optarg << std::endl;
                        databaseOptions.path = optarg;
                        break;
                    case 'h' :
                        system("man documentation/database.man");
                        return 0;
                        break;
                    default :
                        std::cout << "Illegal flag" << std::endl;
                        exit(1);
                }
            } 
        }

        lost::DatabaseBuild(databaseOptions);

        
    } else if (strcmp(argv[1], "pipeline") == 0) {

        static struct option long_options[] =
        {
            {"png",             required_argument, 0, 'a'},
            {"focal-length",    required_argument, 0, 'b'},
            {"pixel-size",      required_argument, 0, 'c'},
            {"fov",             required_argument, 0, 'd'},
            {"centroid-dummy",  optional_argument, 0, 'e'},
            {"centroid-cog",    no_argument, 0, 'f'},
            {"centroid-iwcog",  no_argument, 0, 'g'},
            {"centroid-mag-filter",  required_argument, 0, 'h'},
            {"database",        required_argument, 0, 'i'},
            {"id-dummy",        no_argument, 0, 'j'},
            {"id-gv",           required_argument, 0, 'k'},
            {"id-pyramid",      no_argument, 0, 'l'},
            {"py-tolerance",    required_argument, 0, 'm'},
            {"false-stars",     required_argument, 0, 'n'},
            {"max-mismatch-prob",  required_argument, 0, 'o'},
            {"attitude-dqm",    no_argument, 0, 'p'},
            {"plot",            required_argument, 0, 'q'},
            {"generate",        optional_argument, 0, 'r'},
            {"horizontal-res",  required_argument, 0, 's'},
            {"vertical-res",  required_argument, 0, 't'},
            {"horizontal-fov",  required_argument, 0, 'u'},
            {"ref-brightness-mag",  no_argument, 0, 'v'},
            {"spread-stddev",  required_argument, 0, 'w'},
            {"noise-stddev",  required_argument, 0, 'x'},
            {"boresight-right-asc",  required_argument, 0, 'y'},
            {"boresight-dec",  required_argument, 0, 'z'},
            {"boresight-roll",  required_argument, 0, '{'},
            {"help",            no_argument, 0, '}'},
            {0, 0, 0, 0}
        };

        // default values/flags
        std::map<std::string,std::string> parsedValues = {
            {"png",""},
            {"focal-length", "defaultTBD"},
            {"pixel-size","defaultTBD"},
            {"fov", "defaultTBD"},
            {"centroid-algo", "defaultTBD"},
            {"centroid-mag-filter", "defaultTBD"},
            {"database","defaultTBD"},
            {"id-algo", "defaultTBD"},
            {"attitude-dqm","defaultTBD"},
            {"plot","defaultTBD"},
            {"generate","1"},
            {"horizontal-res","defaultTBD"},
            {"vertical-res","defaultTBD"},
            {"horizontal-fov","defaultTBD"},
            {"ref-brightness-mag","defaultTBD"},
            {"spread-stddev","defaultTBD"},
            {"noise-stddev","defaultTBD"},
            {"boresight-right-asc","defaultTBD"},
            {"boresight-dec","defaultTBD"},
            {"boresight-roll","defaultTBD"}
        };

        int index;
        int option;

        while (optind < argc) {
            if ((option = getopt_long(argc, argv, "a:b:c:d:e::fgh:i:jk:lm:n:o:pq:r::s:t:u:vw:x:y:z:{:", long_options, &index)) != -1) {
                switch (option) {
                    case 'a' :
                        std::cout << "You made the file path " << optarg << std::endl;
                        parsedValues["png"] = optarg;
                        break;
                    case 'b' :
                        std::cout << "You set the focal length to " << optarg << std::endl;
                        parsedValues["focal-length"] = optarg;
                        break;
                    case 'c' :
                        std::cout << "You set the pixel size to " << optarg << std::endl;
                        parsedValues["pixel-size"] = optarg;
                        break;
                    case 'd' :
                        std::cout << "You set the fov to " << optarg << std::endl;
                        parsedValues["fov"] = optarg;
                        break;
                    case 'e' : 
                    {
                        std::string ns = "default";
                        if (optarg != NULL) {
                            ns = optarg;    // put default with the other defaults
                        } 
                        std::cout << "You set the centroid algo to dummy with " << ns << " stars." << std::endl;
                        parsedValues["centroid-algo"] = "dummy," + ns;
                        break; 
                    }
                    case 'f' :
                        std::cout << "You set the centroid algo to cog " << std::endl;
                        parsedValues["centroid-algo"] = "cog";
                        break;
                    case 'g' :
                        std::cout << "You set the centroid algo to iwcog " << std::endl;
                        parsedValues["centroid-algo"] = "iwcog";
                        break;
                    case 'h' :
                        std::cout << "You set the centroid mag filter to to " << optarg << std::endl;
                        parsedValues["centroid-mag-filter"] = optarg;
                        break;
                    case 'i' :
                        std::cout << "You set the database to " << optarg << std::endl;
                        parsedValues["database"] = optarg;
                        break;
                    case 'j' :
                        std::cout << "You set the id algo to dummy" << std::endl;
                        parsedValues["id-algo"] = "dummy";
                        break;
                    case 'k' :
                        std::cout << "You set the id algo to geometric voting" << std::endl;
                        parsedValues["id-algo"] = "gv";
                        break;                    
                    case 'l' :
                        std::cout << "You set the id algo to pyramid" << std::endl;
                        parsedValues["id-algo"] = "pyramid";
                        break;
                    case 'm' :
                        std::cout << "You set the id algo pyramid tolerance to " << optarg << std::endl;
                        parsedValues["id-algo"] += ",tol=";
                        parsedValues["id-algo"] += optarg;
                        break;
                    case 'n' :
                        std::cout << "You set the id algo false stars to " << optarg << std::endl;
                        parsedValues["id-algo"] += ",fs=";
                        parsedValues["id-algo"] += optarg;
                        break;
                    case 'o' :
                        std::cout << "You set the id algo max mismatch probability to " << optarg << std::endl;
                        parsedValues["id-algo"] += ",prob=";
                        parsedValues["id-algo"] += optarg;
                        break;
                    case 'p' :
                        parsedValues["attitude-dqm"] = "true";
                        break;
                    case 'q' :
                        std::cout << "You set the plotted output path to " << optarg << std::endl;
                        parsedValues["plot"] = optarg;
                        break;
                    case 'r' :
                        std::cout << "Generating images! " << optarg << std::endl;
                        if (optarg != NULL) parsedValues["generate"] = optarg;
                        break;
                    case 's' :
                        std::cout << "You set the horizontal res to " << optarg << std::endl;
                        parsedValues["horizontal-res"] = optarg;
                        break;
                    case 't' :
                        std::cout << "You set the vertical res to " << optarg << std::endl;
                        parsedValues["vertical-res"] = optarg;
                        break;
                    case 'u' :
                        std::cout << "You set the horizontal fov to " << optarg << std::endl;
                        parsedValues["horizontal-fov"] = optarg;
                        break;
                    case 'v' :
                        std::cout << "You have a ref brightness magnitude " ;
                        parsedValues["ref-brightness-mag"] = "true";                        
                        break;
                    case 'w' :
                        std::cout << "You set the spread stddev to " << optarg << std::endl;
                        parsedValues["spread-stddev"] = optarg;
                        break;
                    case 'x' :
                        std::cout << "You set the noise stddev to " << optarg << std::endl;
                        parsedValues["noise-stddev"] = optarg;
                        break;
                    case 'y' :
                        std::cout << "You set the boresight right asc to " << optarg << std::endl;
                        parsedValues["boresight-right-asc"] = optarg;
                        break;
                    case 'z' :
                        std::cout << "You set the boresight declination to " << optarg << std::endl;
                        parsedValues["boresight-dec"] = optarg;
                        break;
                    case '{' :
                        std::cout << "You set the boresight roll to " << optarg << std::endl;
                        parsedValues["boresight-roll"] = optarg;
                        break;
                    case '}': 
                        system("man documentation/pipeline.man");
                        return 0;
                        break;
                    default :
                        std::cout << "Illegal flag" << std::endl;
                        exit(1);
                }
            } 
        }

        lost::PipelineRun(parsedValues);

    } else {
        std::cout << "Unrecognized command" << std::endl;
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