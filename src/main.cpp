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
    if (values.path == "stdout") { //TODO make sure I understand this correctly? maybe use isatty function or something?
        std::cout << "Warning: output contains binary contents. Not printed to terminal." << std::endl;
    } else {
        PromptedOutputStream pos = PromptedOutputStream(values.path);
        pos.Stream().write((char *)builder.Buffer(), builder.BufferLength());
    }
    
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

static void PipelineRun(PipelineOptions values) {
    PipelineInputList input = GetPipelineInput(values);
    Pipeline pipeline = SetPipeline(values);
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

    if (argc == 1 || argc == 2) {
        std::cout << "Usage: ./lost database or ./lost pipeline" << std::endl << "Use --help flag on those commands for further help" << std::endl; 
        return 0;
    } 

    std::string command (argv[1]);
    if (command == "database") {

        enum DatabaseEnum {mag, stars, kvector, kvectorMinDistance, kvectorMaxDistance, 
            kvectorDistanceBins, help, output};

        static struct option long_options[] =
        {
            {"mag",     required_argument, 0, mag},
            {"stars",   required_argument, 0, stars},
            {"kvector", no_argument,       0, kvector},
            {"kvector-min-distance", required_argument, 0, kvectorMinDistance},
            {"kvector-max-distance",  required_argument, 0, kvectorMaxDistance},
            {"kvector-distance-bins",  required_argument, 0, kvectorDistanceBins},
            {"help",  no_argument, 0, help},
            {"output", required_argument, 0, output},
            {0, 0, 0, 0}
        };

        lost::DatabaseOptions databaseOptions;
        int index;
        int option;

        while (optind < argc) {
            if ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
                switch (option) {
                    case mag :
                        std::cout << "You raised the maginude to " << optarg << std::endl;
                        databaseOptions.maxMagnitude = atoi(optarg);
                        break;
                    case stars :
                        std::cout << "You lowered the stars to " << optarg << std::endl;
                        databaseOptions.maxStars = atoi(optarg);
                        break;
                    case kvector :
                        std::cout << "You picked the kvector version! " << std::endl;
                        databaseOptions.databaseBuilder = "kvector";
                        break;
                    case kvectorMinDistance :
                        databaseOptions.kvectorMinDistance = atof(optarg);
                        std::cout << "You set the min distance to " << optarg << std::endl;
                        break;
                    case kvectorMaxDistance :
                        std::cout << "You set the max distance to " << optarg << std::endl;
                        databaseOptions.kvectorMaxDistance = atof(optarg);
                        break;
                    case kvectorDistanceBins :
                        std::cout << "You set the number of bins to " << optarg << std::endl;
                        databaseOptions.kvectorDistanceBins = atol(optarg);
                        break;
                    case output :
                        std::cout << "You set the path to " << optarg << std::endl;
                        databaseOptions.path = optarg;
                        break;
                    case help :
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

    } else if (command == "pipeline") {

        enum PipelineEnum {png, focalLength, pixelSize, fov, centroidAlgo, centroidDummyStars, centroidMagFilter, database, idAlgo,
            gvTolerance, pyTolerance, falseStars, maxMismatchProb, attitudeAlgo, plot, generate, horizontalRes, verticalRes, refBrightnessMag,
            spreadStddev, noiseStddev, boresightRightAsc, boresightDec,boresightRoll, help};

        static struct option long_options[] =
        {
            {"png",             required_argument, 0, png},
            {"focal-length",    required_argument, 0, focalLength},
            {"pixel-size",      required_argument, 0, pixelSize},
            {"fov",             required_argument, 0, fov},
            {"centroid-algo", required_argument, 0, centroidAlgo},
            {"centroid-dummy-stars", required_argument, 0, centroidDummyStars},
            {"centroid-mag-filter",  required_argument, 0, centroidMagFilter},
            {"database",        required_argument, 0, database},
            {"id-algo",        required_argument, 0, idAlgo},
            {"gv-tolerance",    required_argument, 0, gvTolerance},
            {"py-tolerance",    required_argument, 0, pyTolerance},
            {"false-stars",     required_argument, 0, falseStars},
            {"max-mismatch-prob",  required_argument, 0, maxMismatchProb},
            {"attitude-algo",    required_argument, 0, attitudeAlgo},
            {"plot",            required_argument, 0, plot},
            {"generate",        optional_argument, 0, generate},
            {"horizontal-res",  required_argument, 0, horizontalRes},
            {"vertical-res",  required_argument, 0, verticalRes},
            {"ref-brightness-mag",  no_argument, 0, refBrightnessMag},
            {"spread-stddev",  required_argument, 0, spreadStddev},
            {"noise-stddev",  required_argument, 0, noiseStddev},
            {"boresight-right-asc",  required_argument, 0, boresightRightAsc},
            {"boresight-dec",  required_argument, 0, boresightDec},
            {"boresight-roll",  required_argument, 0, boresightRoll},
            {"help",            no_argument, 0, 25},
            {0, 0, 0, 0}
        };

        lost::PipelineOptions pipelineOptions;
        int index;
        int option;
        
        while (optind < argc) {
            if ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
                switch (option) {
                    case png :
                        std::cout << "You made the file path " << optarg << std::endl;
                        pipelineOptions.png = optarg;
                        break;
                    case focalLength :
                        std::cout << "You set the focal length to " << optarg << std::endl;
                        pipelineOptions.focalLength = atof(optarg);
                        break;
                    case pixelSize :
                        std::cout << "You set the pixel size to " << optarg << std::endl;
                        pipelineOptions.pixelSize = atof(optarg);
                        break;
                    case fov :
                        std::cout << "You set the fov to " << optarg << std::endl;
                        pipelineOptions.fov = atof(optarg);
                        break;
                    case centroidAlgo :
                    {
                        std::string algo (optarg);
                        if (algo != "dummy" && algo != "cog" && algo != "iwcog") {
                            std::cout << "Unrecognized centroid algorithm!" << std::endl;
                        } else {
                            std::cout << algo << "!" << std::endl;
                            pipelineOptions.centroidAlgo = algo;
                        }                        
                        break;
                    }
                    case centroidDummyStars : 
                        pipelineOptions.dummyCentroidNumStars = atoi(optarg);
                        break;
                    case centroidMagFilter :
                        std::cout << "You set the centroid mag filter to to " << optarg << std::endl;
                        pipelineOptions.centroidMagFilter = atoi(optarg);
                        break;
                    case database :
                        std::cout << "You set the database to " << optarg << std::endl;
                        pipelineOptions.database = optarg;
                        break;
                    case idAlgo :
                    {
                        std::string algo (optarg);
                        if (algo != "dummy" && algo != "gv" && algo != "pyramid") {
                            std::cout << "Unrecognized id algorithm!" << std::endl;
                        } else {
                            std::cout << algo << "!" << std::endl;
                            pipelineOptions.idAlgo = algo;
                        }                        
                        break;
                    }
                    case gvTolerance :
                        pipelineOptions.gvTolerance = atof(optarg);
                        break;
                    case pyTolerance :
                        std::cout << "You set the id algo pyramid tolerance to " << optarg << std::endl;
                        pipelineOptions.pyTolerance = atof(optarg);
                        break;
                    case falseStars :
                        std::cout << "You set the id algo false stars to " << optarg << std::endl;
                        pipelineOptions.pyFalseStars = atoi(optarg);
                        break;
                    case maxMismatchProb :
                        std::cout << "You set the id algo max mismatch probability to " << optarg << std::endl;
                        pipelineOptions.pyMismatchProb = atof(optarg);
                        break;
                    case attitudeAlgo :
                    {
                        std::string algo (optarg);
                        if (algo != "dqm") {
                            std::cout << "Unrecognized id algorithm!" << std::endl;
                        } else {
                            std::cout << algo << "!" << std::endl;
                            pipelineOptions.attitudeAlgo = algo;
                        }                        
                        break; 
                    }
                    case plot :
                        std::cout << "You set the plotted output path to " << optarg << std::endl;
                        pipelineOptions.plot = optarg;
                        break;
                    case generate :
                        std::cout << "Generating images! " << optarg << std::endl;
                        if (optarg != NULL) pipelineOptions.generate = atoi(optarg);
                        break;
                    case horizontalRes :
                        std::cout << "You set the horizontal res to " << optarg << std::endl;
                        pipelineOptions.horizontalRes = atoi(optarg);
                        break;
                    case verticalRes :
                        std::cout << "You set the vertical res to " << optarg << std::endl;
                        pipelineOptions.verticalRes = atoi(optarg);
                        break;
                    case refBrightnessMag :
                        std::cout << "You have a ref brightness magnitude " ;
                        pipelineOptions.referenceBrightness = atoi(optarg);                        
                        break;
                    case spreadStddev :
                        std::cout << "You set the spread stddev to " << optarg << std::endl;
                        pipelineOptions.brightnessDeviation = atof(optarg);
                        break;
                    case noiseStddev :
                        std::cout << "You set the noise stddev to " << optarg << std::endl;
                        pipelineOptions.noiseDeviation = atof(optarg);
                        break;
                    case boresightRightAsc :
                        std::cout << "You set the boresight right asc to " << optarg << std::endl;
                        pipelineOptions.ra = atof(optarg);
                        break;
                    case boresightDec :
                        std::cout << "You set the boresight declination to " << optarg << std::endl;
                        pipelineOptions.dec = atof(optarg);
                        break;
                    case boresightRoll :
                        std::cout << "You set the boresight roll to " << optarg << std::endl;
                        pipelineOptions.roll = atof(optarg);
                        break;
                    case help : 
                        system("man documentation/pipeline.man");
                        return 0;
                        break;
                    default :
                        std::cout << "Illegal flag" << std::endl;
                        exit(1);
                }
            } 
        }

        lost::PipelineRun(pipelineOptions);

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