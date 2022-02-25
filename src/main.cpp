#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <getopt.h>
#include <cstring>
#include <map>
#include <unistd.h>

#include "databases.hpp"
#include "centroiders.hpp"
#include "io.hpp"
#include "man-database.h"
#include "man-pipeline.h"

namespace lost {

static void DatabaseBuild(const DatabaseOptions &values) {
    Catalog narrowedCatalog = NarrowCatalog(CatalogRead(),values.maxMagnitude,values.maxStars);
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseBuilder builder;
    // TODO: allow magnitude and weird
    unsigned char *catalogBuffer = builder.AddSubDatabase(kCatalogMagicValue,SerializeLengthCatalog(narrowedCatalog, false, true));
    SerializeCatalog(narrowedCatalog, false, true, catalogBuffer);

    GenerateDatabases(builder, narrowedCatalog, values);

    std::cerr << "Generated database with " << builder.BufferLength() << " bytes" << std::endl;

    if (isatty(fileno(stdout)) && values.path == "stdout") {
       std::cout << "Warning: output contains binary contents. Not printed to terminal." << std::endl;
    } else {
        PromptedOutputStream pos = PromptedOutputStream(values.path);
        pos.Stream().write((char *)builder.Buffer(), builder.BufferLength());
    }

}

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

}

int main(int argc, char **argv) {

    if (argc == 1) {
        std::cout << "Usage: ./lost database or ./lost pipeline" << std::endl << "Use --help flag on those commands for further help" << std::endl; 
        return 0;
    } 

    std::string command (argv[1]);
    if (command == "--database") {

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

        while (optind < argc && argc != 2) {
            if (optind < 2) optind++;

            if ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
                switch (option) {
                    case mag :
                        std::cout << "You raised the magmitude to " << optarg << std::endl;
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
                        // system("man documentation/database.man");
                        std::cout << documentation_database_man << std::endl;
                        return 0;
                        break;
                    default :
                        std::cout << "Illegal flag" << std::endl;
                        exit(1);
                }
            } 
        }

        lost::DatabaseBuild(databaseOptions);

    } else if (command == "--pipeline") {

        if (argc == 2) {
            std::cout << "Run pipeline --help for further help" << std::endl; 
            return 0;
        }

        enum PipelineEnum {png, focalLength, pixelSize, fov, centroidAlgo, centroidDummyStars, centroidMagFilter, database, idAlgo,
            gvTolerance, pyTolerance, falseStars, maxMismatchProb, attitudeAlgo, plot, generate, horizontalRes, verticalRes, refBrightnessMag,
            spreadStddev, noiseStddev, boresightRightAsc, boresightDec,boresightRoll, help, threshold, plotRawInput, plotInput, plotOutput, 
            printCentroids, compareCentroids, compareStars, printAttitude, compareAttitude};

        static struct option long_options[] =
        {
            {"png",             required_argument, 0, png},
            {"focal-length",    required_argument, 0, focalLength},
            {"pixel-size",      required_argument, 0, pixelSize},
            {"fov",             required_argument, 0, fov},
            {"centroid-algo", required_argument, 0, centroidAlgo},
            {"centroid-dummy-stars", required_argument, 0, centroidDummyStars},
            {"centroid-mag-filter",  required_argument, 0, centroidMagFilter},
            {"id-algo",        required_argument, 0, idAlgo},
            {"gv-tolerance",    required_argument, 0, gvTolerance},
            {"py-tolerance",    required_argument, 0, pyTolerance},
            {"false-stars",     required_argument, 0, falseStars},
            {"max-mismatch-prob",  required_argument, 0, maxMismatchProb},
            {"attitude-algo",    required_argument, 0, attitudeAlgo},
            {"generate",        optional_argument, 0, generate},
            {"horizontal-res",  required_argument, 0, horizontalRes},
            {"vertical-res",  required_argument, 0, verticalRes},
            {"ref-brightness-mag",  no_argument, 0, refBrightnessMag},
            {"spread-stddev",  required_argument, 0, spreadStddev},
            {"noise-stddev",  required_argument, 0, noiseStddev},
            {"boresight-right-asc",  required_argument, 0, boresightRightAsc},
            {"boresight-dec",  required_argument, 0, boresightDec},
            {"boresight-roll",  required_argument, 0, boresightRoll},
            {"threshold",       required_argument, 0, threshold},
            {"plot-raw-input", optional_argument, 0, plotRawInput},
            {"plot-input", optional_argument, 0, plotInput},
            {"plot-output", optional_argument, 0, plotOutput},
            {"print-centroids", optional_argument, 0, printCentroids},
            {"compare-centroids", optional_argument, 0, compareCentroids},
            {"compare-stars", optional_argument, 0, compareStars},
            {"print-attitude", optional_argument, 0, printAttitude},
            {"compareAttitude", optional_argument, 0, compareAttitude},
            {"help",            no_argument, 0, help},
            {0, 0, 0, 0}
        };

        lost::PipelineOptions pipelineOptions;
        int index;
        int option;
        
        while (optind < argc) {
            if (optind < 2) optind++;

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
                        std::cout << "You set the centroid mag filter to " << optarg << std::endl;
                        pipelineOptions.centroidMagFilter = atoi(optarg);
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
                    case generate :
                        if (optarg) {
                            pipelineOptions.generate = atoi(optarg);
                            std::cout << "Generating images! " << optarg << std::endl;
                        } else {
                            std::cout <<"If this hangs, make sure you are using --generate=[number] for setting a number" << std::endl;
                            std::cout << "Generating images! " << std::endl;
                        }
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
                    case threshold : 
                        std::cout << "You set the threshold to " << optarg << std::endl;
                        pipelineOptions.threshold = atof(optarg);
                        break;
                    case plotRawInput :
                        if (optarg) {
                            pipelineOptions.plotRawInput = optarg;
                            std::cout <<"If this hangs, make sure you are using --plot-raw-input=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.plotRawInput = "stdout";
                        }
                        break;
                    case plotInput :
                        if (optarg) {
                            pipelineOptions.plotInput = optarg;
                            std::cout << optarg << std::endl;
                            std::cout <<"If this hangs, make sure you are using --plot-input=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.plotInput = "stdout";
                        }
                        break;
                    case plotOutput :
                        if (optarg) {
                            pipelineOptions.plotOutput = optarg;
                            std::cout <<"If this hangs, make sure you are using --plot-output=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.plotOutput = "stdout";
                        }
                        break;
                    case printCentroids :
                        if (optarg) {
                            pipelineOptions.printCentroids = optarg;
                            std::cout <<"If this hangs, make sure you are using --printCentroids=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.printCentroids = "stdout";
                        }
                        break;
                    case compareCentroids :
                        if (optarg) {
                            pipelineOptions.compareCentroids = optarg;
                            std::cout <<"If this hangs, make sure you are using --compare-centroids=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.compareCentroids = "stdout";
                        }
                        break;
                    case compareStars :
                        if (optarg) {
                            pipelineOptions.compareStars = optarg;
                            std::cout <<"If this hangs, make sure you are using --compare-stars=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.compareStars = "stdout";
                        }
                        break;
                    case printAttitude :
                        if (optarg) {
                            pipelineOptions.printAttitude = optarg;
                            std::cout <<"If this hangs, make sure you are using --print-attitude=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.printAttitude = "stdout";
                        }
                        break;
                    case compareAttitude :
                        if (optarg) {
                            pipelineOptions.compareAttitude = optarg;
                            std::cout <<"If this hangs, make sure you are using --compare-attitude=[string] for setting a path" << std::endl;
                        } else {
                            pipelineOptions.compareAttitude = "stdout";
                        }
                        break;
                    case help : 
                        //system("man documentation/pipeline.man");
                        std::cout << documentation_pipeline_man << std::endl;
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