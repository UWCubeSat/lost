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
    Catalog narrowedCatalog = NarrowCatalog(CatalogRead(), (int)(values.maxMagnitude*100), values.maxStars);
    std::cerr << "Narrowed catalog has " << narrowedCatalog.size() << " stars." << std::endl;

    MultiDatabaseBuilder builder;
    // TODO: allow magnitude and weird
    unsigned char *catalogBuffer = builder.AddSubDatabase(kCatalogMagicValue, SerializeLengthCatalog(narrowedCatalog, false, true));
    SerializeCatalog(narrowedCatalog, false, true, catalogBuffer);

    GenerateDatabases(builder, narrowedCatalog, values);

    std::cerr << "Generated database with " << builder.BufferLength() << " bytes" << std::endl;

    if (isatty(fileno(stdout)) && values.path == "stdout") {
       std::cerr << "Warning: output contains binary contents. Not printed to terminal." << std::endl;
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

}

// https://stackoverflow.com/a/69177115
#define LOST_OPTIONAL_OPTARG()                                   \
    ((optarg == NULL && optind < argc && argv[optind][0] != '-') \
     ? (bool) (optarg = argv[optind++])                          \
     : (optarg != NULL))

int main(int argc, char **argv) {

    if (argc == 1) {
        std::cout << "Usage: ./lost database or ./lost pipeline" << std::endl << "Use --help flag on those commands for further help" << std::endl; 
        return 0;
    } 

    std::string command (argv[1]);
    optind = 2;

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

        while ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
            switch (option) {
            case mag :
                databaseOptions.maxMagnitude = atof(optarg);
                break;
            case stars :
                databaseOptions.maxStars = atoi(optarg);
                break;
            case kvector :
                databaseOptions.kVectorEnabled = true;
                break;
            case kvectorMinDistance :
                databaseOptions.kVectorMinDistance = atof(optarg);
                break;
            case kvectorMaxDistance :
                databaseOptions.kVectorMaxDistance = atof(optarg);
                break;
            case kvectorDistanceBins :
                databaseOptions.kVectorDistanceBins = atol(optarg);
                break;
            case output :
                databaseOptions.path = optarg;
                break;
            case help :
                std::cout << documentation_database_txt << std::endl;
                return 0;
                break;
            default :
                std::cout << "Illegal flag" << std::endl;
                exit(1);
            }
        }


        lost::DatabaseBuild(databaseOptions);

    } else if (command == "pipeline") {

        enum PipelineEnum {png, focalLength, pixelSize, fov, centroidAlgo, centroidDummyStars, centroidMagFilter, database, idAlgo,
            gvTolerance, pyTolerance, falseStars, maxMismatchProb, attitudeAlgo, plot, generate, horizontalRes, verticalRes, refBrightnessMag,
            spreadStddev, noiseStddev, boresightRightAsc, boresightDec,boresightRoll, help, centroidCompareThreshold, attitudeCompareThreshold, 
            plotRawInput, plotInput, plotOutput, printCentroids, compareCentroids, compareStars, printAttitude, compareAttitude};

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
            {"ref-brightness-mag",  required_argument, 0, refBrightnessMag},
            {"spread-stddev",  required_argument, 0, spreadStddev},
            {"noise-stddev",  required_argument, 0, noiseStddev},
            {"boresight-right-asc",  required_argument, 0, boresightRightAsc},
            {"boresight-dec",  required_argument, 0, boresightDec},
            {"boresight-roll",  required_argument, 0, boresightRoll},
            {"plot-raw-input", optional_argument, 0, plotRawInput},
            {"plot-input", optional_argument, 0, plotInput},
            {"plot-output", optional_argument, 0, plotOutput},
            {"print-centroids", optional_argument, 0, printCentroids},
            {"compare-centroids", optional_argument, 0, compareCentroids},
            {"compare-stars", optional_argument, 0, compareStars},
            {"print-attitude", optional_argument, 0, printAttitude},
            {"compare-attitude", optional_argument, 0, compareAttitude},
            {"centroid-compare-threshold",       required_argument, 0, centroidCompareThreshold},
            {"attitude-compare-threshold",       required_argument, 0, attitudeCompareThreshold},
            {"database", required_argument, 0, database},
            {"help",            no_argument, 0, help},
            {0, 0, 0, 0}
        };

        lost::PipelineOptions pipelineOptions;
        int index;
        int option;
        
        while ((option = getopt_long(argc, argv, "", long_options, &index)) != -1) {
            switch (option) {
            case png :
                pipelineOptions.png = optarg;
                break;
            case focalLength :
                pipelineOptions.focalLength = atof(optarg);
                break;
            case pixelSize :
                pipelineOptions.pixelSize = atof(optarg);
                break;
            case fov :
                pipelineOptions.fov = atof(optarg);
                break;
            case centroidAlgo :
                pipelineOptions.centroidAlgo = std::string (optarg);                     
                break;
            case centroidDummyStars : 
                pipelineOptions.dummyCentroidNumStars = atoi(optarg);
                break;
            case centroidMagFilter :
                pipelineOptions.centroidMagFilter = atoi(optarg);
                break;
            case idAlgo :
                pipelineOptions.idAlgo = std::string (optarg);                      
                break;
            case gvTolerance :
                pipelineOptions.gvTolerance = atof(optarg);
                break;
            case pyTolerance :
                pipelineOptions.pyTolerance = atof(optarg);
                break;
            case falseStars :
                pipelineOptions.pyFalseStars = atoi(optarg);
                break;
            case maxMismatchProb :
                pipelineOptions.pyMismatchProb = atof(optarg);
                break;
            case attitudeAlgo :
                pipelineOptions.attitudeAlgo = std::string (optarg);                    
                break; 
            case generate :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.generate = atoi(optarg);
                } else {
                    pipelineOptions.generate = 1;
                }
                break;
            case horizontalRes :
                pipelineOptions.horizontalRes = atoi(optarg);
                break;
            case verticalRes :
                pipelineOptions.verticalRes = atoi(optarg);
                break;
            case refBrightnessMag :
                pipelineOptions.referenceBrightness = atoi(optarg);                        
                break;
            case spreadStddev :
                pipelineOptions.brightnessDeviation = atof(optarg);
                break;
            case noiseStddev :
                pipelineOptions.noiseDeviation = atof(optarg);
                break;
            case boresightRightAsc :
                pipelineOptions.ra = atof(optarg);
                break;
            case boresightDec :
                pipelineOptions.dec = atof(optarg);
                break;
            case boresightRoll :
                pipelineOptions.roll = atof(optarg);
                break;
            case attitudeCompareThreshold : 
                pipelineOptions.attitudeCompareThreshold = atof(optarg);
                break;
            case centroidCompareThreshold:
                pipelineOptions.centroidCompareThreshold = atof(optarg);
                break;
            case plotRawInput :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.plotRawInput = optarg;
                } else {
                    pipelineOptions.plotRawInput = "stdout";
                }
                break;
            case plotInput :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.plotInput = optarg;
                } else {
                    pipelineOptions.plotInput = "stdout";
                }
                break;
            case plotOutput :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.plotOutput = optarg;
                } else {
                    pipelineOptions.plotOutput = "stdout";
                }
                break;
            case printCentroids :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.printCentroids = optarg;
                } else {
                    pipelineOptions.printCentroids = "stdout";
                }
                break;
            case compareCentroids :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.compareCentroids = optarg;
                } else {
                    pipelineOptions.compareCentroids = "stdout";
                }
                break;
            case compareStars :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.compareStars = optarg;
                } else {
                    pipelineOptions.compareStars = "stdout";
                }
                break;
            case printAttitude :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.printAttitude = optarg;
                } else {
                    pipelineOptions.printAttitude = "stdout";
                }
                break;
            case compareAttitude :
                if (LOST_OPTIONAL_OPTARG()) {
                    pipelineOptions.compareAttitude = optarg;
                } else {
                    pipelineOptions.compareAttitude = "stdout";
                }
                break;
            case database: 
                pipelineOptions.database = optarg;
                break;
            case help : 
                std::cout << documentation_pipeline_txt << std::endl;
                return 0;
                break;
            default :
                std::cout << "Illegal flag" << std::endl;
                exit(1);
            }
        }


        lost::PipelineRun(pipelineOptions);

    } else {
        std::cout << "Usage: ./lost database or ./lost pipeline" << std::endl << "Use --help flag on those commands for further help" << std::endl; 
    }
    return 0;
}
