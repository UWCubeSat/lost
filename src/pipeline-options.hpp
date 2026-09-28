// this file uses the "X" pattern.

// Arguments to LOST_CLI_OPTION:
// 1. String used as the command line option.
// 2. Type of the option value.
// 3. Property name
// 4. Default value
// 5. Code to convert optarg into the value.
// 6. The default value if the option is specified with no argument, or kNoDefaultArgument

// To properly align these fields, I recommend using an editor plugin. In Vim, try `vim-lion`; in
// Emacs, try `evil-lion`. With your cursor inside any of the blocks, type `glip,` to aLign the
// Inside of the current Paragraph to comma.

#include <string>

#include "scalar.hpp"

// CAMERA
LOST_CLI_OPTION("png"          , std::string  , png         , "" , optarg       , kNoDefaultArgument)
LOST_CLI_OPTION("focal-length" , scalar      , focalLength , 0  , atof(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("pixel-size"   , scalar      , pixelSize   , -1 , atof(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("fov"          , scalar      , fov         , 20 , atof(optarg) , kNoDefaultArgument)

// PIPELINE STAGES
LOST_CLI_OPTION("centroid-algo"            , std::string, centroidAlgo                  , ""  , optarg                  , "cog")
LOST_CLI_OPTION("centroid-dummy-stars"     , int        , centroidDummyNumStars         , 5   , atoi(optarg)            , kNoDefaultArgument)
LOST_CLI_OPTION("centroid-mag-filter"      , scalar    , centroidMagFilter             , -1  , STR_TO_SCALAR(optarg)  , 5)
LOST_CLI_OPTION("centroid-filter-brightest", int        , centroidFilterBrightest       , -1  , atoi(optarg)            , 10)
LOST_CLI_OPTION("database"                 , std::string, databasePath                  , ""  , optarg                  , kNoDefaultArgument)
LOST_CLI_OPTION("star-id-algo"             , std::string, idAlgo                        , ""  , optarg                  , "pyramid")
LOST_CLI_OPTION("angular-tolerance"        , scalar    , angularTolerance              , .04 , STR_TO_SCALAR(optarg)  , kNoDefaultArgument)
LOST_CLI_OPTION("false-stars-estimate"     , int        , estimatedNumFalseStars        , 500 , atoi(optarg)            , kNoDefaultArgument)
LOST_CLI_OPTION("max-mismatch-probability" , scalar    , maxMismatchProb               , .001, STR_TO_SCALAR(optarg)  , kNoDefaultArgument)
LOST_CLI_OPTION("attitude-algo"            , std::string, attitudeAlgo                  , ""  , optarg                  , "dqm")

// OUTPUT COMPARISON
LOST_CLI_OPTION("centroid-compare-threshold", scalar    , centroidCompareThreshold, 2 , STR_TO_SCALAR(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("attitude-compare-threshold", scalar    , attitudeCompareThreshold, 1 , STR_TO_SCALAR(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("plot-raw-input"            , std::string, plotRawInput            , "", optarg                 , "-")
LOST_CLI_OPTION("plot-input"                , std::string, plotInput               , "", optarg                 , "-")
LOST_CLI_OPTION("plot-expected"             , std::string, plotExpected            , "", optarg                 , "-")
LOST_CLI_OPTION("plot-centroid-indices"     , std::string, plotCentroidIndices     , "", optarg                 , "-")
LOST_CLI_OPTION("plot-output"               , std::string, plotOutput              , "", optarg                 , "-")
LOST_CLI_OPTION("print-expected-centroids"  , std::string, printExpectedCentroids  , "", optarg                 , "-")
LOST_CLI_OPTION("print-input-centroids"     , std::string, printInputCentroids     , "", optarg                 , "-")
LOST_CLI_OPTION("print-actual-centroids"    , std::string, printActualCentroids    , "", optarg                 , "-")
LOST_CLI_OPTION("print-attitude"            , std::string, printAttitude           , "", optarg                 , "-")
LOST_CLI_OPTION("print-expected-attitude"   , std::string, printExpectedAttitude   , "", optarg                 , "-")
LOST_CLI_OPTION("print-speed"               , std::string, printSpeed              , "", optarg                 , "-")
LOST_CLI_OPTION("compare-centroids"         , std::string, compareCentroids        , "", optarg                 , "-")
LOST_CLI_OPTION("compare-star-ids"          , std::string, compareStarIds          , "", optarg                 , "-")
LOST_CLI_OPTION("compare-attitudes"         , std::string, compareAttitudes        , "", optarg                 , "-")

// IMAGE GENERATION
LOST_CLI_OPTION("generate"                    , int     , generate                  , 0     , atoi(optarg)    , 1)
LOST_CLI_OPTION("generate-x-resolution"       , int     , generateXRes              , 1024  , atoi(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-y-resolution"       , int     , generateYRes              , 1024  , atoi(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-centroids-only"     , bool    , generateCentroidsOnly     , false , atobool(optarg) , true)
LOST_CLI_OPTION("generate-zero-mag-photons"   , scalar , generateZeroMagPhotons    , 20000 , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-saturation-photons" , scalar , generateSaturationPhotons , 150   , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-spread-stddev"      , scalar , generateSpreadStdDev      , 1     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-shot-noise"         , bool    , generateShotNoise         , true  , atobool(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("generate-dark-current"       , scalar , generateDarkCurrent       , 0.1   , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-read-noise-stddev"  , scalar , generateReadNoiseStdDev   , .05   , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-ra"                 , scalar , generateRa                , 88    , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-de"                 , scalar , generateDe                , 7     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-roll"               , scalar , generateRoll              , 0     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-random-attitudes"   , bool    , generateRandomAttitudes   , false , atobool(optarg) , true)
LOST_CLI_OPTION("generate-blur-ra"            , scalar , generateBlurRa            , 0     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-blur-de"            , scalar , generateBlurDe            , 0     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-blur-roll"          , scalar , generateBlurRoll          , 0     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-exposure"           , scalar , generateExposure          , 0.2   , STR_TO_SCALAR(optarg)    , 0.1)
LOST_CLI_OPTION("generate-readout-time"       , scalar , generateReadoutTime       , 0     , STR_TO_SCALAR(optarg)    , 0.01)
LOST_CLI_OPTION("generate-oversampling"       , int     , generateOversampling      , 4     , atoi(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-false-stars"        , int     , generateNumFalseStars     , 0     , atoi(optarg)    , 50)
LOST_CLI_OPTION("generate-false-min-mag"      , scalar , generateFalseMinMag       , 8     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-false-max-mag"      , scalar , generateFalseMaxMag       , 1     , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-perturb-centroids"  , scalar , generatePerturbationStddev, 0     , STR_TO_SCALAR(optarg)    , 0.2)
LOST_CLI_OPTION("generate-cutoff-mag"         , scalar , generateCutoffMag         , 6.0   , STR_TO_SCALAR(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-seed"               , int     , generateSeed              , 394859, atoi(optarg)    , kNoDefaultArgument)
LOST_CLI_OPTION("generate-time-based-seed"    , bool    , timeSeed                  , false , atobool(optarg) , true)
