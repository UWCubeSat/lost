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

// CAMERA
LOST_CLI_OPTION("png"          , std::string, png         , "" , optarg       , kNoDefaultArgument)
LOST_CLI_OPTION("focal-length" , float      , focalLength , 0  , atof(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("pixel-size"   , float      , pixelSize   , -1 , atof(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("fov"          , float      , fov         , 20 , atof(optarg) , kNoDefaultArgument)

// PIPELINE STAGES
LOST_CLI_OPTION("centroid-algo"            , std::string, centroidAlgo                 , ""   , optarg       , "cog")
LOST_CLI_OPTION("centroid-dummy-stars"     , int        , centroidDummyNumStars        , 5    , atoi(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("centroid-mag-filter"      , float      , centroidMagFilter            , -1   , atof(optarg) , 5)
LOST_CLI_OPTION("database"                 , std::string, databasePath                 , ""   , optarg       , kNoDefaultArgument)
LOST_CLI_OPTION("star-id-algo"             , std::string, idAlgo                       , ""   , optarg       , "pyramid")
LOST_CLI_OPTION("angular-tolerance"        , float      , angularTolerance             , .04  , atof(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("false-stars-estimate"     , int        , estimatedNumFalseStars       , 500  , atoi(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("bayes-soft-threshold"     , float      , bayesSoftConfidenceThreshold , 0.999, atof(optarg), kNoDefaultArgument)
LOST_CLI_OPTION("bayes-hard-threshold"     , float      , bayesHardConfidenceThreshold , 0.99 , atof(optarg), kNoDefaultArgument)
LOST_CLI_OPTION("bayes-ignored-probability", float      , bayesIgnoredProbability      , 0.001, atof(optarg), kNoDefaultArgument)
LOST_CLI_OPTION("max-mismatch-probability" , float      , maxMismatchProb              , .001 , atof(optarg) , kNoDefaultArgument)
LOST_CLI_OPTION("attitude-algo"            , std::string, attitudeAlgo                 , ""   , optarg       , "dqm")

// OUTPUT COMPARISON
LOST_CLI_OPTION("centroid-compare-threshold", float      , centroidCompareThreshold, 1 , atof(optarg), kNoDefaultArgument)
LOST_CLI_OPTION("attitude-compare-threshold", float      , attitudeCompareThreshold, 1 , atof(optarg), kNoDefaultArgument)
LOST_CLI_OPTION("plot-raw-input"            , std::string, plotRawInput            , "", optarg      , "-")
LOST_CLI_OPTION("plot-input"                , std::string, plotInput               , "", optarg      , "-")
LOST_CLI_OPTION("plot-output"               , std::string, plotOutput              , "", optarg      , "-")
LOST_CLI_OPTION("print-centroids"           , std::string, printCentroids          , "", optarg      , "-")
LOST_CLI_OPTION("print-attitude"            , std::string, printAttitude           , "", optarg      , "-")
LOST_CLI_OPTION("compare-centroids"         , std::string, compareCentroids        , "", optarg      , "-")
LOST_CLI_OPTION("compare-star-ids"          , std::string, compareStarIds          , "", optarg      , "-")
LOST_CLI_OPTION("compare-attitudes"         , std::string, compareAttitudes        , "", optarg      , "-")

// IMAGE GENERATION
LOST_CLI_OPTION("generate"                     , int   , generate               , 0     , atoi(optarg)   , 1)
LOST_CLI_OPTION("generate-x-resolution"        , int   , generateXRes           , 1024  , atoi(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-y-resolution"        , int   , generateYRes           , 1024  , atoi(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-reference-brightness", float , generateRefBrightness  , 100   , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-spread-stddev"       , float , generateSpreadStdDev   , 1     , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-shot-noise"          , bool  , generateShotNoise      , true  , atobool(optarg), kNoDefaultArgument)
LOST_CLI_OPTION("generate-dark-current"        , float , generateDarkCurrent    , 0.1   , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-read-noise-stddev"   , float , generateReadNoiseStdDev, .05   , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-ra"                  , float , generateRa             , 88    , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-de"                  , float , generateDe             , 7     , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-roll"                , float , generateRoll           , 0     , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-blur-ra"             , float , generateBlurRa         , 4     , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-blur-de"             , float , generateBlurDe         , 1     , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-blur-roll"           , float , generateBlurRoll       , 20    , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-sensitivity"         , float , generateSensitivity    , 0.01  , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-exposure"            , float , generateExposure       , 0     , atof(optarg)   , 0.1)
LOST_CLI_OPTION("generate-readout-time"        , float , generateReadoutTime    , 0     , atof(optarg)   , 0.01)
LOST_CLI_OPTION("generate-oversampling"        , int   , generateOversampling   , 4     , atoi(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-false-stars"         , int   , generateNumFalseStars  , 0     , atoi(optarg)   , 50)
LOST_CLI_OPTION("generate-false-min-mag"       , float , generateFalseMinMag    , 8     , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-false-max-mag"       , float , generateFalseMaxMag    , 1     , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("generate-seed"                , int   , generateSeed           , 394859, atoi(optarg)   , kNoDefaultArgument)
