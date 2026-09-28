// see pipeline-options.hpp for more information

#include <string>
#include "scalar.hpp"

LOST_CLI_OPTION("min-mag"                , scalar      , minMag                , 100   , STR_TO_SCALAR(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("max-stars"              , int        , maxStars                , 10000 , atoi(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("min-separation"         , scalar      , minSeparation         , 0.08  , STR_TO_SCALAR(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector"                , bool       , kvector                 , false , atobool(optarg), true)
LOST_CLI_OPTION("kvector-min-distance"   , scalar      , kvectorMinDistance    , 0.5   , STR_TO_SCALAR(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-max-distance"   , scalar      , kvectorMaxDistance    , 15    , STR_TO_SCALAR(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-distance-bins"  , long       , kvectorNumDistanceBins  , 10000 , atol(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("swap-integer-endianness", bool       , swapIntegerEndianness   , false , atobool(optarg), true)
LOST_CLI_OPTION("swap-scalar-endianness", bool       , swapScalarEndianness   , false , atobool(optarg), true)
LOST_CLI_OPTION("output"                 , std::string, outputPath              , "-"   , optarg         , kNoDefaultArgument)
