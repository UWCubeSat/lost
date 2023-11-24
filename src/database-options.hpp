// see pipeline-options.hpp for more information

#include <string>
#include "decimal.hpp"

LOST_CLI_OPTION("min-mag"                , decimal      , minMag                , 100   , STR_TO_DECIMAL(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("max-stars"              , int        , maxStars              , 10000 , atoi(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("min-separation"         , decimal      , minSeparation         , 0.08  , STR_TO_DECIMAL(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector"                , bool       , kvector               , false , atobool(optarg), true)
LOST_CLI_OPTION("kvector-min-distance"   , decimal      , kvectorMinDistance    , 0.5   , STR_TO_DECIMAL(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-max-distance"   , decimal      , kvectorMaxDistance    , 15    , STR_TO_DECIMAL(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-distance-bins"  , long       , kvectorNumDistanceBins, 10000 , atol(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("swap-integer-endianness", bool       , swapIntegerEndianness , false , atobool(optarg), true)
LOST_CLI_OPTION("swap-float-endianness"  , bool       , swapFloatEndianness   , false , atobool(optarg), true)
LOST_CLI_OPTION("output"                 , std::string, outputPath            , "-"   , optarg         , kNoDefaultArgument)
