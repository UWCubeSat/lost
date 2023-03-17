// see pipeline-options.hpp for more information

#include <string>

LOST_CLI_OPTION("min-mag"              , float      , minMag                , 100   , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("max-stars"            , int        , maxStars              , 10000 , atoi(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("min-separation"       , float      , minSeparation         , 0.05  , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector"              , bool       , kvector               , false , atobool(optarg), true)
LOST_CLI_OPTION("kvector-min-distance" , float      , kvectorMinDistance    , 0.5   , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-max-distance" , float      , kvectorMaxDistance    , 15    , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-distance-bins", long       , kvectorNumDistanceBins, 10000 , atol(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("tetra"                , bool       , tetra                 , false , atobool(optarg), true)
LOST_CLI_OPTION("tetra-max-angle"      , float      , tetraMaxAngle         , 12    , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("output"               , std::string, outputPath            , "-"   , optarg         , kNoDefaultArgument)
