// see pipeline-options.hpp for more information

LOST_CLI_OPTION("min-mag"              , float      , minMag                , 0    , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("max-stars"            , int        , maxStars              , 0    , atoi(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector"              , bool       , kvector               , false, atobool(optarg), true)
LOST_CLI_OPTION("kvector-min-distance" , float      , kvectorMinDistance    , 0.5  , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-max-distance" , float      , kvectorMaxDistance    , 15   , atof(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("kvector-distance-bins", long       , kvectorNumDistanceBins, 3000 , atol(optarg)   , kNoDefaultArgument)
LOST_CLI_OPTION("output"               , std::string, outputPath            , "-"  , optarg         , kNoDefaultArgument)
