#pragma once

// Compile Time Float/Double Type
#ifdef LOST_DATABASE_DOUBLE
    typedef double decimal;
    #define STR_TO_DECIMAL(x) std::stod(x)
#else
    typedef float decimal;
    #define STR_TO_DECIMAL(x) std::stof(x)
#endif
