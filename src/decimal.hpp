#ifndef DECIMAL_H
#define DECIMAL_H

#include <cmath>
#include <string>

#ifdef LOST_FLOAT_MODE
    typedef float decimal;
    #define STR_TO_DECIMAL(x) std::stof(x)
#else
    typedef double decimal;
    #define STR_TO_DECIMAL(x) std::stod(x)
#endif

// This should only be used sparingly.
// It's better to verbosely typecast sometimes. Only use these to prevent promotions.
// The reason why this isn't used everywhere instead of the wrapped macros is
// because the code becomes hard to read when there are multiple layers of typecasting.
// With this method, we might have more preprocessing to do BUT the code remains readable
// as the methods remain relatively the same.
#define DECIMAL(x) (decimal) x

// Math Constants wrapped with Decimal typecast
#define DECIMAL_M_E             (decimal) M_E           /* e */
#define DECIMAL_M_LOG2E         (decimal) M_LOG2E       /* log_2 e */
#define DECIMAL_M_LOG10E        (decimal) M_LOG10E      /* log_10 e */
#define DECIMAL_M_LN2           (decimal) M_LN2         /* log_e 2 */
#define DECIMAL_M_LN10          (decimal) M_LN10        /* log_e 10 */
#define DECIMAL_M_PI            (decimal) M_PI          /* pi */
#define DECIMAL_M_PI_2          (decimal) M_PI_2        /* pi/2 */
#define DECIMAL_M_PI_4          (decimal) M_PI_4        /* pi/4 */
#define DECIMAL_M_1_PI          (decimal) M_1_PI        /* 1/pi */
#define DECIMAL_M_2_PI          (decimal) M_2_PI        /* 2/pi */
#define DECIMAL_M_2_SQRTPI      (decimal) M_2_SQRTPI    /* 2/sqrt(pi) */
#define DECIMAL_M_SQRT2         (decimal) M_SQRT2       /* sqrt(2) */
#define DECIMAL_M_SQRT1_2       (decimal) M_SQRT1_2     /* 1/sqrt(2) */

// Math Functions wrapped with Decimal typecast
#define DECIMAL_POW(base,power) (decimal) std::pow(base, power)
#define DECIMAL_SQRT(x)         (decimal) std::sqrt(x)
#define DECIMAL_LOG(x)          (decimal) std::log(x)
#define DECIMAL_EXP(x)          (decimal) std::exp(x)
#define DECIMAL_ERF(x)          (decimal) std::erf(x)

// Rouding methods wrapped with Decimal typecast
#define DECIMAL_ROUND(x)        (decimal) std::round(x)
#define DECIMAL_CEIL(x)         (decimal) std::ceil(x)
#define DECIMAL_FLOOR(x)        (decimal) std::floor(x)
#define DECIMAL_ABS(x)          (decimal) std::abs(x)

// Trig Methods wrapped with Decimal typecast
#define DECIMAL_SIN(x)          (decimal) std::sin(x)
#define DECIMAL_COS(x)          (decimal) std::cos(x)
#define DECIMAL_TAN(x)          (decimal) std::tan(x)
#define DECIMAL_ASIN(x)         (decimal) std::asin(x)
#define DECIMAL_ACOS(x)         (decimal) std::acos(x)
#define DECIMAL_ATAN(x)         (decimal) std::atan(x)
#define DECIMAL_ATAN2(x,y)      (decimal) std::atan2(x,y)

// Float methods wrapped with Decimal typecast
#define DECIMAL_FMA(x,y,z)      (decimal) std::fma(x,y,z)
#define DECIMAL_HYPOT(x,y)      (decimal) std::hypot(x,y)

#endif // decimal.hpp
