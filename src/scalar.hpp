#ifndef SCALAR_H
#define SCALAR_H

#include <cmath>
#include <string>

#ifdef LOST_FLOAT_MODE
    typedef float scalar;
    #define STR_TO_SCALAR(x) std::stof(x)
#else
    typedef double scalar;
    #define STR_TO_SCALAR(x) std::stod(x)
#endif

// This should only be used sparingly.
// It's better to verbosely typecast sometimes. Only use these to prevent promotions.
// The reason why this isn't used everywhere instead of the wrapped macros is
// because the code becomes hard to read when there are multiple layers of typecasting.
// With this method, we might have more preprocessing to do BUT the code remains readable
// as the methods remain relatively the same.
#define SCALAR(x) ((scalar) x)

// Math Constants wrapped with Scalar typecast
#define SCALAR_M_E             ((scalar) M_E)           /* e */
#define SCALAR_M_LOG2E         ((scalar) M_LOG2E)       /* log_2 e */
#define SCALAR_M_LOG10E        ((scalar) M_LOG10E)      /* log_10 e */
#define SCALAR_M_LN2           ((scalar) M_LN2)         /* log_e 2 */
#define SCALAR_M_LN10          ((scalar) M_LN10)        /* log_e 10 */
#define SCALAR_M_PI            ((scalar) M_PI)          /* pi */
#define SCALAR_M_PI_2          ((scalar) M_PI_2)        /* pi/2 */
#define SCALAR_M_PI_4          ((scalar) M_PI_4)        /* pi/4 */
#define SCALAR_M_1_PI          ((scalar) M_1_PI)        /* 1/pi */
#define SCALAR_M_2_PI          ((scalar) M_2_PI)        /* 2/pi */
#define SCALAR_M_2_SQRTPI      ((scalar) M_2_SQRTPI)    /* 2/sqrt(pi) */
#define SCALAR_M_SQRT2         ((scalar) M_SQRT2)       /* sqrt(2) */
#define SCALAR_M_SQRT1_2       ((scalar) M_SQRT1_2)     /* 1/sqrt(2) */

// Math Functions wrapped with Scalar typecast
#define SCALAR_POW(base,power) ((scalar) std::pow(base, power))
#define SCALAR_SQRT(x)         ((scalar) std::sqrt(x))
#define SCALAR_LOG(x)          ((scalar) std::log(x))
#define SCALAR_EXP(x)          ((scalar) std::exp(x))
#define SCALAR_ERF(x)          ((scalar) std::erf(x))

// Rounding methods wrapped with Scalar typecast)
#define SCALAR_ROUND(x)        ((scalar) std::round(x))
#define SCALAR_CEIL(x)         ((scalar) std::ceil(x))
#define SCALAR_FLOOR(x)        ((scalar) std::floor(x))
#define SCALAR_ABS(x)          ((scalar) std::abs(x))

// Trig Methods wrapped with Scalar typecast)
#define SCALAR_SIN(x)          ((scalar) std::sin(x))
#define SCALAR_COS(x)          ((scalar) std::cos(x))
#define SCALAR_TAN(x)          ((scalar) std::tan(x))
#define SCALAR_ASIN(x)         ((scalar) std::asin(x))
#define SCALAR_ACOS(x)         ((scalar) std::acos(x))
#define SCALAR_ATAN(x)         ((scalar) std::atan(x))
#define SCALAR_ATAN2(x,y)      ((scalar) std::atan2(x,y))

// Float methods wrapped with Scalar typecast)
#define SCALAR_FMA(x,y,z)      ((scalar) std::fma(x,y,z))
#define SCALAR_HYPOT(x,y)      ((scalar) std::hypot(x,y))

#endif // scalar.hpp
