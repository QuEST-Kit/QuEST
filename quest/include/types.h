/** @file
 * Definitions of API and backend numerical types.
 */

#ifndef TYPES_H
#define TYPES_H

#include "quest/include/modes.h"
#include "quest/include/precision.h"



/*
 * REAL TYPE ALIASES
 */

typedef FLOAT_TYPE qreal;
typedef INDEX_TYPE qindex;



/*
 * COMPLEX TYPE ALIAS
 */

// when C++ parses this header during backend or C++ user-code compilation...
#ifdef __cplusplus

    // resolve qcomp as the standard C++ complex type
    #include <complex>
    typedef std::complex<FLOAT_TYPE> qcomp;

// when C parses this header, during compilation of C user code...
#else

    // pretend that the API's qcomp is the C complex type
    #include <complex.h>

    // which is either MSVC's custom C complex...
    #ifdef _MSC_VER

        #if (FLOAT_PRECISION == 1)
            typedef _Fcomplex qcomp;

        #elif (FLOAT_PRECISION == 2)
            typedef _Dcomplex qcomp;

        #elif (FLOAT_PRECISION == 4)
            typedef _Lcomplex qcomp;

        #endif

    // or that used by GNU & Clang
    #else
        typedef FLOAT_TYPE _Complex qcomp;

    #endif

#endif



/*
 * COMPLEX TYPE INSTANTIATION
 */

#ifdef __cplusplus

    // qcomp() C++ instantiation is already enabled

#else

    // TODO: below we define macros qcomp(re,im) to spoof
    // C++'s initializer, but this is an antipattern! It
    // causes code like qcomp(*)[] (a valid type) to be
    // attemptedly substituted with the macro at 
    // pre-processing which fails and throws an error.
    // This might plague user code, and also forces us to
    // use smelly type aliasing in matrices.h. Fix this!

    #ifdef _MSC_VER

        // enable qcomp() C instantiation
        #define qcomp(re,im) = (qcomp) {(re), (im)}

    #else

        // enable qcomp() C instantiation
        #define qcomp(re,im) ( (qreal) (re) + I*((qreal) (im)) )

    #endif

#endif



/*
 * COMPLEX TYPE LITERALS
 */

// MSVC C literals are literally impossible, and
// C11 literals are already defined in complex header

#ifdef __cplusplus

    // enables C++14 literals like "3.5i" which are specifically double precision
    using namespace std::complex_literals;

    // the above literals are irksome for C++ users changing precision; for example, qcomp x = 3.5i
    // will not compile in single precision, because complex<double> cannot be assigned to a
    // complex<float>. To spare C++ users from having to change their literals when recompiling,
    // we define variable-precision custom literals, enabling qcomp = 3.5_i. Note custom literal
    // args always have max-precision (hence here are "long double" and "long long int")

    static inline qcomp operator ""_i(long double y) {
        return qcomp(0, (qreal) y); // discarding y precision
    }

    static inline qcomp operator ""_i(unsigned long long int y) {
        return qcomp(0, (qreal) y); // discarding y precision
    }

#endif



/*
 * COMPLEX TYPE ARITHMETIC OVERLOADS
 */

// C11 arithmetic is already defined in complex header, and beautifully
// permits mixing of parameterised types and precisions

#ifdef __cplusplus

    // <complex> defines overloads between complex and same-precision floats,
    // so e.g. complex<double> + double work fine. However, for differing
    // types and precisions, we must ourselves define operators which cast
    // the non-complex type to a qreal (to match qcomp). We must do this for
    // all operators the user might wish to use (e.g. + - * /) upon all
    // natural literal precisions (typically double, e.g. "3.5"), plus all
    // assignment forms (x += 2), and when their argument order is reversed.
    // Furthermore, the user might do arithmetic on complex literals which are
    // not the same precision as qcomp, so compilation will fail depending
    // on the setting of PRECISION. To avoid this, we'll define overloads
    // between all type/precision permutations, always returning qcomp.
    // Via the unholy macros below, we create 312 overloads; since this will
    // no doubt break somebody's build/integration, users can disable this
    // attempt at precision-agnostic arithmetic via DEFINE_ARITHMETIC_OVERLOADS=0

    #ifndef DEFINE_ARITHMETIC_OVERLOADS
    #define DEFINE_ARITHMETIC_OVERLOADS 1
    #endif

    #if DEFINE_ARITHMETIC_OVERLOADS

    // define operator between complex<compprec> and realtype, in either order, casting to qcomp
    #define DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL(op, compprec, realtype) \
        static inline qcomp operator op (const std::complex<compprec>& a, const realtype& b) { \
            return \
                qcomp((qreal) real(a), (qreal) imag(a)) \
                op \
                qcomp((qreal) b, 0); \
        } \
        static inline qcomp operator op (const realtype& b, const std::complex<compprec>& a) { \
            return a op b; \
        } \
        static inline std::complex<compprec>& operator op##= (std::complex<compprec>& a, const realtype& b) { \
            a = a op b; \
            return a; \
        }

    // defines + - * / between complex<compprec> and the given real type
    #define DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(compprec, realtype) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL(+, compprec, realtype) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL(-, compprec, realtype) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL(*, compprec, realtype) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL(/, compprec, realtype)

    // defines + - * / between complex<compprec> and all (relevant) integer types
    #define DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_INTEGER(compprec) \
        DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(compprec, int) \
        DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(compprec, long int) \
        DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(compprec, long long int) \
        DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(compprec, unsigned) \
        DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(compprec, long unsigned) \
        DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(compprec, long long unsigned)

    // define arithmetic between any-precision-decimal complex and any-precision integer
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_INTEGER(float)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_INTEGER(double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_INTEGER(long double)

    // define arithmetic between any-precision decimal complex and any-differing-precision decimal real
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(float,       double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(float,       long double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(double,      long double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(double,      float)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(long double, float)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL(long double, double)

    // defines arithmetic between differing-precision decimal complex, casting to qcomp, except assignment
    #define DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(op, precA, precB) \
        static inline qcomp operator op (const std::complex<precA>& a, const std::complex<precB>& b) { \
            return \
                qcomp((qreal) real(a), (qreal) imag(a)) \
                op \
                qcomp((qreal) real(b), (qreal) imag(b)); \
        } \
        static inline qcomp operator op (const std::complex<precB>& b, const std::complex<precA>& a) { \
            return a op b; \
        }

    // defines + - * / between the given differing-precision complex, casting to qcomp, except assignment
    #define DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(+, precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(-, precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(*, precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(/, precA, precB)

    // define arithmetic between all differing-precision decimal complex, casting to qcomp, except assignment
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(float,  double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(float,  long double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(double, long double)

    // keep these unholy macros away from the users
    #undef DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL
    #undef DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX
    #undef DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL
    #undef DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_INTEGER
    #undef DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX

    #endif // DEFINE_ARITHMETIC_OVERLOADS

#endif



/*
 * CONVENIENCE FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif

    void reportQcomp(qcomp num);

// end de-mangler
#ifdef __cplusplus
}
#endif



#endif // TYPES_H