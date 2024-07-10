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

    // requires C++14
    using namespace std::complex_literals;

#endif



/*
 * COMPLEX TYPE ARITHMETIC OVERLOADS
 */

// C11 arithmetic is already defined in complex header, and beautifully
// permits mixing of parameterised types and precisions

#ifdef __cplusplus

    // C++14 defines overloads between complex and same-precision floats,
    // so e.g. complex<double> + double work fine. However, for differing
    // types and precisions, we must ourselves define operators which cast
    // the non-complex type to a qreal (to match qcomp). We must do this for
    // all operators the user might wish to use (e.g. + - * /), and their
    // assignment forms (x += 2), and when their argument order is reversed.
    // To avoid 108 unique definitions, we employ these unholy macros.

    #define DEFINE_OPERATOR_FOR_TYPE(op, type) \
        inline qcomp operator op (const qcomp& a, const type& b) { \
            return a op qcomp((qreal) b,0); \
        } \
        inline qcomp operator op (const type& b, const qcomp& a) { \
            return a op b; \
        } \
        inline qcomp& operator op##= (qcomp& a, const type& b) { \
            a = a op b; \
            return a; \
        }

    #define DEFINE_COMPLEX_ARITHMETIC_FOR(type) \
        DEFINE_OPERATOR_FOR_TYPE(+, type) \
        DEFINE_OPERATOR_FOR_TYPE(-, type) \
        DEFINE_OPERATOR_FOR_TYPE(*, type) \
        DEFINE_OPERATOR_FOR_TYPE(/, type)

    DEFINE_COMPLEX_ARITHMETIC_FOR(int)
    DEFINE_COMPLEX_ARITHMETIC_FOR(long int)
    DEFINE_COMPLEX_ARITHMETIC_FOR(long long int)
    DEFINE_COMPLEX_ARITHMETIC_FOR(unsigned)
    DEFINE_COMPLEX_ARITHMETIC_FOR(long unsigned)
    DEFINE_COMPLEX_ARITHMETIC_FOR(long long unsigned)

    // we must not redefine the exising overloads between complex<T> and T
    #if (FLOAT_PRECISION != 1)
        DEFINE_COMPLEX_ARITHMETIC_FOR(float)
    #endif
    #if (FLOAT_PRECISION != 2)
        DEFINE_COMPLEX_ARITHMETIC_FOR(double)
    #endif
    #if (FLOAT_PRECISION != 4)
        DEFINE_COMPLEX_ARITHMETIC_FOR(long double)
    #endif

    // keep these unholy macros away from the users
    #undef DEFINE_OPERATOR_FOR_TYPE
    #undef DEFINE_COMPLEX_ARITHMETIC_FOR

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