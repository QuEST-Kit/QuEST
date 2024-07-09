/** @file
 * Definitions of API and backend numerical types.
 */

#ifndef TYPES_H
#define TYPES_H

#include "modes.h"
#include "precision.h"



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
 * COMPLEX TYPE OVERLOADS
 */

#ifdef __cplusplus

    // enable C++ literals (requires C++14)
    using namespace std::complex_literals;

    // qcomp() C++ instantiation is already enabled

#else

    #ifdef _MSC_VER

        // MSVC C literals are literally impossible

        // enable qcomp() C instantiation
        #define qcomp(re,im) = (qcomp) {(re), (im)}

    #else

        // C literals are already enabled (requires C99)

        // enable qcomp() C instantiation
        #define qcomp(re,im) ( (qreal) (re) + I*((qreal) (im)) )

    #endif

#endif



#endif // TYPES_H