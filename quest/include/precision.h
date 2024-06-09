/** @file
 * User-overridable numerical precisions of the API and backend
 */

#ifndef PRECISION_H
#define PRECISION_H

#include "quest/include/modes.h"



/*
 * STATE-INDEXING TYPE
 */

// can be (for example) int, long, long long, unsigned, long unsigned, long long unsigned
#define INDEX_TYPE long long unsigned



/*
 * RE-CONFIGURABLE FLOATING-POINT PRECISION
 */

// assume double precision as default
#ifndef FLOAT_PRECISION
    #define FLOAT_PRECISION 2
#endif

// validate precision is 1 (float), 2 (double) or 4 (long double)
#if ! (FLOAT_PRECISION == 1 || FLOAT_PRECISION == 2 || FLOAT_PRECISION == 4)
    #error "FLOAT_PRECISION must be 1 (float), 2 (double) or 4 (long double)"
#endif 

// infer floating-point type from precision
#if FLOAT_PRECISION == 1
    #define FLOAT_TYPE float
#elif FLOAT_PRECISION == 2
    #define FLOAT_TYPE double
#elif FLOAT_PRECISION == 4
    #define FLOAT_TYPE long double
#endif



/*
 * CHECK PRECISION TYPES ARE COMPATIBLE WITH DEPLOYMENT
 */

#if ENABLE_GPU_ACCELERATION && (FLOAT_PRECISION == 4)
    #error "A quad floating-point precision (FLOAT_PRECISION=4, i.e. long double) is not supported by GPU deployment"
#endif

// Windows MSVC OpenMP doesn't permit operator overloading of the qcomp type,
// as is necessary when performing multithreaded reductions of amplitudes.
// We could support MSVC by separately reducing the real and imaginary components,
// but Bill Gates would have to wrestle me into submission.
#if ENABLE_MULTITHREADING && defined(_MSC_VER)
    #error "Cannot use multi-threading on Windows"
#endif



/*
 * MACROS FOR PRINTING MULTI-WORD MACROS
 */

#define GET_STR_INTERNAL(x) #x
#define GET_STR(x) GET_STR_INTERNAL(x)



/*
 * RE-CONFIGURABLE VALIDATION PRECISION
 */

#ifndef VALIDATION_EPSILON

    #if (FLOAT_PRECISION == 1)
        #define VALIDATION_EPSILON 1E-5

    #elif (FLOAT_PRECISION == 2)
        #define VALIDATION_EPSILON 1E-13

    #elif (FLOAT_PRECISION == 4)
        #define VALIDATION_EPSILON 1E-14

    #endif

#endif



#endif // PRECISION_H