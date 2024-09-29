/** @file
 * User-overridable numerical precisions of the API and backend
 */

#ifndef PRECISION_H
#define PRECISION_H

#include "quest/include/modes.h"



/*
 * STATE-INDEXING TYPE
 */

// can be (for example) int, long, long long, unsigned, long unsigned, long long unsigned.
// We make sure that the backend never relies upon being able to represent negative 
// indices (e.g. as flags) since that would require a strictly signed type. Note this precision
// determines qindex which is user-facing so using unsigned types opens the users to the
// risks of underflowing. Since we never store large collections of this type, there is little 
// benefit in shrinking the type size and facing the associated precision risks. Similarly,
// there is little benefit in making it larger since a 'long long int' can represent 62 qubits,
// which is already well beyond simulability, requiring 64 EiB total at double precision.
#define INDEX_TYPE long long int



/*
 * PAULI STRING INDEXING TYPE
 */

// should never be changed; it is unsigned due to its use in extensive bitwise processing
// (no overflow risks since the API does not use this type), and its precision constrains
// the number of Paulis which can be specified in a PauliStr. Specifically, PauliStr stores
// two PAULI_MASK_TYPE instances, each of which are interpreted as half the digits of a 
// base-4 numeral encoding the Pauli string. A single 64-bit 'long long unsigned' can ergo
// specify only 32 qubits, whereas two can specify more qubits (64) than we can simulate.
// This type is defined purely to avoid littering the source with explicit typing.
#define PAULI_MASK_TYPE long long unsigned int



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

#if COMPILE_CUDA && (FLOAT_PRECISION == 4)
    #error "A quad floating-point precision (FLOAT_PRECISION=4, i.e. long double) is not supported by GPU deployment"
#endif

// Windows MSVC OpenMP doesn't permit operator overloading of the qcomp type,
// as is necessary when performing multithreaded reductions of amplitudes.
// We could support MSVC by separately reducing the real and imaginary components,
// but Bill Gates would have to wrestle me into submission.
#if COMPILE_OPENMP && defined(_MSC_VER)
    #error "Cannot use OpenMP multi-threading on Windows"
#endif



/*
 * RE-CONFIGURABLE VALIDATION PRECISION
 */

#ifndef DEAULT_VALIDATION_EPSILON

    #if (FLOAT_PRECISION == 1)
        #define DEAULT_VALIDATION_EPSILON 1E-5

    #elif (FLOAT_PRECISION == 2)
        #define DEAULT_VALIDATION_EPSILON 1E-13

    #elif (FLOAT_PRECISION == 4)
        #define DEAULT_VALIDATION_EPSILON 1E-14

    #endif

#endif



#endif // PRECISION_H