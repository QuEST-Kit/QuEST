/** @file
 * User-overridable numerical precision of
 * both the QuEST API and backends
 * 
 * @author Tyson Jones
 * @author Milos Prokop (patched trig overloads in v3)
 * 
 * @defgroup precision Precision
 * @ingroup api
 * @brief Macros for controlling QuEST's numerical precision.
 * @{
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
// Still, we use a #define, rather than a typedef, so that the value can be compile-time overridden.

/// @neverdoced
#define INDEX_TYPE long long int

// spoofing above macro as const to doc
#if 0

    /// @notyetdoced
    /// @macrodoc
    typedef long long int INDEX_TYPE;

#endif



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

/// @neverdoced
#define PAULI_MASK_TYPE long long unsigned int

// spoofing above macro as typedef to doc
#if 0

    /// @notyetdoced
    /// @macrodoc
    typedef long long unsigned int PAULI_MASK_TYPE;

#endif



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

// spoofing above macros as typedefs and consts to doc
#if 0

    /// @notyetdoced
    /// @macrodoc
    const int FLOAT_PRECISION = 2;

    /// @notyetdoced
    /// @macrodoc
    typedef double int FLOAT_TYPE;

#endif



/*
 * CHECK PRECISION TYPES ARE COMPATIBLE WITH DEPLOYMENT
 */

#if COMPILE_CUDA && (FLOAT_PRECISION == 4)
    #error "A quad floating-point precision (FLOAT_PRECISION=4, i.e. long double) is not supported by GPU deployment"
#endif



/*
 * RE-CONFIGURABLE DEFAULT VALIDATION PRECISION
 *
 * which is compile-time overridable by pre-defining DEFAULT_VALIDATION_EPSILON (e.g. 
 * in user code before importing QuEST, or passed as a preprocessor constant by the
 * compiler using argument -D), and runtime overridable using setValidationEpsilon()
 */

#ifndef DEFAULT_VALIDATION_EPSILON

    #if FLOAT_PRECISION == 1
        #define DEFAULT_VALIDATION_EPSILON 1E-5

    #elif FLOAT_PRECISION == 2
        #define DEFAULT_VALIDATION_EPSILON 1E-12

    #elif FLOAT_PRECISION == 4
        #define DEFAULT_VALIDATION_EPSILON 1E-13

    #endif

#endif

// spoofing above macros as typedefs and consts to doc
#if 0

    /// @notyetdoced
    /// @macrodoc
    const qreal DEFAULT_VALIDATION_EPSILON = 1E-12;

#endif



/*
 * PRECISION-AGNOSTIC CONVENIENCE MACROS
 */

#if FLOAT_PRECISION == 1
    #define QREAL_FORMAT_SPECIFIER "%.8g"

#elif FLOAT_PRECISION == 2
    #define QREAL_FORMAT_SPECIFIER "%.14g"

#elif FLOAT_PRECISION == 4
    #define QREAL_FORMAT_SPECIFIER "%.17Lg"
    
#endif

// spoofing above macros as typedefs and consts to doc
#if 0

    /// @notyetdoced
    /// @macrodoc
    const char* QREAL_FORMAT_SPECIFIER = "%.14g";

#endif



#endif // PRECISION_H

/** @} */ // (end file-wide doxygen defgroup)
