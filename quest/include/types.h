/** @file
 * Definitions of the variable-precision API 
 * scalar types (qindex, qreal and qcomp) and
 * convenience arithmetic operator overloads.
 * Note both the user code and backend parse
 * this file and may disagree on the underlying
 * qcomp type which is no issue; we ensure the
 * frontend-backend interface is agnostic to it.
 * 
 * @author Tyson Jones
 * @author Ali Rezaei (aided in design)
 *
 * @defgroup types Types
 * @ingroup api
 * @brief Macros for precision-agnostic real and complex arithmetic.
 * @{
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
 *
 * which is already elegantly possible in C++ using qcomp(re,im),
 * but which cannot be extended nicely to C using a macro qcomp(),
 * because it will break compilation of function signatures which
 * accept funciton pointers which return qcomp. For example:
 * myfunc( qcomp (*callback)(int,int) )
 * So instead, we define a new, ugly 'getQcomp()' helper. Aw! :(
 * We define it here in the header, inlining directly into user
 * code, to avoid C & C++ qcomp interoperability issues.
 */

/// @notyetdoced
static inline qcomp getQcomp(qreal re, qreal im) {

    #if defined(__cplusplus)
        return qcomp(re, im);

    #elif defined(_MSC_VER)
        return (qcomp) {re, im};

    #else
        return re + I*im;
        
    #endif
}



/*
 * COMPLEX TYPE LITERALS
 */

// MSVC C literals are literally impossible, and
// C11 literals are already defined in complex header

/// @cond EXCLUDE_FROM_DOXYGEN

#ifdef __cplusplus

    // enables C++14 literals like "3.5i" which are specifically double precision
    using namespace std::complex_literals;

    // the above literals are irksome for C++ users changing precision; for example, qcomp x = 3.5i
    // will not compile in single precision, because complex<double> cannot be assigned to a
    // complex<float>. To spare C++ users from having to change their literals when recompiling,
    // we define variable-precision custom literals, enabling qcomp = 3.5_i. Note custom literal
    // args always have max-precision (hence here are "long double" and "long long int")

    static inline qcomp operator ""_i(long double y) {
        return qcomp(0, static_cast<qreal>(y)); // discarding y precision
    }

    static inline qcomp operator ""_i(unsigned long long int y) {
        return qcomp(0, static_cast<qreal>(y)); // discarding y precision
    }

#endif

/// @endcond // EXCLUDE_FROM_DOXYGEN



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

    // spoofing above macro as const to doc
    #if 0

        /// @notyetdoced
        /// @macrodoc
        const int DEFINE_ARITHMETIC_OVERLOADS = 1;

    #endif


    #if DEFINE_ARITHMETIC_OVERLOADS

    /// @cond EXCLUDE_FROM_DOXYGEN

    // shortcuts for below overload definitions
    #define COMP_TO_QCOMP(a) \
        qcomp( \
            static_cast<qreal>(std::real(a)), \
            static_cast<qreal>(std::imag(a)))

    #define REAL_TO_QCOMP(b) \
        qcomp(static_cast<qreal>(b), 0)

    // define operator between complex<compprec> and realtype, in either order, casting to qcomp
    #define DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL(op, compprec, realtype) \
        static inline qcomp operator op (const std::complex<compprec>& a, const realtype& b) { \
            return \
                COMP_TO_QCOMP(a) \
                op \
                REAL_TO_QCOMP(b); \
        } \
        static inline qcomp operator op (const realtype& b, const std::complex<compprec>& a) { \
            return \
                REAL_TO_QCOMP(b) \
                op \
                COMP_TO_QCOMP(a); \
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
                COMP_TO_QCOMP(a) \
                op \
                COMP_TO_QCOMP(b); \
        }

    // defines + - * / between the given differing-precision complex, casting to qcomp, except assignment
    #define DEFINE_SINGLE_DIRECTION_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(+, precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(-, precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(*, precA, precB) \
        DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX(/, precA, precB)

    #define DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(precA, precB) \
        DEFINE_SINGLE_DIRECTION_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(precA, precB) \
        DEFINE_SINGLE_DIRECTION_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(precB, precA)

    // define arithmetic between all differing-precision decimal complex, casting to qcomp, except assignment
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(float,  double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(float,  long double)
    DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX(double, long double)

    // keep these unholy macros away from the users
    #undef COMP_TO_QCOMP
    #undef REAL_TO_QCOMP
    #undef DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_REAL
    #undef DEFINE_OPERATOR_BETWEEN_COMPLEX_AND_COMPLEX
    #undef DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_REAL
    #undef DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_INTEGER
    #undef DEFINE_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX
    #undef DEFINE_SINGLE_DIRECTION_ARITHMETIC_BETWEEN_COMPLEX_AND_COMPLEX

    /// @endcond // EXCLUDE_FROM_DOXYGEN

    #endif // DEFINE_ARITHMETIC_OVERLOADS

#endif



/*
 * REPORTERS
 */

// both C and C++ declare "reportScalar" accepting
// char array and a qcomp; but each language then
// separately defines convenience overloads for
// C++ strings and using C11 generics to overload

#ifdef __cplusplus

    #include <string>


    /// @notyetdoced
    /// @notyettested
    extern "C" void reportStr(const char* str);


    /// @notyetdoced
    /// @notyettested
    /// @cpponly
    /// @see reportStr()
    void reportStr(std::string str);


    /// @notyetdoced
    /// @notyettested
    extern "C" void reportScalar(const char* label, qcomp num);


    /// @notyetdoced
    /// @notyettested
    void reportScalar(const char* label, qreal num);


    /// @notyetdoced
    /// @notyettested
    /// @cpponly
    /// @see reportScalar()
    void reportScalar(std::string label, qcomp num);


    /// @notyetdoced
    /// @notyettested
    /// @cpponly
    void reportScalar(std::string label, qreal num);

#else

    /// @notyetdoced
    /// @notyettested
    void reportStr(const char* str);


    /// @notyetdoced
    /// @notyettested
    void reportScalar(const char* label, qcomp num);


    /// @private
    void _reportScalar_real(const char* label, qreal num);


    // no need to be doc'd since signatures identical to C++ above
    /// @neverdoced
    #define reportScalar(label, num) \
        _Generic((num), \
            qcomp   : reportScalar,       \
            qreal   : _reportScalar_real, \
            default : _reportScalar_real  \
        )(label, num)

#endif



#endif // TYPES_H

/** @} */ // (end file-wide doxygen defgroup)
