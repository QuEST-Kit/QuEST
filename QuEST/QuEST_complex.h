/*
 * Adapted from the header originally written by Randy Meyers and Dr. Thomas Plum, accessed at
 * http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/2003/0303/cuj0303meyers/index.htm
 */

// Compatibility file for C99 and C++ complex.  This header
// can be included by either C99 or ANSI C++ programs to
// allow complex arithmetic to be written in a common subset.
// Note that C overloads for both the real and complex math
// functions are available after this header has been
// included.

#ifndef QUEST_COMPLEX_H_
#define QUEST_COMPLEX_H_



/*
 * creating precision-specific complex aliases
 */

// C++ uses complex<T>
#ifdef __cplusplus

#include <cmath>
#include <complex>

using namespace std;

typedef complex<float> float_complex;
typedef complex<double> double_complex;
typedef complex<long double> long_double_complex;

// enable C-style funcs in C++
#define creal(x) real(x)
#define cimag(x) imag(x)
#define carg(x) arg(x)
#define cabs(x) abs(x)

#else

// C uses complex type
#include <tgmath.h> // includes <math.h> and <complex.h>

typedef float complex float_complex;
typedef double complex double_complex;
typedef long double complex long_double_complex;

#define float_complex(r,i) ((float)(r) + ((float)(i))*I)
#define double_complex(r,i) ((double)(r) + ((double)(i))*I)
#define long_double_complex(r,i) ((long double)(r) + ((long double)(i))*I)

#endif  // #ifdef __cplusplus



/*
 * creating a single precision-agnostic type
 */

# if QuEST_PREC==1
    # define qcomp float_complex
# elif QuEST_PREC==2
    # define qcomp double_complex
# elif QuEST_PREC==4
    # define qcomp long_double_complex
# endif



/*
 * creating converters to/from QuEST's internal type
 */

# define toComplex(scalar) ((Complex) {.real = creal(scalar), .imag = cimag(scalar)})
# define fromComplex(comp) qcomp(comp.real, comp.imag)



#endif  // #ifndef QUEST_COMPLEX_H_