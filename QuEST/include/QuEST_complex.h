// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Specifies a precision-agnostic type qcomp, which resolves to a complex<T> in C++ and a
 * complex T in C (that provided by complex.h), and which supports operator overloading
 * for easy complex algebra. This allows users to calculate with a natural complex type before 
 * passinfg instances to the QuEST API as a \ref Complex through \ref toComplex and \ref fromComplex
 *
 * Adapted from the header originally written by Randy Meyers and Dr. Thomas Plum, accessed at
 * http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/2003/0303/cuj0303meyers/index.htm
 * Original header doc:
 * Compatibility file for C99 and C++ complex.  This header can be included by either C99 or 
 * ANSI C++ programs to allow complex arithmetic to be written in a common subset. Note that C 
 * overloads for both the real and complex math functions are available after this header has been
 * included.
 *
 * @authors Randy Meyers and Dr. Thomas Plum
 * @author Tyson Jones
 */

#ifndef QUEST_COMPLEX_H
#define QUEST_COMPLEX_H


/*
 * creating precision-specific complex aliases
 */

// hide these from doxygen
// \cond HIDDEN_SYMBOLS   

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

// \endcond


/*
 * creating a single precision-agnostic type
 */

// this horrible hack is needed for doxygen doc
#define qcomp
#undef qcomp

#if QuEST_PREC==1
#define qcomp float_complex
#elif QuEST_PREC==2
#define qcomp double_complex
#elif QuEST_PREC==4
#define qcomp long_double_complex
#endif


/*
 * creating converters to/from QuEST's internal type
 */

#define toComplex(scalar) ((Complex) {.real = creal(scalar), .imag = cimag(scalar)})
#define fromComplex(comp) qcomp(comp.real, comp.imag)


/*
 * creating doc
 */

/** @def qcomp
 *
 * A precision-agnostic operator-overloaded complex number type.  
 * This is a complex analog of \ref qreal and is of single, double or quad
 * precision depending on the value of \ref QuEST_PREC.
 * It resolves to the native complex type provided by <complex.h> 
 * for both C99 and C++11, so can be used with operators.
 * It can be constructed with \p qcomp(real, imag).
 * 
 * For example, in C, 
 * \code 
 * qcomp x = 2 + 3i;
 * x -= 3.2*x; 
 * \endcode
 * and in C++,
 * \code 
 * qcomp x = qcomp(2, 3);
 * x -= 3*x; 
 * \endcode
 *
 * Assuming \p QuEST_PREC=4, qcomp will be 'complex long double' in C and 
 * 'complex<long double>' in C++.
 *
 * Can be converted to/from \ref Complex, the struct accepted by the QuEST
 * interface, using \ref toComplex and \ref fromComplex.
 *
 * @ingroup type
 * @authors Randy Meyers and Dr. Thomas Plum (created C & C++ agnosticism)
 * @author Tyson Jones (created precision agnosticism)
 */

/** @def toComplex(qcomp)
 *
 * Creates a Complex struct, which can be passed to the QuEST API, from a qcomp
 * @ingroup type
 * @author Tyson Jones
 */

/** @def fromComplex(Complex)
 *
 * Converts a Complex struct to a qcomp native type
 * @ingroup type
 * @author Tyson Jones
 */



#endif  // #ifndef QUEST_COMPLEX_H
