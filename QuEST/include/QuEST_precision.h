// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Sets the QuEST constants which depend on the variable (during compilation)
 * state-vector precision. Using single, double or quad floating point precision
 * has consequences for operation precision, total memory requirements,
 * network bandwidth and ultimately runtime.
 *
 * @author Ania Brown
 * @author Tyson Jones (doc)
 */

# include <math.h>

# ifndef QUEST_PRECISION_H
# define QUEST_PRECISION_H


// set default double precision if not set during compilation
# ifndef QuEST_PREC
# define QuEST_PREC 2
# endif


/*
 * Single precision, which uses 4 bytes per amplitude component
 */
# if QuEST_PREC==1
    # define qreal float
    // \cond HIDDEN_SYMBOLS   
    # define MPI_QuEST_REAL MPI_FLOAT
    # define MPI_MAX_AMPS_IN_MSG (1LL<<29) // must be 2^int
    # define REAL_STRING_FORMAT "%.8f"
    # define REAL_QASM_FORMAT "%.8g"
    # define REAL_EPS 1e-5
    # define REAL_SPECIFIER "%f"
    # define absReal(X) fabs(X) // not fabsf(X) - better to return doubles where possible
    // \endcond
/*
 * Double precision, which uses 8 bytes per amplitude component
 */
# elif QuEST_PREC==2
    # define qreal double
    // \cond HIDDEN_SYMBOLS   
    # define MPI_QuEST_REAL MPI_DOUBLE
    # define MPI_MAX_AMPS_IN_MSG (1LL<<28) // must be 2^int
    # define REAL_STRING_FORMAT "%.14f"
    # define REAL_QASM_FORMAT "%.14g"
    # define REAL_EPS 1e-13
    # define REAL_SPECIFIER "%lf"
    # define absReal(X) fabs(X)
    // \endcond
/*
 * Quad precision, which uses 16 bytes per amplitude component.
 * This is not compatible with most GPUs.
 */
# elif QuEST_PREC==4
    # define qreal long double
    // \cond HIDDEN_SYMBOLS   
    # define MPI_QuEST_REAL MPI_LONG_DOUBLE
    # define MPI_MAX_AMPS_IN_MSG (1LL<<27) // must be 2^int
    # define REAL_STRING_FORMAT "%.17Lf"
    # define REAL_QASM_FORMAT "%.17Lg"
    # define REAL_EPS 1e-14
    # define REAL_SPECIFIER "%Lf"
    # define absReal(X) fabsl(X)
    // \endcond
# endif


// the maximum number of qubit registers which can be passed to functions like applyMultiArbitraryPhaseOverrides()
# define MAX_NUM_REGS_APPLY_ARBITRARY_PHASE 100


/** @def QuEST_PREC 
 *
 * Sets the precision of \ref qreal and \ref qcomp floating-point numbers, and 
 * hence the numerical precision of \p Qureg.
 *
 * #QuEST_PREC can be set as a preprocessor macro during compilation, or by 
 * editing its definition in QuEST_precision.h. \n
 * The possible values are:
 * #QuEST_PREC  | qreal         | sizeof(qreal)
 * ------------ |------         | -------------
 * 1            | `float`       | 4 bytes
 * 2            | `double`      | 8 bytes
 * 4            | `long double` | 16 bytes
 *
 * > Note that quad precision (#QuEST_PREC <b>= 4 </b>) is not compatible with most GPUs.
 *
 * @see 
 * - createQureg() and createDensityQureg() for the total memory costs of 
 *   creating registers under these precisions.
 *
 * @ingroup type
 * @author Ania Brown
 * @author Tyson Jones (doc)
 */

/** @def qreal
 *
 * A precision-agnostic floating point number, as determined by \ref QuEST_PREC.
 * Is a single, double or quad precision float when \ref QuEST_PREC is 1, 2 or 4 respectively.
 *
 * @ingroup type
 * @author Ania Brown
 * @author Tyson Jones (doc)
 */

# endif // QUEST_PRECISION_H
