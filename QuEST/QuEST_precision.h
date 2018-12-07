// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details

/** @file
 * Sets the QuEST constants which depend on the variable (during compilation)
 * state-vector precision. Using single, double or quad floating point precision
 * has consequences for operation precision, total memory requirements,
 * network bandwidth and ultimately runtime.
 */

# include <math.h>

# ifndef QUEST_PRECISION_H
# define QUEST_PRECISION_H

// set default double precision if not set during compilation
# ifndef QuEST_PREC
# define QuEST_PREC 2
# endif

# if QuEST_PREC==1
    // SINGLE PRECISION
    // 4 bytes per amplitude component
    # define qreal float
    # define MPI_QuEST_REAL MPI_FLOAT
    # define MPI_MAX_AMPS_IN_MSG (1LL<<29) // must be 2^int
    # define REAL_STRING_FORMAT "%.8f"
    # define REAL_QASM_FORMAT "%.8g"
    # define REAL_EPS 1e-5
    # define absReal(X) fabs(X) // fabsf(X) - better to return doubles where possible
# elif QuEST_PREC==2
    // DOUBLE PRECISION
    // 8 bytes per amplitude component
    # define qreal double
    # define MPI_QuEST_REAL MPI_DOUBLE
    # define MPI_MAX_AMPS_IN_MSG (1LL<<28) // must be 2^int
    # define REAL_STRING_FORMAT "%.14f"
    # define REAL_QASM_FORMAT "%.14g"
    # define REAL_EPS 1e-13
    # define absReal(X) fabs(X)
# elif QuEST_PREC==4
    // QUAD PRECISION
    // 16 bytes / 80-bit precision for most implementations
    // not compatible with most GPUs
    # define qreal long double
    # define MPI_QuEST_REAL MPI_LONG_DOUBLE
    # define MPI_MAX_AMPS_IN_MSG (1LL<<27) // must be 2^int
    # define REAL_STRING_FORMAT "%.17Lf"
    # define REAL_QASM_FORMAT "%.17Lg"
    # define REAL_EPS 1e-14
    # define absReal(X) fabsl(X)
# endif


# endif // QUEST_PRECISION_H
