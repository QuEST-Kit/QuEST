// Distributed under MIT licence. See https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt for details

# include <math.h>

# ifndef QuEST_PRECISION
# define QuEST_PRECISION 

// set precision here if not set during compilation
# ifndef QuEST_PREC
# define QuEST_PREC 2
# endif

# if QuEST_PREC==1
    // SINGLE PRECISION
    # define REAL float
    # define MPI_QuEST_REAL MPI_FLOAT
    # define REAL_STRING_FORMAT "%.8f"
    # define REAL_QASM_FORMAT "%.8g"
    # define REAL_EPS 1e-5
    # define absReal(X) fabs(X) // fabsf(X) - better to return doubles where possible
# elif QuEST_PREC==4
    // QUAD PRECISION
    // 80-bit precision for most implementations
    // not compatible with most GPUs
    # define REAL long double
    # define MPI_QuEST_REAL MPI_LONG_DOUBLE
    # define REAL_STRING_FORMAT "%.17Lf"
    # define REAL_QASM_FORMAT "%.17Lg"
    # define REAL_EPS 1e-14
    # define absReal(X) fabsl(X)
# else
    // DOUBLE PRECISION
    # define REAL double
    # define MPI_QuEST_REAL MPI_DOUBLE
    # define REAL_STRING_FORMAT "%.14f"
    # define REAL_QASM_FORMAT "%.14g"
    # define REAL_EPS 1e-13
    # define absReal(X) fabs(X)
# endif


# endif
