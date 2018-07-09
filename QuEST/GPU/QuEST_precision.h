// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details 

# ifndef QuEST_PRECISION
# define QuEST_PRECISION

// *** EDIT PRECISION HERE
// OPTIONS: 1, 2, 4
# define QuEST_PREC 2

# if QuEST_PREC==1
    // SINGLE PRECISION
    # define REAL float
    # define MPI_QuEST_REAL MPI_FLOAT
    # define REAL_STRING_FORMAT "%.8f"
    # define REAL_EPS 1e-5
# elif QuEST_PREC==4
    // QUAD PRECISION
    // 80-bit precision for most implementations
    # define REAL long double
    # define MPI_QuEST_REAL MPI_LONG_DOUBLE
    # define REAL_STRING_FORMAT "%.17Lf"
    # define REAL_EPS 1e-14
# else
    // DOUBLE PRECISION
    # define REAL double
    # define MPI_QuEST_REAL MPI_DOUBLE
    # define REAL_STRING_FORMAT "%.14f"
    # define REAL_EPS 1e-13
# endif


# endif
