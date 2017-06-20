# ifndef PRECISION
# define PRECISION

// *** EDIT PRECISION HERE
// OPTIONS: 1, 2, 4
# define P 2

# if P==1
	// SINGLE PRECISION
	# define REAL float
	# define MPI_QuEST_REAL MPI_FLOAT
	# define REAL_STRING_FORMAT "%.8f"
# elif P==4
	// QUAD PRECISION
	// 80-bit precision for most implementations
	# define REAL long double
	# define MPI_QuEST_REAL MPI_LONG_DOUBLE
	# define REAL_STRING_FORMAT "%.17Lf"
# else
	// DOUBLE PRECISION
	# define REAL double
	# define MPI_QuEST_REAL MPI_DOUBLE
	# define REAL_STRING_FORMAT "%.14f"
# endif


# endif
