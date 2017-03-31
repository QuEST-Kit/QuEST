/** @file
Structs and specifications for functions that can be used from any environment (local, MPI)
*/

/** Represents an array of complex numbers grouped into an array of real components and an array of coressponding complex components.
*/
typedef struct ComplexArray
{
	double *real; 
	double *imag;
} ComplexArray;

/** Represents one complex number.
*/
typedef struct Complex
{
	double real;
	double imag;
} Complex;

/** Represents a system of qubits.
Qubits are zero-based and the the first qubit is the rightmost
*/
typedef struct MultiQubit
{
	//! Probablilty amplitudes for the multi qubit state
	ComplexArray stateVec; 
	//! Temporary storage for a chunk of the state vector received from another process in the MPI version
	ComplexArray pairStateVec;
	//! Number of qubits in the state
	int numQubits;
	//! Number of probability amplitudes held in stateVec by this process
	//! In the non-MPI version, this is the total number of amplitudes
	long long int numAmps;
	//! The position of the chunk of the state vector held by this process in the full state vector
	int chunkId;
	//! Number of chunks the state vector is broken up into -- the number of MPI processes used
	int numChunks;
} MultiQubit;

/** Information about the environment the program is running in.
In practice, this holds info about MPI ranks and helps to hide MPI initialization code
*/
typedef struct QUESTEnv
{
	int rank;
	int numRanks;
} QUESTEnv;


void createMultiQubit(MultiQubit *multiQubit, int numQubits, QUESTEnv env);

void destroyMultiQubit(MultiQubit multiQubit, QUESTEnv env);

void reportState(MultiQubit multiQubit);

void initStateVec(MultiQubit *multiQubit);

void rotateQubitLocal (MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta);

void rotateQubitDistributed (MultiQubit multiQubit, const int rotQubit,
		Complex rot1, Complex rot2,
                ComplexArray stateVecUp,
                ComplexArray stateVecLo,
                ComplexArray stateVecOut);

double findProbabilityOfZeroLocal (MultiQubit multiQubit,
                const int measureQubit);

double findProbabilityOfZeroDistributed (MultiQubit multiQubit,
                const int measureQubit);

int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber);

void controlPhaseGate (const int numQubits, const int idQubit1, const int idQubit2,
                       double *restrict stateVecReal, double *restrict stateVecImag);

void quadCPhaseGate (const int numQubits, const int idQubit1, const int idQubit2, 
		const int idQubit3, const int idQubit4, double *restrict stateVecReal, 
		double *restrict stateVecImag);

double measureInZero (const int numQubits,
                              const int measureQubit,
                              double *restrict stateVecReal,
                              double *restrict stateVecImag);

double filterOut111 (const int numQubits, const int idQubit1, const int idQubit2, const int idQubit3,
                              double *restrict stateVecReal,
                              double *restrict stateVecImag);

double probOfFilterOut111 (const int numQubits, const int idQubit1, const int idQubit2, const int idQubit3,
                              double *restrict stateVecReal,
                              double *restrict stateVecImag);
