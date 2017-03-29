/** @file
Structs and specifications for functions that can be used from any environment (local, MPI)
*/

typedef struct ComplexArray
{
	double *real;
	double *imag;
} ComplexArray;

typedef struct Complex
{
	double real;
	double imag;
} Complex;

typedef struct MultiQubit
{
	ComplexArray stateVec, pairStateVec;
	int numQubits;
	long long int numAmps;
	int chunkId, numChunks;
} MultiQubit;

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

int chunkIsUpper(int chunkId, int chunkSize, int rotQubit);

void getRotAngle(int chunkIsUpper, Complex *rot1, Complex *rot2, Complex alpha, Complex beta);

int getChunkPairId(int chunkIsUpper, int chunkId, int chunkSize, int rotQubit);

int halfMatrixBlockFitsInChunk(int chunkSize, int rotQubit);

void rotateQubitDistributed (MultiQubit multiQubit, const int rotQubit,
		Complex rot1, Complex rot2,
                ComplexArray stateVecUp,
                ComplexArray stateVecLo,
                ComplexArray stateVecOut);

int isChunkToSkipInFindPZero(int chunkId, int chunkSize, int measureQubit);

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
