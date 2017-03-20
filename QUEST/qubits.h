void initStateVec (const int numQubits,
                   double *restrict stateVecReal,
                   double *restrict stateVecImag);

void rotateQubitLocal (const long long int numTasks, const int numQubits, const int rotQubit,
                double alphaReal, double alphaImag,
                double betaReal,  double betaImag,
                double *restrict stateVecReal, double *restrict stateVecImag);

int chunkIsUpper(int chunkId, int chunkSize, int rotQubit);

void getAlphaBeta(int chunkIsUpper, double *rot1Real, double *rot1Imag, double *rot2Real, double *rot2Imag,
                        double aReal, double aImag, double bReal, double bImag);

int getChunkPairId(int chunkIsUpper, int chunkId, int chunkSize, int rotQubit);

int halfMatrixBlockFitsInChunk(int chunkSize, int rotQubit);

void rotateQubitDistributed (const long long int numTasks, const int numQubits, const int rotQubit,
                double rot1Real, double rot1Imag,
                double rot2Real,  double rot2Imag,
                double *stateVecRealUp, double *stateVecImagUp,
                double *stateVecRealLo, double *stateVecImagLo,
                double *stateVecRealOut, double *stateVecImagOut);

int isChunkToSkipInFindPZero(int chunkId, int chunkSize, int measureQubit);

double findProbabilityOfZeroLocal (const long long int numTasks, const int numQubits,
                const int measureQubit,
                double *restrict stateVecReal,
                double *restrict stateVecImag);

double findProbabilityOfZeroDistributed (const long long int numTasks, const int numQubits,
                const int measureQubit,
                double *restrict stateVecReal,
                double *restrict stateVecImag);

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
