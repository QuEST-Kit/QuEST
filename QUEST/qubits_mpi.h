double calcTotalProbability(int rank, long int numAmpsPerRank, int numRanks, double *stateVecReal, double *stateVecImag);

void rotateQubit(const long int numAmpsPerRank, const int numQubits, const int rotQubit,
                double aRe, double aIm, double bRe,  double beIm,
                double *restrict stateVecReal, double *restrict stateVecImag,
                double *restrict stateVecRealPair, double *restrict stateVecImagPair, int rank);

double findProbabilityOfZero(int rank, const long int numAmpsPerRank, const int numQubits,
                const int measureQubit,
                double *restrict stateVecReal,
                double *restrict stateVecImag);
