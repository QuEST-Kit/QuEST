void initStateVec (const int numQubits,
                   double *restrict stateVecReal,
                   double *restrict stateVecImag);

void rotateQubit (const int numQubits, const int rotQubit,
                  double alphaReal, double alphaImag,
                  double betaReal,  double betaImag,
                  double *restrict stateVecReal, double *restrict stateVecImag);

double findProbabilityOfZero (const int numQubits,
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
