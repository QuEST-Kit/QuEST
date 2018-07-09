// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

# ifndef QuEST_INTERNAL
# define QuEST_INTERNAL

/** @file
 * Internal functions used to implement the public facing API in qubits.h. Do not call these functions
 * directly. 
 */

#include "QuEST_precision.h"

extern const char* errorCodes[];

#ifdef __cplusplus
extern "C" {
#endif

void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type);

/** Measure the probability
of a specified qubit being in the zero state.     

@param[in] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/
REAL findProbabilityOfZero(MultiQubit multiQubit, const int measureQubit);


/** Update the state vector to be consistent with measuring measureQubit=0.
Measure in Zero performs an irreversible change to the state vector: it updates the vector according
to the event that a zero have been measured on the qubit indicated by measureQubit (where 
his label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
then renormalising based on the total probability of measuring measureQubit=0. It then returns the 
probability of making this measurement. 

@param[in,out] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/
REAL measureInZero(MultiQubit multiQubit, const int measureQubit);

// Validation

int validateMatrixIsUnitary(ComplexMatrix2 u);

int validateAlphaBeta(Complex alpha, Complex beta);

int validateUnitVector(REAL ux, REAL uy, REAL uz);


// Error reporting

void exitWithError(int errorCode, const char *func);

void QuESTAssert(int isValid, int errorCode, const char *func);

unsigned long int hashString(char *str);

#ifdef __cplusplus
}
#endif

# endif
