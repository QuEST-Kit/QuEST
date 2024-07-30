/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include "quest/include/types.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include <vector>
#include <string>

using std::vector;
using std::string;



/*
 *  MATRIX ATTRIBUTE FLAGS
 */

const int validate_UNITARITY_UNKNOWN_FLAG = -1;



/*
 * VALIDATION TOGGLE
 */

void validate_enable();

void validate_disable();

bool validate_isEnabled();



/*
 * ENVIRONMENT CREATION
 */

void validate_envNeverInit(bool isQuESTInit, bool isQuESTFinal, const char* caller);

void validate_newEnvDeploymentMode(int isDistrib, int isGpuAccel, int isMultithread, const char* caller);

void validate_newEnvDistributedBetweenPower2Nodes(const char* caller);

void validate_gpuIsCuQuantumCompatible(const char* caller);



/*
 * EXISTING ENVIRONMENT
 */

void validate_envIsInit(const char* caller);



/*
 * DEBUG UTILITIES
 */

void validate_numReportedItems(qindex num, const char* caller);



/*
 * QUREG CREATION
 */

void validate_newQuregParams(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, int numCpuThreads, QuESTEnv env, const char* caller);

void validate_newQuregNotBothMultithreadedAndGpuAccel(int useGpu, int numThreads, const char* caller);

void validate_newQuregAllocs(Qureg qureg, const char* caller);



/*
 * EXISTING QUREG
 */

void validate_quregFields(Qureg qureg, const char* caller);



/*
 * MATRIX CREATION
 */

void validate_newCompMatrParams(int numQubits, const char* caller);
void validate_newDiagMatrParams(int numQubits, const char* caller);
void validate_newFullStateDiagMatrParams(int numQubits, int useDistrib, const char* caller);

void validate_newMatrixAllocs(CompMatr matr, qindex numBytes, const char* caller);
void validate_newMatrixAllocs(DiagMatr matr, qindex numBytes, const char* caller);
void validate_newMatrixAllocs(FullStateDiagMatr matr, qindex numBytes, const char* caller);



/*
 * EXISTING MATRIX
 */

void validate_matrixFields(CompMatr1 matr, const char* caller);
void validate_matrixFields(CompMatr2 matr, const char* caller);
void validate_matrixFields(CompMatr  matr, const char* caller);
void validate_matrixFields(DiagMatr1 matr, const char* caller);
void validate_matrixFields(DiagMatr2 matr, const char* caller);
void validate_matrixFields(DiagMatr  matr, const char* caller);
void validate_matrixFields(FullStateDiagMatr matr, const char* caller);

void validate_matrixNumNewElems(int numQubits, vector<vector<qcomp>> elems, const char* caller);
void validate_matrixNumNewElems(int numQubits, vector<qcomp> elems, const char* caller);

void validate_matrixNewElemsDontContainUnsyncFlag(qcomp firstElem, const char* caller);

void validate_fullStateDiagMatrNewElems(FullStateDiagMatr matr, qindex startInd, qindex numElems, const char* caller);

void validate_matrixIsSynced(CompMatr matr, const char* caller);
void validate_matrixIsSynced(DiagMatr matr, const char* caller);
void validate_matrixIsSynced(FullStateDiagMatr matr, const char* caller);

void validate_matrixIsUnitary(CompMatr1 matr, const char* caller);
void validate_matrixIsUnitary(CompMatr2 matr, const char* caller);
void validate_matrixIsUnitary(CompMatr  matr, const char* caller);
void validate_matrixIsUnitary(DiagMatr1 matr, const char* caller);
void validate_matrixIsUnitary(DiagMatr2 matr, const char* caller);
void validate_matrixIsUnitary(DiagMatr  matr, const char* caller);
void validate_matrixIsUnitary(FullStateDiagMatr matr, const char* caller);

void validate_matrixIsCompatibleWithQureg(FullStateDiagMatr matr, Qureg qureg, const char* caller);



/*
 * QUREG INITIALISATIONS
 */

void validate_initClassicalStateIndex(Qureg qureg, qindex ind, const char* caller);



/*
 * PAULI STRINGS AND SUMS CREATION
 */

void validate_newPauliStrParams(const char* paulis, int* indices, int numPaulis, int maxNumPaulis, const char* caller);

void validate_newPauliStrNumChars(int numPaulis, int numIndices, const char* caller); // called by C++ only

void validate_newPauliStrNumPaulis(int numPaulis, int maxNumPaulis, const char* caller); // called by C++ only

void validate_newPauliStrSumParams(qindex numTerms, const char* caller);

void validate_newPauliStrSumMatchingListLens(qindex numStrs, qindex numCoeffs, const char* caller); // called by C++ only

void validate_newPauliStrSumAllocs(PauliStrSum sum, qindex numBytesStrings, qindex numBytesCoeffs, const char* caller);


/*
 * PAULI STRING SUM PARSING
 */

void validate_parsedPauliStrSumLineIsInterpretable(bool isInterpretable, string line, qindex lineIndex, const char* caller);

void validate_parsedPauliStrSumCoeffIsValid(bool isCoeffValid, string line, qindex lineIndex, const char* caller);

void validate_parsedPauliStrSumLineHasConsistentNumPaulis(int numPaulis, int numLinePaulis, string line, qindex lineIndex, const char* caller);

void validate_parsedStringIsNotEmpty(bool stringIsNotEmpty, const char* caller);



/*
 * EXISTING PAULI STRING SUMS
 */

void validate_pauliStrSumFields(PauliStrSum sum, const char* caller);



/*
 * OPERATOR PARAMETERS
 */

void validate_target(Qureg qureg, int target, const char* caller);



/*
 * FILE IO
 */

void validate_canReadFile(string fn, const char* caller);



#endif // VALIDATION_HPP