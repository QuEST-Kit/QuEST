/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include "environment.h"
#include "qureg.h"
#include "structures.h"

#include <vector>



/*
 * ENVIRONMENT CREATION
 */

void validate_envNeverInit(bool isQuESTInit, bool isQuESTFinal, const char* caller);

void validate_newEnvDeploymentMode(int isDistrib, int isGpuAccel, int isMultithread, const char* caller);

void validate_newEnvDistributedBetweenPower2Nodes(int numNodes, const char* caller);

void validate_gpuIsCuQuantumCompatible(const char* caller);



/*
 * EXISTING ENVIRONMENT
 */

void validate_envInit(const char* caller);



/*
 * QUREG CREATION
 */

void validate_newQuregParams(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, int numCpuThreads, QuESTEnv env, const char* caller);

void validate_newOrExistingQuregAllocs(Qureg qureg, bool isNewQureg, const char* caller);

void validate_newQuregNotBothMultithreadedAndGpuAccel(int useGpu, int numThreads, const char* caller);



/*
 * EXISTING QUREG
 */

void validate_quregInit(Qureg qureg, const char* caller);



/*
 * MATRIX CREATION
 */

void validate_newMatrixNumQubits(int numQubits, const char* caller);

void validate_newOrExistingMatrixAllocs(CompMatrN matr, bool isNewMatr, const char* caller);



/*
 * EXISTING MATRIX
 */

void validate_matrixInit(CompMatr1 matr, const char* caller);
void validate_matrixInit(CompMatr2 matr, const char* caller);
void validate_matrixInit(CompMatrN matr, const char* caller);

void validate_numMatrixElems(int numQubits, std::vector<std::vector<qcomp>> elems, const char* caller);

void validate_matrixElemsDontContainUnsyncFlag(qcomp firstElem, const char* caller);

void validate_matrixIsSynced(CompMatrN matr, const char* caller);

void validate_matrixIsUnitary(CompMatr1 matr, const char* caller);
void validate_matrixIsUnitary(CompMatr2 matr, const char* caller);
void validate_matrixIsUnitary(CompMatrN matr, const char* caller);



/*
 * QUREG INITIALISATIONS
 */

void validate_initClassicalStateIndex(Qureg qureg, qindex ind, const char* caller);



/*
 * OPERATOR PARAMETERS
 */

void validate_target(Qureg qureg, int target, const char* caller);



#endif // VALIDATION_HPP