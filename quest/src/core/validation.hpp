/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 * 
 * @author Tyson Jones
 */

#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include "quest/include/types.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"

#include <vector>
#include <string>

using std::vector;
using std::string;



/*
 *  PLACEHOLDER STRUCT FIELDS
 */

const int validate_STRUCT_PROPERTY_UNKNOWN_FLAG = -1;



/*
 * VALIDATION ERROR HANDLER
 */

void validateconfig_setErrorHandler(void (*callback)(const char* func, const char* msg));



/*
 * VALIDATION TOGGLE
 */

void validateconfig_enable();

void validateconfig_disable();

bool validateconfig_isEnabled();



/*
 * VALIDATION PRECISION
 */

void validateconfig_setEpsilon(qreal eps);

void validateconfig_setEpsilonToDefault();

qreal validateconfig_getEpsilon();



/*
 * ENVIRONMENT CREATION
 */

void validate_envNeverInit(bool isQuESTInit, bool isQuESTFinal, const char* caller);

void validate_newEnvDeploymentMode(int isDistrib, int isGpuAccel, int isMultithread, const char* caller);

void validate_newEnvDistributedBetweenPower2Nodes(const char* caller);

void validate_newEnvNodesEachHaveUniqueGpu(const char* caller);

void validate_gpuIsCuQuantumCompatible(const char* caller);



/*
 * EXISTING ENVIRONMENT
 */

void validate_envIsInit(const char* caller);



/*
 * DEBUG UTILITIES
 */

void validate_randomSeeds(unsigned* seeds, int numSeeds, const char* caller);

void validate_newEpsilonValue(qreal eps, const char* caller);

void validate_newMaxNumReportedScalars(qindex numRows, qindex numCols, const char* caller);

void validate_newMaxNumReportedSigFigs(int numSigFigs, const char* caller);

void validate_newNumReportedNewlines(int numNewlines, const char* caller);

void validate_numReportedNewlinesAboveZero(const char* caller);

void validate_numPauliChars(const char* paulis, const char* caller);

void validate_reportedPauliStrStyleFlag(int flag, const char* caller);



/*
 * QUREG CREATION
 */

void validate_newQuregParams(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, int numCpuThreads, QuESTEnv env, const char* caller);

void validate_newQuregAllocs(Qureg qureg, const char* caller);



/*
 * EXISTING QUREG
 */

void validate_quregFields(Qureg qureg, const char* caller);

void validate_quregIsStateVector(Qureg qureg, const char* caller);

void validate_quregIsDensityMatrix(Qureg qureg, const char* caller);



/*
 * MATRIX CREATION
 */

void validate_newCompMatrParams(int numQubits, const char* caller);
void validate_newDiagMatrParams(int numQubits, const char* caller);
void validate_newFullStateDiagMatrParams(int numQubits, int useDistrib, int useGpu, int useMultithread, const char* caller);

void validate_newMatrixAllocs(CompMatr matr, const char* caller);
void validate_newMatrixAllocs(DiagMatr matr, const char* caller);
void validate_newMatrixAllocs(FullStateDiagMatr matr, const char* caller);



/*
 * MATRIX INITIALISATION
 */

void validate_matrixNumNewElems(int numQubits, vector<vector<qcomp>> elems, const char* caller);
void validate_matrixNumNewElems(int numQubits, vector<qcomp> elems, const char* caller);

void validate_matrixNewElemsPtrNotNull(qcomp* elems, const char* caller);
void validate_matrixNewElemsPtrNotNull(qcomp** elems, qindex numRows, const char* caller);

void validate_fullStateDiagMatrNewElems(FullStateDiagMatr matr, qindex startInd, qindex numElems, const char* caller);

void validate_matrixNumQubitsMatchesParam(int numMatrQubits, int numSetterQubits, const char* caller);
void validate_declaredNumElemsMatchesVectorLength(qindex numElems, qindex vecLength, const char* caller);

void validate_multiVarFuncQubits(int numMatrQubits, int* numQubitsPerVar, int numVars, const char* caller);
void validate_funcVarSignedFlag(int areSigned, const char* caller);

void validate_matrixRowsAllSameSize(vector<vector<qcomp>> matrix, const char* caller);



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

void validate_unitaryExponentIsReal(qcomp exponent, const char* caller);

void validate_matrixIsHermitian(CompMatr1 matr, const char* caller);
void validate_matrixIsHermitian(CompMatr2 matr, const char* caller);
void validate_matrixIsHermitian(CompMatr  matr, const char* caller);
void validate_matrixIsHermitian(DiagMatr1 matr, const char* caller);
void validate_matrixIsHermitian(DiagMatr2 matr, const char* caller);
void validate_matrixIsHermitian(DiagMatr  matr, const char* caller);
void validate_matrixIsHermitian(FullStateDiagMatr matr, const char* caller);

void validate_matrixExpIsHermitian(DiagMatr          matr, qreal exponent, const char* caller);
void validate_matrixExpIsHermitian(FullStateDiagMatr matr, qreal exponent, const char* caller);

void validate_matrixExpIsNonDiverging(DiagMatr          matr, qcomp exponent, const char* caller);
void validate_matrixExpIsNonDiverging(FullStateDiagMatr matr, qcomp exponent, const char* caller);

void validate_matrixDimMatchesTargets(CompMatr1 matr, int numTargs, const char* caller);
void validate_matrixDimMatchesTargets(CompMatr2 matr, int numTargs, const char* caller);
void validate_matrixDimMatchesTargets(CompMatr  matr, int numTargs, const char* caller);
void validate_matrixDimMatchesTargets(DiagMatr1 matr, int numTargs, const char* caller);
void validate_matrixDimMatchesTargets(DiagMatr2 matr, int numTargs, const char* caller);
void validate_matrixDimMatchesTargets(DiagMatr  matr, int numTargs, const char* caller);

void validate_matrixAndQuregAreCompatible(FullStateDiagMatr matr, Qureg qureg, bool expecOnly, const char* caller);



/*
 * SUPEROPERATOR CREATION
 */

void validate_newSuperOpParams(int numQubits, const char* caller);

void validate_newSuperOpAllocs(SuperOp op, const char* caller);

void validate_newInlineSuperOpDimMatchesVectors(int numDeclaredQubits, vector<vector<qcomp>> matrix, const char* caller);



/*
 * SUPEROPERATOR INITIALISATION
 */

void validate_superOpNewMatrixDims(SuperOp op, vector<vector<qcomp>> matrix, const char* caller);

void validate_superOpFieldsMatchPassedParams(SuperOp op, int numQb, const char* caller);



/*
 * EXISTING SUPEROPERATOR
 */

void validate_superOpFields(SuperOp op, const char* caller);

void validate_superOpIsSynced(SuperOp op, const char* caller);

void validate_superOpDimMatchesTargs(SuperOp op, int numTargets, const char* caller);



/*
 * KRAUS MAP CREATION
 */

void validate_newKrausMapParams(int numQubits, int numMatrices, const char* caller);

void validate_newKrausMapAllocs(KrausMap map, const char* caller);

void validate_newInlineKrausMapDimMatchesVectors(int numQubits, int numOperators, vector<vector<vector<qcomp>>> matrices, const char* caller);



/*
 * KRAUS MAP INITIALISATION
 */

void validate_krausMapNewMatrixDims(KrausMap map, vector<vector<vector<qcomp>>> matrices, const char* caller);

void validate_krausMapFieldsMatchPassedParams(KrausMap map, int numQb, int numOps, const char* caller);



/*
 * EXISTING KRAUS MAP
 */

void validate_krausMapFields(KrausMap map, const char* caller);

void validate_krausMapIsSynced(KrausMap map, const char* caller);

void validate_krausMapIsCPTP(KrausMap map, const char* caller);

void validate_krausMapMatchesTargets(KrausMap map, int numTargets, const char* caller);



/*
 * PAULI STRINGS AND SUMS CREATION
 */

void validate_newPauliStrParams(const char* paulis, int* indices, int numPaulis, int maxNumPaulis, const char* caller);
void validate_newPauliStrParams(int*        paulis, int* indices, int numPaulis, int maxNumPaulis, const char* caller);

void validate_newPauliStrNumChars(int numPaulis, int numIndices, const char* caller); // used by C++ only API
void validate_newPauliStrNumPaulis(int numPaulis, int maxNumPaulis, const char* caller); // used by C++ only API

void validate_newPauliStrSumParams(qindex numTerms, const char* caller);

void validate_newPauliStrSumMatchingListLens(qindex numStrs, qindex numCoeffs, const char* caller); // used by C++ only API

void validate_newPauliStrSumAllocs(PauliStrSum sum, qindex numBytesStrings, qindex numBytesCoeffs, const char* caller);



/*
 * PAULI STRING SUM PARSING
 */

void validate_parsedPauliStrSumLineIsInterpretable(bool isInterpretable, string line, qindex lineIndex, const char* caller);

void validate_parsedPauliStrSumCoeffIsValid(bool isCoeffValid, string line, qindex lineIndex, const char* caller);

void validate_parsedPauliStrSumLineHasConsistentNumPaulis(int numPaulis, int numLinePaulis, string line, qindex lineIndex, const char* caller);

void validate_parsedStringIsNotEmpty(bool stringIsNotEmpty, const char* caller);



/*
 * EXISTING PAULI STRING
 */

void validate_pauliStrTargets(Qureg qureg, PauliStr str, const char* caller);

void validate_controlAndPauliStrTargets(Qureg qureg, int ctrl, PauliStr str, const char* caller);

void validate_controlsAndPauliStrTargets(Qureg qureg, int* ctrls, int numCtrls, PauliStr str, const char* caller);



/*
 * EXISTING PAULI STRING SUMS
 */

void validate_pauliStrSumFields(PauliStrSum sum, const char* caller);

void validate_pauliStrSumIsHermitian(PauliStrSum sum, const char* caller);

void validate_pauliStrSumTargets(PauliStrSum sum, Qureg qureg, const char* caller);

void validate_pauliStrSumCanInitMatrix(FullStateDiagMatr matr, PauliStrSum sum, const char* caller);



/*
 * BASIS STATE INDICES
 */

void validate_basisStateIndex(Qureg qureg, qindex ind, const char* caller);
void validate_basisStateRowCol(Qureg qureg, qindex row, qindex col, const char* caller);

void validate_basisStateIndices(Qureg qureg, qindex startInd, qindex numInds, const char* caller);
void validate_basisStateRowCols(Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols, const char* caller);

void validate_localAmpIndices(Qureg qureg, qindex localStartInd, qindex numInds, const char* caller);



/*
 * QUBIT INDICES
 */

void validate_target(Qureg qureg, int target, const char* caller);
void validate_targets(Qureg qureg, int* targets, int numTargets, const char* caller);

void validate_controls(Qureg qureg, int* ctrls, int numCtrls, const char* caller);
void validate_controlStates(int* states, int numCtrls, const char* caller);
void validate_controlsMatchStates(int numCtrls, int numStates, const char* caller); // C++ only

void validate_controlAndTarget(Qureg qureg, int ctrl, int targ, const char* caller);
void validate_controlsAndTarget(Qureg qureg, int* ctrls, int numCtrls, int targ, const char* caller);
void validate_controlAndTargets(Qureg qureg, int ctrl, int* targs, int numTargs, const char* caller);
void validate_controlsAndTargets(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, const char* caller);

void validate_twoTargets(Qureg qureg, int target1, int target2, const char* caller);
void validate_controlAndTwoTargets(Qureg qureg, int ctrl, int targ1, int targ2, const char* caller);
void validate_controlsAndTwoTargets(Qureg qureg, int* ctrls, int numCtrls, int targ1, int targ2, const char* caller);



/*
 * MEASUREMENT PARAMETERS
 */

void validate_measurementOutcomeIsValid(int outcome, const char* caller);

void validate_measurementOutcomesAreValid(int* outcomes, int numOutcomes, const char* caller);

void validate_measurementOutcomeProbNotZero(int outcome, qreal prob, const char* caller);

void validate_measurementOutcomesProbNotZero(int* outcomes, int numQubits, qreal prob, const char* caller);

void validate_measurementOutcomesFitInGpuMem(Qureg qureg, int numQubits, const char* caller);

void validate_measurementProbsAreNormalised(vector<qreal> probs, const char* caller);

void validate_measurementOutcomesMatchTargets(int numQubits, int numOutcomes, const char* caller); // C++ interface



/*
 * MISC GATE PARAMETERS
 */

void validate_rotationAxisNotZeroVector(qreal x, qreal y, qreal z, const char* caller);

void validate_mixedAmpsFitInNode(Qureg qureg, int numTargets, const char* caller);

void validate_trotterParams(Qureg qureg, int order, int reps, const char* caller);



/*
 * DECOHERENCE 
 */

void validate_probability(qreal prob, const char* caller);

void validate_oneQubitDepashingProb(qreal prob, const char* caller);
void validate_twoQubitDepashingProb(qreal prob, const char* caller);

void validate_oneQubitDepolarisingProb(qreal prob, const char* caller);
void validate_twoQubitDepolarisingProb(qreal prob, const char* caller);

void validate_oneQubitDampingProb(qreal prob, const char* caller);

void validate_oneQubitPauliChannelProbs(qreal pX, qreal pY, qreal pZ, const char* caller);



/*
 * QUREG COMBINATION
 */

void validate_quregCanBeWorkspace(Qureg quregA, Qureg quregB, const char* caller);

void validate_quregsCanBeMixed(Qureg quregOut, Qureg quregIn, const char* caller);

void validate_quregsCanBeSuperposed(Qureg qureg1, Qureg qureg2, Qureg qureg3, const char* caller);

void validate_quregCanBeInitialisedToPureState(Qureg qureg, Qureg pure, const char* caller);

void validate_quregsCanBeCloned(Qureg quregA, Qureg quregB, const char* caller);

void validate_quregsCanBeProducted(Qureg quregA, Qureg quregB, const char* caller);

void validate_throwErrorBecauseCalcFidOfDensMatrNotYetImplemented(const char* caller);

void validate_fidelityIsReal(qcomp fid, const char* caller);

void validate_buresDistanceInnerProdIsNormalised(qreal mag, const char* caller); // may modify mag

void validate_purifiedDistanceIsNormalised(qcomp fid, const char* caller); // may modify mag



/*
 * QUREG MODIFICATION
 */

void validate_quregRenormProbIsNotZero(qreal prob, const char* caller);

void validate_numInitRandomPureStates(qindex numPureStates,  const char* caller);



/*
 * EXPECTATION VALUES
 */

void validate_expecPauliStrValueIsReal(qcomp value, bool isDensMatr, const char* caller);

void validate_expecPauliStrSumValueIsReal(qcomp value, bool isDensMatr, const char* caller);

void validate_densMatrExpecDiagMatrValueIsReal(qcomp value, qcomp exponent, const char* caller);



/*
 * PARTIAL TRACE
 */


void validate_quregCanBeReduced(Qureg qureg, int numTraceQubits, const char* caller);

void validate_quregCanBeSetToReducedDensMatr(Qureg out, Qureg in, int numTraceQubits, const char* caller);



/*
 * FILE IO
 */

void validate_canReadFile(string fn, const char* caller);



/*
 * TEMPORARY ALLOCATIONS
 */

void validate_tempAllocSucceeded(bool succeeded, qindex numElems, qindex numBytesPerElem, const char* caller);



#endif // VALIDATION_HPP