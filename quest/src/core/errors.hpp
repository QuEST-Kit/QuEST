/** @file
 * Defensively designed functions for checking internal preconditions, 
 * and raising internal errors. These primarily check that that
 * hardware accelerators are behaving as expected, and that runtime
 * deployment is consistent with the compiled deployment modes.
 * 
 * @author Tyson Jones
 */

#ifndef ERRORS_HPP
#define ERRORS_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include <string>

using std::string;



/*
 * VALIDATION ERRORS
 */

void error_validationMessageVarWasIllFormed(string msg, string illFormedVar);

void error_validationMessageVarNotSubstituted(string msg, string var);

void error_validationMessageContainedUnsubstitutedVars(string msg);

void error_validationEncounteredUnsupportedDistributedDenseMatrix();

void error_validationListUniquenessCheckExceededMaskSize();



/*
 * ENVIRONMENT ERRORS
 */

void error_allocOfQuESTEnvFailed();



/*
 * MEMORY ERRORS
 */

void error_memSizeQueriedButWouldOverflow();



/*
 * COMMUNICATION ERRORS
 */

void error_commNotInit();

void error_commAlreadyInit();

void error_commButEnvNotDistributed();

void error_commButQuregNotDistributed();

void error_commOutOfBounds();

void error_commWithSameRank();

void error_commGivenInconsistentNumSubArraysANodes();

void error_commNumMessagesExceedTagMax();

void assert_commBoundsAreValid(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps);

void assert_commPayloadIsPowerOf2(qindex numAmps);

void assert_commQuregIsDistributed(Qureg qureg);

void assert_commFullStateDiagMatrIsDistributed(FullStateDiagMatr matr);

void assert_pairRankIsDistinct(Qureg qureg, int pairRank);

void assert_bufferSendRecvDoesNotOverlap(qindex sendInd, qindex recvInd, qindex numAmps);

void assert_receiverCanFitSendersEntireState(Qureg receiver, Qureg sender);

void assert_receiverCanFitSendersEntireElems(Qureg receiver, FullStateDiagMatr sender);



/*
 * LOCALISER ERRORS
 */

void error_localiserNumCtrlStatesInconsistentWithNumCtrls();

void error_localiserGivenPauliTensorOrGadgetWithoutXOrY();

void error_localiserPassedStateVecToChannelComCheck();

void error_localiserGivenDistribMatrixAndLocalQureg();

void error_localiserFailedToAllocTempMemory();

void error_localiserGivenPauliStrWithoutXorY();

void error_localiserGivenNonUnityGlobalFactorToZTensor();

void assert_localiserSuccessfullyAllocatedTempMemory(qcomp* ptr, bool isGpu);

void assert_localiserGivenStateVec(Qureg qureg);

void assert_localiserGivenDensMatr(Qureg qureg);

void assert_localiserPartialTraceGivenCompatibleQuregs(Qureg inQureg, Qureg outQureg, int numTargs);

void error_calcFidStateVecDistribWhileDensMatrLocal();

void assert_localiserDistribQuregSpooferGivenValidQuregs(Qureg local, Qureg distrib);



/*
 * ACCELERATOR ERRORS
 */

void assert_highPauliStrSumMaskIsZero(PauliStrSum sum);

void assert_numQubitsMatchesQubitStatesAndTemplateParam(int numQubits, int numQubitStates, int templateParam, string label="qubit");

void assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(int numCtrls, int numCtrlStates, int templateParam);

void assert_numTargsMatchesTemplateParam(int numTargs, int templateParam);

void assert_exponentMatchesTemplateParam(qcomp exponent, bool hasPower);
void assert_exponentMatchesTemplateParam(qcomp exponent, bool hasPower, bool useRealPow);

void assert_mixedQuregIsDensityMatrix(Qureg qureg);

void assert_mixedQuregIsStatevector(Qureg qureg);

void assert_mixedQuregIsDistributed(Qureg qureg);

void assert_mixedQuregIsLocal(Qureg qureg);

void assert_mixedQuregsAreBothOrNeitherDistributed(Qureg a, Qureg b);

void error_mixQuregsAreLocalDensMatrAndDistribStatevec();

void assert_fullStateDiagMatrIsLocal(FullStateDiagMatr matr);

void assert_fullStateDiagMatrIsDistributed(FullStateDiagMatr matr);

void assert_acceleratorQuregIsDistributed(Qureg qureg);

void assert_quregAndFullStateDiagMatrAreBothOrNeitherDistrib(Qureg qureg, FullStateDiagMatr matr);

void assert_calcFidStateVecIsLocal(Qureg qureg);

void assert_calcFidTempGpuAllocSucceeded(qcomp* ptr);

void assert_calcExpecDiagTempGpuAllocSucceeded(qcomp* ptr);

void assert_innerProductedSameDimQuregsHaveSameGpuAccel(Qureg quregA, Qureg quregB);

void assert_partialTraceQuregsAreIdenticallyDeployed(Qureg inQureg, Qureg outQureg);



/*
 * BUFFER PACKING ERRORS
 */

void error_noCtrlsGivenToBufferPacker();

void assert_bufferPackerGivenIncreasingQubits(int qubit1, int qubit2, int qubit3);



/*
 * BACKEND PRECONDITION ERRORS
 */

void assert_quregAndFullStateDiagMatrHaveSameDistrib(Qureg qureg, FullStateDiagMatr matr);

void assert_quregDistribAndFullStateDiagMatrLocal(Qureg qureg, FullStateDiagMatr matr);

void assert_superposedQuregDimsAndDeploysMatch(Qureg facOut, Qureg in1, Qureg in2);



/*
 * CPU ERRORS
 */

void error_cpuThreadsQueriedButEnvNotMultithreaded();



/*
 * GPU ERRORS
 */

void error_gpuQueriedButGpuNotCompiled();

void error_gpuSyncedButGpuNotCompiled();

void error_gpuAllocButGpuNotCompiled();

void error_gpuDeallocButGpuNotCompiled();

void error_gpuCopyButGpuNotCompiled();

void error_gpuSimButGpuNotCompiled();

void error_gpuCacheModifiedButGpuNotCompiled();

void error_gpuCopyButQuregNotGpuAccelerated();

void error_gpuCopyButMatrixNotGpuAccelerated();

void error_gpuMemSyncQueriedButEnvNotGpuAccelerated();

void error_gpuUnexpectedlyInaccessible();

void error_gpuDeadCopyMatrixFunctionCalled();

void assert_gpuIsAccessible();

void assert_gpuHasBeenBound(bool isBound);

void assert_quregIsGpuAccelerated(Qureg qureg);

void assert_mixQuregTempGpuAllocSucceeded(qcomp* gpuPtr);

void assert_quregGpuBufferIsNotGraftedToMatrix(Qureg qureg, FullStateDiagMatr matr);

void assert_applyFullStateDiagMatrTempGpuAllocSucceeded(qcomp* gpuPtr);



/*
 * CUDA ERRORS
 */

void error_cudaCallFailed(const char* msg, const char* func, const char* caller, const char* file, int line);



/*
 * THRUST ERRORS
 */

void error_thrustTempGpuAllocFailed();



/*
 * CUQUANTUM ERRORS
 */

void error_cuQuantumInitOrFinalizedButNotCompiled();

void error_cuQuantumTempCpuAllocFailed();



/*
 * PAULI ERRORS 
 */

void error_pauliStrShiftedByIllegalAmount();



/*
 * UTILITY ERRORS 
 */

void error_nodeUnexpectedlyContainedNoElems();

void error_utilsGetBraIndGivenNonDensMatr();

void error_utilsGetPrefixIndGivenSuffixQubit();

void error_utilsGetPrefixBraIndGivenNonDensMatr();

void error_utilsGetPrefixBraIndGivenSuffixQubit();

void error_utilsIsBraQubitInSuffixGivenNonDensMatr();

void error_utilsGivenGlobalIndexOutsideNode();

void assert_utilsGivenStateVec(Qureg qureg);

void assert_utilsGivenDensMatr(Qureg qureg);

void assert_utilsGivenNonZeroEpsilon(qreal eps);



/*
 * PARSING ERRORS
 */

void error_attemptedToParseComplexFromInvalidString();

void error_attemptedToParseRealFromInvalidString();

void error_attemptedToParseOutOfRangeReal();

void error_attemptedToParsePauliStringFromInvalidString();

void error_attemptedToParseUnrecognisedPauliChar();

void error_couldNotReadFile();



/*
 * RANDOMISER ERRORS
 */

void error_randomiserGivenNonNormalisedProbList();



/*
 * PRINTER ERRORS
 */

void error_printerFailedToAllocTempMemory();

void assert_printerGivenNonNegativeNumNewlines();

void assert_printerGivenPositiveNumNewlines();



#endif // ERRORS_HPP