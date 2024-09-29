/** @file
 * Defensively designed functions for checking internal preconditions, 
 * and raising internal errors. These primarily check that that
 * hardware accelerators are behaving as expected, and that runtime
 * deployment is consistent with the compiled deployment modes.
 */

#ifndef ERRORS_HPP
#define ERRORS_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include <string>

using std::string;



/*
 * DEVELOPMENT ERRORS
 */

void error_functionNotImplemented(const char* caller);



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

void assert_validCommBounds(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps);

void assert_quregIsDistributed(Qureg qureg);

void assert_pairRankIsDistinct(Qureg qureg, int pairRank);

void assert_bufferSendRecvDoesNotOverlap(qindex sendInd, qindex recvInd, qindex numAmps);



/*
 * LOCALISER ERRORS
 */

void error_localiserNumCtrlStatesInconsistentWithNumCtrls();

void error_localiserGivenPauliTensorOrGadgetWithoutXOrY();

void error_localiserPassedStateVecToChannelComCheck();

void assert_localiserGivenDensMatr(Qureg qureg);

void assert_localiserPartialTraceGivenCompatibleQuregs(Qureg inQureg, Qureg outQureg, int numTargs);



/*
 * ACCELERATOR ERRORS
 */

void assert_numQubitsMatchesQubitStatesAndTemplateParam(int numQubits, int numQubitStates, int templateParam, string label="qubit");

void assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(int numCtrls, int numCtrlStates, int templateParam);

void assert_numTargsMatchesTemplateParam(int numTargs, int templateParam);



/*
 * BUFFER PACKING ERRORS
 */

void error_noCtrlsGivenToBufferPacker();

void assert_bufferPackerGivenIncreasingQubits(int qubit1, int qubit2, int qubit3);



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

void assert_gpuIsAccessible();

void assert_quregIsGpuAccelerated(Qureg qureg);



/*
 * CUDA ERRORS
 */

void error_cudaCallFailed(const char* msg, const char* func, const char* caller, const char* file, int line);



/*
 * CUQUANTUM ERRORS
 */

void error_cuQuantumInitOrFinalizedButNotCompiled();



/*
 * UTILITY ERRORS 
 */

void error_nodeUnexpectedlyContainedNoElems();

void error_utilsGetBraIndGivenNonDensMatr();

void error_utilsGetPrefixIndGivenSuffixQubit();

void error_utilsGetPrefixBraIndGivenNonDensMatr();

void error_utilsGetPrefixBraIndGivenSuffixQubit();

void error_utilsIsBraQubitInSuffixGivenNonDensMatr();



/*
 * PARSING ERRORS
 */

void error_attemptedToParseComplexFromInvalidString();

void error_attemptedToParseRealFromInvalidString();

void error_attemptedToParseOutOfRangeReal();

void error_attemptedToParsePauliStringFromInvalidString();

void error_attemptedToParseUnrecognisedPauliChar();

void error_couldNotReadFile();



#endif // ERRORS_HPP