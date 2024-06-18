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



/*
 * VALIDATION ERRORS
 */

void error_validationMessageVarWasIllFormed(std::string msg, std::string illFormedVar);

void error_validationMessageVarNotSubstituted(std::string msg, std::string var);

void error_validationMessageContainedUnsubstitutedVars(std::string msg);



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




/*
 * CPU ERRORS
 */

void error_cpuThreadsQueriedButEnvNotMultithreaded();



/*
 * GPU ERRORS
 */

void error_gpuQueriedButGpuNotCompiled();

void error_gpuAllocButGpuNotCompiled();

void error_gpuDeallocButGpuNotCompiled();

void error_gpuCopyButGpuNotCompiled();

void error_gpuSimButGpuNotCompiled();

void error_gpuCopyButQuregNotGpuAccelerated();

void error_gpuUnexpectedlyInaccessible();

void assert_gpuIsAccessible();

void assert_quregIsGpuAccelerated(Qureg qureg);



/*
 * CUQUANTUM ERRORS
 */

void error_cuQuantumInitOrFinalizedButNotCompiled();



/*
 * UTILITY ERRORS 
 */

void assert_shiftedQuregIsDensMatr(Qureg qureg);




#endif // ERRORS_HPP