/** @file
 * Defensively designed functions for checking internal preconditions, 
 * and raising internal errors. These primarily check that that
 * hardware accelerators are behaving as expected, and that runtime
 * deployment is consistent with the compiled deployment modes.
 */

#include "types.h"
#include "qureg.h"
#include "../comm/comm_config.hpp"
#include "../gpu/gpu_config.hpp"

#include <iostream>
#include <string>



/*
 * INTERNAL ERROR RESPONSE
 */

void raiseInternalError(std::string errorMsg) {

    if (comm_getRank() == 0)
        std::cout 
            << "\n\n"
            << "A fatal internal QuEST error occurred. "
            << errorMsg << " "
            << "Please report this to the developers. QuEST will now exit..."
            << "\n"
            << std::endl;

    exit(EXIT_FAILURE);
}



/*
 * VALIDATION ERRORS
 */

void error_validationMessageVarWasIllFormed(std::string msg, std::string illFormedVar) {

    raiseInternalError("User input validation failed and an error string was attemptedly prepared, but an ill-formed variable was attempted substitution. This variable was \"" + illFormedVar + "\" and was attempted substitution into message:\n" + msg + "\n");
}

void error_validationMessageVarNotSubstituted(std::string msg, std::string var) {

    raiseInternalError("User input validation failed and an error string was attemptedly prepared. However, the internal variable \"" + var + "\" was unable to be found and substituted into the message, which was:\n" + msg + "\n");
}

void error_validationMessageContainedUnsubstitutedVars(std::string msg) {

    raiseInternalError("User input validation failed and an error string was prepared. However, the message contained unexpected (and potentially ill-formed) unsubstituted variables. The message was:\n" + msg + "\n");
}



/*
 * ENVIRONMENT ERRORS
 */

void error_allocOfQuESTEnvFailed() {

    raiseInternalError("Attempted memory allocation for the newly created QuESTEnv unexpectedly failed.");
}



/*
 * MEMORY ERRORS
 */

void error_memSizeQueriedButWouldOverflow() {

    raiseInternalError("Attempted to obtain memory necessary to allocate local Qureg partition but it overflowed size_t despite prior validation.");
}



/*
 * COMMUNICATION ERRORS
 */

void error_commNotInit() {

    raiseInternalError("MPI was queried but the MPI environment had not yet been initialised.");
}

void error_commAlreadyInit() {

    raiseInternalError("The MPI communication environment was attemptedly re-initialised despite the QuEST environment already existing.");
}

void error_commButEnvNotDistributed() {

    raiseInternalError("A function attempted to invoke communication despite QuEST being compiled in non-distributed mode.");
}

void error_commButQuregNotDistributed() {

    raiseInternalError("A function attempted to invoke communication of a Qureg which was not distributed.");
}

void error_commOutOfBounds() {

    raiseInternalError("A function invoked communication which attempted to exchange amplitudes between arrays at invalid bounds.");
}

void error_commWithSameRank() {

    raiseInternalError("A distributed function attempted to communicate to a pair rank equal to its own rank.");
}

void assert_validCommBounds(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps) {

    bool valid = (
        sendInd >= 0 &&
        recvInd >= 0 &&
        numAmps >  0 &&
        sendInd + numAmps <= qureg.numAmpsPerNode &&
        recvInd + numAmps <= qureg.numAmpsPerNode
    );

    if (!valid)
        error_commOutOfBounds();
}

void assert_quregIsDistributed(Qureg qureg) {

    if (!qureg.isDistributed)
        error_commButQuregNotDistributed();
}

void assert_pairRankIsDistinct(Qureg qureg, int pairRank) {

    if (pairRank == qureg.rank)
        error_commWithSameRank();
}



/*
 * CPU ERRORS
 */

void error_cpuThreadsQueriedButEnvNotMultithreaded() {

    raiseInternalError("A function attempted to query CPU thread information but QuEST is not running with multithreading enabled.");
}



/*
 * GPU ERRORS
 */

void error_gpuQueriedButGpuNotCompiled() {

    raiseInternalError("A function attempted to query GPU properties but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuSyncedButGpuNotCompiled() {

    raiseInternalError("A function attempted to synchronise the GPU, but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuAllocButGpuNotCompiled() {

    raiseInternalError("A function (likely Qureg or CompMatrN creation) attempted to allocate GPU memory but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuDeallocButGpuNotCompiled() {

    raiseInternalError("A function (likely Qureg or CompMatrN destruction) attempted to deallocate GPU memory but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuCopyButGpuNotCompiled() {

    raiseInternalError("A function attempted to access GPU memory but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuSimButGpuNotCompiled() {

    raiseInternalError("A function attempted to dispatch a simulation to the GPU but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuCopyButQuregNotGpuAccelerated() {

    raiseInternalError("A function attempted to access GPU memory of a Qureg which is not GPU accelerated.");
}

void error_gpuCopyButCompMatrNotGpuAccelerated() {

    raiseInternalError("A function attempted to access GPU memory of a CompMatrN which is not GPU accelerated.");
}

void error_gpuUnexpectedlyInaccessible() {

    raiseInternalError("A function internally assumed (as a precondition) that QuEST was compiled with GPU-acceleration enabled, and that one was physically accessible, though this was untrue.");
}

void error_gpuMemSyncQueriedButEnvNotGpuAccelerated() {

    raiseInternalError("A function checked whether persistent GPU memory (such as in a CompMatrN) had been synchronised, but the QuEST environment is not GPU accelerated.");  
}

void assert_quregIsGpuAccelerated(Qureg qureg) {

    if (!qureg.isGpuAccelerated)
        error_gpuCopyButQuregNotGpuAccelerated();
}

void assert_gpuIsAccessible() {

    if (!gpu_isGpuCompiled() || !gpu_isGpuAvailable())
        error_gpuUnexpectedlyInaccessible();
}



/*
 * CUDA ERRORS
 */

void error_cudaCallFailed(const char* msg, const char* func, const char* caller, const char* file, int line) {

    // using operator overloads to cast const char[] literals to std::string, to concat with const char*.
    std::string err = "";
    err += "A CUDA (or cuQuantum) API function (\"";
    err += func;
    err += "\", called by \"";
    err += caller;
    err += "()\" at line ";
    err += std::to_string(line);
    err += " of file ";
    err += file;
    err += ") unexpectedly failed with error message: \"";
    err += msg;
    err += "\". ";
    raiseInternalError(err);
}



/*
 * CUQUANTUM ERRORS
 */

void error_cuQuantumInitOrFinalizedButNotCompiled() {

    raiseInternalError("Attempted to initialize or finalize cuQuantum, but cuQuantum was not compiled.");
}



/*
 * UTILITY ERRORS 
 */

void assert_shiftedQuregIsDensMatr(Qureg qureg) {
    if (!qureg.isDensityMatrix)
        raiseInternalError("A function attempted to obtain the shifted Choi indices of a state-vector, as if it were a density matrix.");
}
