/** @file
 * Defensively designed functions for checking internal preconditions 
 * and raising internal errors. These primarily check that that
 * hardware accelerators are behaving as expected, and that runtime
 * deployment is consistent with the compiled deployment modes.
 * 
 * @author Tyson Jones
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "quest/src/core/printer.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <string>

using std::string;



/**
 * @todo
 * This design doesn't pass info useful for debugging,
 * line erroneous caller's line-number or function name.
 * Consider refactoring to be more similar to CUDA_CHECK().
 */




/*
 * INTERNAL ERROR RESPONSE
 */

void raiseInternalError(string errorMsg) {

    print(string("")
        + "\n\n"
        + "A fatal internal QuEST error occurred. "
        + errorMsg + " "
        + "Please report this to the QuEST developers. QuEST will now exit..."
        + "\n"
    );

    exit(EXIT_FAILURE);
}



/*
 * VALIDATION ERRORS
 */

void error_validationMessageVarWasIllFormed(string msg, string illFormedVar) {

    raiseInternalError("User input validation failed and an error message was attemptedly prepared, but an ill-formed variable was attempted substitution. This variable was \"" + illFormedVar + "\" and was attempted substitution into message:\n" + msg + "\n");
}

void error_validationMessageVarNotSubstituted(string msg, string var) {

    raiseInternalError("User input validation failed and an error message was attemptedly prepared. However, the internal variable \"" + var + "\" was unable to be found and substituted into the message, which was:\n" + msg + "\n");
}

void error_validationMessageContainedUnsubstitutedVars(string msg) {

    raiseInternalError("User input validation failed and an error message was prepared. However, the message contained unexpected (and potentially ill-formed) unsubstituted variables. The message was:\n" + msg + "\n");
}

void error_validationEncounteredUnsupportedDistributedDenseMatrix() {

    raiseInternalError("User input validation processed a matrix it believed was both dense and distributed, though no such data structure currently exists.");
}

void error_validationListUniquenessCheckExceededMaskSize() {

    raiseInternalError("User input validation was checking uniqueness of an index list using bitmasks but encountered an index larger than the number of bits in the mask. This should have been caught by prior validation.");
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

    raiseInternalError("Attempted to obtain memory necessary to allocate a distributed object's single-node partition but it overflowed size_t despite prior validation.");
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

void error_commGivenInconsistentNumSubArraysANodes() {

    raiseInternalError("A distributed function was given a different number of per-node subarray lengths than exist nodes.");
}

void error_commNumMessagesExceedTagMax() {

    raiseInternalError("A function attempted to communicate via more messages than permitted (since there would be more uniquely-tagged messages than the tag upperbound).");
}

void assert_commBoundsAreValid(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps) {

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

void assert_commPayloadIsPowerOf2(qindex numAmps) {

    if (!isPowerOf2(numAmps))
        raiseInternalError("A communication function was given a payload which was unexpectedly not a power of 2, as implied by preconditions.");
}

void assert_commQuregIsDistributed(Qureg qureg) {

    if (!qureg.isDistributed)
        error_commButQuregNotDistributed();
}

void assert_commFullStateDiagMatrIsDistributed(FullStateDiagMatr matr) {
    
    if (!matr.isDistributed)
        raiseInternalError("A function attempted to invoke communication of a FullStateDiagMatr which was not distributed.");
}

void assert_pairRankIsDistinct(Qureg qureg, int pairRank) {

    if (pairRank == qureg.rank)
        error_commWithSameRank();
}

void assert_bufferSendRecvDoesNotOverlap(qindex sendInd, qindex recvInd, qindex numAmps) {

    if (sendInd < recvInd + numAmps)
        raiseInternalError("A distributed function attempted to send and receive portions of the buffer which overlapped.");
}

void assert_receiverCanFitSendersEntireState(Qureg receiver, Qureg sender) {

    if (receiver.numAmpsPerNode < sender.numAmps)
        raiseInternalError("A distributed function attempted to broadcast a Qureg's entire state into another Qureg's buffer, which could not fit all global amps.");
}

void assert_receiverCanFitSendersEntireElems(Qureg receiver, FullStateDiagMatr sender) {

    if (receiver.numAmpsPerNode < sender.numElems)
        raiseInternalError("A distributed function attempted to broadcast the entirety of a FullStateDiagMatr's elements into a Qureg's buffer, which is too small to contain all global elements.");
}



/*
 * LOCALISER ERRORS
 */

void error_localiserNumCtrlStatesInconsistentWithNumCtrls() {

    raiseInternalError("An inconsistent number of ctrls and ctrlStates were passed to a function in localiser.cpp.");
}

void error_localiserGivenPauliTensorOrGadgetWithoutXOrY() {

    raiseInternalError("The localiser was asked to simulate a Pauli tensor or gadget which contained no X or Y Paulis, which is a special case reserved for phase gadgets.");
}

void error_localiserPassedStateVecToChannelComCheck() {

    raiseInternalError("The localiser queried whether a channel would invoke communication upon a statevector.");
}

void error_localiserGivenDistribMatrixAndLocalQureg() {

    raiseInternalError("A localiser function was given a distributed FullStateDiagMatr but a non-distributed Qureg, which are incompatible.");
}

void error_localiserFailedToAllocTempMemory() {

    raiseInternalError("A localiser function attempted and failed to allocate temporary memory.");
}

void error_localiserGivenPauliStrWithoutXorY() {

    raiseInternalError("A localiser function was given a PauliStr which unexpectedly contained no X or Y Paulis.");
}

void error_localiserGivenNonUnityGlobalFactorToZTensor() {

    raiseInternalError("A localiser function to apply a PauliStr (as a tensor, not a gadget) was given a PauliStr containing only Z and I, along with a non-unity global factor. This is an illegal combination.");
}

void assert_localiserSuccessfullyAllocatedTempMemory(qcomp* ptr, bool isGpu) {

    if (mem_isAllocated(ptr))
        return;

    string platform = (isGpu)? "GPU" : "CPU";
    raiseInternalError("A localiser function attempted and failed to allocate temporary " + platform + " memory.");
}

void assert_localiserGivenStateVec(Qureg qureg) {

    if (qureg.isDensityMatrix)
        raiseInternalError("The localiser received a density matrix to a function defined strictly for statevectors.");
}

void assert_localiserGivenDensMatr(Qureg qureg) {

    if (!qureg.isDensityMatrix)
        raiseInternalError("The localiser received a statevector to a function defined only upon density matrices.");
}

void assert_localiserPartialTraceGivenCompatibleQuregs(Qureg inQureg, Qureg outQureg, int numTargs) {

    if (!inQureg.isDensityMatrix || !outQureg.isDensityMatrix)
        raiseInternalError("A non-density matrix was given to localiser's partial trace function.");

    if (inQureg.isDistributed != outQureg.isDistributed)
        raiseInternalError("Quregs of different distributions (one distributed, one not) was given to localiser's partial trace function.");

    if (inQureg.isGpuAccelerated != outQureg.isGpuAccelerated)
        raiseInternalError("Quregs of different GPU-accelerations (one accelerated, the other not) was given to localiser's partial trace function.");

    if (inQureg.numQubits - numTargs != outQureg.numQubits)
        raiseInternalError("Inconsistent Qureg sizes and number of traced qubits given to localiser's partial trace function.");
}

void error_calcFidStateVecDistribWhileDensMatrLocal() {

    raiseInternalError("A localiser function attempted to compute the fidelity between a local density matrix and a distributed statevector, which is an illegal combination.");
}

void assert_localiserDistribQuregSpooferGivenValidQuregs(Qureg local, Qureg distrib) {

    if (local.isDistributed || !distrib.isDistributed)
        raiseInternalError("The localiser attempted to spoof a distributed Qureg using a local one, but received invalid Qureg distributions.");
}



/*
 * BACKEND PRECONDITION ERRORS
 */

void assert_highPauliStrSumMaskIsZero(PauliStrSum sum) {

    // this is safe to enumerate sum here, since the calling function 
    // enumerates sum repeatedly (once for each amplitude, albeit in
    // parallel). ergo at absolute worst, this check doubles runtime,
    // with the slowdown exponentially vanishing with increasing #qubits

    for (qindex n=0; n<sum.numTerms; n++)
        if (sum.strings[n].highPaulis != 0)
            raiseInternalError("A CPU or GPU subroutine received a PauliStrSum within which a PauliStr contained a non-identity Pauli in the 'highPaulis' mask, which is illegal for this function.");
}

void assert_numQubitsMatchesQubitStatesAndTemplateParam(int numQubits, int numQubitStates, int templateParam, string label) {

    if (numQubits != numQubitStates)
        raiseInternalError("A CPU or GPU subroutine received an inconsistent number of " + label + "s and " + label + "-states from accelerator.cpp.");

    // template parameter of -1 is always valid (it indicates the routine has not been compile-time optimised)
    if (templateParam == -1)
        return;

    if (templateParam != numQubits)
        raiseInternalError("A CPU or GPU subroutine received a number of qubits inconsistent with its compile-time template parameter, as dispatched by accelerator.cpp.");
}

void assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(int numCtrls, int numCtrlStates, int templateParam) {

    assert_numQubitsMatchesQubitStatesAndTemplateParam(numCtrls, numCtrlStates, templateParam, "control");
}

void assert_numTargsMatchesTemplateParam(int numTargs, int templateParam) {

    // template parameter of -1 is always valid (it indicates the routine has not been compile-time optimised)
    if (templateParam == -1)
        return;

    if (templateParam != numTargs)
        raiseInternalError("A CPU or GPU subroutine received a number of targets inconsistent with its compile-time template parameter, as dispatched by accelerator.cpp.");
}

void assert_exponentMatchesTemplateParam(qcomp exponent, bool hasPower) {

    // require hasPower==false => exponent==1
    if (!hasPower && exponent != qcomp(1,0))
        raiseInternalError("A CPU or GPU subroutine received a matrix exponent that was inconsistent with its compile-time template parameter, as dispatched by accelerator.cpp.");
}

void assert_exponentMatchesTemplateParam(qcomp exponent, bool hasPower, bool useRealPow) {

    assert_exponentMatchesTemplateParam(exponent, hasPower);

    if (!hasPower && useRealPow)
        raiseInternalError("A CPU or GPU subroutine received an invalid combination of HasPower=false and UseRealPow=false template parameters");
}



/*
 * BUFFER PACKING ERRORS
 */

void error_noCtrlsGivenToBufferPacker() {

    raiseInternalError("A function attempted to (superfluously) pack the communication buffer but specified no control qubits, which would lead to the buffer being entirely filled and leave no room to receive amps.");
}

void assert_bufferPackerGivenIncreasingQubits(int qubit1, int qubit2, int qubit3) {

    if (qubit1 >= qubit2 || qubit2 >= qubit3)
        raiseInternalError("A function attempted to pack a buffer using non-increasing qubit indices.");
}



/*
 * ACCELERATOR ERRORS
 */

void assert_mixedQuregIsDensityMatrix(Qureg qureg) {

    if (!qureg.isDensityMatrix)
        raiseInternalError("An internal function invoked by mixQureg() received a statevector where a density matrix was expected.");
}

void assert_mixedQuregIsStatevector(Qureg qureg) {

    if (qureg.isDensityMatrix)
        raiseInternalError("An internal function invoked by mixQureg() received a density matrix where a statevector was expected.");
}

void assert_mixedQuregIsDistributed(Qureg qureg) {

    if (!qureg.isDistributed)
        raiseInternalError("An internal function invoked by mixQureg() received a non-distributed Qureg where a distributed one was expected.");
}

void assert_mixedQuregIsLocal(Qureg qureg) {

    if (qureg.isDistributed)
        raiseInternalError("An internal function invoked by mixQureg() received a distributed Qureg where a non-distributed one was expected.");
}

void assert_mixedQuregsAreBothOrNeitherDistributed(Qureg a, Qureg b) {

    if (a.isDistributed != b.isDistributed)
        raiseInternalError("An internal function invoked by mixQureg() received density-matrix Quregs of inconsistent distribution.");
}

void assert_mixQuregTempGpuAllocSucceeded(qcomp* gpuPtr) {

    if (!mem_isAllocated(gpuPtr))
        raiseInternalError("An internal function invoked by mixQureg() attempted to allocate temporary GPU memory but failed.");
}

void error_mixQuregsAreLocalDensMatrAndDistribStatevec() {

    raiseInternalError("An internal function invoked by mixQureg() received a non-distributed density matrix and a distributed statevector, which is an illegal combination.");
}

void assert_fullStateDiagMatrIsLocal(FullStateDiagMatr matr) {

    if (matr.isDistributed)
        raiseInternalError("An accelerator function received a distributed FullStateDiagMatr where a non-distributed one was expected.");
}

void assert_fullStateDiagMatrIsDistributed(FullStateDiagMatr matr) {

    if (!matr.isDistributed)
        raiseInternalError("An accelerator function received a non-distributed FullStateDiagMatr where a distributed one was expected.");
}

void assert_acceleratorQuregIsDistributed(Qureg qureg) {

    if (!qureg.isDistributed)
        raiseInternalError("An accelerator function invoked received non-distributed Qureg where a distributed one was expected.");
}

void assert_quregAndFullStateDiagMatrAreBothOrNeitherDistrib(Qureg qureg, FullStateDiagMatr matr) {

    if (qureg.isDistributed != matr.isDistributed)
        raiseInternalError("An accelerator function unexpectedly received a Qureg and FullStateDiagMatr with different distributions.");
}

void assert_quregGpuBufferIsNotGraftedToMatrix(Qureg qureg, FullStateDiagMatr matr) {

    // permit both pointers to be null-ptr, of course
    if (!mem_isAllocated(matr.gpuElems))
        return;

    if (matr.gpuElems == qureg.gpuCommBuffer)
        raiseInternalError("An accelerator function received a FullStateDiagMatr with a GPU pointer which was a Qureg's GPU communication buffer, in a setting where the buffer was separately needed.");
}

void assert_applyFullStateDiagMatrTempGpuAllocSucceeded(qcomp* gpuPtr) {

    if (!mem_isAllocated(gpuPtr))
        raiseInternalError("An internal function invoked by applying a FullStateDiagMatr upon a density matrix attempted to allocate temporary GPU memory but failed.");
}

void assert_calcFidStateVecIsLocal(Qureg qureg) {

    if (qureg.isDistributed)
        raiseInternalError("An accelerator function involved with calculating the fidelity between a density matrix and statevector was given an illegally distributed statevector.");
}

void assert_calcFidTempGpuAllocSucceeded(qcomp* ptr) {

    if (!mem_isAllocated(ptr))
        raiseInternalError("An accelerator function involved with calculating the fidelity between a density matrix and statevector failed to allocate temporary GPU memory.");
}

void assert_calcExpecDiagTempGpuAllocSucceeded(qcomp* ptr) {

    if (!mem_isAllocated(ptr))
        raiseInternalError("An accelerator function involved with calculating the expectation value of a FullStateDiagMatr upon a density matrix failed to allocate temporary GPU memory.");
}

void assert_innerProductedSameDimQuregsHaveSameGpuAccel(Qureg quregA, Qureg quregB) {

    if (quregA.isGpuAccelerated != quregB.isGpuAccelerated)
        raiseInternalError("The accelerator was asked to compute the inner product between a GPU-accelerated and non-accelerated Qureg, where were both statevectors or density matrices, which is an illegal combination.");
}

void assert_partialTraceQuregsAreIdenticallyDeployed(Qureg inQureg, Qureg outQureg) {

    bool valid = (
        inQureg.isDistributed == outQureg.isDistributed &&
        inQureg.isGpuAccelerated == outQureg.isGpuAccelerated);

    if (!valid)
        raiseInternalError("An accelerator function involved in computing the partial trace received Quregs with differing distributions or GPU-accelerations.");
}




/*
 * BACKEND PRECONDITION ERRORS
 */

void assert_quregAndFullStateDiagMatrHaveSameDistrib(Qureg qureg, FullStateDiagMatr matr) {

    if (qureg.isDistributed != matr.isDistributed)
        raiseInternalError("A Qureg and FullStateDiagMatr had unexpectedly mismatching distribution statuses.");
}

void assert_quregDistribAndFullStateDiagMatrLocal(Qureg qureg, FullStateDiagMatr matr) {

    if (!qureg.isDistributed)
        raiseInternalError("The Qureg was unexpectedly non-distributed.");
        
    if (matr.isDistributed)
        raiseInternalError("The FullStateDiagMatr was unexpectedly distributed.");
}

void assert_superposedQuregDimsAndDeploysMatch(Qureg facOut, Qureg in1, Qureg in2) {

    if (
        facOut.isDistributed    != in1.isDistributed    || in1.isDistributed    != in2.isDistributed    ||
        facOut.isDensityMatrix  != in1.isDensityMatrix  || in1.isDensityMatrix  != in2.isDensityMatrix  ||
        facOut.isGpuAccelerated != in1.isGpuAccelerated || in1.isGpuAccelerated != in2.isGpuAccelerated ||
        facOut.numQubits        != in1.numQubits        || in1.numQubits        != in2.numQubits
    )
        raiseInternalError("An internal function *_setQuregToSuperposition() received Quregs of mismatching dimensions and/or deployments.");
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

    raiseInternalError("A function (likely Qureg or CompMatr creation) attempted to allocate GPU memory but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuDeallocButGpuNotCompiled() {

    raiseInternalError("A function (likely Qureg or CompMatr destruction) attempted to deallocate GPU memory but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuCopyButGpuNotCompiled() {

    raiseInternalError("A function attempted to access GPU memory but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuSimButGpuNotCompiled() {

    raiseInternalError("A function attempted to dispatch a simulation to the GPU but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuCacheModifiedButGpuNotCompiled() {

    raiseInternalError("A function attempted to allocate or clear the GPU cache memory but QuEST was not compiled with GPU acceleration enabled.");
}

void error_gpuCopyButQuregNotGpuAccelerated() {

    raiseInternalError("A function attempted to access GPU memory of a Qureg which is not GPU accelerated.");
}

void error_gpuCopyButMatrixNotGpuAccelerated() {

    raiseInternalError("A function attempted to access GPU memory of a matrix (e.g. a CompMatr or DiagMatr) which is not GPU accelerated.");
}

void error_gpuUnexpectedlyInaccessible() {

    raiseInternalError("A function internally assumed (as a precondition) that QuEST was compiled with GPU-acceleration enabled, and that one was physically accessible, though this was untrue.");
}

void error_gpuMemSyncQueriedButEnvNotGpuAccelerated() {

    raiseInternalError("A function checked whether persistent GPU memory (such as in a CompMatr) had been synchronised, but the QuEST environment is not GPU accelerated.");  
}

void error_gpuDeadCopyMatrixFunctionCalled() {

    raiseInternalError("The internal GPU function copyMatrixIfGpuCompiled() was called, though is intended as dead-code - matrices needing copying to GPU should be stored as flat row-wise lists.");
}

void assert_quregIsGpuAccelerated(Qureg qureg) {

    if (!qureg.isGpuAccelerated)
        error_gpuCopyButQuregNotGpuAccelerated();
}

void assert_gpuIsAccessible() {

    if (!gpu_isGpuCompiled() || !gpu_isGpuAvailable())
        error_gpuUnexpectedlyInaccessible();
}

void assert_gpuHasBeenBound(bool isBound) {

    if (!isBound)
        raiseInternalError("An internal GPU-querying function was illegally called before local GPUs had been bound to local MPI processes.");
}



/*
 * CUDA ERRORS
 */

void error_cudaCallFailed(const char* msg, const char* func, const char* caller, const char* file, int line) {

    // using operator overloads to cast const char[] literals to std::string, to concat with const char*.
    string err = "";
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
 * THRUST ERRORS
 */


void error_thrustTempGpuAllocFailed() {

    raiseInternalError("Thrust failed to allocate temporary GPU memory.");
}



/*
 * CUQUANTUM ERRORS
 */

void error_cuQuantumInitOrFinalizedButNotCompiled() {

    raiseInternalError("Attempted to initialise or finalise cuQuantum, but cuQuantum was not compiled.");
}

void error_cuQuantumTempCpuAllocFailed() {

    raiseInternalError("Attempted allocation of temporary host-memory for a cuQuantum routine failed.");
}



/*
 * PAULI ERRORS 
 */

void error_pauliStrShiftedByIllegalAmount() {

    raiseInternalError("A PauliStr was attemptedly shifted (likely invoked by its application upon a density matrix) by an illegal amount (e.g. negative, or that exceeding the PauliStr bitmask length).");
}



/*
 * UTILITY ERRORS 
 */

void error_nodeUnexpectedlyContainedNoElems() {

    raiseInternalError("A function queried which distributed elements within a specified range overlapped the node's stored elements, but the node contained no overlap. This situation should have been prior handled."); 
}

void error_utilsGetBraIndGivenNonDensMatr() {

    raiseInternalError("A function attempted to obtain the index of a bra-qubit of a state-vector, as if it were a density matrix.");
}

void error_utilsGetPrefixIndGivenSuffixQubit() {

    raiseInternalError("A function passed a suffix qubit to a utilities function expecting a prefix qubit.");
}

void error_utilsGetPrefixBraIndGivenNonDensMatr() {

    raiseInternalError("A function attempted to obtain the prefix index of a bra-qubit of a state-vector, as if it were a density matrix.");

}

void error_utilsGetPrefixBraIndGivenSuffixQubit() {

    raiseInternalError("A function attmpted to obtain the prefix index of a bra-qubit, but passed a suffix qubit.");
}

void error_utilsIsBraQubitInSuffixGivenNonDensMatr() {

    raiseInternalError("A functiion queried whether a qubit's corresponding bra-qubit was in the suffix substate, but the Qureg was not a density matrix.");
}

void error_utilsGivenGlobalIndexOutsideNode() {

    // this error might NOT be thrown by all nodes, causing non-consensus crash. Eh!
    raiseInternalError("A utility function was asked for the corresponding local index of a global index which did not exist in the calling node.");
}

void assert_utilsGivenStateVec(Qureg qureg) {

    if (qureg.isDensityMatrix)
        raiseInternalError("A utility function was given a density matrix where a statevector was expected.");
}

void assert_utilsGivenDensMatr(Qureg qureg) {

    if (!qureg.isDensityMatrix)
        raiseInternalError("A utility function was given a statevector where a density matrix was expected.");
}

void assert_utilsGivenNonZeroEpsilon(qreal eps) {

    if (eps == 0)
        raiseInternalError("A utility function (e.g. isUnitary) received an epsilon of zero, which should have precluded it being called.");
}



/*
 * PARSING ERRORS
 */

void error_attemptedToParseComplexFromInvalidString() {

    raiseInternalError("A function attempted to parse a string to a qcomp but the string was not validly formatted. This should have been caught by prior user validation.");
}

void error_attemptedToParseRealFromInvalidString() {

    raiseInternalError("A function attempted to parse a string to a qreal but the string was not validly formatted. This should have been caught by prior user validation.");
}

void error_attemptedToParseOutOfRangeReal() {

    raiseInternalError("A function attempted to parse a string to a qreal but the numerical value of the string literal exceeded the range of the qreal. This should have been caught by prior user validation.");
}

void error_attemptedToParsePauliStringFromInvalidString() {

    raiseInternalError("A function attempted to parse a string as a sequence of Pauli operators but the string was not validly formatted. This should have been caught by prior user validation.");
}

void error_attemptedToParseUnrecognisedPauliChar() {

    raiseInternalError("A function attempted to parse an unrecognised character as a Pauli operator. This should have been caught by prior validation.");
}

void error_couldNotReadFile() {

    raiseInternalError("A function failed to open and read a file that previous validation confirmed was readable.");
}



/*
 * RANDOMISER ERRORS
 */

void error_randomiserGivenNonNormalisedProbList() {

    raiseInternalError("The randomiser was asked to sample from a list of probabilities which did not sum to one, within epsilon error.");
}



/*
 * PRINTER ERRORS
 */

void error_printerFailedToAllocTempMemory() {

    raiseInternalError("A printer utility attempted and failed to allocate temporary memory, which likely results from the attemptedly printed object being too large.");
}

void assert_printerGivenNonNegativeNumNewlines() {

    int min = 0;

    if (printer_getNumTrailingNewlines() < min)
        raiseInternalError("A printer utility attempted to print a negative number of newlines, informed by the user-set number of trailing newlines, which should have been caught by prior validation.");
}

void assert_printerGivenPositiveNumNewlines() {

    int min = 1;

    if (printer_getNumTrailingNewlines() < min)
        raiseInternalError("A printer utility attempted to print one fewer than the user-set number of trailing newlines; but that number was zero! This violates prior validation.");
}
