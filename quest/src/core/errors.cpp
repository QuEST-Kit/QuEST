/** @file
 * Defensively designed functions for checking internal preconditions, 
 * and raising internal errors. These primarily check that that
 * hardware accelerators are behaving as expected, and that runtime
 * deployment is consistent with the compiled deployment modes.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "quest/src/core/printer.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <string>

using std::string;



/*
 * INTERNAL ERROR RESPONSE
 */

void raiseInternalError(string errorMsg) {

    print(string("")
        + "\n\n"
        + "A fatal internal QuEST error occurred. "
        + errorMsg + " "
        + "Please report this to the developers. QuEST will now exit..."
        + "\n"
    );

    exit(EXIT_FAILURE);
}



/*
 * DEVELOPMENT ERRORS
 */

void error_functionNotImplemented(const char* caller) {

    string name = caller;
    raiseInternalError("The function '" + name + "' has not yet been implemented.");
}



/*
 * VALIDATION ERRORS
 */

void error_validationMessageVarWasIllFormed(string msg, string illFormedVar) {

    raiseInternalError("User input validation failed and an error string was attemptedly prepared, but an ill-formed variable was attempted substitution. This variable was \"" + illFormedVar + "\" and was attempted substitution into message:\n" + msg + "\n");
}

void error_validationMessageVarNotSubstituted(string msg, string var) {

    raiseInternalError("User input validation failed and an error string was attemptedly prepared. However, the internal variable \"" + var + "\" was unable to be found and substituted into the message, which was:\n" + msg + "\n");
}

void error_validationMessageContainedUnsubstitutedVars(string msg) {

    raiseInternalError("User input validation failed and an error string was prepared. However, the message contained unexpected (and potentially ill-formed) unsubstituted variables. The message was:\n" + msg + "\n");
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

    if (sendInd + numAmps > recvInd)
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



/*
 * ACCELERATOR ERRORS
 */

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
 * BACKEND PRECONDITION ERRORS
 */

void assert_mixedQuregIsDensityMatrix(Qureg qureg) {

    if (!qureg.isDensityMatrix)
        raiseInternalError("An internal function invoked by mixQuregs() received a statevector where a density matrix was expected.");
}

void assert_mixedQuregIsStatevector(Qureg qureg) {

    if (qureg.isDensityMatrix)
        raiseInternalError("An internal function invoked by mixQuregs() received a density matrix where a statevector was expected.");
}

void assert_mixedQuregIsDistributed(Qureg qureg) {

    if (!qureg.isDistributed)
        raiseInternalError("An internal function invoked by mixQuregs() received a non-distributed Qureg where a distributed one was expected.");
}

void assert_mixedQuregIsLocal(Qureg qureg) {

    if (qureg.isDistributed)
        raiseInternalError("An internal function invoked by mixQuregs() received a distributed Qureg where a non-distributed one was expected.");
}

void assert_mixedQuregsAreBothOrNeitherDistributed(Qureg a, Qureg b) {

    if (a.isDistributed != b.isDistributed)
        raiseInternalError("An internal function invoked by mixQuregs() received density-matrix Quregs of inconsistent distribution.");
}

void assert_mixQuregTempGpuAllocSucceeded(qcomp* gpuPtr) {

    if (gpuPtr == nullptr)
        raiseInternalError("An internal function invoked by mixQuregs() attempted to allocate temporary GPU memory but failed.");
}

void error_mixQuregsAreLocalDensMatrAndDistribStatevec() {

    raiseInternalError("An internal function invoked by mixQuregs() received a non-distributed density matrix and a distributed statevector, which is an illegal combination.");
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

    if (matr.gpuElems == qureg.gpuCommBuffer)
        raiseInternalError("An accelerator function received a FullStateDiagMatr with a GPU pointer which was a Qureg's GPU communication buffer, in a setting where the buffer was separately needed.");
}

void assert_applyFullStateDiagMatrTempGpuAllocSucceeded(qcomp* gpuPtr) {

    if (gpuPtr == nullptr)
        raiseInternalError("An internal function invoked by applying a FullStateDiagMatr upon a density matrix attempted to allocate temporary GPU memory but failed.");
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

void error_cuQuantumCompiledButNotCuda() {

    raiseInternalError(
        "Preprocessor COMPILE_CUQUANTUM was set, but COMPILE_CUDA was not. These are not intended to be exclusive fields, so "
        "a function which depended upon both flags could not continue");
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

    raiseInternalError("Attempted to initialize or finalize cuQuantum, but cuQuantum was not compiled.");
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

    raiseInternalError("Randomiser was asked to sample from a list of probabilities which did not sum to one, within epsilon error.");
}
