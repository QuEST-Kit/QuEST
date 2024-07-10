/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#include "quest/include/modes.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"
#include "quest/include/structures.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>



/*
 * INVALID INPUT ERROR MESSAGES
 * which can contain variables with syntax ${VAR1} ${VAR2}, substituted at error-throw with
 * runtime parameters via assertThat(..., {{"${VAR1}",1}, {"${VAR2}",2}}, ...)
 */

namespace report {


    /*
     *  ENVIRONMENT CREATION
     */

    std::string QUEST_ENV_ALREADY_INIT =
        "The QuEST environment has already been initialised. This can only be performed once during program execution.";

    std::string QUEST_ENV_ALREADY_FINAL =
        "The QuEST environment has already been finalised, and can thereafter never be re-initialised since this leads to undefined MPI behaviour.";


    std::string INVALID_OPTION_FOR_ENV_IS_DISTRIB =
        "Argument 'useDistribution' must be 1 or 0 to respectively indicate whether or not to distribute the new environment, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    std::string INVALID_OPTION_FOR_ENV_IS_GPU_ACCEL =
        "Argument 'useGpuAcceleration' must be 1 or 0 to respectively indicate whether or not to GPU-accelerate the new environment, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    std::string INVALID_OPTION_FOR_ENV_IS_MULTITHREAD =
        "Argument 'useMultithreading' must be 1 or 0 to respectively indicate whether or not to enable multithreading in the new environment, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";


    std::string CANNOT_MULTITHREAD_ENV_WITHOUT_OPENMP_COMPILATION =
        "Cannot create a multithreaded environment because the QuEST source was not compiled with OpenMP enabled.";

    std::string CANNOT_DISTRIB_ENV_WITHOUT_MPI_COMMPILATION =
        "Cannot create a distributed environment because the QuEST source was not compiled with MPI enabled.";

    std::string CANNOT_GPU_ACCEL_ENV_WITH_GPU_COMPILATION =
        "Cannot create a GPU-accelerated environment because the QuEST source was not compiled with GPU acceleration.";


    std::string CANNOT_GPU_ACCEL_ENV_WITH_NO_AVAILABLE_GPUS =
        "Cannot create a GPU-accelerated environment because there are no GPUs available.";

    std::string CANNOT_DISTRIB_ENV_BETWEEN_NON_POW_2_NODES =
        "Cannot distribute QuEST between ${NUM_NODES} nodes; must use a power-of-2 number of nodes.";

    std::string CUQUANTUM_DEPLOYED_ON_GPU_WITHOUT_MEM_POOLS =
        "Cannot use cuQuantum since your GPU does not support memory pools. Please recompile with cuQuantum disabled to fall-back to using Thrust and custom kernels.";

    
    /*
     * EXISTING QUESTENV
     */

    std::string QUEST_ENV_NOT_INIT =
        "The QuEST environment is not initialised. Please first call initQuESTEnv() or initCustomQuESTEnv().";


    /*
     * QUREG CREATION
     */

    std::string NON_POSITIVE_NUM_QUBITS_IN_INIT_QUREG =
        "Cannot create Qureg of ${NUM_QUBITS} qubits; must contain one or more qubits.";


    std::string NEW_QUREG_AMPS_WOULD_EXCEED_QINDEX = 
        "Cannot create Qureg of ${NUM_QUBITS} qubits: the statevector would contain more amplitudes (2^${NUM_QUBITS}) than the maximum which can be addressed by the qindex type (2^${MAX_QUBITS}). See reportQuESTEnv().";

    std::string NEW_DENS_QUREG_AMPS_WOULD_EXCEED_QINDEX = 
        "Cannot create density Qureg of ${NUM_QUBITS} qubits: the density matrix would contain more amplitudes (4^${NUM_QUBITS}) than can be addressed by the qindex type (4^${MAX_QUBITS}). See reportQuESTEnv().";

    std::string NEW_QUREG_LOCAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits distributed over ${NUM_NODES} nodes because the necessary local memory (in bytes) of each node would overflow size_t. In this deployment, the maximum number of similar qubits is ${MAX_QUBITS}. See reportQuESTEnv().";


    std::string NEW_DISTRIB_QUREG_HAS_TOO_FEW_AMPS = 
        "Cannot create a distributed Qureg of only ${NUM_QUBITS} qubits between ${NUM_NODES} nodes, because each node would contain fewer than one amplitude of the statevector. The minimum size is ${MIN_QUBITS} qubits; see reportQuESTEnv(). Consider disabling distribution for this Qureg.";

    std::string NEW_DENS_DISTRIB_QUREG_HAS_TOO_FEW_AMPS = 
        "Cannot create a distributed density Qureg of only ${NUM_QUBITS} qubits between ${NUM_NODES} nodes, because each node would contain fewer than a column's worth of amplitudes of the density matrix. The minimum size is ${MIN_QUBITS} qubits; see reportQuESTEnv(). Consider disabling distribution for this Qureg.";


    std::string INVALID_OPTION_FOR_QUREG_IS_DENS_MATR = 
        "Argument 'isDensityMatrix' must be 1 or 0 to respectively indicate whether the Qureg should be instantiated as a potentially-mixed density matrix or a strictly-pure state-vector.";

    std::string INVALID_OPTION_FOR_QUREG_IS_DISTRIB = 
        "Argument 'useDistribution' must be 1 or 0 to respectively indicate whether or not to distribute the new Qureg, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    std::string INVALID_OPTION_FOR_QUREG_IS_GPU_ACCEL = 
        "Argument 'useGpuAcceleration' must be 1 or 0 to respetively indicate whether or not to GPU-accelerate the new Qureg, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    std::string INVALID_OPTION_FOR_QUREG_IS_MULTITHREAD = 
        "Argument 'useMultithreading' must be 1 or 0 to respectively indicate whether or not to use multithreading when processing the new Qureg, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";


    std::string NEW_DISTRIB_QUREG_IN_NON_DISTRIB_ENV =
        "Cannot distribute a Qureg when in a non-distributed QuEST environment.";

    std::string NEW_GPU_QUREG_IN_NON_GPU_ACCEL_ENV =
        "Cannot allocate a Qureg to a GPU when in a non-GPU-accelerated QuEST environment.";

    std::string NEW_MULTITHREAD_QUREG_IN_NON_MULTITHREAD_ENV =
        "Cannot enable multithreaded processing of a Qureg created in a non-multithreaded QuEST environment.";


    std::string NEW_GPU_QUREG_CANNOT_USE_MULTITHREADING = 
        "Cannot simultaneously GPU-accelerate and multithread a Qureg. Please disable multithreading, or set it to ${AUTO_DEPLOYMENT_FLAG} for QuEST to automatically disable it when deploying to GPU.";


    std::string NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CPU_MEM =
        "The non-distributed Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits would be too large (${QCOMP_BYTES} * ${EXP_BASE}^${NUM_QUBITS} bytes) to fit into a single node's RAM (${RAM_SIZE} bytes). See reportQuESTEnv(), and consider using distribution.";

    std::string NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CPU_MEM = 
        "The distributed Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits, together with its commnication buffer, would be too large (2 * ${QCOMP_BYTES} * ${EXP_BASE}^${NUM_QUBITS} bytes) to fit into the combined RAM of all ${NUM_NODES} nodes (${NUM_NODES} * ${RAM_SIZE} bytes). See reportQuESTEnv().";

    std::string NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CURRENT_GPU_MEM =
        "The non-distributed GPU-accelerated Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits would be too large (${QCOMP_BYTES} * ${EXP_BASE}^${NUM_QUBITS} bytes) to fit into a single node's available GPU memory (currently ${MIN_VRAM_AVAIL} bytes free). Consider additionally using distribution, or disabling GPU-acceleration (though this may greatly increase runtime).";

    std::string NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CURRENT_GPU_MEM =
        "The distributed GPU-accelerated Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits, together with its commnication buffer, is too large; one or more of the ${NUM_GPUS} GPUs has insufficient available memory (only ${MIN_VRAM_AVAIL} bytes) to store its Qureg partition (${QCOMP_BYTES} * 2^${LOG2_NUM_AMPS} bytes) bytes. Consider disabling GPU-acceleration.";


    std::string FAILED_NEW_ALLOC_OF_CPU_AMPS = 
        "Allocation of memory to store the CPU amplitudes failed.";

    std::string FAILED_NEW_ALLOC_OF_GPU_AMPS = 
        "Allocation of memory to store the GPU amplitudes failed.";

    std::string FAILED_NEW_ALLOC_OF_CPU_COMM_BUFFER = 
        "Allocation of memory for the distributed CPU communication buffer failed.";

    std::string FAILED_NEW_ALLOC_OF_GPU_COMM_BUFFER = 
        "Allocation of memory for the distributed GPU communication buffer failed.";


    /*
     * EXISTING QUREG
     */

    std::string INVALID_EXISTING_ALLOC_OF_CPU_AMPS = 
        "Invalid Qureg. The CPU memory was seemingly unallocated.";
    
    std::string INVALID_EXISTING_ALLOC_OF_CPU_COMM_BUFFER = 
        "Invalid Qureg. The distributed CPU communication buffer was seemingly unallocated.";

    std::string INVALID_EXISTING_ALLOC_OF_GPU_AMPS = 
        "Invalid Qureg. The GPU memory was seemingly unallocated.";
    
    std::string INVALID_EXISTING_ALLOC_OF_GPU_COMM_BUFFER = 
        "Invalid Qureg. The distributed GPU communication buffer was seemingly unallocated.";

    std::string INVALID_EXISTING_QUREG_FIELDS = 
        "Invalid Qureg; invalid or incompatible fields isDensityMatrix=${DENS_MATR}, numQubits=${NUM_QUBITS}, numAmps=${NUM_AMPS}. It is likely this Qureg was not initialised with createQureg().";


    /*
     * MATRIX CREATION
     */

    std::string NEW_MATRIX_NUM_QUBITS_NOT_POSITIVE = 
        "Cannot create a matrix which acts upon ${NUM_QUBITS} qubits; must target one or more qubits.";

    std::string NEW_MATRIX_LOCAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a matrix which acts upon ${NUM_QUBITS} qubits since the necessary memory size (${QCOMP_BYTES} * 2^${DUB_QUBITS} bytes) would overflow size_t, and be intractably slow to serially process. The maximum size matrix targets ${MAX_QUBITS} qubits.";


    std::string NEW_MATRIX_CPU_ALLOC_FAILED = 
        "Attempted allocation of memory (${NUM_BYTES} bytes in RAM) failed.";

    std::string NEW_MATRIX_GPU_ALLOC_FAILED = 
        "Attempted allocation of GPU memory (${NUM_BYTES} bytes in VRAM) failed.";


    /*
     * EXISTING MATRIX
     */

    std::string INVALID_COMP_MATR_1_FIELDS =
        "Invalid CompMatr1. Targeted ${NUM_QUBITS} qubits (instead of 1) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 2x2). It is likely this matrix was not initialised with getCompMatr1().";

    std::string INVALID_COMP_MATR_2_FIELDS =
        "Invalid CompMatr2. Targeted ${NUM_QUBITS} qubits (instead of 2) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 4x4). It is likely this matrix was not initialised with getCompMatr2().";

    std::string INVALID_COMP_MATR_FIELDS =
        "Invalid CompMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createCompMatr().";

    std::string INVALID_COMP_MATR_CPU_MEM_ALLOC =
        "Invalid CompMatr. One or more rows of the 2D CPU memory (RAM) was seemingly unallocated. It is likely this matrix was not initialised with createCompMatr().";

    std::string INVALID_COMP_MATR_GPU_MEM_ALLOC =
        "Invalid CompMatr. The GPU memory (VRAM) was seemingly unallocated. It is likely this matrix was not initialised with createCompMatr().";


    std::string INVALID_DIAG_MATR_1_FIELDS =
        "Invalid DiagMatr1. Targeted ${NUM_QUBITS} qubits (instead of 1) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 2x2). It is likely this matrix was not initialised with getDiagMatr1().";

    std::string INVALID_DIAG_MATR_2_FIELDS =
        "Invalid DiagMatr2. Targeted ${NUM_QUBITS} qubits (instead of 2) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 4x4). It is likely this matrix was not initialised with getDiagMatr2().";

    std::string INVALID_DIAG_MATR_FIELDS =
        "Invalid DiagMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createDiagMatr().";

    std::string INVALID_DIAG_MATR_CPU_MEM_ALLOC =
        "Invalid DiagMatr. The CPU memory (RAM) was seemingly unallocated. It is likely this matrix was not initialised with createDiagMatr().";

    std::string INVALID_DIAG_MATR_GPU_MEM_ALLOC =
        "Invalid DiagMatr. The GPU memory (VRAM) was seemingly unallocated. It is likely this matrix was not initialised with createDiagMatr().";


    std::string MATRIX_NEW_ELEMS_CONTAINED_GPU_SYNC_FLAG = 
        "The new elements contained a reserved, forbidden value as the first element, used internally to detect that whether GPU memory has not synchronised. The value was intended to be extremely unlikely to be used by users - go buy a lottery ticket! If you insist on using this value in the first element, add a numerically insignificant perturbation.";

    std::string COMP_MATR_NOT_SYNCED_TO_GPU = 
        "The CompMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please first call syncCompMatr() after manually modifying elements, or overwrite all elements with setCompMatr().";

    std::string DIAG_MATR_NOT_SYNCED_TO_GPU = 
        "The DiagMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please first call syncDiagMatr() after manually modifying elements, or overwrite all elements with setDiagMatr().";


    std::string COMP_MATR_NEW_ELEMS_WRONG_NUM_ROWS =
        "Incompatible number of rows (${NUM_GIVEN_ROWS}) of elements given to overwrite a ${NUM_QUBITS}-qubit CompMatr, which expects ${NUM_EXPECTED_ROWS} rows.";

    std::string COMP_MATR_NEW_ELEMS_WRONG_ROW_DIM =
        "One or more rows contained an incompatible number of elements (${NUM_GIVEN_ELEMS}). The ${NUM_QUBITS}-qubit CompMatr expects a square ${EXPECTED_DIM}x${EXPECTED_DIM} matrix.";

    std::string DIAG_MATR_WRONG_NUM_NEW_ELEMS = 
        "Incompatible number of elements (${NUM_GIVEN_ELEMS}) assigned to a ${NUM_QUBITS}-qubit DiagMatr, which expects ${NUM_EXPECTED_ELEMS} elements.";
    

    std::string MATRIX_NOT_UNITARY = 
        "The given matrix was not (approximately) unitary.";


    /*
     * QUREG INITIALISATIONS
     */

    std::string INVALID_STATE_INDEX = 
        "Classical state index ${STATE_IND} is invalid for the given ${NUM_QUBITS} qubit Qureg. Index must be greater than or equal to zero, and cannot equal nor exceed the number of unique classical states (2^${NUM_QUBITS} = ${NUM_STATES}).";

    
    /*
     * OPERATOR PARAMETERS
     */
    
    std::string INVALID_TARGET_QUBIT = 
        "Invalid target qubit (${TARGET}). Must be greater than or equal to zero, and less than the number of qubits in the Qureg (${NUM_QUBITS}).";
}



/*
 * INVALID INPUT RESPONSE BEHAVIOUR
 */

// default C/C++ compatible error response is to simply exit in fail state
void default_invalidQuESTInputError(const char* msg, const char* func) {

    // safe to call even before MPI has been setup
    if (comm_getRank() == 0)
        std::cout 
            << "QuEST encountered a validation error during function '" << func << "':\n"
            << msg << "\nExiting..." 
            << std::endl;

    // force a synch because otherwise non-main nodes may exit before print, and MPI
    // will then attempt to instantly abort all nodes, losing the error message.
    comm_sync();

    exit(EXIT_FAILURE);
}

// enable default error response to be user-overriden as a weak symbol (even in C, and on Windows)
extern "C" {

    #ifndef _WIN32
        #pragma weak invalidQuESTInputError
        void invalidQuESTInputError(const char* msg, const char* func) {
            default_invalidQuESTInputError(msg, func);
        }
    #elif defined(_WIN64)
        #pragma comment(linker, "/alternatename:invalidQuESTInputError=default_invalidQuESTInputError")   
    #else
        #pragma comment(linker, "/alternatename:_invalidQuESTInputError=_default_invalidQuESTInputError")
    #endif

} // end C++ de-mangler



/*
 * UTILITIES
 */

// map like "${X}" -> 5, with max-size signed int values to prevent overflows.
// in C++11, these can be initialised with {{"${X}", 5}, ...}
using tokenSubs = std::map<std::string, long long int>;

std::string getStringWithSubstitutedVars(std::string oldStr, tokenSubs varsAndVals) {

    std::string newStr = oldStr;

    // substitute every var,val pair into newStr
    for (auto varAndVal : varsAndVals) {

        // unpack var -> val
        std::string var = std::get<0>(varAndVal);
        std::string val = std::to_string(std::get<1>(varAndVal));

        // assert var is well-formed 
        if (var[0] != '$' || var[1] != '{' || var.back() != '}' )
            error_validationMessageVarWasIllFormed(newStr, var);

        // assert var appears at least once in string
        if (newStr.find(var) == std::string::npos)
            error_validationMessageVarNotSubstituted(newStr, var);

        // replace every occurrence of var with val
        size_t ind = newStr.find(var);
        while (ind != std::string::npos) {
            newStr = newStr.replace(ind, var.length(), val);
            ind = newStr.find(var, ind);
        }
    }

    // assert there is no $ left in the strings
    if (newStr.find("$") != std::string::npos)
        error_validationMessageContainedUnsubstitutedVars(newStr);

    return newStr;
}

void assertThat(bool valid, std::string msg, const char* func) {

    if (!valid)
        invalidQuESTInputError(msg.c_str(), func);
}

void assertThat(bool valid, std::string msg, tokenSubs vars, const char* func) {

    std::string newMsg = getStringWithSubstitutedVars(msg, vars);
    assertThat(valid, newMsg, func);
}



/*
 * ENVIRONMENT CREATION
 */

void validate_envNeverInit(bool isQuESTInit, bool isQuESTFinal, const char* caller) {

    assertThat(!isQuESTInit, report::QUEST_ENV_ALREADY_INIT, caller);
    assertThat(!isQuESTFinal, report::QUEST_ENV_ALREADY_FINAL, caller);
}

void validate_newEnvDeploymentMode(int isDistrib, int isGpuAccel, int isMultithread, const char* caller) {

    // deployment flags must be boolean or auto
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(isDistrib     == 0 || isDistrib     == 1 || isDistrib     == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_DISTRIB,     vars, caller);
    assertThat(isGpuAccel    == 0 || isGpuAccel    == 1 || isGpuAccel    == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_GPU_ACCEL,   vars, caller);
    assertThat(isMultithread == 0 || isMultithread == 1 || isMultithread == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_MULTITHREAD, vars, caller);

    // if a deployment is forced (i.e. not auto), assert that the backend binaries are compiled correctly
    if (isDistrib == 1)
        assertThat(comm_isMpiCompiled(), report::CANNOT_DISTRIB_ENV_WITHOUT_MPI_COMMPILATION, caller);

    if (isMultithread == 1)
        assertThat(cpu_isOpenmpCompiled(), report::CANNOT_MULTITHREAD_ENV_WITHOUT_OPENMP_COMPILATION, caller);

    if (isGpuAccel == 1) {
        assertThat(gpu_isGpuCompiled(), report::CANNOT_GPU_ACCEL_ENV_WITH_GPU_COMPILATION, caller);

        // additionally require GPUs are runtime discoverable
        assertThat(gpu_isGpuAvailable(), report::CANNOT_GPU_ACCEL_ENV_WITH_NO_AVAILABLE_GPUS, caller);
    }

    // note we do not require that when (isDistrib && isGpuAccel), the MPI is CUDA-aware and that
    // the GPU is MPI-supporting. We have a runtime fallback to route through CPU memory in that situation.

    // note we do not here check whether, when distributed, a power-of-2 number of nodes are used,
    // because that requires we first initialise MPI, which we wish the caller to explicitly perform
}

void validate_newEnvDistributedBetweenPower2Nodes(const char* caller) {

    // note that we do NOT finalize MPI before erroring below, because that would necessitate
    // every node (launched by mpirun) serially print the error message, causing spam.
    // Instead, we permit the evil of every MPI process calling exit() and MPI aborting when
    // encountering the first non-zero exit code.

    int numNodes = comm_getNumNodes(); // callable even when not distributed
    tokenSubs vars = {{"${NUM_NODES}", numNodes}};

    assertThat(isPowerOf2(numNodes), report::CANNOT_DISTRIB_ENV_BETWEEN_NON_POW_2_NODES, vars, caller);
}

void validate_gpuIsCuQuantumCompatible(const char* caller) {

    assertThat(gpu_doesGpuSupportMemPools(), report::CUQUANTUM_DEPLOYED_ON_GPU_WITHOUT_MEM_POOLS, caller);

    // TODO:
    // check other requirements like compute-capability?
}



/*
 * EXISTING ENVIRONMENT
 */

void validate_envInit(const char* caller) {

    assertThat(isQuESTEnvInit(), report::QUEST_ENV_NOT_INIT, caller);
}



/*
 * QUREG CREATION
 */

void assertQuregNonEmpty(int numQubits, const char* caller) {

    assertThat(numQubits >= 1, report::NON_POSITIVE_NUM_QUBITS_IN_INIT_QUREG, {{"${NUM_QUBITS}",numQubits}}, caller);
}

void assertQuregDeployFlagsRecognised(int isDensMatr, int isDistrib, int isGpuAccel, int isMultithread, const char* caller) {

    // qureg type must be boolean
    assertThat(isDensMatr == 0 || isDensMatr == 1, report::INVALID_OPTION_FOR_QUREG_IS_DENS_MATR, caller);

    // deployment flags must be boolean or auto
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(isDistrib     == 0 || isDistrib     == 1 || isDistrib     == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_QUREG_IS_DISTRIB,     vars, caller);
    assertThat(isGpuAccel    == 0 || isGpuAccel    == 1 || isGpuAccel    == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_QUREG_IS_GPU_ACCEL,   vars, caller);
    assertThat(isMultithread == 0 || isMultithread == 1 || isMultithread == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_QUREG_IS_MULTITHREAD, vars, caller);
}

void assertQuregDeploysEnabledByEnv(int isDistrib, int isGpuAccel, int isMultithread, QuESTEnv env, const char* caller) {

    // qureg cannot deploy to backend not already enabled by the environment
    if (!env.isDistributed)
        assertThat(isDistrib     == 0 || isDistrib     == modeflag::USE_AUTO, report::NEW_DISTRIB_QUREG_IN_NON_DISTRIB_ENV, caller);
    if (!env.isGpuAccelerated)
        assertThat(isGpuAccel    == 0 || isGpuAccel    == modeflag::USE_AUTO, report::NEW_GPU_QUREG_IN_NON_GPU_ACCEL_ENV, caller);
    if (!env.isMultithreaded)
        assertThat(isMultithread == 0 || isMultithread == modeflag::USE_AUTO, report::NEW_MULTITHREAD_QUREG_IN_NON_MULTITHREAD_ENV, caller);
}

void assertQuregTotalNumAmpsDontExceedMaxIndex(int numQubits, int isDensMatr, const char* caller) {

    int maxNumQubits = mem_getMaxNumQubitsBeforeIndexOverflow(isDensMatr);

    // make message specific to statevector or density matrix
    std::string msg = (isDensMatr)? report::NEW_DENS_QUREG_AMPS_WOULD_EXCEED_QINDEX : report::NEW_QUREG_AMPS_WOULD_EXCEED_QINDEX;
    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits}, 
        {"${MAX_QUBITS}", maxNumQubits}};

    assertThat(numQubits <= maxNumQubits, msg, vars, caller);
}

void assertQuregLocalMemDoesntExceedMaxSizeof(int numQubits, int isDensMatr, int isDistrib, QuESTEnv env, const char* caller) {

    // assume distributed if possible (via auto-deployer)
    int numQuregNodes = (isDistrib == 0 || ! env.isDistributed)? 1 : env.numNodes;
    int maxNumQubits = (int) mem_getMaxNumQubitsBeforeLocalMemSizeofOverflow(isDensMatr, numQuregNodes);

    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits},
        {"${NUM_NODES}",  numQuregNodes},
        {"${MAX_QUBITS}", maxNumQubits},
        {"${IS_DENS}",    isDensMatr}};

    assertThat(numQubits <= maxNumQubits, report::NEW_QUREG_LOCAL_MEM_WOULD_EXCEED_SIZEOF, vars, caller);
}

void assertQuregNotDistributedOverTooManyNodes(int numQubits, int isDensMatr, int isDistrib, QuESTEnv env, const char* caller) {

    // only validate when distribution is forced (auto-deployer will never over-distribute)
    if (isDistrib != 1)
        return;

    // make message specific to statevector or density matrix
    std::string msg = (isDensMatr)? report::NEW_DENS_DISTRIB_QUREG_HAS_TOO_FEW_AMPS : report::NEW_DISTRIB_QUREG_HAS_TOO_FEW_AMPS;
    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits},
        {"${NUM_NODES}",  env.numNodes},
        {"${MIN_QUBITS}", std::floor(std::log2(env.numNodes))}};

    int minQubits = mem_getMinNumQubitsForDistribution(env.numNodes);
    assertThat(numQubits >= minQubits, msg, vars, caller);
}

void assertQuregFitsInCpuMem(int numQubits, int isDensMatr, int isDistrib, QuESTEnv env, const char* caller) {

    // if we can determine the total CPU memory available...
    try {
        // TODO: 
        //      we should broadcast this over nodes to find global minimum, forming consensus.
        //      would a broadcast in a try/catch be unwise?
        size_t memPerNode = mem_tryGetLocalRamCapacityInBytes();

        // check whether qureg (considering if distributed) fits between node memory(s).
        // note this sets numQuregNodes=1 only when distribution is impossible/switched-off,
        // but not when it would later be automatically disabled. that's fine; the auto-deployer
        // will never disable distribution if local RAM can't store the qureg, so we don't need to
        // validate the auto-deployed-to-non-distributed scenario. we only need to ensure that
        // auto-deploying-to-distribution is permitted by memory capacity.
        int numQuregNodes = (isDistrib == 0 || ! env.isDistributed)? 1 : env.numNodes;
        bool quregFitsInMem = mem_canQuregFitInMemory(numQubits, isDensMatr, numQuregNodes, memPerNode);

        tokenSubs vars = {
            {"${IS_DENS}",     isDensMatr},
            {"${NUM_QUBITS}",  numQubits},
            {"${QCOMP_BYTES}", sizeof(qcomp)},
            {"${EXP_BASE}",    (isDensMatr)? 4 : 2},
            {"${RAM_SIZE}",    memPerNode}};

        // make error message specific to whether qureg is local or distributed
        if (numQuregNodes == 1)
            assertThat(quregFitsInMem, report::NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CPU_MEM, vars, caller);
        else {
            vars["${NUM_NODES}"] = numQuregNodes;
            assertThat(quregFitsInMem, report::NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CPU_MEM, vars, caller);
        }
    } 

    // skip above checks if we cannot know total RAM. if we have too little RAM,
    // subsequent memory alloc validation will likely trigger anyway.
    catch(mem::COULD_NOT_QUERY_RAM e) {};
}

void assertQuregFitsInGpuMem(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, QuESTEnv env, const char* caller) {

    // validate when GPU-acceleration is possible; EVEN if it's auto! Because the auto-deployer will never
    // fall-back to painfully slow CPU if the GPU memory is filled. So if GPU-accel is possible, it must fit qureg.
    if (isGpuAccel == 0 || env.isGpuAccelerated == 0)
        return;

    // we consult the current available local GPU memory (being more strict than is possible for RAM)
    size_t localCurrGpuMem = gpu_getCurrentAvailableMemoryInBytes();

    // TODO:
    //      we should MPI reduce localCurrGpuMem to get minimum between all GPUs, in case of inhomogeneous
    //      GPU load, in order to form node consensus. However, this requires broadcasting
    //      type 'size_t' which needs to be specified as an MPI_TYPE, which is implementation
    //      specific. Don't care for boilerplate like this: stackoverflow.com/questions/40807833

    // check whether qureg (considering if distributed) fits between node GPU memory(s).
    // note this sets numQuregNodes=1 only when distribution is impossible/switched-off,
    // but not when it would later be automatically disabled. that's fine; the auto-deployer
    // will never disable distribution if local GPU memory can't store the qureg, so we don't 
    // need to validate the auto-deployed-to-non-distributed scenario. we only need to ensure 
    // that auto-deploying-to-distribution is permitted by GPU memory capacity.
    int numQuregNodes = (isDistrib == 0 || ! env.isDistributed)? 1 : env.numNodes;
    bool quregFitsInMem = mem_canQuregFitInMemory(numQubits, isDensMatr, numQuregNodes, localCurrGpuMem);

    tokenSubs vars = {
        {"${IS_DENS}",        isDensMatr},
        {"${NUM_QUBITS}",     numQubits},
        {"${QCOMP_BYTES}",    sizeof(qcomp)},
        {"${MIN_VRAM_AVAIL}", localCurrGpuMem}};

    // make error message specific to whether qureg is local or distributed
    if (numQuregNodes == 1) {
        vars["${EXP_BASE}"] = (isDensMatr)? 4 : 2;
        assertThat(quregFitsInMem, report::NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CURRENT_GPU_MEM, vars, caller);

    // when distributed, comm buffers are considered (hence +1 below)
    } else {
        vars["${LOG2_NUM_AMPS}"] = 1 + mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numQuregNodes);
        vars["${NUM_GPUS}"] = numQuregNodes;
        assertThat(quregFitsInMem, report::NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CURRENT_GPU_MEM, vars, caller);
    }
}

void validate_newQuregNotBothMultithreadedAndGpuAccel(int useGpu, int useMultithread, const char* caller) {

    // note either or both of useGpu and useMultithread are permitted to be modeflag::USE_AUTO (=-1)
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(useGpu != 1 || useMultithread != 1, report::NEW_GPU_QUREG_CANNOT_USE_MULTITHREADING, vars, caller);
}

void validate_newQuregParams(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, int isMultithread, QuESTEnv env, const char* caller) {
    assertQuregNonEmpty(numQubits, caller);
    assertQuregDeployFlagsRecognised(isDensMatr, isDistrib, isGpuAccel, isMultithread, caller);
    assertQuregDeploysEnabledByEnv(isDistrib, isGpuAccel, isMultithread, env, caller);
    assertQuregTotalNumAmpsDontExceedMaxIndex(numQubits, isDensMatr, caller);
    assertQuregLocalMemDoesntExceedMaxSizeof(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregNotDistributedOverTooManyNodes(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInCpuMem(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInGpuMem(numQubits, isDensMatr, isDistrib, isGpuAccel, env, caller);

    validate_newQuregNotBothMultithreadedAndGpuAccel(isGpuAccel, isMultithread, caller);
}

void validate_newOrExistingQuregAllocs(Qureg qureg, bool isNewQureg, const char* caller) {

    // determine if all relevant arrays are correctly allocated (in order of report precedence)...
    bool isValid = true;
    std::string errMsg = "";

    // CPU amps should always be allocated
    if (qureg.cpuAmps == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_NEW_ALLOC_OF_CPU_AMPS : 
            report::INVALID_EXISTING_ALLOC_OF_CPU_AMPS;
        isValid = false;
    }

    // CPU communication buffer is only allocated if distributed to >1 nodes
    else if (qureg.isDistributed && qureg.cpuCommBuffer == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_NEW_ALLOC_OF_CPU_COMM_BUFFER :
            report::INVALID_EXISTING_ALLOC_OF_CPU_COMM_BUFFER;
        isValid = false;
    }

    // GPU amps are only allocated in GPU mode
    else if (qureg.isGpuAccelerated && qureg.gpuAmps == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_NEW_ALLOC_OF_GPU_AMPS :
            report::INVALID_EXISTING_ALLOC_OF_GPU_AMPS;
        isValid = false;
    }

    // GPU communication buffer is only allocated in GPU mode, and when distributed to >1 nodes
    else if (qureg.isGpuAccelerated && qureg.isDistributed && qureg.gpuCommBuffer == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_NEW_ALLOC_OF_GPU_COMM_BUFFER :
            report::INVALID_EXISTING_ALLOC_OF_GPU_COMM_BUFFER;
        isValid = false;
    }

    // if the qureg was only just (attemptedly) created, and one or more arrays were not correctly allocated...
    if (isNewQureg && !isValid) {

        // then de-allocate all the successfully allocated arrays to avoid memory leak in subsequent error
        if (qureg.cpuAmps != NULL)
            cpu_deallocAmps(qureg.cpuAmps);
        if (qureg.cpuCommBuffer != NULL)
            cpu_deallocAmps(qureg.cpuCommBuffer);
        if (qureg.gpuAmps != NULL)
            gpu_deallocAmps(qureg.gpuAmps);
        if (qureg.gpuCommBuffer != NULL)
            gpu_deallocAmps(qureg.gpuCommBuffer);
    }

    // throw error or continue
    assertThat(isValid, errMsg, caller);
}



/*
 * EXISTING QUREG
 */

void validate_quregInit(Qureg qureg, const char* caller) {

    // attempt to detect the Qureg was not initialised with createQureg by the 
    // struct fields being randomised, and ergo being dimensionally incompatible
    bool valid = true;
    valid &= (qureg.numQubits > 0);
    valid &= (qureg.isDensityMatrix == 0 || qureg.isDensityMatrix == 1);
    valid &= (qureg.numAmps == powerOf2(((qureg.isDensityMatrix)? 2:1) * qureg.numQubits));

    tokenSubs vars = {
        {"${DENS_MATR}", qureg.isDensityMatrix},
        {"${NUM_QUBITS}", qureg.numQubits},
        {"${NUM_AMPS}", qureg.numAmps}};
        
    assertThat(valid, report::INVALID_EXISTING_QUREG_FIELDS, vars, caller);

    // ensure platform-specific arrays are all not-NULL. note this
    // won't catch when Qureg was un-initialised because most
    // compilers will not automatically set struct pointers to NULL, but
    // it will catch the allocation somehow failing or being overwritten
    bool isNewQureg = false;
    validate_newOrExistingQuregAllocs(qureg, isNewQureg, __func__);
}



/*
 * MATRIX CREATION
 */

void assertNewMatrixNotTooBig(int numQubits, const char* caller) {

    // assert the total memory (in bytes) to store the matrix
    // (i.e. all of its row arrays combined) does not exceed
    // max size_t, so sizeof doesn't overflow. We may never
    // actually need to compute all memory, but the threshold
    // for causing this overflow is already impractically huge,
    // and this is the 'smallest' max-size threshold we can test
    // without querying RAM and VRAM occupancy. We can't do the
    // latter ever since the CompMatr might never be dispatched
    // to the GPU.

    // the total ComplexMatrixN memory equals that of a same-size density matrix
    bool isDensMatr = true;
    int maxQubits = mem_getMaxNumQubitsBeforeIndexOverflow(isDensMatr);

    tokenSubs vars = {
        {"${NUM_QUBITS}",  numQubits},
        {"${MAX_QUBITS}",  maxQubits},
        {"${DUB_QUBITS}",   2*numQubits},
        {"${QCOMP_BYTES}", sizeof(qcomp)}};

    assertThat(numQubits < maxQubits, report::NEW_MATRIX_LOCAL_MEM_WOULD_EXCEED_SIZEOF, vars, caller);
}

void validate_newMatrixNumQubits(int numQubits, const char* caller) {

    assertThat(numQubits >= 1, report::NEW_MATRIX_NUM_QUBITS_NOT_POSITIVE, caller);
    assertNewMatrixNotTooBig(numQubits, caller);
}

void validate_newMatrixAllocs(CompMatr matr, qindex numBytes, const char* caller) {
    tokenSubs vars = {{"${NUM_BYTES}", numBytes}};

    // assert CPU array of rows was malloc'd successfully
    assertThat(matr.cpuElems != NULL, report::NEW_MATRIX_CPU_ALLOC_FAILED, vars, caller);

    // assert each CPU row was calloc'd successfully
    for (qindex r=0; r<matr.numRows; r++)
        assertThat(matr.cpuElems[r] != NULL, report::NEW_MATRIX_CPU_ALLOC_FAILED, vars, caller);
    
    // optionally assert GPU memory was malloc'd successfully
    if (getQuESTEnv().isGpuAccelerated)
        assertThat(matr.gpuElems != NULL, report::NEW_MATRIX_GPU_ALLOC_FAILED, vars, caller);
}

void validate_newMatrixAllocs(DiagMatr matr, qindex numBytes, const char* caller) {
    tokenSubs vars = {{"${NUM_BYTES}", numBytes}};

    // assert CPU array of rows was malloc'd successfully
    assertThat(matr.cpuElems != NULL, report::NEW_MATRIX_CPU_ALLOC_FAILED, vars, caller);

    // optionally assert GPU memory was malloc'd successfully
    if (getQuESTEnv().isGpuAccelerated)
        assertThat(matr.gpuElems != NULL, report::NEW_MATRIX_GPU_ALLOC_FAILED, vars, caller);
}



/*
 * EXISTING MATRICES
 */

// T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr
template <class T>
void assertMatrixFieldsAreValid(T matr, int expectedNumQb, std::string errMsg, const char* caller) {

    qindex dim = getMatrixDim(matr);
    tokenSubs vars = {
        {"${NUM_QUBITS}", matr.numQubits},
        {"${NUM_ROWS}",   dim}};

    // assert correct fixed-size numQubits (caller gaurantees this passes for dynamic-size),
    // where the error message string already contains the expected numQb
    assertThat(matr.numQubits == expectedNumQb, errMsg, vars, caller);

    qindex expectedDim = powerOf2(matr.numQubits);
    assertThat(matr.numQubits >= 1, errMsg, vars, caller);
    assertThat(dim == expectedDim,  errMsg, vars, caller);
}

// T can be CompMatr or DiagMatr (the only matrix structs with pointers)
template <class T>
void assertMatrixAllocsAreValid(T matr, std::string cpuErrMsg, std::string gpuErrMsg, const char* caller) {

    // assert CPU memory is allocated
    assertThat(matr.cpuElems != NULL, cpuErrMsg, caller);

    // we do not check that each CPU memory row (of CompMatr; irrelevant to DiagMatr)
    // is not-NULL because it's really unlikely that inner memory wasn't allocaed but
    // outer was and this wasn't somehow caught by post-creation validation (i.e. the 
    // user manually malloc'd memory after somehow getting around const fields). Further,
    // since this function is called many times (i.e. each time the user passes a matrix
    // to a simulation function like multiQubitUnitary()), it may be inefficient to 
    // serially process each row pointer.

    // optionally assert GPU memory is allocated
    if (getQuESTEnv().isGpuAccelerated)
        assertThat(matr.gpuElems != NULL, gpuErrMsg, caller);
}

void validate_matrixFields(CompMatr1 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 1, report::INVALID_COMP_MATR_1_FIELDS, caller);
}
void validate_matrixFields(CompMatr2 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 2, report::INVALID_COMP_MATR_2_FIELDS, caller);
}
void validate_matrixFields(CompMatr matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_COMP_MATR_FIELDS, caller);
    assertMatrixAllocsAreValid(matr, report::INVALID_COMP_MATR_CPU_MEM_ALLOC, report::INVALID_COMP_MATR_GPU_MEM_ALLOC, caller);
}

void validate_matrixFields(DiagMatr1 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 1, report::INVALID_DIAG_MATR_1_FIELDS, caller);
}
void validate_matrixFields(DiagMatr2 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 2, report::INVALID_DIAG_MATR_2_FIELDS, caller);
}
void validate_matrixFields(DiagMatr matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_DIAG_MATR_FIELDS, caller);
    assertMatrixAllocsAreValid(matr, report::INVALID_DIAG_MATR_CPU_MEM_ALLOC, report::INVALID_DIAG_MATR_GPU_MEM_ALLOC, caller);
}

void validate_matrixNumNewElems(int numQubits, std::vector<std::vector<qcomp>> elems, const char* caller) {

    qindex dim = powerOf2(numQubits);
    tokenSubs vars = {
        {"${NUM_QUBITS}",        numQubits},
        {"${NUM_EXPECTED_ROWS}", dim},
        {"${NUM_GIVEN_ROWS}",    elems.size()}};

    assertThat(elems.size() == dim, report::COMP_MATR_NEW_ELEMS_WRONG_NUM_ROWS, vars, caller);

    for(auto & row : elems) {

        vars = {
            {"${NUM_QUBITS}",      numQubits},
            {"${EXPECTED_DIM}",    dim},
            {"${NUM_GIVEN_ELEMS}", row.size()}};

        assertThat(row.size() == dim, report::COMP_MATR_NEW_ELEMS_WRONG_ROW_DIM, vars, caller);
    }
}
void validate_matrixNumNewElems(int numQubits, std::vector<qcomp> elems, const char* caller) {

    qindex dim = powerOf2(numQubits);
    tokenSubs vars = {
        {"${NUM_QUBITS}",        numQubits},
        {"${NUM_EXPECTED_ELEMS}", dim},
        {"${NUM_GIVEN_ELEMS}",    elems.size()}};

    assertThat(elems.size() == dim, report::DIAG_MATR_WRONG_NUM_NEW_ELEMS, vars, caller);
}

void validate_matrixNewElemsDontContainUnsyncFlag(qcomp firstElem, const char* caller) {

    // we permit the matrix to contain the GPU-mem-unsync flag in CPU-only mode,
    // to avoid astonishing a CPU-only user with a GPU-related error message
    if (!getQuESTEnv().isGpuAccelerated)
        return;

    assertThat(!gpu_doCpuAmpsHaveUnsyncMemFlag(firstElem), report::MATRIX_NEW_ELEMS_CONTAINED_GPU_SYNC_FLAG, caller);
}

// type T can be CompMatr or DiagMatr
template <class T>
void assertMatrixIsSynced(T matr, std::string errMsg, const char* caller) {

    // checks fields (including memory allocations)
    validate_matrixFields(matr, caller);

    // we don't need to perform any sync check in CPU-only mode
    if (matr.gpuElems != NULL)
        return;

    // check if GPU amps have EVER been overwritten; we sadly cannot check the LATEST changes were pushed though
    assertThat(gpu_haveGpuAmpsBeenSynced(matr.gpuElems), errMsg, caller);
}

// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr
template <class T> 
void assertMatrixIsUnitary(T matr, const char* caller) {

    // TODO: 
    //  this function always serially processes the CPU elements,
    //  but given the matrix may have persistent GPU memory, it is
    //  prudent to GPU process it instead. I propose we...
    //      - perform current serial utils check if matr is small (e.g. <=5 qubit)
    //      - perform multithread check if matr is not gpu-accel
    //      - perform gpu check if matr is gpu-accel (gpu mem gauranteed already set)

    // TODO: to make above changes, we need to remove this template; only CompMatr needs accel checks

    // assert CPU amps are unitary
    assertThat(util_isUnitary(matr), report::MATRIX_NOT_UNITARY, caller);
}

void validate_matrixIsSynced(CompMatr matr, const char* caller) {
    assertMatrixIsSynced(matr, report::COMP_MATR_NOT_SYNCED_TO_GPU, caller);
}
void validate_matrixIsSynced(DiagMatr matr, const char* caller) {
    assertMatrixIsSynced(matr, report::DIAG_MATR_NOT_SYNCED_TO_GPU, caller);
}

void validate_matrixIsUnitary(CompMatr1 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(CompMatr2 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(CompMatr matr, const char* caller) {
    validate_matrixIsSynced(matr, caller); // also checks fields
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(DiagMatr1 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(DiagMatr2 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(DiagMatr matr, const char* caller) {
    validate_matrixIsSynced(matr, caller); // also checks fields
    assertMatrixIsUnitary(matr, caller);
}


/*
 * QUREG INITIALISATIONS
 */

void validate_initClassicalStateIndex(Qureg qureg, qindex ind, const char* caller) {

    qindex maxIndExcl = powerOf2(qureg.numQubits);

    tokenSubs vars = {
        {"${STATE_IND}",  ind},
        {"${NUM_QUBITS}", qureg.numQubits},
        {"${NUM_STATES}", maxIndExcl}};

    assertThat(ind >= 0 && ind < maxIndExcl, report::INVALID_STATE_INDEX, vars, caller);
}



/*
 * OPERATOR PARAMETERS
 */

void validate_target(Qureg qureg, int target, const char* caller) {

    tokenSubs vars = {
        {"${TARGET}",  target},
        {"${NUM_QUBITS}", qureg.numQubits}};

    assertThat(target >= 0 && target < qureg.numQubits, report::INVALID_TARGET_QUBIT, vars, caller);
}