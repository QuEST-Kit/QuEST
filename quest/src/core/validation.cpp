/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#include "quest/include/modes.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/comm/communication.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <cstdlib>
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
        "The distributed GPU-accelerated Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS}, together with its commnication buffer, is too large; one or more of the ${NUM_GPUS} GPUs has insufficient available memory (only ${MIN_VRAM_AVAIL} bytes) to store its Qureg partition (${QCOMP_BYTES} * 2^${LOG2_NUM_AMPS} bytes) bytes. Consider disabling GPU-acceleration.";


    std::string FAILED_ALLOC_OF_CPU_AMPS = 
        "Allocation of memory to store the CPU amplitudes failed.";

    std::string FAILED_ALLOC_OF_GPU_AMPS = 
        "Allocation of memory to store the GPU amplitudes failed.";

    std::string FAILED_ALLOC_OF_CPU_COMM_BUFFER = 
        "Allocation of memory for the distributed CPU communication buffer failed.";

    std::string FAILED_ALLOC_OF_GPU_COMM_BUFFER = 
        "Allocation of memory for the distributed GPU communication buffer failed.";


    std::string INVALID_ALLOC_OF_CPU_AMPS = 
        "Invalid Qureg state. The CPU memory was seemingly unallocated.";
    
    std::string INVALID_ALLOC_OF_CPU_COMM_BUFFER = 
        "Invalid Qureg state. The distributed CPU communication buffer was seemingly unallocated.";

    std::string INVALID_ALLOC_OF_GPU_AMPS = 
        "Invalid Qureg state. The GPU memory was seemingly unallocated.";
    
    std::string INVALID_ALLOC_OF_GPU_COMM_BUFFER = 
        "Invalid Qureg state. The distributed GPU communication buffer was seemingly unallocated.";
}



/*
 * INVALID INPUT RESPONSE BEHAVIOUR
 */

// default C/C++ compatible error response is to simply exit in fail state
void default_invalidQuESTInputError(const char* msg, const char* func) {

    // safe to call even before MPI has been setup
    if (comm_getRank() == 0)
        std::cout 
            << std::endl
            << "QuEST encountered a validation error during function '" << func << "':\n"
            << msg << "\nExiting..." 
            << std::endl;

    // force a synch because otherwise non-main nodes may exit before print, and MPI
    // will then attempt to instantly abort all nodes, losing the error message.
    comm_synch();

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

void validate_existingEnv(QuESTEnv env, const char* caller) {

    // TOOD:
    // confirm all the mode settings are correct, etc
}

void validate_envNotYetInit(const char* caller) {

    // TODO:
    // consult a comm/ singleton?
}

void validate_envDeploymentMode(int isDistrib, int isGpuAccel, int isMultithread, const char* caller) {

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

void validate_envDistributedBetweenPower2Nodes(int numNodes, const char* caller) {

    // note that we do NOT finalize MPI before erroring below, because that would necessitate
    // every node (launched by mpirun) serially print the error message, causing spam.
    // Instead, we permit the evil of every MPI process calling exit() and MPI aborting when
    // encountering the first non-zero exit code.

    if (!isPowerOf2(numNodes))
        assertThat(false, report::CANNOT_DISTRIB_ENV_BETWEEN_NON_POW_2_NODES, {{"${NUM_NODES}",numNodes}}, caller);
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

    // cannot store more amplitudes than can be counted by the qindex type (even when distributed)
    qindex maxNumAmps = std::numeric_limits<qindex>::max();
    int maxNumQubits = std::floor(std::log2(maxNumAmps) / (qreal) ((isDensMatr)? 2 : 1));

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

    // cannot store fewer than 1 statevec amp or 1 densmatr column per node
    qindex numAmpsOrCols = powerOf2(numQubits);
    assertThat(numAmpsOrCols >= env.numNodes, msg, vars, caller);
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

void validate_quregNotBothMultithreadedAndGpuAccel(int useGpu, int useMultithread, const char* caller) {

    // note either or both of useGpu and useMultithread are permitted to be modeflag::USE_AUTO (=-1)
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(useGpu != 1 || useMultithread != 1, report::NEW_GPU_QUREG_CANNOT_USE_MULTITHREADING, vars, caller);
}

void validate_quregParams(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, int isMultithread, QuESTEnv env, const char* caller) {
    assertQuregNonEmpty(numQubits, caller);
    assertQuregDeployFlagsRecognised(isDensMatr, isDistrib, isGpuAccel, isMultithread, caller);
    assertQuregDeploysEnabledByEnv(isDistrib, isGpuAccel, isMultithread, env, caller);
    assertQuregTotalNumAmpsDontExceedMaxIndex(numQubits, isDensMatr, caller);
    assertQuregLocalMemDoesntExceedMaxSizeof(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregNotDistributedOverTooManyNodes(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInCpuMem(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInGpuMem(numQubits, isDensMatr, isDistrib, isGpuAccel, env, caller);

    validate_quregNotBothMultithreadedAndGpuAccel(isGpuAccel, isMultithread, caller);
}

void validate_quregAllocs(Qureg qureg, bool isNewQureg, const char* caller) {

    // determine if all relevant arrays are correctly allocated (in order of report precedence)...
    bool isValid = true;
    std::string errMsg = "";

    // CPU amps should always be allocated
    if (qureg.cpuAmps == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_ALLOC_OF_CPU_AMPS : 
            report::INVALID_ALLOC_OF_CPU_AMPS;
        isValid = false;
    }

    // CPU communication buffer is only allocated if distributed to >1 nodes
    else if (qureg.isDistributed && qureg.cpuCommBuffer == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_ALLOC_OF_CPU_COMM_BUFFER :
            report::INVALID_ALLOC_OF_CPU_COMM_BUFFER;
        isValid = false;
    }

    // GPU amps are only allocated in GPU mode
    else if (qureg.isGpuAccelerated && qureg.gpuAmps == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_ALLOC_OF_GPU_AMPS :
            report::INVALID_ALLOC_OF_GPU_AMPS;
        isValid = false;
    }

    // GPU communication buffer is only allocated in GPU mode, and when distributed to >1 nodes
    else if (qureg.isGpuAccelerated && qureg.isDistributed && qureg.gpuCommBuffer == NULL) {
        errMsg = (isNewQureg)?
            report::FAILED_ALLOC_OF_GPU_COMM_BUFFER :
            report::INVALID_ALLOC_OF_GPU_COMM_BUFFER;
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
