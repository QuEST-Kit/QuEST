/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#include "quest/include/modes.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
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
    tokenSubs subs = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(isDistrib     == 0 || isDistrib     == 1 || isDistrib     == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_DISTRIB, subs, caller);
    assertThat(isGpuAccel    == 0 || isGpuAccel    == 1 || isGpuAccel    == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_GPU_ACCEL, subs, caller);
    assertThat(isMultithread == 0 || isMultithread == 1 || isMultithread == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_MULTITHREAD, subs, caller);

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
