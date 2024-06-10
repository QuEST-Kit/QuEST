/** @file
 * API definitions for managing QuESTEnv instances, which
 * themselves control and query the deployment environment. 
 */

#include "quest/include/environment.h"
#include "quest/include/modes.h"

#include "quest/src/core/autodeployer.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/comm/communication.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"



/*
 * PRIVATE QUESTENV CREATION INNER FUNCTIONS
 */


QuESTEnv validateAndCreateCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread, const char* caller) {

    // ensure the chosen deployment is compiled and supported by hardware.
    // note that these error messages will be printed by every node because
    // validation occurs before comm_init() below, so all processes spawned
    // by mpirun believe they are each the main rank. This seems unavoidable.
    validate_envNotYetInit(caller);
    validate_envDeploymentMode(useDistrib, useGpuAccel, useMultithread, caller);

    // overwrite deployments left as modeflag::USE_AUTO
    autodep_chooseQuESTEnvDeployment(useDistrib, useGpuAccel, useMultithread);

    // bind deployment info to QuESTEnv (may be overwritten still below)
    QuESTEnv env;
    env.isDistributed = useDistrib;
    env.isGpuAccelerated = useGpuAccel;
    env.isMultithreaded = useMultithread;

    // assume no distribution, then revise below
    env.rank = 0;
    env.numNodes = 1;

    // initialise distribution (even for a single node), though we may subsequently error and finalize
    if (useDistrib) {
        comm_init();
        env.rank = comm_getRank();
        env.numNodes = comm_getNumNodes();
    }

    // validates numNodes=2^N and otherwise calls comm_end() before throwing error
    validate_envDistributedBetweenPower2Nodes(env.numNodes, caller);

    // TODO:
    // validate 2^N local GPUs

    // in multi-GPU settings, bind each MPI process to one GPU
    if (useDistrib && useGpuAccel)
        gpu_bindLocalGPUsToNodes(env.rank);

    // TODO: setup RNG

    return env;
}



/*
 * PUBLIC FUNCTIONS
 */


// enable invocation by both C and C++ binaries
extern "C" {


QuESTEnv createCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread) {

    return validateAndCreateCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread, __func__);
}


QuESTEnv createQuESTEnv() {

    return validateAndCreateCustomQuESTEnv(modeflag::USE_AUTO, modeflag::USE_AUTO, modeflag::USE_AUTO, __func__);
}


void destroyQuESTEnv(QuESTEnv env) {

    if (env.isDistributed)
        comm_end();
}


// end de-mangler
}
