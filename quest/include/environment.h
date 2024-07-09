/** @file
 * API signatures for managing QuESTEnv instances.
 */

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


/*
 * QuESTEnv is a struct of which there will be a single, immutable
 * main instance, statically instantiated inside environment.cpp,
 * accessible anywhere via a getter, and which is consulted for
 * determining the deployment configuration. Users can obtain a
 * local copy of this struct with getQuESTEnv().
 * 
 * QuESTEnv is not declared const, because it is impossible to 
 * runtime initialise a const global variable like a singleton.
 * A shame, but this poses no risk; the internal singleton is not 
 * obtainable nor modifiable outside of environment.cpp since its 
 * getter returns a copy. 
 */


typedef struct {

    // deployment mode
    const int isMultithreaded;
    const int isGpuAccelerated;
    const int isDistributed;

    // distributed configuration
    const int rank;
    const int numNodes;

    // TODO: RNG seeds

} QuESTEnv;



void initQuESTEnv();

void initCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread);

void finalizeQuESTEnv();

void reportQuESTEnv();

int isQuESTEnvInit();

QuESTEnv getQuESTEnv();



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // ENVIRONMENT_H