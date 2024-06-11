/** @file
 * API signatures for managing QuESTEnv instances.
 */

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



typedef struct QuESTEnv {

    // deployment mode
    int isGpuAccelerated;
    int isDistributed;
    int isMultithreaded;

    // distributed configuration
    int rank;
    int numNodes;

    // TODO: RNG seeds

} QuESTEnv;



QuESTEnv createQuESTEnv();

QuESTEnv createCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread);

void destroyQuESTEnv(QuESTEnv env);

void reportQuESTEnv(QuESTEnv env);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // ENVIRONMENT_H