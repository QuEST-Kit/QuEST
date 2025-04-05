/** @file
 * API signatures for managing the QuEST
 * execution environment.
 * 
 * @author Tyson Jones
 * @author Richard Meister (aided in design)
 * 
 * @defgroup environment Environment
 * @ingroup api
 * @brief Data structures for managing the QuEST execution environment.
 * @{
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
 */

/// @notdoced
typedef struct {

    // deployment mode
    int isMultithreaded;
    int isGpuAccelerated;
    int isDistributed;

    // distributed configuration
    int rank;
    int numNodes;

} QuESTEnv;


/// @notdoced
void initQuESTEnv();

/// @notdoced
void initCustomQuESTEnv(int useDistrib, int useGpuAccel, int useMultithread);

/// @notdoced
void finalizeQuESTEnv();

/// @notdoced
void syncQuESTEnv();

/// @notdoced
/// @nottested
void reportQuESTEnv();

/// @notdoced
int isQuESTEnvInit();

/// @notdoced
QuESTEnv getQuESTEnv();



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // ENVIRONMENT_H

/** @} */ // (end file-wide doxygen defgroup)
