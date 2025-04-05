/** @file
 * API signatures for creating and managing Quregs.
 * 
 * @author Tyson Jones
 * 
 * @defgroup qureg Qureg
 * @ingroup api
 * @brief Data structures for representing quantum states.
 * @{
 */

#ifndef QUREG_H
#define QUREG_H

#include "quest/include/types.h"



/*
 * These signatures are divided into three partitions; those which are
 * natively C and C++ compatible (first partition), then those which are
 * only exposed to C++ (second partition) because they return 'qcomp' 
 * which cannot cross the C++-to-C ABI, then C++only convenience functions.
 * The first partition defines the doc groups, and the latter partition 
 * functions are added into them.
 */



/*
 * C and C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/** 
 * @defgroup qureg_structs Structs
 * @brief Data structures for representing quantum registers.
 * @{
 */


/// @notdoced
typedef struct {

    // deployment configuration
    int isMultithreaded;
    int isGpuAccelerated;
    int isDistributed;

    // distributed configuration
    int rank;
    int numNodes;
    int logNumNodes;

    // dimension
    int isDensityMatrix;
    int numQubits;
    qindex numAmps;
    qindex logNumAmps;

    // distributed load
    qindex numAmpsPerNode;
    qindex logNumAmpsPerNode;
    qindex logNumColsPerNode;

    // amplitudes in CPU and GPU memory
    qcomp* cpuAmps;
    qcomp* gpuAmps;

    // communication buffer in CPU and GPU memory
    qcomp* cpuCommBuffer;
    qcomp* gpuCommBuffer;

} Qureg;

/** @} */



/** 
 * @defgroup qureg_create Constructors
 * @brief Functions for creating statevectors and density matrices.
 * @{
 */


/// @notdoced
Qureg createQureg(int numQubits);


/// @notdoced
Qureg createDensityQureg(int numQubits);


/// @notdoced
Qureg createForcedQureg(int numQubits);


/// @notdoced
Qureg createForcedDensityQureg(int numQubits);


/// @notdoced
Qureg createCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread);


/// @notdoced
Qureg createCloneQureg(Qureg qureg);


/** @} */



/** 
 * @defgroup qureg_destroy Destructors
 * @brief Functions for destroying existing Qureg.
 * @{
 */


/// @notdoced
void destroyQureg(Qureg qureg);


/** @} */



/** 
 * @defgroup qureg_report Reporters
 * @brief Functions for printing Qureg states or reporting their configuration.
 * @{
 */


/// @notdoced
/// @nottested
void reportQuregParams(Qureg qureg);


/// @notdoced
/// @nottested
void reportQureg(Qureg qureg);


/** @} */



/** 
 * @defgroup qureg_sync Synchronisation
 * @brief Functions for copying memory between a Qureg's CPU (RAM) and GPU (VRAM) memory. 
 * @details These functions are only necessary when the user wishes to manually probe or
 *          modify the Qureg amplitudes (rather than use functions like getQuregAmps() and
 *          setQuregAmps()), to ensure that the CPU and GPU copies of the Qureg are identical.
 *          These functions have no effect when running without GPU-acceleration, but remain 
 *          legal and harmless to call, to achieve platform agnosticism.
 * @{
 */


/// @notdoced
/// @nottested
void syncQuregToGpu(Qureg qureg);


/// @notdoced
/// @nottested
void syncQuregFromGpu(Qureg qureg);


/// @notdoced
/// @nottested
void syncSubQuregToGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps);


/// @notdoced
/// @nottested
void syncSubQuregFromGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps);


/** @} */



/** 
 * @defgroup qureg_get Getters
 * @brief Functions for obtaining amplitudes from statevectors or density matrices.
 * @{
 */


/// @notdoced
void getQuregAmps(qcomp* outAmps, Qureg qureg, qindex startInd, qindex numAmps);


/// @notdoced
void getDensityQuregAmps(qcomp** outAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);


/** @} */


// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ ONLY FUNCTIONS
 *
 * which are not directly C-compatible because they pass or
 * return qcomp primitives by-value (rather than by pointer).
 * This is prohibited because the C and C++ ABI does not agree
 * on a complex type, though C's _Complex has the same memory
 * layout as C++'s std::complex<>. To work around this, the 
 * below functions have a C-compatible wrapper defined in
 * wrappers.h which passes/receives the primitives by pointer;
 * a qcomp ptr can be safely passed from the C++ source binary
 * the user's C binary. These functions use the existing doxygen
 * doc groups defined above
 */


/// @ingroup qureg_get
/// @notdoced
qcomp getQuregAmp(Qureg qureg, qindex index);


/// @ingroup qureg_get
/// @notdoced
qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column);



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

#include <vector>


/// @ingroup qureg_get
/// @nottested
/// @notvalidated
/// @notdoced
/// @cpponly
std::vector<qcomp> getQuregAmps(Qureg qureg, qindex startInd, qindex numAmps);


/// @ingroup qureg_get
/// @nottested
/// @notvalidated
/// @notdoced
/// @cpponly
std::vector<std::vector<qcomp>> getDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);


#endif // __cplusplus



/*
 * ORPHAN DOC GROUP
 */

/** 
 * @defgroup qureg_setters Setters
 * @brief See @ref init_amps "Amplitude initialisations".
 */



#endif // QUREG_H

/** @} */ // (end file-wide doxygen defgroup)
