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

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



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



/// @notdoced
void destroyQureg(Qureg qureg);



/// @notdoced
/// @nottested
void reportQuregParams(Qureg qureg);

/// @notdoced
/// @nottested
void reportQureg(Qureg qureg);



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



/// @notdoced
void getQuregAmps(qcomp* outAmps, Qureg qureg, qindex startInd, qindex numAmps);

/// @notdoced
void getDensityQuregAmps(qcomp** outAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);


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
 * the user's C binary.
 */

/// @notdoced
qcomp getQuregAmp(Qureg qureg, qindex index);

/// @notdoced
qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column);



#endif // QUREG_H

/** @} (end doxygen defgroup) */
