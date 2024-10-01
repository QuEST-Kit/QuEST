/** @file
 * API signatures for creating and managing Quregs.
 */

#ifndef QUREG_H
#define QUREG_H

#include "types.h"

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



typedef struct {

    // deployment configuration
    const int isMultithreaded;
    const int isGpuAccelerated;
    const int isDistributed;

    // distributed configuration
    const int rank;
    const int numNodes;
    const int logNumNodes;

    // dimension
    const int isDensityMatrix;
    const int numQubits;
    const qindex numAmps;
    const qindex logNumAmps;

    // distributed load
    const qindex numAmpsPerNode;
    const qindex logNumAmpsPerNode;
    const qindex logNumColsPerNode;

    // amplitudes in CPU and GPU memory
    qcomp* cpuAmps;
    qcomp* gpuAmps;

    // communication buffer in CPU and GPU memory
    qcomp* cpuCommBuffer;
    qcomp* gpuCommBuffer;

} Qureg;



Qureg createQureg(int numQubits);

Qureg createDensityQureg(int numQubits);

Qureg createCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread);

Qureg createCloneQureg(Qureg qureg);

void destroyQureg(Qureg qureg);

void reportQuregParams(Qureg qureg);
void reportQureg(Qureg qureg);
void reportQuregToFile(Qureg qureg, char* fn);

void syncQuregToGpu  (Qureg qureg);
void syncQuregFromGpu(Qureg qureg);

void syncSubQuregToGpu  (Qureg qureg, qindex startInd, qindex numAmps);
void syncSubQuregFromGpu(Qureg qureg, qindex startInd, qindex numAmps);

void syncSubDensityQuregToGpu  (Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);
void syncSubDensityQuregFromGpu(Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);


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

#ifdef __cplusplus

qcomp getQuregAmp(Qureg qureg, qindex index);

qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column);

#endif // __cplusplus



#endif // QUREG_H