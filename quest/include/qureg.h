/** @file
 * API signatures for creating and managing Quregs.
 */

#ifndef QUREG_H
#define QUREG_H

#include "quest/include/types.h"

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
    const qindex numAmpsPerNode;
    const qindex logNumAmpsPerNode;

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

void destroyQureg(Qureg qureg);

void reportQureg(Qureg qureg);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // QUREG_H