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



typedef struct Qureg
{
    // deployment configuration
    int isGpuAccelerated;
    int isDistributed;
    int isMultithreaded;

    // distributed configuration
    int rank;
    int numNodes;
    int logNumNodes;

    // dimension
    int isDensityMatrix;
    int numQubits;
    qindex numAmps;
    qindex numAmpsPerNode;
    qindex logNumAmpsPerNode;

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