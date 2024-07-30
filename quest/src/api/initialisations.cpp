/** @file
 * API definitions for initialisaing Qureg states.
 *
 * Note many of these can be replaced with bespoke,
 * faster OpenMP and CUDA implementations although
 * that this may prove a premature optimisation
 * since initialisations are expectedly called 
 * relatively very few times in simulation. 
 * Think carefully!
 */

#include "qureg.h"

#include "../core/validation.hpp"
#include "../core/bitwise.hpp"
#include "../src/gpu/gpu_config.hpp"

#include <algorithm>
#include <cmath>


// enable invocation by both C and C++ binaries
extern "C" {


void initBlankState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // fill CPU memory with 0
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(0,0));

    // overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
}


void initZeroState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // fill CPU memory with 0
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(0,0));

    // set the first global element to 1 (valid for both statevecs and density matrices)
    if (qureg.rank == 0)
        qureg.cpuAmps[0] = 0;

    // overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
}


void initPlusState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // |+>    = sum_i 1/sqrt(2^N) |i>  where 2^N = numAmps
    // |+><+| = sum_ij 1/2^N |i><j|    where 2^N = sqrt(numAmps)
    qreal elem = 1.0 / sqrt(qureg.numAmps);
        
    // fill CPU memory with elem
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(elem,0));

    // overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
}


void initClassicalState(Qureg qureg, qindex ind) {
    validate_quregFields(qureg, __func__);
    validate_initClassicalStateIndex(qureg, ind, __func__);

    // fill CPU memory with 0
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(0,0));

    // on-diagonal |ind><ind| = ||(2^N+1)ind >>
    if (qureg.isDensityMatrix)
        ind *= (1 + powerOf2(qureg.numQubits));

    // only one node needs to write a 1 elem
    if (qureg.rank == ind / qureg.numAmpsPerNode) // floors
        qureg.cpuAmps[ind % qureg.numAmpsPerNode] = 1;

    // all nodes overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
}


} // end de-mangler