/** @file
 * CUDA GPU-accelerated signatures of the subroutines called by accelerator.cpp.
 */

#ifndef GPU_SUBROUTINES_HPP
#define GPU_SUBROUTINES_HPP

#include "qureg.h"
#include "structures.h"



void gpu_statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matrix);

void gpu_statevec_oneTargetGate_subB(Qureg qureg, qcomp fac0, qcomp fac1);



#endif // GPU_SUBROUTINES_HPP