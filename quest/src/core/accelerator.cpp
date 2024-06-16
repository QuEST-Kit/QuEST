/** @file
 * Internal functions for choosing which accelerator backend
 * (CPU or GPU) to call to effect local simulation subroutines 
 * upon Quregs. The data needed for these subroutines must 
 * already be localised into the appropriate memory (RAM vs VRAM)
 * and location (qureg's amplitudes or buffer space), as is
 * performed by localiser.cpp. These subroutines are ergo
 * embarrassingly parallel.
 */

#include "quest/include/qureg.h"
#include "quest/include/structures.h"

#include "quest/src/cpu/cpu_subroutines.hpp"
#include "quest/src/gpu/gpu_subroutines.hpp"



void statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matrix) {

    if (qureg.isGpuAccelerated)
        gpu_statevec_oneTargetGate_subA(qureg, target, matrix);
    else
        cpu_statevec_oneTargetGate_subA(qureg, target, matrix);

}
void statevec_oneTargetGate_subB(Qureg qureg, qcomp fac0, qcomp fac1) {

    if (qureg.isGpuAccelerated)
        gpu_statevec_oneTargetGate_subB(qureg, fac0, fac1);
    else
        cpu_statevec_oneTargetGate_subB(qureg, fac0, fac1);
}