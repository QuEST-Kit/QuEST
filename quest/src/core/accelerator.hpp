/** @file
 * Internal functions for choosing which accelerator backend
 * (CPU or GPU) to call to effect local simulation subroutines 
 * upon Quregs. The data needed for these subroutines must 
 * already be localised into the appropriate memory (RAM vs VRAM)
 * and location (qureg's amplitudes or buffer space), as is
 * performed by localiser.cpp. These subroutines are ergo
 * embarrassingly parallel.
 */

#ifndef ACCELERATOR_HPP
#define ACCELERATOR_HPP

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"



void statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matrix);

void statevec_oneTargetGate_subB(Qureg qureg, qcomp fac0, qcomp fac1);



#endif // ACCELERATOR_HPP