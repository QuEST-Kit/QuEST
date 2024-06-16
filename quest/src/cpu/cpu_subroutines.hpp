/** @file
 * CPU signatures of the subroutines called by accelerator.cpp. 
 */

#ifndef CPU_SUBROUTINES_HPP
#define CPU_SUBROUTINES_HPP

#include "quest/include/qureg.h"
#include "quest/include/structures.h"



void cpu_statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matrix);

void cpu_statevec_oneTargetGate_subB(Qureg qureg, qcomp fac0, qcomp fac1);



#endif // CPU_SUBROUTINES_HPP