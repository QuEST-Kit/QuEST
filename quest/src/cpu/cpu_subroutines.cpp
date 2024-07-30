/** @file
 * CPU OpenMP-accelerated definitions of the subroutines called by
 * accelerator.cpp. 
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "../core/bitwise.hpp"



void cpu_statevec_oneTargetGate_subA(Qureg qureg, int target, CompMatr1 matrix) {

    qindex numIts = qureg.numAmpsPerNode / 2;

    qcomp m00 = matrix.elems[0][0];
    qcomp m01 = matrix.elems[0][1];
    qcomp m10 = matrix.elems[1][0];
    qcomp m11 = matrix.elems[1][1];
    
    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex j=0; j<numIts; j++) {

        qindex i0 = insertBit(j, target, 0);
        qindex i1 = flipBit(i0, target);
        
        qcomp amp0 = qureg.cpuAmps[i0];
        qcomp amp1 = qureg.cpuAmps[i1];
        
        qureg.cpuAmps[i0] = m00*amp0 + m01*amp1;
        qureg.cpuAmps[i1] = m10*amp0 + m11*amp1;
    }
}

void cpu_statevec_oneTargetGate_subB(Qureg qureg, qcomp fac0, qcomp fac1) {

    qindex numIts = qureg.numAmpsPerNode;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex i=0; i<numIts; i++)
        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*qureg.cpuCommBuffer[i];
}