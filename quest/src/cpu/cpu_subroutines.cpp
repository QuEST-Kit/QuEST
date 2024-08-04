/** @file
 * CPU OpenMP-accelerated definitions of the subroutines called by
 * accelerator.cpp. Some of these definitions are templated, defining
 * multiple versions optimised (at compile-time) for handling different
 * numbers of control qubits.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/indexer.hpp"

#include <vector>

using std::vector;
using namespace index_flags;



/*
 * COMMUNICATION BUFFER PACKING
 */


template <CtrlFlag ctrlFlag>
qindex cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    CtrlIndParams params = getParamsInformingIndsWhereCtrlsAreActive(qureg.numAmpsPerNode, ctrls, ctrlStates);

    #pragma omp parallel for
    for (qindex n=0; n<params.numInds; n++) {

        // pack the first segment of the buffer with all active amps
        qindex i = getNthIndWhereCtrlsAreActive<ctrlFlag>(n, params);
        qureg.cpuCommBuffer[n] = qureg.cpuAmps[i];
    }

    // return the total number of packed amps to caller
    return params.numInds;
}



/*
 * ANY-CTRL ONE-TARG MATRIX
 *
 * which are templated and require explicit instantiation below
 */


template <CtrlFlag ctrlFlag>
void cpu_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    CtrlTargIndParams params = getParamsInformingIndsWhereCtrlsAreActiveAndTargIsOne(qureg.numAmpsPerNode, ctrls, ctrlStates, targ);

    #pragma omp parallel for
    for (qindex n=0; n<params.numInds; n++) {

        // each iteration locates and modifies two amplitudes
        qindex i1 = getNthIndWhereCtrlsAreActiveAndTargIsOne<ctrlFlag>(n, params);
        qindex i0 = flipBit(i1, targ);

        // note the two amplitudes are likely strided and not adjacent (separated by 2^t)
        qcomp amp0 = qureg.cpuAmps[i0];
        qcomp amp1 = qureg.cpuAmps[i1];

        qureg.cpuAmps[i0] = matr.elems[0][0]*amp0 + matr.elems[0][1]*amp1;
        qureg.cpuAmps[i1] = matr.elems[1][0]*amp0 + matr.elems[1][1]*amp1;
    }
}


template <CtrlFlag ctrlFlag>
void cpu_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {

    CtrlIndParams params = getParamsInformingIndsWhereCtrlsAreActive(qureg.numAmpsPerNode, ctrls, ctrlStates);

    #pragma omp parallel for
    for (qindex n=0; n<params.numInds; n++) {

        // each iteration scales one amplitude
        qindex i = getNthIndWhereCtrlsAreActive<ctrlFlag>(n, params);
        int bit = getBit(i, targ);
        
        qureg.cpuAmps[i] *= matr.elems[bit];
    }
}


template <CtrlFlag ctrlFlag>
void cpu_statevec_anyCtrlOneTargDenseMatrix_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    CtrlIndParams params = getParamsInformingIndsWhereCtrlsAreActive(qureg.numAmpsPerNode, ctrls, ctrlStates);

    #pragma omp parallel for
    for (qindex n=0; n<params.numInds; n++) {

        // each iteration modifies one local amp
        qindex i = getNthIndWhereCtrlsAreActive<ctrlFlag>(n, params);

        // the received buffer amplitudes begin at 0+numInds
        qindex l = n + params.numInds;
        qcomp amp = qureg.cpuCommBuffer[l];

        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*amp;
    }
}


INSTANTIATE_TEMPLATED_FUNC_WITH_ALL_CTRL_FLAGS(
    qindex, cpu_statevec_packAmpsIntoBuffer, 
    (Qureg qureg, vector<int> ctrls, vector<int> ctrlStates))

INSTANTIATE_TEMPLATED_FUNC_WITH_ALL_CTRL_FLAGS(
    void, cpu_statevec_anyCtrlOneTargMatrix_subA, 
    (Qureg, vector<int>, vector<int>, int, CompMatr1))

INSTANTIATE_TEMPLATED_FUNC_WITH_ALL_CTRL_FLAGS(
    void, cpu_statevec_anyCtrlOneTargMatrix_subA, 
    (Qureg, vector<int>, vector<int>, int, DiagMatr1))

INSTANTIATE_TEMPLATED_FUNC_WITH_ALL_CTRL_FLAGS(
    void, cpu_statevec_anyCtrlOneTargDenseMatrix_subB, 
    (Qureg, vector<int>, vector<int>, qcomp, qcomp))
