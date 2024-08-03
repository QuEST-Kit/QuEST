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
#include "quest/src/core/errors.hpp"

#include <vector>

using std::vector;
using namespace index_flags;



/*
 * ANY-CTRL ONE-TARG MATRIX TEMPLATES
 */


template <CtrlFlag ctrlFlag>
void cpu_statevector_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr)
{
    // note ctrls is modified and pointed to by params, so much persist for function call
    CtrlTargIndParams params = getParamsInformingIndsWhereCtrlsAreActiveAndTargIsOne(ctrls, ctrlStates, targ);

    // each iteration will modify 2 amplitudes (targ=0,1), and each ctrl halves needed modifications
    qindex numIts = qureg.numAmpsPerNode / powerOf2(1 + ctrls.size());

    #pragma omp parallel for
    for (qindex n=0; n<numIts; n++) {

        qindex i1 = getNthIndWhereCtrlsAreActiveAndTargIsOne<ctrlFlag>(n, targ, params);
        qindex i0 = flipBit(i1, targ);

        // modify both amplitudes, even though they are likely strided and not adjacent
        qcomp amp0 = qureg.cpuAmps[i0];
        qcomp amp1 = qureg.cpuAmps[i1];

        qureg.cpuAmps[i0] = matr.elems[0][0]*amp0 + matr.elems[0][1]*amp1;
        qureg.cpuAmps[i1] = matr.elems[1][0]*amp0 + matr.elems[1][1]*amp1;
    }
}


template <CtrlFlag ctrlFlag>
void cpu_statevector_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr)
{
    // note ctrls is modified and pointed to by params, so much persist for function call
    CtrlIndParams params = getParamsInformingIndsWhereCtrlsAreActive(ctrls, ctrlStates);

    // each iteration will modify 1 amplitude, and each ctrl halves needed modifications
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    #pragma omp parallel for
    for (qindex n=0; n<numIts; n++) {

        qindex i = getNthIndWhereCtrlsAreActive<ctrlFlag>(n, params);
        int bit = getBit(i, targ);

        qureg.cpuAmps[i] *= matr.elems[bit];
    }
}



/*
 *  INSTANTIATING TEMPLATES
 */


INSTANTIATE_TEMPLATED_VOID_FUNC_WITH_ALL_CTRL_FLAGS(
    cpu_statevector_anyCtrlOneTargMatrix_subA, 
    (Qureg, vector<int>, vector<int>, int, CompMatr1))

INSTANTIATE_TEMPLATED_VOID_FUNC_WITH_ALL_CTRL_FLAGS(
    cpu_statevector_anyCtrlOneTargMatrix_subA, 
    (Qureg, vector<int>, vector<int>, int, DiagMatr1))
