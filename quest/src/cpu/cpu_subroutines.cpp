/** @file
 * CPU OpenMP-accelerated definitions of the subroutines called by
 * accelerator.cpp. Some of these definitions are templated, defining
 * multiple versions optimised (at compile-time) for handling different
 * numbers of control qubits; such functions are proceeded by macro
 * INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(), to force the compilation
 * of their needed versions within this translation unit for later linkage.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/accelerator.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <vector>

using std::vector;



/*
 * COMMUNICATION BUFFER PACKING
 */


template <int NumCtrls>
void cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    // amplitudes are packed at an offset into the buffer
    qindex offset = getSubBufferSendInd(qureg);

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    
    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        qureg.cpuCommBuffer[offset + n] = qureg.cpuAmps[i];
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) )



/*
 * ONE-TARGET DENSE MATRIX
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations, and each iteration modifies two amplitudes
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size() + 1);

    auto sortedQubits   = util_getSorted(ctrls, {targ});
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ}, {0});

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    int numQubitBits = numCtrlBits + 1;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i0 = nth local index where ctrl bits are in specified states and targ is 0
        qindex i0 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);
        qindex i1 = flipBit(i0, targ);

        // note the two amplitudes are likely strided and not adjacent (separated by 2^t)
        qcomp amp0 = qureg.cpuAmps[i0];
        qcomp amp1 = qureg.cpuAmps[i1];

        qureg.cpuAmps[i0] = matr.elems[0][0]*amp0 + matr.elems[0][1]*amp1;
        qureg.cpuAmps[i1] = matr.elems[1][0]*amp0 + matr.elems[1][1]*amp1;
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, int, CompMatr1) )


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    
    qindex offset = getSubBufferRecvInd();

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*qureg.cpuCommBuffer[offset+n];
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) )



/*
 * MANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs>
void cpu_statevec_anyCtrlAnyTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {
    
    // TODO
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, vector<int>, CompMatr) )



/*
 * DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs>
void cpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {
    
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

    // each control qubit halves the needed iterations, each of which will modify 1 amplitude
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    // use template params to compile-time unroll loops in insertBits() and getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, targs.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // j = nth local index where ctrls are active (in the specified states)
        qindex j = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // i = global index corresponding to j
        qindex i = concatenateBits(qureg.rank, j, qureg.logNumAmpsPerNode);

        // t = value of targeted bits, which may be in the prefix substate
        qindex t = getValueOfBits(i, targs.data(), numTargBits);

        qureg.cpuAmps[i] *= matr.cpuElems[t];
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr) )
