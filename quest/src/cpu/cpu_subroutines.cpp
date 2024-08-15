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

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    
    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex j = insertBits(n, sortedCtrls.data(), numCtrlBits, 0);
        qindex i = activateBits(j, ctrlStateMask);

        qureg.cpuCommBuffer[n] = qureg.cpuAmps[i];
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
        qindex j  = insertBits(n, sortedQubits.data(), numQubitBits, 0);
        qindex i0 = activateBits(j, qubitStateMask);
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

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex j = insertBits(n, sortedCtrls.data(), numCtrlBits, 0);
        qindex i = activateBits(j, ctrlStateMask);

        // l = index of nth received buffer amp
        qindex l = n + numIts;
        qcomp amp = qureg.cpuCommBuffer[l];

        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*amp;
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) )



/*
 * MANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs>
void cpu_statevec_anyCtrlAnyTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {
    
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

    // we tested a variant of this function where a mask for each ctrl-targ state is calculated
    // upfront (of which there are numTargAmps many), replacing all setBits() calls with
    // activateBits(), which lacks the runtime loop and does not need compile-time unrolling.
    // curiously, that worsened performance in all regimes!

    // each control qubit halves iterations, each of which modifies 2^(targs.size()) amplitudes
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size() + targs.size());

    // prepare a mask which yields ctrls in specified state, and targs in all-zero
    auto sortedQubits   = util_getSorted(ctrls, targs);
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, targs, vector<int>(targs.size(),0));

    // attempt to use compile-time variables to automatically optimise/unroll dependent loops
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, targs.size());

    // compiler will infer these at compile-time if possible
    int numQubitBits = numCtrlBits + numTargBits;
    qindex numTargAmps = powerOf2(numTargBits);

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // a private cache needed by each thread to update each iteration's amplitudes (may be compile-time sized)
        vector<qcomp> cache(numTargAmps);

        // i0 = nth local index where ctrls are active and targs are all zero
        qindex j  = insertBits(n, sortedQubits.data(), numQubitBits, 0); // loop may be unrolled
        qindex i0 = activateBits(j, qubitStateMask);

        // collect and cache all to-be-modified amps (loop might be unrolled)
        for (qindex k=0; k<numTargAmps; k++) {

            // i = nth local index where ctrls are active and targs form value k
            qindex i = setBits(i0, targs.data(), numTargBits, k); // loop may be unrolled
            cache[k] = qureg.cpuAmps[i];
        }

        // modify each amplitude (loop might be unrolled)
        for (qindex l=0; l<numTargAmps; l++) {

            // i = nth local index where ctrls are active and targs form value l
            qindex i = setBits(i0, targs.data(), numTargBits, l); // loop may be unrolled
            qureg.cpuAmps[i] = 0;
        
            // loop may be unrolled
            for (qindex k=0; k<numTargAmps; k++)
                qureg.cpuAmps[i] += matr.cpuElems[l][k] * cache[k];
        }
    }
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
        qindex k = insertBits(n, sortedCtrls.data(), numCtrlBits, 0);
        qindex j = activateBits(k, ctrlStateMask);

        // i = global index corresponding to j
        qindex i = concatenateBits(qureg.rank, j, qureg.logNumAmpsPerNode);

        // t = value of targeted bits, which may be in the prefix substate
        qindex t = getValueOfBits(i, targs.data(), numTargBits);

        qureg.cpuAmps[i] *= matr.cpuElems[t];
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr) )
