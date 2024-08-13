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

    // each control qubit halves the needed iterations
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    auto qubits = util_getSorted(ctrls);
    qindex mask = util_getBitMask(ctrls, ctrlStates);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrls, NumCtrls, ctrls.size());

    #pragma omp parallel for
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex j = insertBits(n, qubits.data(), numCtrls, 0);
        qindex i = activateBits(j, mask);

        qureg.cpuCommBuffer[n] = qureg.cpuAmps[i];
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) )



/*
 * MATRICES
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    // each control qubit halves the needed iterations, and each iteration modifies two amplitudes
    qindex numInds = qureg.numAmpsPerNode / powerOf2(ctrls.size() + 1);

    auto qubits = util_getSorted(ctrls, {targ});
    qindex mask = util_getBitMask(ctrls, ctrlStates, {targ}, {0});

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrls, NumCtrls, ctrls.size());

    #pragma omp parallel for
    for (qindex n=0; n<numInds; n++) {

        // i0 = nth local index where ctrl bits are in specified states and targ is 0
        qindex j  = insertBits(n, qubits.data(), numCtrls+1, 0);
        qindex i0 = activateBits(j, mask);
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

    // each control qubit halves the needed iterations
    qindex numIts = qureg.numAmpsPerNode / powerOf2(NumCtrls);

    auto qubits = util_getSorted(ctrls);
    qindex mask = util_getBitMask(ctrls, ctrlStates);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrls, NumCtrls, ctrls.size());

    #pragma omp parallel for
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex j = insertBits(n, qubits.data(), numCtrls, 0);
        qindex i = activateBits(j, mask);

        // l = index of nth received buffer amp
        qindex l = n + numIts;
        qcomp amp = qureg.cpuCommBuffer[l];

        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*amp;
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) )





template <int NumCtrls, int NumTargs>
void cpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {

    // each control qubit halves the needed iterations, each of which will modify 1 amplitude
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    auto sortedCtrls = util_getSorted(ctrls);
    qindex ctrlMask  = util_getBitMask(ctrls, ctrlStates);

    // use template params to compile-time unroll loops in insertBits() and getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrls, NumCtrls, ctrls.size());
    SET_VAR_AT_COMPILE_TIME(int, numTargs, NumTargs, targs.size());

    #pragma omp parallel for
    for (qindex n=0; n<numIts; n++) {

        // j = nth local index where ctrls are in the specified states
        qindex k = insertBits(n, sortedCtrls.data(), numCtrls, 0);
        qindex j = activateBits(k, ctrlMask);

        // i = global index corresponding to j
        qindex i = concatenateBits(qureg.rank, j, qureg.logNumAmpsPerNode);

        // t = value of targeted bits, which may be in the prefix substate
        qindex t = getValueOfBits(i, targs.data(), numTargs);

        qureg.cpuAmps[i] *= matr.cpuElems[t];
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr) )









// TODO:
// HMMM OH GOD!!!! THE NUMBER OF COMBINATIONS IS INSANE!!!!
// be careful not to let it blow up! Obiously we wanna keep NumTargs small anyway, because
// else the 'amp' stack overflows


// template <int NumCtrls, int NumTargs>
// void cpu_statevec_anyCtrlAnyTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {


//     // TODO:
//     // you need to make this work when targs > NumTargs somehow!!


//     // compile-time convenience variables
//     constexpr int NumQubits = NumCtrls + NumTargs;
//     constexpr qindex NumTargAmps = powerOf2(NumTargs);

//     auto [qubits, mask] = getSortedQubitsAndMask(ctrls, ctrlStates, targs, vector<int>(NumTargs,0));

//     qindex numInds = qureg.numAmpsPerNode / powerOf2(qubits.size());

//     #pragma omp parallel for
//     for (qindex n=0; n<numInds; n++) {

//         qindex k = insertBits(n, qubits.data(), NumQubits, 0);
//         qindex j = activateBits(k, mask);

//         qcomp cache[NumTargAmps];

//         // compiler might unroll this loop
//         for (qindex l=0; l<NumTargAmps; l++) {
//             qindex i = setBits(j, targs.data(), NumTargs, l);
//             cache[l] = qureg.cpuAmps[i];
//         }

//         for (qindex l=0; l<NumTargAmps; l++) {
//             qindex i = setBits(j, targs.data(), NumTargs, l);

//             qureg.cpuAmps[i] = 0;
//             for (qindex m=0; m<NumTargAmps; m++)
//                 qureg.cpuAmps[i] += matr.cpuElems[l][m] * cache[m];
//         }

//     }

// }

// accel_INSTANTIATE_FUNC_OPTIMISED_FOR_ANY_CTRLS_AND_SOME_TARGS(
//     void, cpu_statevec_anyCtrlAnyTargDenseMatr_subA, 
//     Qureg, vector<int>, vector<int>, vector<int>, CompMatr
// );








/*

template <int NumCtrls, int NumTargs>
void masks_NEW_cpu_statevec_anyCtrlAnyTargCompMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {

    // compile-time convenience variables
    constexpr int NumQubits = NumCtrls + NumTargs;
    constexpr qindex NumTargAmps = powerOf2(NumTargs);

    auto [qubits, qubitStateMask] = getSortedQubitsAndMask(ctrls, ctrlStates, targs, vector<int>(NumTargs,0));

    qindex numInds = qureg.numAmpsPerNode / powerOf2(qubits.size());


    // TODO:
    // HMMM we could precompute all 2^n targ masks (ctrl-states + value of targs)
    // then we'd just call activeBits(masks[t]) instead of calling setBits


    qindex targStateMask[NumTargAmps];
    for (qindex m=0; m<NumTargAmps; m++)
        targStateMask[m] = std::get<0>( getSortedQubitsAndMask(ctrls, ctrlS)


    #pragma omp parallel for
    for (qindex n=0; n<numInds; n++) {

        qindex k  = insertBits(n, qubits.data(), NumQubits, 0);
        qindex j = activateBits(k, qubitStateMask);

        qcomp amps[numTargAmps];

        for (qindex t=0; t<powerOf2(NumTargs); t++) {
            qindex i = setBits(j, targs.data(), NumTargs, t);
            amps[t] = qureg.cpuAmps[i]
        }


        // TODO:

    }



}
*/





