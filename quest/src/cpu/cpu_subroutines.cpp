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



// TODO: hmm

#include "quest/src/core/accelerator.hpp"



/*
 * COMMUNICATION BUFFER PACKING
 */



// TODO:
// HMMM having inner templated functions means we avoid GPU subroutines being duplicated when
// they don't need to be (because they might not call the inner templated kernel) but...
// - it's only a FEW subroutines that ever use templating, not the whole backend!
// - that means we must duplicate the logic/call for choosing the compile-time func
// - we would move dispatch logic out of accelerator, even though compile-time choosing is an
//   act of acceleration
// - kernel invocations in GPU subroutines would then look like:
//   MACRO(
//            kernel_statevec_anyCtrlOneTargDenseMatr_subA,
//             <<<numBlocks,NUM_THREADS_PER_BLOCK>>> (amps, params, targ, m00,m01,m10,m11) );
//
//   i.e. the CUDA <<< syntax would be in the macro call
// 
// SO I think I will make accelerator perform the compile-time flags, and it can also sort
// out the control states etc





template <int NumCtrls>
void cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {
    
    // TODO:
    // this is ctrl specific right now, rather than generically about qubits, because it is
    // accepting optionally-none ctrlStates which we automatically populate to be 1. We should
    // do this automatic defaulting/population at a higher bottleneck (localiser? accelerator?), 
    // NOT here, then we can make this generic to qubits

    auto [qubits, stateMask] = getSortedQubitsAndMask(ctrls, ctrlStates, {}, {});
    qindex numInds = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    SET_VAR_AT_COMPILE_TIME(int, numCtrls, NumCtrls, ctrls.size());

    #pragma omp parallel for
    for (qindex n=0; n<numInds; n++) {

        qindex j = insertBits(n, qubits.data(), numCtrls, 0);
        qindex i = activateBits(j, stateMask);

        qureg.cpuCommBuffer[n] = qureg.cpuAmps[i];
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) );




/*
 * ANY-CTRL ONE-TARG MATRIX
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    // TODO:
    // we don't yet know if this is any faster than templated AnyTarg - we have to test.
    // if it's no faster, delete it


    // TODO: 
    //  in new templated inner function, must assert NumCtrls == ctrls.size()

    // compile-time
    constexpr int NumQubits = NumCtrls + 1;

    auto [qubits, stateMask] = getSortedQubitsAndMask(ctrls, ctrlStates, {targ}, {0});
    qindex numInds = qureg.numAmpsPerNode / powerOf2(qubits.size());

    #pragma omp parallel for
    for (qindex n=0; n<numInds; n++) {

        qindex j  = insertBits(n, qubits.data(), NumQubits, 0);
        qindex i0 = activateBits(j, stateMask);
        qindex i1 = flipBit(i0, targ);

        // note the two amplitudes are likely strided and not adjacent (separated by 2^t)
        qcomp amp0 = qureg.cpuAmps[i0];
        qcomp amp1 = qureg.cpuAmps[i1];

        qureg.cpuAmps[i0] = matr.elems[0][0]*amp0 + matr.elems[0][1]*amp1;
        qureg.cpuAmps[i1] = matr.elems[1][0]*amp0 + matr.elems[1][1]*amp1;
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, int, CompMatr1) );



template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    auto [qubits, stateMask] = getSortedQubitsAndMask(ctrls, ctrlStates, {}, {});
    qindex numIts = qureg.numAmpsPerNode / powerOf2(NumCtrls);


    #pragma omp parallel for
    for (qindex n=0; n<numIts; n++) {

        qindex j = insertBits(n, qubits.data(), NumCtrls, 0);
        qindex i = activateBits(j, stateMask);

        // the received buffer amplitudes begin at 0+numInds
        qindex l = n + numIts;
        qcomp amp = qureg.cpuCommBuffer[l];

        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*amp;
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) );





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





// template <int NumCtrls, int NumTargs>
// void cpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {

//     // TODO:
//     // explain: targ may exceed logNumAmpsPerNode


//     // TODO:
//     //      need to handle when targs exceeds prcompiled-max

//     auto [sortedCtrls, ctrlMask] = getSortedQubitsAndMask(ctrls, ctrlStates, {}, {});
//     qindex numIts = qureg.numAmpsPerNode / powerOf2(NumCtrls);

//     #pragma omp parallel for
//     for (qindex n=0; n<numIts; n++) {

//         qindex k = insertBits(n, sortedCtrls.data(), NumCtrls, 0);
//         qindex j = activateBits(k, ctrlMask);
//         qindex i = concatenateBits(qureg.rank, j, qureg.logNumAmpsPerNode);
//         qindex t = getValueOfBits(i, targs.data(), NumTargs);

//         qureg.cpuAmps[i] *= matr.cpuElems[t];
//     }
// }

// accel_INSTANTIATE_FUNC_OPTIMISED_FOR_ANY_CTRLS_AND_SOME_TARGS(
//     void, cpu_statevec_anyCtrlAnyTargDiagMatr_sub, 
//     Qureg, vector<int>, vector<int>, vector<int>, DiagMatr
// )

