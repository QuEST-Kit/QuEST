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
#include "quest/src/comm/comm_indices.hpp"
#include "quest/src/cpu/cpu_subroutines.hpp"

#include <vector>
#include <algorithm>

using std::vector;



/*
 * COMMUNICATION BUFFER PACKING
 */


template <int NumQubits>
qindex cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubitInds, vector<int> qubitStates) {

    assert_numQubitsMatchesQubitStatesAndTemplateParam(qubitInds.size(), qubitStates.size(), NumQubits);

    // each control qubit halves the needed iterations
    qindex numIts = qureg.numAmpsPerNode / powerOf2(qubitInds.size());

    // amplitudes are packed at an offset into the buffer
    qindex offset = getSubBufferSendInd(qureg);

    auto sortedQubitInds = util_getSorted(qubitInds);
    auto qubitStateMask  = util_getBitMask(qubitInds, qubitStates);
    
    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, qubitInds.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where qubits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedQubitInds.data(), numBits, qubitStateMask);

        // pack the potentially-strided amplitudes into a contiguous sub-buffer
        qureg.cpuCommBuffer[offset + n] = qureg.cpuAmps[i];
    }

    // return the number of packed amps
    return numIts;
}


qindex cpu_statevec_packPairSummedAmpsIntoBuffer(Qureg qureg, int qubit1, int qubit2, int qubit3, int bit2) {
    
    assert_bufferPackerGivenIncreasingQubits(qubit1, qubit2, qubit3);

    // pack eighth of buffer with pre-summed amp pairs
    qindex numIts = qureg.numAmpsPerNode / 8;

    // amplitudes are packed at an offset into the buffer
    qindex offset = getSubBufferSendInd(qureg);

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i000 = nth local index where all qubits are 0
        qindex i000 = insertThreeZeroBits(n, qubit3, qubit2, qubit1);
        qindex i0b0 = setBit(i000, qubit2, bit2);
        qindex i1b1 = flipTwoBits(i0b0, qubit3, qubit1);

        qureg.cpuCommBuffer[offset + n] = qureg.cpuAmps[i0b0] + qureg.cpuAmps[i1b1];
    }

    // return the number of packed amps
    return numIts;
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( qindex, cpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) )



/* 
 * SWAPS
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the number of iterations, each of which modifies 2 amplitudes, and skips 2
    qindex numIts = qureg.numAmpsPerNode / powerOf2(2 + ctrls.size());

    auto sortedQubits   = util_getSorted(ctrls, {targ2, targ1});
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ2, targ1}, {0, 1});
    
    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    int numQubitBits = numCtrlBits + 2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i01 = nth local index where ctrls are active, targ2=0 and targ1=1
        qindex i01 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);
        qindex i10 = flipTwoBits(i01, targ2, targ1);

        std::swap(qureg.cpuAmps[i01], qureg.cpuAmps[i10]);
    }
}


template <int NumCtrls>
void cpu_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the number of received amplitudes
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    
    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrls are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // j = index of nth received amplitude from pair rank in buffer
        qindex j = n + offset;

        // unpack the continuous sub-buffer among the strided local amplitudes
        qureg.cpuAmps[i] = qureg.cpuCommBuffer[j];
    }
}


template <int NumCtrls>
void cpu_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the number of iterations, each of which modifies one of the two target qubit states
    qindex numIts = qureg.numAmpsPerNode / powerOf2(1 + ctrls.size());

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    auto sortedQubits   = util_getSorted(ctrls, {targ});
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ}, {targState});
    
    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    int numQubitBits = numCtrlBits + 1;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrls and targ are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);
    
        // j = index of nth received amplitude from pair rank in buffer
        qindex j = n + offset;

        // unpack the continuous sub-buffer among the strided local amplitudes
        qureg.cpuAmps[i] = qureg.cpuCommBuffer[j];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlSwap_subA, (Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlSwap_subB, (Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlSwap_subC, (Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState) )



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


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations, and each iteration modifies one amplitude
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    
    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // j = index of nth received amplitude from pair rank in buffer
        qindex j = n + offset;

        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*qureg.cpuCommBuffer[j];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, int, CompMatr1) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) )



/*
 * TWO-TARGET DENSE MATRIX
 */


template <int NumCtrls> 
void cpu_statevec_anyCtrlTwoTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, CompMatr2 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations, and each iteration modifies four amplitudes
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size() + 2);

    auto sortedQubits   = util_getSorted(ctrls, {targ1, targ2});
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, {targ1, targ2}, {0, 0});

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    int numQubitBits = numCtrlBits + 2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i0 = nth local index where ctrl bits are in specified states and both targs are 0
        qindex i00 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);
        qindex i01 = flipBit(i00, targ1);
        qindex i10 = flipBit(i00, targ2);
        qindex i11 = flipBit(i01, targ2);

        // note amplitudes are not necessarily adjacent, nor uniformly spaced
        qcomp amp00 = qureg.cpuAmps[i00];
        qcomp amp01 = qureg.cpuAmps[i01];
        qcomp amp10 = qureg.cpuAmps[i10];
        qcomp amp11 = qureg.cpuAmps[i11];

        // amps[i_n] = sum_j elems[n][j] amp[i_n]
        qureg.cpuAmps[i00] = matr.elems[0][0]*amp00 + matr.elems[0][1]*amp01 + matr.elems[0][2]*amp10 + matr.elems[0][3]*amp11;
        qureg.cpuAmps[i01] = matr.elems[1][0]*amp00 + matr.elems[1][1]*amp01 + matr.elems[1][2]*amp10 + matr.elems[1][3]*amp11;
        qureg.cpuAmps[i10] = matr.elems[2][0]*amp00 + matr.elems[2][1]*amp01 + matr.elems[2][2]*amp10 + matr.elems[2][3]*amp11;
        qureg.cpuAmps[i11] = matr.elems[3][0]*amp00 + matr.elems[3][1]*amp01 + matr.elems[3][2]*amp10 + matr.elems[3][3]*amp11;
    }
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlTwoTargDenseMatr_sub, (Qureg, vector<int>, vector<int>, int, int, CompMatr2) )



/*
 * MANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj>
void cpu_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {
    
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

    // TODO:
    // this function allocates powerOf2(targs.size())-sized caches for each thread, sometimes in
    // heap. At the ~max non-distributed double CompMatr of 16 qubits = 64 GiB, this is 1 MiB 
    // per thread; for a conceivable 100 thread execution, this is 100 MiB being alloc/dealloced
    // at every call. It is debatable whether this justifies pre-allocating persistent cache space
    // (one for each thread, to avoid false sharing), similar to GPU's AnyTargDenseMatr, though
    // for an order of magnitude fewer threads, and using non-coalesced memory. Certainly making
    // persistent heap caches is inadvisable when the cache fits in the stack (currently automated 
    // using std::vector). Perhaps we should keep the current re-allocs, constrain that this 
    // function is only called for few-targets (e.g. <= qureg.numQubits - 5), and define another
    // function for almost-all target matrices which uses persistent heap memory, wherein the 
    // optimal parallelisation scheme is anyway different.

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

    // create an explicit parallel region to avoid re-initialisation of vectors every iteration
    #pragma omp parallel if(qureg.isMultithreaded)
    {
        // create a private cache for every thread (might be compile-time sized, and in heap or stack)
        vector<qcomp> cache(numTargAmps);

        #pragma omp for
        for (qindex n=0; n<numIts; n++) {

            // i0 = nth local index where ctrls are active and targs are all zero
            qindex i0 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);

            // collect and cache all to-be-modified amps (loop might be unrolled)
            for (qindex j=0; j<numTargAmps; j++) {

                // i = nth local index where ctrls are active and targs form value j
                qindex i = setBits(i0, targs.data(), numTargBits, j); // loop may be unrolled
                cache[j] = qureg.cpuAmps[i];
            }

            // modify each amplitude (loop might be unrolled)
            for (qindex k=0; k<numTargAmps; k++) {

                // i = nth local index where ctrls are active and targs form value k
                qindex i = setBits(i0, targs.data(), numTargBits, k); // loop may be unrolled
                qureg.cpuAmps[i] = 0;
            
                // loop may be unrolled
                for (qindex j=0; j<numTargAmps; j++) {

                    // optionally conjugate matrix elems on the fly to avoid pre-modifying heap structure
                    qcomp elem = matr.cpuElems[k][j];
                    if constexpr (ApplyConj)
                        elem = conj(elem);

                    qureg.cpuAmps[i] += elem * cache[j];
                }
            }
        }
    }
}


INSTANTIATE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDenseMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, CompMatr) )



/*
 * ONE-TARG DIAGONAL MATRIX
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations, each of which will modify 1 amplitude
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    // use template params to compile-time unroll loops in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // j = nth local index where ctrls are active (in the specified states)
        qindex j = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // i = global index corresponding to j
        qindex i = concatenateBits(qureg.rank, j, qureg.logNumAmpsPerNode);

        int b = getBit(i, targ);
        qureg.cpuAmps[i] *= matr.elems[b];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, int, DiagMatr1) )



/*
 * TWO-TARG DIAGONAL MATRIX
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlTwoTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, DiagMatr2 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations, each of which will modify 1 amplitude
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    // use template params to compile-time unroll loops in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // j = nth local index where ctrls are active (in the specified states)
        qindex j = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // i = global index corresponding to j
        qindex i = concatenateBits(qureg.rank, j, qureg.logNumAmpsPerNode);

        int k = getTwoBits(i, targ2, targ1);
        qureg.cpuAmps[i] *= matr.elems[k];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlTwoTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, int, int, DiagMatr2) )



/*
 * ANY-TARG DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj>
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

        // optionally conjugate matrix elems on the fly to avoid pre-modifying heap structure
        qcomp elem = matr.cpuElems[t];
        if constexpr (ApplyConj)
            elem = conj(elem);

        qureg.cpuAmps[i] *= elem;
    }
}


INSTANTIATE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr) )



/*
 * PAULI TENSOR AND GADGET
 */


template <int NumCtrls, int NumTargs>
void cpu_statevector_anyCtrlPauliTensorOrGadget_subA(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> suffixTargsXY, 
    qindex suffixMaskXY, qindex allMaskYZ, qcomp powI, qcomp thisAmpFac, qcomp otherAmpFac
) {
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(suffixTargsXY.size(), NumTargs);
    
    // TODO:
    //  should we attempt to OpenMP parallelise the inner loop when there are many paulis?
    //  can we achieve this using something like if(numOuterIts < nthreads) ?

    // each outer iteration handles all assignments of the target qubits and each ctrl halves the outer iterations
    qindex numOuterIts = qureg.numAmpsPerNode / powerOf2(suffixTargsXY.size() + ctrls.size());

    // prepare a mask which yields ctrls in specified state, and targs in all-zero
    auto sortedQubits   = util_getSorted(ctrls, suffixTargsXY);
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, suffixTargsXY, vector<int>(suffixTargsXY.size(),0));

    // use template params to compile-time unroll loops in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, suffixTargsXY.size());

    // compiler will infer these at compile-time if possible
    int numQubitBits = numCtrlBits + numTargBits;
    qindex numTargAmps = powerOf2(numTargBits);

    // each inner iteration modifies 2 amplitudes (may be compile-time sized) 
    qindex numInnerIts = numTargAmps / 2;
    
    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numOuterIts; n++) {

        // i0 = nth local index where ctrls are active and targs are all zero
        qindex i0 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);

        // loop may be unrolled
        for (qindex m=0; m<numInnerIts; m++) {

            // iA = nth local index where targs have value m, iB = (last - nth) such index
            qindex iA = setBits(i0, suffixTargsXY.data(), suffixTargsXY.size(), m);
            qindex iB = flipBits(iA, suffixMaskXY);

            // jA = global index corresponding to iA
            qindex jA = concatenateBits(qureg.rank, iA, qureg.logNumAmpsPerNode);
            qindex jB = concatenateBits(qureg.rank, iB, qureg.logNumAmpsPerNode);

            // determine whether to multiply amps by +-1 or +-i
            qcomp pmPowA = powI * (1 - 2 * getBitMaskParity(jA & allMaskYZ));
            qcomp pmPowB = powI * (1 - 2 * getBitMaskParity(jB & allMaskYZ));

            qcomp ampA = qureg.cpuAmps[iA];
            qcomp ampB = qureg.cpuAmps[iB];

            // mix or swap scaled amp pair
            qureg.cpuAmps[iA] = (thisAmpFac * ampA) + (otherAmpFac * pmPowB * ampB);
            qureg.cpuAmps[iB] = (thisAmpFac * ampB) + (otherAmpFac * pmPowA * ampA);
        }
    }
}


template <int NumCtrls>
void cpu_statevector_anyCtrlPauliTensorOrGadget_subB(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates,
    qindex suffixMaskXY, qindex bufferMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
) {
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());
    
    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // j = buffer index of amp to be mixed with i
        qindex j = flipBits(n, bufferMaskXY) + offset;

        // k = global index of amp at buffer index j
        qindex k = concatenateBits(qureg.rank, flipBits(i, suffixMaskXY), qureg.logNumAmpsPerNode);

        // determine whether to multiply buffer amp by +-1 or +-i
        int negParity = getBitMaskParity(k & allMaskYZ);
        qcomp pmPowI = (1 - 2*negParity) * powI;

        qureg.cpuAmps[i] *= fac0;
        qureg.cpuAmps[i] += fac1 * pmPowI * qureg.cpuCommBuffer[j];
    }
}


template <int NumCtrls>
void cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, 
    qcomp fac0, qcomp fac1
) {
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // each control qubit halves the needed iterations, each of which modifies 1 amp
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    qcomp facs[] = {fac0, fac1};
    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    auto targMask      = util_getBitMask(targs);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // j = global index corresponding to i
        qindex j = concatenateBits(qureg.rank, i, qureg.numAmpsPerNode);

        // apply phase to amp depending on parity of targets in global index 
        int p = getBitMaskParity(j & targMask);
        qureg.cpuAmps[j] *= facs[p];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevector_anyCtrlPauliTensorOrGadget_subA, (Qureg, vector<int>, vector<int>, vector<int>, qindex, qindex, qcomp, qcomp, qcomp) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevector_anyCtrlPauliTensorOrGadget_subB, (Qureg, vector<int>, vector<int>, qindex, qindex, qindex, qcomp, qcomp, qcomp) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub, (Qureg, vector<int>, vector<int>, vector<int>, qcomp, qcomp) )



/*
 * QUREG COMBINATION
 */


void cpu_densmatr_mixQureg_subA(qreal outProb, Qureg outQureg, qreal inProb, Qureg inDensMatr) {

    qindex numIts = outQureg.numAmpsPerNode;
    qcomp* out = outQureg.cpuAmps;
    qcomp* in = inDensMatr.cpuAmps;
    
    #pragma omp parallel for if(outQureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++)
        out[n] = (outProb * out[n]) + (inProb * in[n]);
}


void cpu_densmatr_mixQureg_subB(qreal outProb, Qureg outQureg, qreal inProb, Qureg inStateVec) {

    qindex numIts = outQureg.numAmpsPerNode;
    qindex dim = inStateVec.numAmps;
    qcomp* out = outQureg.cpuAmps;
    qcomp* in = inStateVec.cpuAmps;
    
    #pragma omp parallel for if(outQureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // (i,j) = row & column of outQureg corresponding to n
        qindex i = n % dim;
        qindex j = n / dim; // floors

        out[n] = (outProb * out[n]) + (inProb * in[i] * conj(in[j]));
    }
}


void cpu_densmatr_mixQureg_subC(qreal outProb, Qureg outQureg, qreal inProb) {

    // received inQureg's entire statevector amplitudes into every node's buffer
    qindex numIts = outQureg.numAmpsPerNode;
    qindex dim = powerOf2(outQureg.numQubits);
    qcomp* out = outQureg.cpuAmps;
    qcomp* in = outQureg.cpuCommBuffer;

    #pragma omp parallel for if(outQureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // m = global index of local index n
        qindex m = concatenateBits(outQureg.rank, n, outQureg.logNumAmpsPerNode);

        // (i,j) = global row & column of outQureg corresponding to n
        qindex i = m % dim;
        qindex j = m / dim; // floors

        out[n] = (outProb * out[n]) + (inProb * in[i] * conj(in[j]));
    }
}



/*
 * ONE-QUBIT DEPHASING
 */


void cpu_densmatr_oneQubitDephasing_subA(Qureg qureg, int ketQubit, qreal prob) {

    // half of all local amps are scaled, and each iteration modifies two
    qindex numIts = qureg.numAmpsPerNode / 4;

    // loop constants
    qreal fac = util_getOneQubitDephasingFactor(prob);
    int braQubit = util_getBraQubit(ketQubit, qureg);

    // TODO:
    // this enumeration order is suboptimal and seems unnecessary in this simple two
    // bit scenario, where we are modifying but not at all mixing two strided and
    // potentially very-distant amplitudes. It is of course trivial to split this
    // into two separate loops accessing monotonically increasing indices, although
    // we then pay double the caching costs when ketQubit is low-index. We can also
    // turn it into two nested loops to force monotonically increasing index access,
    // but then the parallelisation is not optimal when ketQubit is high-index. Experiment
    // with what's fastest and replace below or delete this comment!

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i01 = nth local index of |*0*><*1*|
        qindex i01 = insertTwoBits(n, braQubit, 0, ketQubit, 1);
        qindex i10 = insertTwoBits(n, braQubit, 1, ketQubit, 0);

        qureg.cpuAmps[i01] *= fac;
        qureg.cpuAmps[i10] *= fac;
    }
}


void cpu_densmatr_oneQubitDephasing_subB(Qureg qureg, int ketQubit, qreal prob) {

    // half of all local amps are scaled
    qindex numIts = qureg.numAmpsPerNode / 2;
    
    // loop constants
    qreal fac = util_getOneQubitDephasingFactor(prob);
    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    
    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where bra-qubit differs from ket-qubit
        qindex i = insertBit(n, ketQubit, ! braBit);

        qureg.cpuAmps[i] *= fac;
    }
}



/*
 * TWO-QUBIT DEPHASING
 */



void cpu_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubitA, int qubitB, qreal prob) {

    // TODO: 
    // test whether use of subB has identical performance, or whether changing i=n below
    // non-negligibly accelerates the routine; if so, make a templated inner func.

    // the rank-agnostic version is identical to the subB algorithm below, because the
    // queried bits of the global index i below will always be in the suffix substate.
    // We still define separate _subA and _subB routines because they differ in cuQuantum.
    cpu_densmatr_twoQubitDephasing_subB(qureg, qubitA, qubitB, prob);
}


void cpu_densmatr_twoQubitDephasing_subB(Qureg qureg, int ketQubitA, int ketQubitB, qreal prob) {

    // 75% of amps are updated, but we just enumerate all for simplicity
    qindex numIts = qureg.numAmpsPerNode;

    // loop constants
    qreal term = util_getTwoQubitDephasingTerm(prob);

    int braQubitA = util_getBraQubit(ketQubitA, qureg);
    int braQubitB = util_getBraQubit(ketQubitB, qureg);

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = global index of nth local amp
        qindex i = concatenateBits(qureg.rank, n, qureg.logNumAmpsPerNode);

        int bitA = getBit(i, ketQubitA) ^ getBit(i, braQubitA);
        int bitB = getBit(i, ketQubitB) ^ getBit(i, braQubitB);

        // determine whether or not to modify this amplitude...
        int flag = bitA | bitB;

        // by multiplying by 1 or (1 + term)
        qureg.cpuAmps[n] *= 1 + (flag * term);
    }
}



/*
 * ONE-QUBIT DEPOLARISING
 */


void cpu_densmatr_oneQubitDepolarising_subA(Qureg qureg, int ketQubit, qreal prob) {

    // all amps are modified, and each iteration modifies 4
    qindex numIts = qureg.numAmpsPerNode / 4;

    // for brevity
    qcomp* amps = qureg.cpuAmps;

    int braQubit = util_getBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitDepolarisingFactors(prob);

    auto facAA = factors.c1;
    auto facBB = factors.c2;
    auto facAB = factors.c3;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i00 = nth local index where both qubits are 0
        qindex i00 = insertTwoBits(n, braQubit, 0, ketQubit, 0);
        qindex i01 = flipBit(i00, ketQubit);
        qindex i10 = flipBit(i00, braQubit);
        qindex i11 = flipBit(i01, braQubit);

        // modify 4 amps, mixing a pair, and scaling the other
        qcomp amp00 = amps[i00];
        qcomp amp11 = amps[i11];
        amps[i00] = (facAA * amp00) + (facBB * amp11);
        amps[i01] *= facAB;
        amps[i10] *= facAB;
        amps[i11] = (facAA * amp11) + (facBB * amp00);
    }
}


void cpu_densmatr_oneQubitDepolarising_subB(Qureg qureg, int ketQubit, qreal prob) {

    // all amps are modified, and each iteration modifies 2
    qindex numIts = qureg.numAmpsPerNode / 2;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitDepolarisingFactors(prob);

    auto facAA = factors.c1;
    auto facBB = factors.c2;
    auto facAB = factors.c3;

    // TODO:
    // each iteration below modifies 2 independent amps without mixing,
    // which we can trivially split into two loops which may improve
    // per-iteration caching performance; test if this outweights the 
    // cost of re-iteration

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // iAA = nth local index where ket qubit agrees with bra qubit
        qindex iAA = insertBit(n, ketQubit, braBit);

        // jBB = buffer index of amp iAA where ket and bra qubits are flipped
        qindex jBB = n + offset;

        // iAB = nth local index where ket qubit disagrees with bra qubit
        qindex iAB = insertBit(n, ketQubit, ! braBit);

        qureg.cpuAmps[iAA] *= facAA;
        qureg.cpuAmps[iAA] += facBB * qureg.cpuCommBuffer[jBB];
        qureg.cpuAmps[iAB] *= facAB;
    }
}



/*
 * TWO-QUBIT DEPOLARISING
 */


void cpu_densmatr_twoQubitDepolarising_subA(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // all amps are scaled (although 1/16 of them will be unchanged)
    qindex numIts  = qureg.numAmpsPerNode;

    // bra-qubits corresponding to ket-qubits
    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braQb2 = util_getBraQubit(ketQb2, qureg);

    auto c3 = util_getTwoQubitDepolarisingFactors(prob).c3;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // determine whether to modify amp
        int flag1 = !(getBit(n, ketQb1) ^ getBit(n, braQb1));
        int flag2 = !(getBit(n, ketQb2) ^ getBit(n, braQb2));
        int mod   = !(flag1 & flag2);

        // multiply 15/16 of all amps by (1 + c3)
        qureg.cpuAmps[n] *= 1 + c3 * mod;
    }
}


void cpu_densmatr_twoQubitDepolarising_subB(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // one quarter of amps will be modified, and four are mixed each iteration
    qindex numIts = qureg.numAmpsPerNode / 16;

    // bra-qubits corresponding to ket-qubits
    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braQb2 = util_getBraQubit(ketQb2, qureg);

    auto factors = util_getTwoQubitDepolarisingFactors(prob);
    auto c1 = factors.c1;
    auto c2 = factors.c2;

    // for brevity
    qcomp* amps = qureg.cpuAmps;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {
    
        // i0000 = nth local index where all bra = ket = 00
        qindex i0000 = insertFourZeroBits(n, braQb2, braQb1, ketQb2, ketQb1);
        qindex i0101 = flipTwoBits(i0000, braQb1, ketQb1);
        qindex i1010 = flipTwoBits(i0000, braQb2, ketQb2);
        qindex i1111 = flipTwoBits(i0101, braQb2, ketQb2);
        
        // mix 1/16 of all amps in groups of 4
        qcomp term = amps[i0000] + amps[i0101] + amps[i1010] + amps[i1111];

        amps[i0000] = c1*amps[i0000] + c2*term;
        amps[i0101] = c1*amps[i0101] + c2*term;
        amps[i1010] = c1*amps[i1010] + c2*term;
        amps[i1111] = c1*amps[i1111] + c2*term;
    }
}


void cpu_densmatr_twoQubitDepolarising_subC(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // scale 25% of amps but iterate all
    qindex numIts = qureg.numAmpsPerNode;

    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    auto c3 = util_getTwoQubitDepolarisingFactors(prob).c3;

    // TODO:
    // are we really inefficiently enumerating all amps and applying a non-unity
    // factor to only 25%?! Is this because we do not know braBit2 and ergo 
    // cannot be sure a direct enumeration is accessing indicies in a monotonically
    // increasing order? Can that really outweigh a 3x slowdown?! Test and fix!

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // decide whether or not to modify nth local
        bool flag1 = getBit(n, ketQb1) == getBit(n, braQb1); 
        bool flag2 = getBit(n, ketQb2) == braBit2;
        bool mod   = !(flag1 & flag2);

        // scale amp by 1 or (1 + c3)
        qureg.cpuAmps[n] *= 1 + c3 * mod;
    }
}


void cpu_densmatr_twoQubitDepolarising_subD(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // 25% of local amps are modified, two in each iteration
    qindex numIts = qureg.numAmpsPerNode / 8;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    auto factors = util_getTwoQubitDepolarisingFactors(prob);
    auto c1 = factors.c1;
    auto c2 = factors.c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i000 = nth local index where all suffix bits are 0
        qindex i000 = insertThreeZeroBits(n, braQb1, ketQb2, ketQb1);
        qindex i0b0 = setBit(i000, ketQb2, braBit2);
        qindex i1b1 = flipTwoBits(i0b0, braQb1, ketQb1);

        // j = nth received amp in buffer
        qindex j = n + offset;

        // mix pair of amps using buffer
        qcomp amp0b0 = qureg.cpuAmps[i0b0];
        qcomp amp1b1 = qureg.cpuAmps[i1b1];

        qureg.cpuAmps[i0b0] = c1*amp0b0 + c2*(amp1b1 + qureg.cpuCommBuffer[j]);
        qureg.cpuAmps[i1b1] = c1*amp1b1 + c2*(amp0b0 + qureg.cpuCommBuffer[j]);
    }
}


void cpu_densmatr_twoQubitDepolarising_subE(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // scale 25% of amps but iterate all
    qindex numIts = qureg.numAmpsPerNode;

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    qreal c3 = util_getTwoQubitDepolarisingFactors(prob).c3;

    // TODO:
    // are we really inefficiently enumerating all amps and applying a non-unity
    // factor to only 25%?! Is this because we do not know braBit2 and ergo 
    // cannot be sure a direct enumeration is accessing indicies in a monotonically
    // increasing order? Can that really outweigh a 3x slowdown?! Test and fix!

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // choose whether to modify amp
        bool flag1 = getBit(n, ketQb1) == braBit1; 
        bool flag2 = getBit(n, ketQb2) == braBit2;
        bool mod   = !(flag1 & flag2);
        
        // multiply amp by 1 or (1 + c3)
        qureg.cpuAmps[n] *=  1 + c3 * mod;
    }
}


void cpu_densmatr_twoQubitDepolarising_subF(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // modify 25% of local amps, one per iteration
    qindex numIts = qureg.numAmpsPerNode / 4;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    auto factors = util_getTwoQubitDepolarisingFactors(prob);
    auto c1 = factors.c1;
    auto c2 = factors.c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where suffix ket qubits equal prefix bra qubits
        qindex i = insertTwoBits(n, ketQb2, braBit2, ketQb1, braBit1);

        // j = nth received amp in buffer
        qindex j = n + offset;

        // mix local amp with received buffer amp
        qureg.cpuAmps[i] *= c1;
        qureg.cpuAmps[i] += c2 * qureg.cpuCommBuffer[j];
    }
}


void cpu_densmatr_twoQubitDepolarising_subG(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // modify 25% of local amps, one per iteration
    qindex numIts = qureg.numAmpsPerNode / 4;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    auto factors = util_getTwoQubitDepolarisingFactors(prob);
    auto c1 = factors.c1;
    auto c2 = factors.c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where suffix ket qubits equal prefix bra qubits
        qindex i = insertTwoBits(n, ketQb2, braBit2, ketQb1, braBit1);

        // j = nth received amp in buffer
        qindex j = n + offset;

        // overwrite local amp with buffer amp
        qureg.cpuAmps[i] = (c2 / c1) * qureg.cpuCommBuffer[j];
    }
}



/*
 * PAULI CHANNEL
 */


void cpu_densmatr_oneQubitPauliChannel_subA(Qureg qureg, int ketQubit, qreal pI, qreal pX, qreal pY, qreal pZ) {

    // all amps are modified, and each iteration modifies 4
    qindex numIts = qureg.numAmpsPerNode / 4;

    int braQubit = util_getBraQubit(ketQubit, qureg);

    auto factors = util_getOneQubitPauliChannelFactors(pI, pX, pY, pZ);
    auto facAA = factors.c1;
    auto facBB = factors.c2;
    auto facAB = factors.c3;
    auto facBA = factors.c4;

    // for brevity
    qcomp* amps = qureg.cpuAmps;

    // TODO:
    // each iteration modifies 4 amps in two separable mixed pairs, which may
    // lead to sub-optimal caching. Iterating twice and modifying a single pair
    // might lead to better performance, though note the stride from i00 to i11
    // will always be adverserially large. Test this!

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i00 = nth local index where both qubits are 0
        qindex i00 = insertTwoBits(n, braQubit, 0, ketQubit, 0);
        qindex i01 = flipBit(i00, ketQubit);
        qindex i10 = flipBit(i00, braQubit);
        qindex i11 = flipBit(i01, braQubit);

        // modify 4 amps in 2 separable pairs
        qcomp amp00 = amps[i00];
        qcomp amp01 = amps[i01];
        qcomp amp10 = amps[i10];
        qcomp amp11 = amps[i11];

        amps[i00] = (facAA * amp00) + (facBB * amp11);
        amps[i01] = (facAB * amp01) + (facBA * amp10);
        amps[i10] = (facAB * amp10) + (facBA * amp01);
        amps[i11] = (facAA * amp11) + (facBB * amp00);
    }
}


void cpu_densmatr_oneQubitPauliChannel_subB(Qureg qureg, int ketQubit, qreal pI, qreal pX, qreal pY, qreal pZ) {

    // all amps are modified, and each iteration modifies 2
    qindex numIts = qureg.numAmpsPerNode / 2;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    int braInd = util_getPrefixBraInd(ketQubit, qureg);
    int braBit = getBit(qureg.rank, braInd);

    auto factors = util_getOneQubitPauliChannelFactors(pI, pX, pY, pZ);
    auto facAA = factors.c1;
    auto facBB = factors.c2;
    auto facAB = factors.c3;
    auto facBA = factors.c4;

    // TODO:
    // each iteration below modifies 2 independent amps without mixing,
    // which we can trivially split into two loops which may improve
    // per-iteration caching performance; test if this outweights the 
    // cost of re-iteration

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // iAA = nth local index where ket qubit agrees with bra, i.e. |.A.><.A.|
        qindex iAA = insertBit(n, ketQubit, braBit);
        qindex iAB = flipBit(iAA, ketQubit);

        // jBB = buffer index of amp to be mixed with iAA's amp, i.e. |.B.><.B.|
        qindex jBB = iAB + offset;
        qindex jBA = iAA + offset;

        // mix each local amp with a received buffer amp
        qureg.cpuAmps[iAA] *= facAA;
        qureg.cpuAmps[iAA] += facBB * qureg.cpuCommBuffer[jBB];

        qureg.cpuAmps[iAB] *= facAB;
        qureg.cpuAmps[iAB] += facBA * qureg.cpuCommBuffer[jBA];
    }
}



/*
 * AMPLITUDE DAMPING CHANNEL
 */


void cpu_densmatr_oneQubitDamping_subA(Qureg qureg, int ketQubit, qreal prob) {

    // each iteration modifies 4 amps
    qindex numIts = qureg.numAmpsPerNode / 4;

    int braQubit = util_getBraQubit(ketQubit, qureg);

    auto factors = util_getOneQubitDampingFactors(prob);
    auto c1 = factors.c1;
    auto c2 = factors.c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i00 = nth local index where bra and ket qubits are 0
        qindex i00 = insertTwoBits(n, braQubit, 0, ketQubit, 0);
        qindex i01 = flipBit(i00, ketQubit);
        qindex i10 = flipBit(i00, braQubit);
        qindex i11 = flipBit(i01, braQubit);
        
        // mix both-zero amp with both-one amp (but not vice versa)
        qureg.cpuAmps[i00] += prob * qureg.cpuAmps[i11];

        // scale other amps
        qureg.cpuAmps[i01] *= c1;
        qureg.cpuAmps[i10] *= c1;
        qureg.cpuAmps[i11] *= c2;
    }
}


void cpu_densmatr_oneQubitDamping_subB(Qureg qureg, int qubit, qreal prob) {

    // half of all local amps are scaled
    qindex numIts = qureg.numAmpsPerNode / 2;

    auto c2 = util_getOneQubitDampingFactors(prob).c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where qubit=1
        qindex i= insertBit(n, qubit, 1);
        qureg.cpuAmps[i] *= c2;
    }
}


void cpu_densmatr_oneQubitDamping_subC(Qureg qureg, int ketQubit, qreal prob) {

    // half of all local amps are scaled
    qindex numIts = qureg.numAmpsPerNode / 2;

    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    auto c1 = util_getOneQubitDampingFactors(prob).c1;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ket differs from bra
        qindex i = insertBit(n, ketQubit, ! braBit);
        qureg.cpuAmps[i] *= c1;
    }
}


void cpu_densmatr_oneQubitDamping_subD(Qureg qureg, int qubit, qreal prob) {

    // half of all local amps are combined with buffer
    qindex numIts = qureg.numAmpsPerNode / 2;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ket is 0
        qindex i = insertBit(n, qubit, 0);

        // j = nth received buffer indes
        qindex j = n + offset;

        qureg.cpuAmps[i] += prob * qureg.cpuCommBuffer[j];
    }
}



/*
 * PARTIAL TRACE
 */


template <int NumTargs>
void cpu_densmatr_partialTrace_sub(Qureg inQureg, Qureg outQureg, vector<int> targs, vector<int> pairTargs) {

    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

    // each outer iteration sets one element of outQureg
    qindex numOuterIts = outQureg.numAmpsPerNode;

    // targs and allTargs are sorted, but pairTargs is arbitrarily ordered (though corresponding targs)
    auto allTargs = util_getSorted(targs, pairTargs);

    // use template param to compile-time unroll below loops
    SET_VAR_AT_COMPILE_TIME(int, numTargPairs, NumTargs, targs.size());
    
    // may be inferred at compile-time
    int numAllTargs = 2*numTargPairs;
    qindex numInnerIts = powerOf2(numTargPairs);

    // TODO:
    // note our parallelisation of only the outer-loop assumes that the number of 
    // amps in outQureg equals or exceeds the number of threads. Ergo tracing out 
    // all but very few qubits will leave threads idle; when only a single qubit
    // remains, the below code would be serial. In that scenario, we should
    // parallelise the inner loop, or preclude this scenario in validation.

    // consult inQureg for multithreading, because total iters = inQureg dim
    #pragma omp parallel for if(inQureg.isMultithreaded)
    for (qindex n=0; n<numOuterIts; n++) {

        // k = nth local index of inQureg where all targs and pairs are zero
        qindex k = insertBits(n, allTargs.data(), numAllTargs, 0); // loop may be unrolled

        // each outQureg amp results from summing 2^targs inQureg amps
        qcomp outAmp = 0;

        // loop may be unrolled
        for (qindex j=0; j<numInnerIts; j++) {

            // i = nth local index of inQureg where targs=j and pairTargs=j
            qindex i = k;
            i = setBits(i, targs    .data(), numTargPairs, j); // loop may be unrolled
            i = setBits(i, pairTargs.data(), numTargPairs, j); // loop may be unrolled

            outAmp += inQureg.cpuAmps[i];
        }

        outQureg.cpuAmps[n] = outAmp;
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, cpu_densmatr_partialTrace_sub, (Qureg, Qureg, vector<int>, vector<int>) )
