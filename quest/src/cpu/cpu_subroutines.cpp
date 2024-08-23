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
void cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubitInds, vector<int> qubitStates) {

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
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) )



/* 
 * SWAPS
 */


template <int NumCtrls> void cpu_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

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


template <int NumCtrls> void cpu_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

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


template <int NumCtrls> void cpu_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState) {

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

        // j = index of nth received amplitude from pair rank in buffer
        qindex j = n + offset;

        qureg.cpuAmps[i] = fac0*qureg.cpuAmps[i] + fac1*qureg.cpuCommBuffer[j];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, int, CompMatr1) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) )



/*
 * MANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs>
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
                for (qindex j=0; j<numTargAmps; j++)
                    qureg.cpuAmps[i] += matr.cpuElems[k][j] * cache[j];
            }
        }
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDenseMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, CompMatr) )



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
 * DECOHERENCE
 */


void cpu_densmatr_oneQubitDephasing_subA(Qureg qureg, int ketQubit, qreal prob) {

    // half of all local amps are scaled, and each iteration modifies two
    qindex numIts = qureg.numAmpsPerNode / 4;

    // loop constants
    qreal fac = 1 - 2*prob;
    int braQubit = ketQubit + qureg.numQubits;

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
    qreal fac = 1 - 2*prob;
    int braBit = getBit(qureg.rank, ketQubit - qureg.logNumColsPerNode);
    
    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where bra-qubit differs from ket-qubit
        qindex i = insertBit(n, ketQubit, ! braBit);

        qureg.cpuAmps[i] *= fac;
    }
}


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
    qreal term = - 4 * prob / 3;
    int braQubitA = ketQubitA + qureg.numQubits;
    int braQubitB = ketQubitB = qureg.numQubits;

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
