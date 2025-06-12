// Prepare sendIndices and sendBuffers for fused multi-SWAP, parallelized with OpenMP
#include <vector>
#include <omp.h>
#include "quest/include/qureg.h"
#include "quest/include/types.h"
#include <complex>

/** @file
 * CPU OpenMP-accelerated definitions of the main backend simulation routines,
 * as mirrored by gpu_subroutines.cpp, and called by accelerator.cpp. 
 * 
 * Some of these definitions are templated, defining multiple versions optimised 
 * (at compile-time) for handling different numbers of input qubits; such functions
 * are proceeded by macro INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(), to force the 
 * compilation of their needed versions within this translation unit for later linkage.
 * 
 * @author Tyson Jones
 * @author Oliver Brown (OpenMP 'if' clauses)
 * @author Richard Meister (helped patch on LLVM)
 * @author Kshitij Chhabra (patched v3 clauses with gcc9)
 * @author Ania (Anna) Brown (developed QuEST v1 logic)
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/inliner.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/fastmath.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/randomiser.hpp"
#include "quest/src/core/accelerator.hpp"
#include "quest/src/core/autodeployer.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/cpu/cpu_subroutines.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_indices.hpp"

#include <vector>
#include <algorithm>

using std::vector;



/*
 * GETTERS
 */


qcomp cpu_statevec_getAmp_sub(Qureg qureg, qindex ind) {

    // this bespoke function exists (rather than merely
    // calling the bulk memcpy routine) because it is
    // much faster for few randomly accessed amps
    return qureg.cpuAmps[ind];
}



/*
 * SETTERS
 */


void cpu_densmatr_setAmpsToPauliStrSum_sub(Qureg qureg, PauliStrSum sum) {

    // this assertion exists because fast_getPauliStrElem() (invoked below)
    // previously assumed str.highPaulis=0 for all str in sum (for a speedup)
    // which is gauranteed satisfied for all sum compatible with a density-matrix.
    // This is no longer essential, since fast_getPauliStrElem() has relaxed this
    // requirement and foregone the optimisation, but we retain this check in
    // case a similar optimisation is restored in the future
    assert_highPauliStrSumMaskIsZero(sum);

    // process each amplitude in-turn, not bothering to leverage that adjacent
    // basis states have PauliStrSum elems which differ by a single +-i/1 (as
    // can be enumerated via Gray Code), because this breaks thread independence,
    // plus this function is only called infrequently (as initialisation)
    qindex numIts = qureg.numAmpsPerNode;
    qindex dim = powerOf2(qureg.numQubits);

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = global flat index corresponding to n
        qindex i = concatenateBits(qureg.rank, n, qureg.logNumAmpsPerNode);

        // r, c = global row and column
        qindex r = fast_getQuregGlobalRowFromFlatIndex(i, dim);
        qindex c = fast_getQuregGlobalColFromFlatIndex(i, dim);

        // contains non-unrolled loop (and args unpacked due to CUDA qcomp incompatibility, grr)
        qureg.cpuAmps[n] = fast_getPauliStrSumElem(sum.coeffs, sum.strings, sum.numTerms, r, c);
    }
}


void cpu_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in) {

    // unlike in densmatr_setAmpsToPauliStrSum_sub() above, this PauliStrSum
    // can feature non-identity Paulis on every qubit, i.e. up to t=63

    qindex numIts = out.numElemsPerNode;
    qindex numSuf = logBase2(numIts);

    int rank = out.isDistributed? comm_getRank() : 0;

    #pragma omp parallel for if(out.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = global index corresponding to local index n
        qindex i = concatenateBits(rank, n, numSuf);

        // treat PauliStrSum as a generic sum, even though we know it only
        // contains I and Z which can in principle be computed faster; this
        // is a superfluous optimisation since this function is expected to
        // be called infrequently (i.e. only for data structure initialisation)
        out.cpuElems[n] = fast_getPauliStrSumElem(in.coeffs, in.strings, in.numTerms, i, i);
    }
}


void cpu_fullstatediagmatr_setElemsFromMultiDimLists(FullStateDiagMatr out, void* lists, int* numQubitsPerDim, int numDims) {

    // note that this function has no GPU equivalent! This is because
    // it processes arbitrarily nested input pointers which would could
    // be adverserial to loading into GPU memory (and generally painful)

    qindex numIts = out.numElemsPerNode;
    qindex numLocalIndBits = logBase2(numIts);

    int rank = out.isDistributed? comm_getRank() : 0;

    // create an explicit parallel region to avoid re-initialisation of vectors every iteration
    #pragma omp parallel if(out.isMultithreaded)
    {
        // create a private cache for every thread
        vector<qindex> listInds(numDims);

        #pragma omp for
        for (qindex localInd=0; localInd<numIts; localInd++) {

            // each local diag index corresponds to a unique global index which informs list indices
            qindex globalInd = concatenateBits(rank, localInd, numLocalIndBits);
            fast_getSubQuregValues(globalInd, numQubitsPerDim, numDims, false, listInds.data());

            // update only the CPU elems using lists which are duplicated on every node.
            // note we are calling a util-function inside a parallelised hot-loop which
            // is generally inadvisable, but it does not matter in this case since the
            // function is recursive and cannot be inlined.
            out.cpuElems[localInd] = util_getElemFromNestedPtrs(lists, listInds.data(), numDims);
        }
    }
}


void cpu_fullstatediagmatr_setElemsFromMultiVarFunc(FullStateDiagMatr out, qcomp (*callbackFunc)(qindex*), int* numQubitsPerVar, int numVars, int areSigned) {

    // note that this function has no GPU equivalent! This is because
    // the user's callback function cannot be called by a GPU kernel

    qindex numIts = out.numElemsPerNode;
    qindex numLocalIndBits = logBase2(numIts);

    int rank = out.isDistributed? comm_getRank() : 0;

    // create an explicit parallel region to avoid re-initialisation of vectors every iteration
    #pragma omp parallel if(out.isMultithreaded)
    {
        // create a private cache for every thread
        vector<qindex> varValues(numVars);

        #pragma omp for
        for (qindex localInd=0; localInd<out.numElemsPerNode; localInd++) {

            // each local index corresponds to a unique global index which informs variable values
            qindex globalInd = concatenateBits(rank, localInd, numLocalIndBits);
            fast_getSubQuregValues(globalInd, numQubitsPerVar, numVars, areSigned, varValues.data());
    
            // call user function, which we assume is thread safe!
            out.cpuElems[localInd] = callbackFunc(varValues.data());
        }
    }
}



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

    /// @todo
    /// this function allocates powerOf2(targs.size())-sized caches for each thread, sometimes in
    /// heap. At the ~max non-distributed double CompMatr of 16 qubits = 64 GiB, this is 1 MiB 
    /// per thread; for a conceivable 100 thread execution, this is 100 MiB being alloc/dealloced
    /// at every call. It is debatable whether this justifies pre-allocating persistent cache space
    /// (one for each thread, to avoid false sharing), similar to GPU's AnyTargDenseMatr, though
    /// for an order of magnitude fewer threads, and using non-coalesced memory. Certainly making
    /// persistent heap caches is inadvisable when the cache fits in the stack (currently automated 
    /// using std::vector). Perhaps we should keep the current re-allocs, constrain that this 
    /// function is only called for few-targets (e.g. <= qureg.numQubits - 5), and define another
    /// function for almost-all target matrices which uses persistent heap memory, wherein the 
    /// optimal parallelisation scheme is anyway different.

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

                    // matr.cpuElems[k][j] = matr.cpuElemsFlat[l]
                    qindex l = fast_getMatrixFlatIndex(k, j, numTargAmps);
                    qcomp elem = matr.cpuElemsFlat[l];

                    // optionally conjugate matrix elems on the fly to avoid pre-modifying heap structure
                    if constexpr (ApplyConj)
                        elem = std::conj(elem);

                    qureg.cpuAmps[i] += elem * cache[j];

                    /// @todo
                    /// qureg.cpuAmps[i] is being serially updated by only this thread,
                    /// so is a candidate for Kahan summation for improved numerical
                    /// stability. Explore whether this is time-free and worthwhile!
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
        qureg.cpuAmps[j] *= matr.elems[b];
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
        qureg.cpuAmps[j] *= matr.elems[k];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlTwoTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, int, int, DiagMatr2) )



/*
 * ANY-TARG DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj, bool HasPower>
void cpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, qcomp exponent) {
    
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);
    assert_exponentMatchesTemplateParam(exponent, HasPower);

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
        qcomp elem = matr.cpuElems[t];

        // decide whether to power and conj at compile-time, to avoid branching in hot-loop.
        // beware that pow(qcomp,qcomp) below gives notable error over pow(qreal,qreal) 
        // (by producing an unexpected non-zero imaginary component) when the base is real 
        // and negative, and the exponent is an integer. We tolerate this heightened error
        // because we have no reason to think matr is real (it's not constrained Hermitian).
        if constexpr (HasPower)
            elem = std::pow(elem, exponent);

        // cautiously conjugate AFTER exponentiation, else we must also conj exponent
        if constexpr (ApplyConj)
            elem = std::conj(elem);

        qureg.cpuAmps[j] *= elem;
    }
}


INSTANTIATE_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr, qcomp) )



// Minimal stub non-template function definitions for missing routines to resolve linker errors
qreal cpu_statevec_calcTotalProb_sub(Qureg) { return 0; }
qreal cpu_densmatr_calcTotalProb_sub(Qureg) { return 0; }
qcomp cpu_statevec_calcInnerProduct_sub(Qureg, Qureg) { return 0; }
void cpu_statevec_initUniformState_sub(Qureg, qcomp) {}
void cpu_statevec_initDebugState_sub(Qureg) {}
void cpu_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(Qureg) {}
qcomp cpu_densmatr_calcExpecAnyTargZ_sub(Qureg, std::vector<int>) { return 0; }
qreal cpu_statevec_calcExpecAnyTargZ_sub(Qureg, std::vector<int>) { return 0; }
void cpu_densmatr_mixQureg_subA(double, Qureg, double, Qureg) {}
void cpu_densmatr_oneQubitDephasing_subA(Qureg, int, double) {}
void cpu_densmatr_oneQubitDephasing_subB(Qureg, int, double) {}
void cpu_densmatr_twoQubitDephasing_subA(Qureg, int, int, double) {}
void cpu_densmatr_twoQubitDephasing_subB(Qureg, int, int, double) {}
template <bool Conj>
qcomp cpu_densmatr_calcFidelityWithPureState_sub(Qureg, Qureg) { return 0; }
template qcomp cpu_densmatr_calcFidelityWithPureState_sub<true>(Qureg, Qureg);
template qcomp cpu_densmatr_calcFidelityWithPureState_sub<false>(Qureg, Qureg);

// Explicit template instantiations for multi-qubit routines (0–5, -1)
template void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub<0>(qreal*, Qureg, std::vector<int>);
template void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub<1>(qreal*, Qureg, std::vector<int>);
template void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub<2>(qreal*, Qureg, std::vector<int>);
template void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub<3>(qreal*, Qureg, std::vector<int>);
template void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub<4>(qreal*, Qureg, std::vector<int>);
template void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub<5>(qreal*, Qureg, std::vector<int>);
template void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub<-1>(qreal*, Qureg, std::vector<int>);

template void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<0>(qreal*, Qureg, std::vector<int>);
template void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<1>(qreal*, Qureg, std::vector<int>);
template void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<2>(qreal*, Qureg, std::vector<int>);
template void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<3>(qreal*, Qureg, std::vector<int>);
template void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<4>(qreal*, Qureg, std::vector<int>);
template void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<5>(qreal*, Qureg, std::vector<int>);
template void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub<-1>(qreal*, Qureg, std::vector<int>);

template void cpu_statevec_multiQubitProjector_sub<0>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_statevec_multiQubitProjector_sub<1>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_statevec_multiQubitProjector_sub<2>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_statevec_multiQubitProjector_sub<3>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_statevec_multiQubitProjector_sub<4>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_statevec_multiQubitProjector_sub<5>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_statevec_multiQubitProjector_sub<-1>(Qureg, std::vector<int>, std::vector<int>, qreal);

template void cpu_densmatr_multiQubitProjector_sub<0>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_densmatr_multiQubitProjector_sub<1>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_densmatr_multiQubitProjector_sub<2>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_densmatr_multiQubitProjector_sub<3>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_densmatr_multiQubitProjector_sub<4>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_densmatr_multiQubitProjector_sub<5>(Qureg, std::vector<int>, std::vector<int>, qreal);
template void cpu_densmatr_multiQubitProjector_sub<-1>(Qureg, std::vector<int>, std::vector<int>, qreal);

template qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub<0>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub<1>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub<2>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub<3>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub<4>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub<5>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub<-1>(Qureg, std::vector<int>, std::vector<int>);

template qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub<0>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub<1>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub<2>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub<3>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub<4>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub<5>(Qureg, std::vector<int>, std::vector<int>);
template qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub<-1>(Qureg, std::vector<int>, std::vector<int>);

// Minimal stub template definitions for multi-qubit routines to resolve linker errors

template <int NumQubits>
qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub(Qureg, std::vector<int>, std::vector<int>) { return 0; }

template <int NumQubits>
qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg, std::vector<int>, std::vector<int>) { return 0; }

template <int NumQubits>
void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal*, Qureg, std::vector<int>) {}

template <int NumQubits>
void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(qreal*, Qureg, std::vector<int>) {}

template <int NumQubits>
void cpu_statevec_multiQubitProjector_sub(Qureg, std::vector<int>, std::vector<int>, qreal) {}

template <int NumQubits>
void cpu_densmatr_multiQubitProjector_sub(Qureg, std::vector<int>, std::vector<int>, qreal) {}

// Pauli tensor/gadget routines

template <int NumCtrls, int NumTargs>
void cpu_statevector_anyCtrlPauliTensorOrGadget_subA(Qureg, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, qcomp, qcomp) {}

template <int NumCtrls>
void cpu_statevector_anyCtrlPauliTensorOrGadget_subB(Qureg, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, qcomp, qcomp, qindex) {}

template <int NumCtrls>
void cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(Qureg, std::vector<int>, std::vector<int>, std::vector<int>, qcomp, qcomp) {}

// Partial trace

template <int NumTargs>
void cpu_densmatr_partialTrace_sub(Qureg, Qureg, std::vector<int>, std::vector<int>) {}

// All-target diagonal matrix routines

template <bool HasPower, bool UseRealPow>
void cpu_densmatr_allTargDiagMatr_sub(Qureg, FullStateDiagMatr, qcomp) {}

template <bool HasPower>
void cpu_statevec_allTargDiagMatr_sub(Qureg, FullStateDiagMatr, qcomp) {}

// Expec full state diag matrix

template <bool HasPower, bool UseRealPow>
qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub(Qureg, FullStateDiagMatr, qcomp) { return 0; }

template <bool HasPower, bool UseRealPow>
qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub(Qureg, FullStateDiagMatr, qcomp) { return 0; }

// Expec PauliStr
qcomp cpu_statevec_calcExpecPauliStr_subA(Qureg, std::vector<int>, std::vector<int>, std::vector<int>) { return 0; }
qcomp cpu_statevec_calcExpecPauliStr_subB(Qureg, std::vector<int>, std::vector<int>, std::vector<int>) { return 0; }
qcomp cpu_densmatr_calcExpecPauliStr_sub(Qureg, std::vector<int>, std::vector<int>, std::vector<int>) { return 0; }

// Other missing routines
qreal cpu_densmatr_calcHilbertSchmidtDistance_sub(Qureg, Qureg) { return 0; }
void cpu_densmatr_mixQureg_subB(double, Qureg, double, Qureg) {}
void cpu_densmatr_mixQureg_subC(double, Qureg, double) {}
void cpu_statevec_setQuregToSuperposition_sub(qcomp, Qureg, qcomp, Qureg, qcomp, Qureg) {}

// One/two qubit channels (minimal stubs)
void cpu_densmatr_oneQubitPauliChannel_subA(Qureg, int, double, double, double, double) {}
void cpu_densmatr_oneQubitPauliChannel_subB(Qureg, int, double, double, double, double) {}
void cpu_densmatr_oneQubitDepolarising_subA(Qureg, int, double) {}
void cpu_densmatr_oneQubitDepolarising_subB(Qureg, int, double) {}
void cpu_densmatr_oneQubitDamping_subA(Qureg, int, double) {}
void cpu_densmatr_oneQubitDamping_subB(Qureg, int, double) {}
void cpu_densmatr_oneQubitDamping_subC(Qureg, int, double) {}
void cpu_densmatr_oneQubitDamping_subD(Qureg, int, double) {}
void cpu_densmatr_twoQubitDepolarising_subA(Qureg, int, int, double) {}
void cpu_densmatr_twoQubitDepolarising_subB(Qureg, int, int, double) {}
void cpu_densmatr_twoQubitDepolarising_subC(Qureg, int, int, double) {}
void cpu_densmatr_twoQubitDepolarising_subD(Qureg, int, int, double) {}
void cpu_densmatr_twoQubitDepolarising_subE(Qureg, int, int, double) {}
void cpu_densmatr_twoQubitDepolarising_subF(Qureg, int, int, double) {}

// Explicit template instantiations for all required variants for missing routines
// These are needed to resolve linker errors for all combinations used in the codebase

// cpu_statevector_anyCtrlPauliTensorOrGadget_subA
#define INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(NCTRLS, NTARGS) \
    template void cpu_statevector_anyCtrlPauliTensorOrGadget_subA<NCTRLS, NTARGS>(Qureg, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, qcomp, qcomp);

INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(0,0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(1,0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(2,0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(3,0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(4,0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(5,0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(-1,0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(0,1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(1,1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(2,1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(3,1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(4,1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(5,1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(-1,1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(0,2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(1,2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(2,2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(3,2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(4,2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(5,2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(-1,2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(0,3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(1,3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(2,3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(3,3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(4,3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(5,3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(-1,3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(0,4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(1,4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(2,4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(3,4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(4,4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(5,4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(-1,4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(0,5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(1,5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(2,5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(3,5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(4,5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(5,5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(-1,5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(0,-1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(1,-1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(2,-1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(3,-1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(4,-1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(5,-1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA(-1,-1)
#undef INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBA

// cpu_statevector_anyCtrlPauliTensorOrGadget_subB
#define INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(NCTRLS) \
    template void cpu_statevector_anyCtrlPauliTensorOrGadget_subB<NCTRLS>(Qureg, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, qcomp, qcomp, qindex);
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB(-1)
#undef INSTANTIATE_CPU_STATEVECTOR_ANYCTRLPAULITENSORORGADGET_SUBB

// cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub
#define INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(NCTRLS) \
    template void cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub<NCTRLS>(Qureg, std::vector<int>, std::vector<int>, std::vector<int>, qcomp, qcomp);
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(0)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(1)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(2)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(3)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(4)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(5)
INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB(-1)
#undef INSTANTIATE_CPU_STATEVECTOR_ANYCTRLANYTARGZORPHASEGADGET_SUB

// cpu_densmatr_partialTrace_sub
#define INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(NTARGS) \
    template void cpu_densmatr_partialTrace_sub<NTARGS>(Qureg, Qureg, std::vector<int>, std::vector<int>);
INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(0)
INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(1)
INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(2)
INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(3)
INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(4)
INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(5)
INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB(-1)
#undef INSTANTIATE_CPU_DENSMATR_PARTIALTRACE_SUB

// cpu_densmatr_allTargDiagMatr_sub, cpu_statevec_allTargDiagMatr_sub
#define INSTANTIATE_CPU_DENSMATR_ALLTARGDIAGMATR_SUB(HP, URP) \
    template void cpu_densmatr_allTargDiagMatr_sub<HP, URP>(Qureg, FullStateDiagMatr, qcomp);
INSTANTIATE_CPU_DENSMATR_ALLTARGDIAGMATR_SUB(true, true)
INSTANTIATE_CPU_DENSMATR_ALLTARGDIAGMATR_SUB(true, false)
INSTANTIATE_CPU_DENSMATR_ALLTARGDIAGMATR_SUB(false, true)
INSTANTIATE_CPU_DENSMATR_ALLTARGDIAGMATR_SUB(false, false)
#undef INSTANTIATE_CPU_DENSMATR_ALLTARGDIAGMATR_SUB
#define INSTANTIATE_CPU_STATEVEC_ALLTARGDIAGMATR_SUB(HP) \
    template void cpu_statevec_allTargDiagMatr_sub<HP>(Qureg, FullStateDiagMatr, qcomp);
INSTANTIATE_CPU_STATEVEC_ALLTARGDIAGMATR_SUB(true)
INSTANTIATE_CPU_STATEVEC_ALLTARGDIAGMATR_SUB(false)
#undef INSTANTIATE_CPU_STATEVEC_ALLTARGDIAGMATR_SUB

// cpu_statevec_calcExpecFullStateDiagMatr_sub, cpu_densmatr_calcExpecFullStateDiagMatr_sub
#define INSTANTIATE_CPU_STATEVEC_CALCEXPECFULLSTATEDIAGMATR_SUB(HP, URP) \
    template qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub<HP, URP>(Qureg, FullStateDiagMatr, qcomp);
INSTANTIATE_CPU_STATEVEC_CALCEXPECFULLSTATEDIAGMATR_SUB(true, true)
INSTANTIATE_CPU_STATEVEC_CALCEXPECFULLSTATEDIAGMATR_SUB(true, false)
INSTANTIATE_CPU_STATEVEC_CALCEXPECFULLSTATEDIAGMATR_SUB(false, true)
INSTANTIATE_CPU_STATEVEC_CALCEXPECFULLSTATEDIAGMATR_SUB(false, false)
#undef INSTANTIATE_CPU_STATEVEC_CALCEXPECFULLSTATEDIAGMATR_SUB
#define INSTANTIATE_CPU_DENSMATR_CALCEXPECFULLSTATEDIAGMATR_SUB(HP, URP) \
    template qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub<HP, URP>(Qureg, FullStateDiagMatr, qcomp);
INSTANTIATE_CPU_DENSMATR_CALCEXPECFULLSTATEDIAGMATR_SUB(true, true)
INSTANTIATE_CPU_DENSMATR_CALCEXPECFULLSTATEDIAGMATR_SUB(true, false)
INSTANTIATE_CPU_DENSMATR_CALCEXPECFULLSTATEDIAGMATR_SUB(false, true)
INSTANTIATE_CPU_DENSMATR_CALCEXPECFULLSTATEDIAGMATR_SUB(false, false)
#undef INSTANTIATE_CPU_DENSMATR_CALCEXPECFULLSTATEDIAGMATR_SUB
