/** @file
 * CPU OpenMP-accelerated definitions of the main backend simulation routines,
 * as mirrored by gpu_subroutines.cpp, and called by accelerator.cpp. 
 * 
 * These 'hot-loop' functions use cpu_qcomp arithmetic operators, in lieu of qcomp
 * (i.e. std::complex) which has compiler-specific performance pitfalls. BEWARE
 * that passing a cpu_qcomp by-value to a function inside an OpenMP parallel region
 * can cause MSVC to crash during compilation, so be sure to pass by reference!
 * 
 * Some of these definitions are templated, defining multiple versions optimised 
 * (at compile-time) for handling different numbers of input qubits; such functions
 * are proceeded by macro INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(), to force the 
 * compilation of their needed versions within this translation unit for later linkage.
 * 
 * @author Tyson Jones
 * @author Oliver Brown (OpenMP 'if' clauses)
 * @author Luc Jaulmes (optimised initUniformState)
 * @author Richard Meister (helped patch on LLVM)
 * @author Kshitij Chhabra (patched v3 clauses with gcc9)
 * @author Ania (Anna) Brown (developed QuEST v1 logic)
 */

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
#include "quest/src/cpu/cpu_qcomp.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/cpu/cpu_subroutines.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_indices.hpp"

#include <vector>
#include <algorithm>

using std::vector;



/*
 * OPENMP QCOMP REDUCTION
 */


// As of Apr 2026, custom reductions are not supported by MSVC,
// even when using the LLVM OpenMP runtime (grr!). So all qcomp
// reductions within this file reduce the real and imaginary
// components (each, a qreal) separately



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

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* coeffs = getCpuQcompPtr(sum.coeffs);

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
        amps[n] = fast_getPauliStrSumElem(coeffs, sum.strings, sum.numTerms, r, c);
    }
}


void cpu_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* inCoeffs = getCpuQcompPtr(in.coeffs);
    cpu_qcomp* outElems = getCpuQcompPtr(out.cpuElems);

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
        outElems[n] = fast_getPauliStrSumElem(inCoeffs, in.strings, in.numTerms, i, i);
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
            // function is recursive and cannot be inlined. We are furthermore copying
            // qcomp directly, rather than cpu_qcomp, because no arithmetic is performed
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
            out.cpuElems[localInd] = callbackFunc(varValues.data()); // no cpu_qcomp casting needed
        }
    }
}



/*
 * COMMUNICATION BUFFER PACKING
 */


template <int NumQubits>
qindex cpu_statevec_packAmpsIntoBuffer(Qureg qureg, SmallView qubitInds, SmallView qubitStates) {

    assert_numQubitsMatchesQubitStatesAndTemplateParam(qubitInds.size(), qubitStates.size(), NumQubits);

    // use cpu_qcomp (in lieu of qcomp) even though no arithmetic happens below - just for consistency!
    cpu_qcomp* amps   = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

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
        buffer[offset + n] = amps[i];
    }

    // return the number of packed amps
    return numIts;
}


qindex cpu_statevec_packPairSummedAmpsIntoBuffer(Qureg qureg, int qubit1, int qubit2, int qubit3, int bit2) {
    
    assert_bufferPackerGivenIncreasingQubits(qubit1, qubit2, qubit3);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps   = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

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

        buffer[offset + n] = amps[i0b0] + amps[i1b1];
    }

    // return the number of packed amps
    return numIts;
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( qindex, cpu_statevec_packAmpsIntoBuffer, (Qureg, SmallView, SmallView) )



/* 
 * SWAPS
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlSwap_subA(Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ1, int targ2) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp (in lieu of qcomp) even though no arithmetic happens below - just for consistency!
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

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

        std::swap(amps[i01], amps[i10]);
    }
}


template <int NumCtrls>
void cpu_statevec_anyCtrlSwap_subB(Qureg qureg, SmallView ctrls, SmallView ctrlStates) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps   = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

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
        amps[i] = buffer[j];
    }
}


template <int NumCtrls>
void cpu_statevec_anyCtrlSwap_subC(Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ, int targState) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps   = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

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
        amps[i] = buffer[j];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlSwap_subA, (Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ1, int targ2) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlSwap_subB, (Qureg qureg, SmallView ctrls, SmallView ctrlStates) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlSwap_subC, (Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ, int targState) )



/*
 * ONE-TARGET DENSE MATRIX
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ, CompMatr1 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    auto elems = getCpuQcompsMatrix<2>(matr.elems); // MSVC requires explicit template param, bah!

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
        cpu_qcomp amp0 = amps[i0];
        cpu_qcomp amp1 = amps[i1];

        amps[i0] = elems[0][0]*amp0 + elems[0][1]*amp1;
        amps[i1] = elems[1][0]*amp0 + elems[1][1]*amp1;
    }
}


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, SmallView ctrls, SmallView ctrlStates, qcomp fac0, qcomp fac1) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps   = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);
    cpu_qcomp f0 = getCpuQcomp(fac0);
    cpu_qcomp f1 = getCpuQcomp(fac1);

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

        amps[i] = f0*amps[i] + f1*buffer[j];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, SmallView, SmallView, int, CompMatr1) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, SmallView, SmallView, qcomp, qcomp) )



/*
 * TWO-TARGET DENSE MATRIX
 */


template <int NumCtrls> 
void cpu_statevec_anyCtrlTwoTargDenseMatr_sub(Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ1, int targ2, CompMatr2 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    auto elems = getCpuQcompsMatrix<4>(matr.elems); // MSVC requires explicit template param, bah!

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
        cpu_qcomp amp00 = amps[i00];
        cpu_qcomp amp01 = amps[i01];
        cpu_qcomp amp10 = amps[i10];
        cpu_qcomp amp11 = amps[i11];

        // amps[i_n] = sum_j matr.elems[n][j] amp[i_n]
        amps[i00] = elems[0][0]*amp00 + elems[0][1]*amp01 + elems[0][2]*amp10 + elems[0][3]*amp11;
        amps[i01] = elems[1][0]*amp00 + elems[1][1]*amp01 + elems[1][2]*amp10 + elems[1][3]*amp11;
        amps[i10] = elems[2][0]*amp00 + elems[2][1]*amp01 + elems[2][2]*amp10 + elems[2][3]*amp11;
        amps[i11] = elems[3][0]*amp00 + elems[3][1]*amp01 + elems[3][2]*amp10 + elems[3][3]*amp11;
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlTwoTargDenseMatr_sub, (Qureg, SmallView, SmallView, int, int, CompMatr2) )



/*
 * MANY-TARGET DENSE MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj, bool ApplyTransp>
void cpu_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, SmallView ctrls, SmallView ctrlStates, SmallView targs, CompMatr matr) {
    
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.cpuElemsFlat);

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
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, targs, util_getConstantList(0,targs.size()));

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
        vector<cpu_qcomp> cache(numTargAmps);

        #pragma omp for
        for (qindex n=0; n<numIts; n++) {

            // i0 = nth local index where ctrls are active and targs are all zero
            qindex i0 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);

            // collect and cache all to-be-modified amps (loop might be unrolled)
            for (qindex j=0; j<numTargAmps; j++) {

                // i = nth local index where ctrls are active and targs form value j
                qindex i = setBits(i0, targs.data(), numTargBits, j); // loop may be unrolled
                cache[j] = amps[i];
            }

            // modify each amplitude (loop might be unrolled)
            for (qindex k=0; k<numTargAmps; k++) {

                // i = nth local index where ctrls are active and targs form value k
                qindex i = setBits(i0, targs.data(), numTargBits, k); // loop may be unrolled
                amps[i] = getCpuQcomp(0, 0);
            
                // loop may be unrolled
                for (qindex j=0; j<numTargAmps; j++) {

                    // matr.cpuElemsFlat[l] = matr.cpuElems[k][j] OR matr.cpuElems[j][k]
                    qindex l;
                    if constexpr (ApplyTransp)
                        l = fast_getMatrixFlatIndex(j, k, numTargAmps);
                    else
                        l = fast_getMatrixFlatIndex(k, j, numTargAmps);
                    cpu_qcomp elem = elems[l];

                    // optionally conjugate matrix elems on the fly to avoid pre-modifying heap structure
                    if constexpr (ApplyConj)
                        elem = conj(elem);

                    amps[i] += elem * cache[j];

                    /// @todo
                    /// qureg.cpuAmps[i] is being serially updated by only this thread,
                    /// so is a candidate for Kahan summation for improved numerical
                    /// stability. Explore whether this is time-free and worthwhile!
                    ///
                    /// BEWARE that Kahan summation may be incompatible with
                    /// the commutator tricks used in base_qcomp's (ancestor
                    /// of cpu_qcomp) arithmetic operator overloads. Check
                    /// base_qcomp.hpp before implementing compensation.
                }
            }
        }
    }
}


INSTANTIATE_TWO_BOOL_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDenseMatr_sub, (Qureg, SmallView, SmallView, SmallView, CompMatr) )



/*
 * ONE-TARG DIAGONAL MATRIX
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ, DiagMatr1 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.elems);

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
        amps[j] *= elems[b];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlOneTargDiagMatr_sub, (Qureg, SmallView, SmallView, int, DiagMatr1) )



/*
 * TWO-TARG DIAGONAL MATRIX
 */


template <int NumCtrls>
void cpu_statevec_anyCtrlTwoTargDiagMatr_sub(Qureg qureg, SmallView ctrls, SmallView ctrlStates, int targ1, int targ2, DiagMatr2 matr) {

    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.elems);

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
        amps[j] *= elems[k];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevec_anyCtrlTwoTargDiagMatr_sub, (Qureg, SmallView, SmallView, int, int, DiagMatr2) )



/*
 * ANY-TARG DIAGONAL MATRIX
 */


template <int NumCtrls, int NumTargs, bool ApplyConj, bool HasPower>
void cpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, SmallView ctrls, SmallView ctrlStates, SmallView targs, DiagMatr matr, qcomp exponent) {
    
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);
    assert_exponentMatchesTemplateParam(exponent, HasPower);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.cpuElems);
    cpu_qcomp expo   = getCpuQcomp(exponent);
    (void) expo; // silence when unused

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
        cpu_qcomp elem = elems[t];

        // decide whether to power and conj at compile-time, to avoid branching in hot-loop.
        // beware that pow(qcomp,qcomp) below gives notable error over pow(qreal,qreal) 
        // (by producing an unexpected non-zero imaginary component) when the base is real 
        // and negative, and the exponent is an integer. We tolerate this heightened error
        // because we have no reason to think matr is real (it's not constrained Hermitian).
        if constexpr (HasPower)
            elem = pow(elem, expo);

        // cautiously conjugate AFTER exponentiation, else we must also conj exponent
        if constexpr (ApplyConj)
            elem = conj(elem);

        amps[j] *= elem;
    }
}


INSTANTIATE_TWO_BOOL_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, SmallView, SmallView, SmallView, DiagMatr, qcomp) )


/// @todo
/// there is currently no density matrix version of anyCtrlAnyTargDiagMatr_sub();
/// instead, operations.cpp invokes the statevector version twice as it does for
/// dense matrices. This re-enumeration of the state however can be avoided since
/// the matrix is diagonal, as done below for cpu_densmatr_allTargDiagMatr_sub()



/*
 * ALL-TARGS DIAGONAL MATRIX
 */


template <bool HasPower>
void cpu_statevec_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    assert_quregAndFullStateDiagMatrHaveSameDistrib(qureg, matr);
    assert_exponentMatchesTemplateParam(exponent, HasPower);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.cpuElems);
    cpu_qcomp expo   = getCpuQcomp(exponent);
    (void) expo; // silence when unused

    // every iteration modifies one amp, using one element
    qindex numIts = qureg.numAmpsPerNode;

    #pragma omp parallel for if(qureg.isMultithreaded||qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        cpu_qcomp elem = elems[n];

        // compile-time decide if applying power to avoid in-loop branching.
        // beware that pow(qcomp,qcomp) below gives notable error over pow(qreal,qreal) 
        // (by producing an unexpected non-zero imaginary component) when the base is real 
        // and negative, and the exponent is an integer. We tolerate this heightened error
        // because we have no reason to think matr is real (it's not constrained Hermitian).
        if constexpr (HasPower)
            elem = pow(elem, expo);

        amps[n] *= elem;
    }
}


template <bool HasPower, bool ApplyLeft, bool ApplyRight, bool ConjRight>
void cpu_densmatr_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    // unlike other functions, this function handles all scenarios of...
    // - matr -> matr qureg conj(matr)
    // - matr -> matr qureg
    // - matr ->      qureg matr
    // and all of the above where matr is raised to a power. This is an
    // optimisation permitted by diagonality of matr, avoiding superfluous
    // re-enumeration of the state otherwise invoked by operations.cpp

    assert_exponentMatchesTemplateParam(exponent, HasPower);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.cpuElems);
    cpu_qcomp expo   = getCpuQcomp(exponent);
    (void) expo; // silence when unused

    // every iteration modifies one qureg amp, using one or two matr elements
    qindex numIts = qureg.numAmpsPerNode;

    #pragma omp parallel for if(qureg.isMultithreaded||matr.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // the nth local amplitude will be multiplied by fac
        cpu_qcomp fac = getCpuQcomp(1,0);

        // update fac to effect rho -> (matr * rho) or (matr^exponent * rho)
        if constexpr (ApplyLeft) {

            // i = global row of nth local amp
            qindex i = fast_getQuregGlobalRowFromFlatIndex(n, matr.numElems);
            cpu_qcomp term = elems[i];

            // compile-time decide if applying power to avoid in-loop branching...
            // (beware that complex pow() is numerically unstable as detailed below)
            if constexpr (HasPower)
                term = pow(term, expo);

            fac = term;
        }

        // update fac to additional include rho -> (rho * matr) or 
        // (rho * conj(matr)), or the same exponentiated
        if constexpr (ApplyRight) {

            // m = global index corresponding to n
            qindex m = concatenateBits(qureg.rank, n, qureg.logNumAmpsPerNode);

            // j = global column corresponding to n
            qindex j = fast_getQuregGlobalColFromFlatIndex(m, matr.numElems);
            cpu_qcomp term = elems[j];

            // beware that pow(qcomp,qcomp) below gives notable error over pow(qreal,qreal) 
            // (by producing an unexpected non-zero imaginary component) when the base is real 
            // and negative, and the exponent is an integer. We tolerate this heightened error
            // because we have no reason to think matr is real (it's not constrained Hermitian).
            if constexpr (HasPower)
                term = pow(term, expo);

            // conj strictly after pow, to effect conj(matr^exponent)
            if constexpr (ConjRight)
                term = conj(term);

            fac *= term;
        }

        amps[n] *= fac;
    }
}


template void cpu_statevec_allTargDiagMatr_sub<true> (Qureg, FullStateDiagMatr, qcomp);
template void cpu_statevec_allTargDiagMatr_sub<false>(Qureg, FullStateDiagMatr, qcomp);

template void cpu_densmatr_allTargDiagMatr_sub<false, true,  true,  true>  (Qureg, FullStateDiagMatr, qcomp); // matr qureg conj(matr)
template void cpu_densmatr_allTargDiagMatr_sub<false, true,  false, false> (Qureg, FullStateDiagMatr, qcomp); // matr qureg
template void cpu_densmatr_allTargDiagMatr_sub<false, false, true,  false> (Qureg, FullStateDiagMatr, qcomp); //      qureg matr
template void cpu_densmatr_allTargDiagMatr_sub<true,  true,  true,  true>  (Qureg, FullStateDiagMatr, qcomp); // matr^P qureg conj(matr^P)
template void cpu_densmatr_allTargDiagMatr_sub<true,  true,  false, false> (Qureg, FullStateDiagMatr, qcomp); // matr^P qureg
template void cpu_densmatr_allTargDiagMatr_sub<true,  false, true,  false> (Qureg, FullStateDiagMatr, qcomp); //      qureg matr^P



/*
 * PAULI TENSOR AND GADGET
 */


template <int NumTargs>
INLINE void applyPauliUponAmpPair(
    cpu_qcomp* amps, qindex& v, qindex& i0, int* indXY, int& numXY, 
    qindex& maskXY, qindex& maskYZ, cpu_qcomp& ampFac, cpu_qcomp& pairAmpFac
) {
    // This is a subroutine of cpu_statevector_anyCtrlPauliTensorOrGadget_subA() below
    // called in a hot-loop (hence it is here inlined) which exists because the caller
    // chooses one of two possible OpenMP parallelisation granularities. All args are
    // pass-by-reference for performance, and because passing the cpu_qcomp types by-
    // value causes a stack overflow during compilation with MSVC with OpenMP enabled;
    // but only at double and quad precision (single is fine), and only when Catch2 is
    // also being compiled (through the tests)... Hours of my life forever lost!

    // remind compiler when NumTargs is compile-time to unroll loop in setBits()
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, numXY);

    // iA = nth local index where targs have value m, iB = (last - nth) such index
    qindex iA = setBits(i0, indXY, numTargBits, v);
    qindex iB = flipBits(iA, maskXY);

    // sign of amps due to Y and Z (excludes Y's i factor which is inside pairAmpFac)
    int signA = fast_getPlusOrMinusMaskedBitParity(iA, maskYZ);
    int signB = fast_getPlusOrMinusMaskedBitParity(iB, maskYZ);

    // mix or swap scaled amp pair (where pairAmpFac includes Y's i factor)
    cpu_qcomp ampA = amps[iA];
    cpu_qcomp ampB = amps[iB];
    amps[iA] = (ampFac * ampA) + (pairAmpFac * signB * ampB);
    amps[iB] = (ampFac * ampB) + (pairAmpFac * signA * ampA);
}


template <int NumCtrls, int NumTargs>
void cpu_statevector_anyCtrlPauliTensorOrGadget_subA(
    Qureg qureg, SmallView ctrls, SmallView ctrlStates, 
    SmallView x, SmallView y, SmallView z, qcomp ampFac, qcomp pairAmpFac
) {
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);
    assert_numTargsMatchesTemplateParam(x.size() + y.size(), NumTargs);

    // we will scale pairAmp below by i^numY, so that each amp need only choose the +-1 sign
    pairAmpFac *= util_getPowerOfI(y.size());

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp f0 = getCpuQcomp(ampFac);
    cpu_qcomp f1 = getCpuQcomp(pairAmpFac);
    
    // only X and Y count as targets
    auto sortedTargsXY = util_getSorted(util_getConcatenated(x, y));

    // prepare a mask which yields ctrls in specified state, and X-Y targs in all-zero
    auto sortedQubits   = util_getSorted(ctrls, sortedTargsXY);
    auto qubitStateMask = util_getBitMask(ctrls, ctrlStates, sortedTargsXY, util_getConstantList(0, sortedTargsXY.size()));

    // prepare masks for extracting Pauli parities
    auto maskXY = util_getBitMask(sortedTargsXY);
    auto maskYZ = util_getBitMask(util_getConcatenated(y, z));

    // use template params to compile-time unroll loops in insertBits() and inner-loop below
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());
    SET_VAR_AT_COMPILE_TIME(int, numTargBits, NumTargs, sortedTargsXY.size());
    int numQubitBits = numCtrlBits + numTargBits;

    // each outer iteration handles all assignments of the target qubits, and each ctrl halves the outer iterations
    qindex numOuterIts = qureg.numAmpsPerNode / powerOf2(numCtrlBits + numTargBits);

    // each inner iteration modifies 2 amplitudes (may be compile-time sized) 
    qindex numInnerIts = powerOf2(numTargBits) / 2; // divides evenly

    // we must choose whether to parallelise the outer or inner loop, since
    // their relative sizes change exponentially depending on numTargAmps.
    // For example, when x+y=O(1), the outer loop is O(2^N) and is parallelised 
    // while the inner loop is O(1). In contrast, when x+y=O(N), the outer loop 
    // is O(1) and left serial while the inner loop is O(2^N) and parallelised.
    // While outer-loop parallelisation can be simply decided with a runtime 
    // if() clause in the pragma, the inner-loop parallelisation must be chosen
    // at compile-time to avoid a huge per-out-loop slowdown when x+y~1. 
    // We will opt to parallelise the outer loop, leaving the inner loop serial,
    // whenever each thread has at least 1 iteration for itself. And of course
    // we serialise both inner and outer loops when qureg multithreading is off.

    if (!qureg.isMultithreaded || numOuterIts >= cpu_getAvailableNumThreads()) {
    
        // parallel
        #pragma omp parallel for if(qureg.isMultithreaded)
        for (qindex n=0; n<numOuterIts; n++) {

            // i0 = nth local index where ctrls are active and targs are all zero
            qindex i0 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);

            // serial
            for (qindex v=0; v<numInnerIts; v++)
                applyPauliUponAmpPair<NumTargs>(
                    amps, v, i0, sortedTargsXY.data(), numTargBits, maskXY, maskYZ, f0, f1);
        }

    } else {

        // serial
        for (qindex n=0; n<numOuterIts; n++) {

            // i0 = nth local index where ctrls are active and targs are all zero
            qindex i0 = insertBitsWithMaskedValues(n, sortedQubits.data(), numQubitBits, qubitStateMask);

            // parallel
            #pragma omp parallel for
            for (qindex v=0; v<numInnerIts; v++)
                applyPauliUponAmpPair<NumTargs>(
                    amps, v, i0, sortedTargsXY.data(), numTargBits, maskXY, maskYZ, f0, f1);
        }
    }
}


template <int NumCtrls>
void cpu_statevector_anyCtrlPauliTensorOrGadget_subB(
    Qureg qureg, SmallView ctrls, SmallView ctrlStates,
    SmallView x, SmallView y, SmallView z, qcomp ampFac, qcomp pairAmpFac, qindex bufferMaskXY
) {
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // we will scale pairAmp by i^numY, so that each amp need only choose the +-1 sign
    pairAmpFac *= util_getPowerOfI(y.size());

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps   = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);
    cpu_qcomp f0 = getCpuQcomp(ampFac);
    cpu_qcomp f1 = getCpuQcomp(pairAmpFac);
    
    // each control qubit halves the needed iterations
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    auto maskXY = util_getBitMask(util_getConcatenated(x, y));
    auto maskYZ = util_getBitMask(util_getConcatenated(y, z));

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // j = buffer index of amp to be mixed with i
        qindex j = flipBits(n, bufferMaskXY) + offset;

        // k = local index of j-th buffer amplitude in its original node
        qindex k = flipBits(i, maskXY);
        int sign = fast_getPlusOrMinusMaskedBitParity(k, maskYZ);

        amps[i] *= f0;
        amps[i] += f1 * sign * buffer[j]; // f1=pairAmpFac includes Y's i factors
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, cpu_statevector_anyCtrlPauliTensorOrGadget_subA, (Qureg, SmallView, SmallView, SmallView, SmallView, SmallView, qcomp, qcomp) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevector_anyCtrlPauliTensorOrGadget_subB, (Qureg, SmallView, SmallView, SmallView, SmallView, SmallView, qcomp, qcomp, qindex) )



/*
 * PHASE TENSOR AND GADGET
 */


template <int NumCtrls>
void cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    Qureg qureg, SmallView ctrls, SmallView ctrlStates, SmallView targs, 
    qcomp fac0, qcomp fac1
) {
    assert_numCtrlsMatchesNumCtrlStatesAndTemplateParam(ctrls.size(), ctrlStates.size(), NumCtrls);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp facs[] = {getCpuQcomp(fac0), getCpuQcomp(fac1)};

    // each control qubit halves the needed iterations, each of which modifies 1 amp
    qindex numIts = qureg.numAmpsPerNode / powerOf2(ctrls.size());

    auto sortedCtrls   = util_getSorted(ctrls);
    auto ctrlStateMask = util_getBitMask(ctrls, ctrlStates);
    auto targMask      = util_getBitMask(targs);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numCtrlBits, NumCtrls, ctrls.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ctrl bits are in specified states
        qindex i = insertBitsWithMaskedValues(n, sortedCtrls.data(), numCtrlBits, ctrlStateMask);

        // apply phase to amp depending on parity of targets
        int p = getBitMaskParity(i & targMask);
        amps[i] *= facs[p];
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub, (Qureg, SmallView, SmallView, SmallView, qcomp, qcomp) )



/*
 * QUREG COMBINATION
 */


template <int NumQuregs>
void cpu_statevec_setQuregToWeightedSum_sub(Qureg outQureg, vector<qcomp> coeffs, vector<Qureg> inQuregs) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* outAmps = getCpuQcompPtr(outQureg.cpuAmps);
    cpu_qcomp* inFacs = getCpuQcompPtr(coeffs.data());

    qindex numIts = outQureg.numAmpsPerNode;

    // use template param to compile-time unroll inner loop below
    SET_VAR_AT_COMPILE_TIME(int, numQuregs, NumQuregs, inQuregs.size());

    #pragma omp parallel for if(outQureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // unrolled when inQuregs.size() <= 5
        cpu_qcomp amp = getCpuQcomp(0, 0);
        for (int q=0; q<numQuregs; q++)
            amp += inFacs[q] * getCpuQcompPtr(inQuregs[q].cpuAmps)[n];

        // must not modify cpuAmps[n] before computing the amp since
        // outQureg can legally appear among inQuregs
        outAmps[n] = amp;
    }
}


void cpu_densmatr_mixQureg_subA(qreal outProb, Qureg outQureg, qreal inProb, Qureg inDensMatr) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* outAmps = getCpuQcompPtr(outQureg.cpuAmps);
    cpu_qcomp* inAmps = getCpuQcompPtr(inDensMatr.cpuAmps);
    
    qindex numIts = outQureg.numAmpsPerNode;

    #pragma omp parallel for if(outQureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++)
        outAmps[n] = (outProb * outAmps[n]) + (inProb * inAmps[n]);
}


void cpu_densmatr_mixQureg_subB(qreal outProb, Qureg outQureg, qreal inProb, Qureg inStateVec) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* outAmps = getCpuQcompPtr(outQureg.cpuAmps);
    cpu_qcomp* inAmps = getCpuQcompPtr(inStateVec.cpuAmps);
    
    qindex numIts = outQureg.numAmpsPerNode;
    qindex dim = inStateVec.numAmps;

    #pragma omp parallel for if(outQureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // (i,j) = row & column of outQureg corresponding to n
        qindex i = fast_getQuregGlobalRowFromFlatIndex(n, dim);
        qindex j = fast_getQuregGlobalColFromFlatIndex(n, dim);

        outAmps[n] = (outProb * outAmps[n]) + (inProb * inAmps[i] * conj(inAmps[j]));
    }
}


void cpu_densmatr_mixQureg_subC(qreal outProb, Qureg outQureg, qreal inProb) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* outAmps = getCpuQcompPtr(outQureg.cpuAmps);
    cpu_qcomp* inAmps = getCpuQcompPtr(outQureg.cpuCommBuffer);

    // received inQureg's entire statevector amplitudes into every node's buffer
    qindex numIts = outQureg.numAmpsPerNode;
    qindex dim = powerOf2(outQureg.numQubits);

    #pragma omp parallel for if(outQureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // m = global index of local index n
        qindex m = concatenateBits(outQureg.rank, n, outQureg.logNumAmpsPerNode);

        // (i,j) = global row & column of outQureg corresponding to n
        qindex i = fast_getQuregGlobalRowFromFlatIndex(m, dim);
        qindex j = fast_getQuregGlobalColFromFlatIndex(m, dim);

        outAmps[n] = (outProb * outAmps[n]) + (inProb * inAmps[i] * conj(inAmps[j]));
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_QUREGS( void, cpu_statevec_setQuregToWeightedSum_sub, (Qureg, vector<qcomp>, vector<Qureg>) )



/*
 * ONE-QUBIT DEPHASING
 */


void cpu_densmatr_oneQubitDephasing_subA(Qureg qureg, int ketQubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // half of all local amps are scaled, and each iteration modifies two
    qindex numIts = qureg.numAmpsPerNode / 4;

    // loop constants
    qreal fac = util_getOneQubitDephasingFactor(prob);
    int braQubit = util_getBraQubit(ketQubit, qureg);

    /// @todo
    /// this enumeration order is suboptimal and seems unnecessary in this simple two
    /// bit scenario, where we are modifying but not at all mixing two strided and
    /// potentially very-distant amplitudes. It is of course trivial to split this
    /// into two separate loops accessing monotonically increasing indices, although
    /// we then pay double the caching costs when ketQubit is low-index. We can also
    /// turn it into two nested loops to force monotonically increasing index access,
    /// but then the parallelisation is not optimal when ketQubit is high-index. Experiment
    /// with what's fastest and replace below or delete this comment!

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i01 = nth local index of |*0*><*1*|
        qindex i01 = insertTwoBits(n, braQubit, 0, ketQubit, 1);
        qindex i10 = insertTwoBits(n, braQubit, 1, ketQubit, 0);

        amps[i01] *= fac;
        amps[i10] *= fac;
    }
}


void cpu_densmatr_oneQubitDephasing_subB(Qureg qureg, int ketQubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    
    // half of all local amps are scaled
    qindex numIts = qureg.numAmpsPerNode / 2;

    // loop constants
    qreal fac = util_getOneQubitDephasingFactor(prob);
    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    
    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where bra-qubit differs from ket-qubit
        qindex i = insertBit(n, ketQubit, ! braBit);

        amps[i] *= fac;
    }
}



/*
 * TWO-QUBIT DEPHASING
 */


void cpu_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubitA, int qubitB, qreal prob) {

    /// @todo 
    /// test whether use of subB has identical performance, or whether changing i=n below
    /// non-negligibly accelerates the routine; if so, make a templated inner func.

    // the rank-agnostic version is identical to the subB algorithm below, because the
    // queried bits of the global index i below will always be in the suffix substate.
    // We still define separate _subA and _subB routines because they differ in cuQuantum.
    cpu_densmatr_twoQubitDephasing_subB(qureg, qubitA, qubitB, prob);
}


void cpu_densmatr_twoQubitDephasing_subB(Qureg qureg, int ketQubitA, int ketQubitB, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

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
        amps[n] *= 1 + (flag * term);
    }
}



/*
 * ONE-QUBIT DEPOLARISING
 */


void cpu_densmatr_oneQubitDepolarising_subA(Qureg qureg, int ketQubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // all amps are modified, and each iteration modifies 4
    qindex numIts = qureg.numAmpsPerNode / 4;

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
        cpu_qcomp amp00 = amps[i00];
        cpu_qcomp amp11 = amps[i11];
        amps[i00] = (facAA * amp00) + (facBB * amp11);
        amps[i01] *= facAB;
        amps[i10] *= facAB;
        amps[i11] = (facAA * amp11) + (facBB * amp00);
    }
}


void cpu_densmatr_oneQubitDepolarising_subB(Qureg qureg, int ketQubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

    // all amps are modified, and each iteration modifies 2
    qindex numIts = qureg.numAmpsPerNode / 2;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    auto factors = util_getOneQubitDepolarisingFactors(prob);
    auto facAA = factors.c1;
    auto facBB = factors.c2;
    auto facAB = factors.c3;

    /// @todo
    /// each iteration below modifies 2 independent amps without mixing,
    /// which we can trivially split into two loops which may improve
    /// per-iteration caching performance; test if this outweights the 
    /// cost of re-iteration

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // iAA = nth local index where ket qubit agrees with bra qubit
        qindex iAA = insertBit(n, ketQubit, braBit);

        // jBB = buffer index of amp iAA where ket and bra qubits are flipped
        qindex jBB = n + offset;

        // iAB = nth local index where ket qubit disagrees with bra qubit
        qindex iAB = insertBit(n, ketQubit, ! braBit);

        amps[iAA] *= facAA;
        amps[iAA] += facBB * buffer[jBB];
        amps[iAB] *= facAB;
    }
}



/*
 * TWO-QUBIT DEPOLARISING
 */


void cpu_densmatr_twoQubitDepolarising_subA(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // all amps are scaled (although 1/16 of them will be unchanged)
    qindex numIts  = qureg.numAmpsPerNode;

    // bra-qubits corresponding to ket-qubits
    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braQb2 = util_getBraQubit(ketQb2, qureg);

    auto c3 = util_getTwoQubitDepolarisingFactors(prob).c3;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // determine whether to modify amp
        bool flag1 = getBit(n, ketQb1) == getBit(n, braQb1);
        bool flag2 = getBit(n, ketQb2) == getBit(n, braQb2);
        bool mod = !(flag1 & flag2);

        // multiply 15/16-th of all amps by (1 + c3)
        amps[n] *= 1 + c3 * mod;
    }
}


void cpu_densmatr_twoQubitDepolarising_subB(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // one quarter of amps will be modified, and four are mixed each iteration
    qindex numIts = qureg.numAmpsPerNode / 16;

    // bra-qubits corresponding to ket-qubits
    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braQb2 = util_getBraQubit(ketQb2, qureg);

    // factors used in amp -> c1*amp + c2(sum of other three amps)
    auto factors = util_getTwoQubitDepolarisingFactors(prob);
    auto c1 = factors.c1;
    auto c2 = factors.c2;

    // below, we compute term = (sum of all four amps) for brevity, so adjust c1
    c1 -= c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {
    
        // i0000 = nth local index where all bra = ket = 00
        qindex i0000 = insertFourZeroBits(n, braQb2, braQb1, ketQb2, ketQb1);
        qindex i0101 = flipTwoBits(i0000, braQb1, ketQb1);
        qindex i1010 = flipTwoBits(i0000, braQb2, ketQb2);
        qindex i1111 = flipTwoBits(i0101, braQb2, ketQb2);
        
        // mix 1/16 of all amps in groups of 4
        cpu_qcomp term = amps[i0000] + amps[i0101] + amps[i1010] + amps[i1111];

        amps[i0000] = c1*amps[i0000] + c2*term;
        amps[i0101] = c1*amps[i0101] + c2*term;
        amps[i1010] = c1*amps[i1010] + c2*term;
        amps[i1111] = c1*amps[i1111] + c2*term;
    }
}


void cpu_densmatr_twoQubitDepolarising_subC(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // scale 25% of amps but iterate all
    qindex numIts = qureg.numAmpsPerNode;

    int braQb1 = util_getBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    auto c3 = util_getTwoQubitDepolarisingFactors(prob).c3;

    /// @todo
    /// are we really inefficiently enumerating all amps and applying a non-unity
    /// factor to only 25%?! Is this because we do not know braBit2 and ergo 
    /// cannot be sure a direct enumeration is accessing indicies in a monotonically
    /// increasing order? Can that really outweigh a 3x slowdown?! Test and fix!

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // decide whether or not to modify nth local
        bool flag1 = getBit(n, ketQb1) == getBit(n, braQb1); 
        bool flag2 = getBit(n, ketQb2) == braBit2;
        bool mod = !(flag1 & flag2);

        // scale amp by 1 or (1 + c3)
        amps[n] *= 1 + c3 * mod;
    }
}


void cpu_densmatr_twoQubitDepolarising_subD(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

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
        cpu_qcomp amp0b0 = amps[i0b0];
        cpu_qcomp amp1b1 = amps[i1b1];

        amps[i0b0] = c1*amp0b0 + c2*(amp1b1 + buffer[j]);
        amps[i1b1] = c1*amp1b1 + c2*(amp0b0 + buffer[j]);
    }
}


void cpu_densmatr_twoQubitDepolarising_subE(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // all amplitudes are scaled; 25% by c1 and 75% by 1 + c3
    qindex numIts = qureg.numAmpsPerNode;

    auto factors = util_getTwoQubitDepolarisingFactors(prob);
    qreal fac0 = 1 + factors.c3;
    qreal fac1 = factors.c1 - fac0;

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {
        
        // choose factor by which to sacle amp
        bool same1 = getBit(n, ketQb1) == braBit1; 
        bool same2 = getBit(n, ketQb2) == braBit2;
        bool flag = same1 & same2;

        // scale amp by c1 or (1+c3)
        amps[n] *= fac1 * flag + fac0;;
    }
}


void cpu_densmatr_twoQubitDepolarising_subF(Qureg qureg, int ketQb1, int ketQb2, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

    // modify 25% of local amps, one per iteration
    qindex numIts = qureg.numAmpsPerNode / 4;

    // received amplitudes may begin at an arbitrary offset in the buffer
    qindex offset = getBufferRecvInd();

    int braBit1 = util_getRankBitOfBraQubit(ketQb1, qureg);
    int braBit2 = util_getRankBitOfBraQubit(ketQb2, qureg);

    auto c2 = util_getTwoQubitDepolarisingFactors(prob).c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where suffix ket qubits equal prefix bra qubits
        qindex i = insertTwoBits(n, ketQb2, braBit2, ketQb1, braBit1);

        // j = nth received amp in buffer
        qindex j = n + offset;

        // mix local amp with received buffer amp
        amps[i] += c2 * buffer[j];
    }
}



/*
 * PAULI CHANNEL
 */


void cpu_densmatr_oneQubitPauliChannel_subA(Qureg qureg, int ketQubit, qreal pI, qreal pX, qreal pY, qreal pZ) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // all amps are modified, and each iteration modifies 4
    qindex numIts = qureg.numAmpsPerNode / 4;

    int braQubit = util_getBraQubit(ketQubit, qureg);

    auto factors = util_getOneQubitPauliChannelFactors(pI, pX, pY, pZ);
    auto facAA = factors.c1;
    auto facBB = factors.c2;
    auto facAB = factors.c3;
    auto facBA = factors.c4;

    /// @todo
    /// each iteration modifies 4 amps in two separable mixed pairs, which may
    /// lead to sub-optimal caching. Iterating twice and modifying a single pair
    /// might lead to better performance, though note the stride from i00 to i11
    /// will always be adverserially large. Test this!

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i00 = nth local index where both qubits are 0
        qindex i00 = insertTwoBits(n, braQubit, 0, ketQubit, 0);
        qindex i01 = flipBit(i00, ketQubit);
        qindex i10 = flipBit(i00, braQubit);
        qindex i11 = flipBit(i01, braQubit);

        // modify 4 amps in 2 separable pairs
        cpu_qcomp amp00 = amps[i00];
        cpu_qcomp amp01 = amps[i01];
        cpu_qcomp amp10 = amps[i10];
        cpu_qcomp amp11 = amps[i11];

        amps[i00] = (facAA * amp00) + (facBB * amp11);
        amps[i01] = (facAB * amp01) + (facBA * amp10);
        amps[i10] = (facAB * amp10) + (facBA * amp01);
        amps[i11] = (facAA * amp11) + (facBB * amp00);
    }
}


void cpu_densmatr_oneQubitPauliChannel_subB(Qureg qureg, int ketQubit, qreal pI, qreal pX, qreal pY, qreal pZ) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

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

    /// @todo
    /// each iteration below modifies 2 independent amps without mixing,
    /// which we can trivially split into two loops which may improve
    /// per-iteration caching performance; test if this outweights the 
    /// cost of re-iteration

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // iAA = nth local index where ket qubit agrees with bra, i.e. |.A.><.A.|
        qindex iAA = insertBit(n, ketQubit, braBit);
        qindex iAB = flipBit(iAA, ketQubit);

        // jBB = buffer index of amp to be mixed with iAA's amp, i.e. |.B.><.B.|
        qindex jBB = iAB + offset;
        qindex jBA = iAA + offset;

        // mix each local amp with a received buffer amp
        amps[iAA] *= facAA;
        amps[iAA] += facBB * buffer[jBB];

        amps[iAB] *= facAB;
        amps[iAB] += facBA * buffer[jBA];
    }
}



/*
 * AMPLITUDE DAMPING CHANNEL
 */


void cpu_densmatr_oneQubitDamping_subA(Qureg qureg, int ketQubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

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
        amps[i00] += prob * amps[i11];

        // scale other amps
        amps[i01] *= c1;
        amps[i10] *= c1;
        amps[i11] *= c2;
    }
}


void cpu_densmatr_oneQubitDamping_subB(Qureg qureg, int qubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // half of all local amps are scaled
    qindex numIts = qureg.numAmpsPerNode / 2;

    auto c2 = util_getOneQubitDampingFactors(prob).c2;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where qubit=1
        qindex i= insertBit(n, qubit, 1);
        amps[i] *= c2;
    }
}


void cpu_densmatr_oneQubitDamping_subC(Qureg qureg, int ketQubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // half of all local amps are scaled
    qindex numIts = qureg.numAmpsPerNode / 2;

    int braBit = util_getRankBitOfBraQubit(ketQubit, qureg);
    auto c1 = util_getOneQubitDampingFactors(prob).c1;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where ket differs from bra
        qindex i = insertBit(n, ketQubit, ! braBit);
        amps[i] *= c1;
    }
}


void cpu_densmatr_oneQubitDamping_subD(Qureg qureg, int qubit, qreal prob) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

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

        amps[i] += prob * buffer[j];
    }
}



/*
 * PARTIAL TRACE
 */


template <int NumTargs>
void cpu_densmatr_partialTrace_sub(Qureg inQureg, Qureg outQureg, SmallView targs, SmallView pairTargs) {

    assert_numTargsMatchesTemplateParam(targs.size(), NumTargs);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* inAmps = getCpuQcompPtr(inQureg.cpuAmps);
    cpu_qcomp* outAmps = getCpuQcompPtr(outQureg.cpuAmps);

    // each outer iteration sets one element of outQureg
    qindex numOuterIts = outQureg.numAmpsPerNode;

    // targs and allTargs are sorted, but pairTargs is arbitrarily ordered (though corresponding targs)
    auto allTargsSorted = util_getSorted(targs, pairTargs);

    // use template param to compile-time unroll below loops
    SET_VAR_AT_COMPILE_TIME(int, numTargPairs, NumTargs, targs.size());
    
    // may be inferred at compile-time
    int numAllTargs = 2*numTargPairs;
    qindex numInnerIts = powerOf2(numTargPairs);

    /// @todo
    /// note our parallelisation of only the outer-loop assumes that the number of 
    /// amps in outQureg equals or exceeds the number of threads. Ergo tracing out 
    /// all but very few qubits will leave threads idle; when only a single qubit
    /// remains, the below code would be serial. In that scenario, we should
    /// parallelise the inner loop, or preclude this scenario in validation.

    // consult inQureg for multithreading, because total iters = inQureg dim
    #pragma omp parallel for if(inQureg.isMultithreaded)
    for (qindex n=0; n<numOuterIts; n++) {

        // k = nth local index of inQureg where all targs and pairs are zero
        qindex k = insertBits(n, allTargsSorted.data(), numAllTargs, 0); // loop may be unrolled

        // each outQureg amp results from summing 2^targs inQureg amps
        cpu_qcomp outAmp = getCpuQcomp(0,0);

        // loop may be unrolled
        for (qindex j=0; j<numInnerIts; j++) {

            // i = nth local index of inQureg where targs=j and pairTargs=j
            qindex i = k;
            i = setBits(i, targs    .data(), numTargPairs, j); // loop may be unrolled
            i = setBits(i, pairTargs.data(), numTargPairs, j); // loop may be unrolled

            outAmp += inAmps[i];
        }

        outAmps[n] = outAmp;
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, cpu_densmatr_partialTrace_sub, (Qureg, Qureg, SmallView, SmallView) )



/*
 * PROBABILITIES
 */


qreal cpu_statevec_calcTotalProb_sub(Qureg qureg) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    /// @todo
    /// check whether OpenMP is performing a numerically stable
    /// reduction, e.g. via 'parallel summation', to avoid the
    /// otherwise catastrophic cancellatin errors. If not, we
    /// should implement Kahan summation (parallelise by 
    /// having each thread Kahan-sum independently before a
    /// final serial combination). This invokes several times
    /// as many arithmetic operations (4x?) but we are anyway
    /// memory-bandwidth bound
    ///
    /// BEWARE that Kahan summation may be incompatible with
    /// the commutator tricks used in base_qcomp's (ancestor
    /// of cpu_qcomp) arithmetic operator overloads. Check
    /// base_qcomp.hpp before implementing compensation.

    qreal prob = 0;

    // every amp, iterated independently, contributes to the probability
    qindex numIts = qureg.numAmpsPerNode;

    #pragma omp parallel for reduction(+:prob) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++)
        prob += norm(amps[n]);

    return prob;
}


qreal cpu_densmatr_calcTotalProb_sub(Qureg qureg) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    /// @todo
    /// check whether OpenMP is performing a numerically stable
    /// reduction, e.g. via 'parallel summation', to avoid the
    /// otherwise catastrophic cancellatin errors. If not, we
    /// should implement Kahan summation (parallelise by 
    /// having each thread Kahan-sum independently before a
    /// final serial combination). This invokes several times
    /// as many arithmetic operations (4x?) but we are anyway
    /// memory-bandwidth bound.
    ///
    /// BEWARE that Kahan summation may be incompatible with
    /// the commutator tricks used in base_qcomp's (ancestor
    /// of cpu_qcomp) arithmetic operator overloads. Check
    /// base_qcomp.hpp before implementing compensation.

    qreal prob = 0;

    // iterate each column, of which one amp (the diagonal) contributes
    qindex numIts = powerOf2(qureg.logNumColsPerNode);
    qindex numAmpsPerCol = powerOf2(qureg.numQubits);
    qindex firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);

    #pragma omp parallel for reduction(+:prob) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = local index of nth local diagonal element
        qindex i = fast_getQuregLocalIndexOfDiagonalAmp(n, firstDiagInd, numAmpsPerCol);
        prob += real(amps[i]);
    }

    return prob;
}


template <int NumQubits>
qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, SmallView qubits, SmallView outcomes) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    qreal prob = 0;

    // each iteration visits one amp per 2^qubits.size() amps
    // (>=1 since all qubits are in suffix, so qubits.size() <= suffix size) 
    qindex numIts = qureg.numAmpsPerNode / powerOf2(qubits.size());

    auto sortedQubits = util_getSorted(qubits); // all in suffix
    auto qubitStateMask = util_getBitMask(qubits, outcomes);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, qubits.size());

    #pragma omp parallel for reduction(+:prob) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = nth local index where qubits are in the specified outcome state
        qindex i = insertBitsWithMaskedValues(n, sortedQubits.data(), numBits, qubitStateMask);

        prob += norm(amps[i]);
    }

    return prob;
}


template <int NumQubits>
qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg qureg, SmallView qubits, SmallView outcomes) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // note that qubits are only ket qubits for which the corresponding bra-qubit is in the suffix;
    // this function is not invoked upon nodes where prefix bra-qubits do not correspond to given outcomes

    qreal prob = 0;

    // each iteration visits one relevant diagonal amp (= one column)
    qindex numIts = powerOf2(qureg.logNumColsPerNode - qubits.size());
    qindex numAmpsPerCol = powerOf2(qureg.numQubits);
    qindex firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);

    auto sortedQubits = util_getSorted(qubits); // all in suffix, with corresponding bra's all in suffix
    auto qubitStateMask = util_getBitMask(qubits, outcomes);

    // use template param to compile-time unroll loop in insertBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, qubits.size());

    #pragma omp parallel for reduction(+:prob) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = local statevector index of nth local basis state with a contributing diagonal
        qindex i = insertBitsWithMaskedValues(n, sortedQubits.data(), numBits, qubitStateMask); // may be unrolled at compile-time

        // j = local, flat, density-matrix index of diagonal amp corresponding to state i
        qindex j = fast_getQuregLocalIndexOfDiagonalAmp(i, firstDiagInd, numAmpsPerCol);

        prob += real(amps[j]);
    }

    return prob;
}


template <int NumQubits>
void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, SmallView qubits) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // every amp contributes to a statevector prob
    qindex numIts = qureg.numAmpsPerNode;

    // use template param to compile-time unroll loop in getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, qubits.size());
    qindex numOutcomes = powerOf2(numBits);

    // decide whether to parallelise below amp-clearing, since outProbs ~ dim of a qureg
    bool parallelise = numBits > MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_MULTITHREADING;
    (void)parallelise; // suppress unused warning when not-compiling openmp)

    // clear amps (may be compile-time unrolled, or parallelised)
    #pragma omp parallel for if(parallelise)
    for (int i=0; i<numOutcomes; i++)
        outProbs[i] = 0;
    
    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        qreal prob = norm(amps[n]);

        // i = global index corresponding to n
        qindex i = concatenateBits(qureg.rank, n, qureg.logNumAmpsPerNode);

        // j = outcome index corresponding to prob
        qindex j = getValueOfBits(i, qubits.data(), numBits); // loop therein may be unrolled

        #pragma omp atomic
        outProbs[j] += prob;
    }
}


template <int NumQubits>
void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, SmallView qubits) {

    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);
    
    // iterate every column, each contributing one element (the diagonal)
    qindex numIts = powerOf2(qureg.logNumColsPerNode);
    qindex numAmpsPerCol = powerOf2(qureg.numQubits);
    qindex firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);

    // use template param to compile-time unroll loop in getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, qubits.size());
    qindex numOutcomes = powerOf2(numBits);

    // decide whether to parallelise below amp-clearing, since outProbs ~ dim of a qureg
    bool parallelise = numBits > MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_MULTITHREADING;
    (void)parallelise; // suppress unused warning when not-compiling openmp)
    
    // clear amps; be compile-time unrolled, and/or parallelised (independent of qureg)
    #pragma omp parallel for if(parallelise)
    for (int i=0; i<numOutcomes; i++)
        outProbs[i] = 0;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = local index of nth local diagonal element
        qindex i = fast_getQuregLocalIndexOfDiagonalAmp(n, firstDiagInd, numAmpsPerCol);
        qreal prob = real(amps[i]);

        // j = global index of i
        qindex j = concatenateBits(qureg.rank, i, qureg.logNumAmpsPerNode);

        // k = outcome index corresponding to basis state j
        qindex k = getValueOfBits(j, qubits.data(), numBits); // loop therein may be unrolled

        #pragma omp atomic
        outProbs[k] += prob;
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( qreal, cpu_statevec_calcProbOfMultiQubitOutcome_sub, (Qureg, SmallView, SmallView) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( qreal, cpu_densmatr_calcProbOfMultiQubitOutcome_sub, (Qureg, SmallView, SmallView) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub, (qreal* outProbs, Qureg, SmallView) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub, (qreal* outProbs, Qureg, SmallView) )



/*
 * INNER PRODUCTS
 */


qcomp cpu_statevec_calcInnerProduct_sub(Qureg quregA, Qureg quregB) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* ampsA = getCpuQcompPtr(quregA.cpuAmps);
    cpu_qcomp* ampsB = getCpuQcompPtr(quregB.cpuAmps);

    // separately reduce real and imag components to make MSVC happy
    qreal prodRe = 0;
    qreal prodIm = 0;

    // every local amp contributes to the reduction
    qindex numIts = quregA.numAmpsPerNode;

    #pragma omp parallel for reduction(+:prodRe,prodIm) if(quregA.isMultithreaded||quregB.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {
        cpu_qcomp term = conj(ampsA[n]) * ampsB[n];

        prodRe += real(term);
        prodIm += imag(term);
    }

    return qcomp(prodRe, prodIm);
}


qreal cpu_densmatr_calcHilbertSchmidtDistance_sub(Qureg quregA, Qureg quregB) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* ampsA = getCpuQcompPtr(quregA.cpuAmps);
    cpu_qcomp* ampsB = getCpuQcompPtr(quregB.cpuAmps);

    qreal dist = 0;

    // every local amp contributes to the reduction
    qindex numIts = quregA.numAmpsPerNode;

    #pragma omp parallel for reduction(+:dist) if(quregA.isMultithreaded||quregA.isMultithreaded)
    for (qindex n=0; n<numIts; n++)
        dist += norm(ampsA[n] - ampsB[n]); // |A-B|^2

    return dist; // do not sqrt yet
}


template <bool Conj>
qcomp cpu_densmatr_calcFidelityWithPureState_sub(Qureg rho, Qureg psi) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* rhoAmps = getCpuQcompPtr(rho.cpuAmps);
    cpu_qcomp* psiAmps = getCpuQcompPtr(psi.cpuAmps);

    // separately reduce real and imag components to make MSVC happy
    qreal fidRe = 0;
    qreal fidIm = 0;

    // every local density matrix amp contributes to the reduction
    qindex numIts = rho.numAmpsPerNode;

    #pragma omp parallel for reduction(+:fidRe,fidIm) if(rho.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = global index of nth local amp of rho
        qindex i = concatenateBits(rho.rank, n, rho.logNumAmpsPerNode);

        // r, c = global row and column indices corresponding to i
        qindex r = getBitsRightOfIndex(i, rho.numQubits);
        qindex c = getBitsLeftOfIndex(i, rho.numQubits-1);

        // collect amps involved in this term
        cpu_qcomp rhoAmp = rhoAmps[n];
        cpu_qcomp rowAmp = psiAmps[r];
        cpu_qcomp colAmp = psiAmps[c]; // likely to be last iteration's amp in cache

        // compute term of <psi|rho^dagger|psi> or <psi|rho|psi>
        if constexpr (Conj) {
            rhoAmp = conj(rhoAmp);
            colAmp = conj(colAmp);
        } else
            rowAmp = conj(rowAmp);

        cpu_qcomp term = rhoAmp * rowAmp * colAmp;
        fidRe += real(term);
        fidIm += imag(term);
    }

    return qcomp(fidRe, fidIm);
}


template qcomp cpu_densmatr_calcFidelityWithPureState_sub<true >(Qureg, Qureg);
template qcomp cpu_densmatr_calcFidelityWithPureState_sub<false>(Qureg, Qureg);



/*
 * PAULI EXPECTATION VALUES
 */


qreal cpu_statevec_calcExpecAnyTargZ_sub(Qureg qureg, SmallView targs) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // this is the only expec-val routine gauranteed to be real,
    // regardless of state normalisation and numerical errors
    qreal value = 0;

    // each iteration contributes one term to the sum
    qindex numIts = qureg.numAmpsPerNode;
    qindex targMask = util_getBitMask(targs);

    #pragma omp parallel for reduction(+:value) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        int sign = fast_getPlusOrMinusMaskedBitParity(n, targMask);
        value += sign * norm(amps[n]);
    }

    return value;
}


qcomp cpu_densmatr_calcExpecAnyTargZ_sub(Qureg qureg, SmallView targs) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // separately reduce real and imag components to make MSVC happy
    qreal valueRe = 0;
    qreal valueIm = 0;

    // each column contributes one amp to sum
    qindex numIts = powerOf2(qureg.logNumColsPerNode);
    qindex numAmpsPerCol = powerOf2(qureg.numQubits);
    qindex firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);

    qindex targMask = util_getBitMask(targs);

    #pragma omp parallel for reduction(+:valueRe,valueIm) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = local index of nth local diagonal element
        qindex i = fast_getQuregLocalIndexOfDiagonalAmp(n, firstDiagInd, numAmpsPerCol);

        // r = global row of nth local diagonal, which determines amp sign
        qindex r = n + firstDiagInd;
        int sign = fast_getPlusOrMinusMaskedBitParity(r, targMask);
        cpu_qcomp term = sign * amps[i];

        valueRe += real(term);
        valueIm += imag(term);
    }

    return qcomp(valueRe, valueIm);
}


qcomp cpu_statevec_calcExpecPauliStr_subA(Qureg qureg, SmallView x, SmallView y, SmallView z) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // separately reduce real and imag components to make MSVC happy
    qreal valueRe = 0;
    qreal valueIm = 0;

    // all local amps appear twice, and each iteration contributes two amps
    qindex numIts = qureg.numAmpsPerNode;

    qindex maskXY = util_getBitMask(util_getConcatenated(x, y));
    qindex maskYZ = util_getBitMask(util_getConcatenated(y, z));

    #pragma omp parallel for reduction(+:valueRe,valueIm) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // j = local index of amp which combines with nth local amp
        qindex j = flipBits(n, maskXY);

        // sign = +-1 induced by Y and Z (excludes Y i factors)
        int sign = fast_getPlusOrMinusMaskedBitParity(j, maskYZ);
        cpu_qcomp term = sign * conj(amps[n]) * amps[j];

        valueRe += real(term);
        valueIm += imag(term);
    }

    // scale by i^numY (because sign above exlcuded i)
    return qcomp(valueRe,valueIm) * util_getPowerOfI(y.size());
}


qcomp cpu_statevec_calcExpecPauliStr_subB(Qureg qureg, SmallView x, SmallView y, SmallView z) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps   = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* buffer = getCpuQcompPtr(qureg.cpuCommBuffer);

    /// @todo
    /// this is identical to the subA() version above, except that
    /// qureg.cpuAmps[j] becomes qureg.cpuCommBuffer[j]. We could
    /// ergo replace subA() with an invocation of subB(), binding
    /// the buffer to the amps ptr. Would this affect/interfere
    /// with memory movement optimisations? I doubt so, but check
    /// and if not, perform the replacement to reduce code-dupe!

    // separately reduce real and imag components to make MSVC happy
    qreal valueRe = 0;
    qreal valueIm = 0;

    // all local amps contribute to the sum
    qindex numIts = qureg.numAmpsPerNode;

    qindex maskXY = util_getBitMask(util_getConcatenated(x, y));
    qindex maskYZ = util_getBitMask(util_getConcatenated(y, z));

    #pragma omp parallel for reduction(+:valueRe,valueIm) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // j = buffer index of amp to be multiplied with nth local amp
        qindex j = flipBits(n, maskXY);

        // sign = +-1 induced by Y and Z (excludes Y i factors)
        int sign = fast_getPlusOrMinusMaskedBitParity(j, maskYZ);
        cpu_qcomp term = sign * conj(amps[n]) * buffer[j];

        valueRe += real(term);
        valueIm += imag(term);
    }

    // scale by i^numY (because sign above exlcuded i)
    return qcomp(valueRe,valueIm) * util_getPowerOfI(y.size());
}


qcomp cpu_densmatr_calcExpecPauliStr_sub(Qureg qureg, SmallView x, SmallView y, SmallView z) {

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // separately reduce real and imag components to make MSVC happy
    qreal valueRe = 0;
    qreal valueIm = 0;

    // each column contributes one amp to sum
    qindex numIts = powerOf2(qureg.logNumColsPerNode);
    qindex numAmpsPerCol = powerOf2(qureg.numQubits);
    qindex firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);

    // these masks indicate global paulis (i.e. not just suffix)
    qindex maskXY = util_getBitMask(util_getConcatenated(x, y));
    qindex maskYZ = util_getBitMask(util_getConcatenated(y, z));

    #pragma omp parallel for reduction(+:valueRe,valueIm) if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // r = global row of nth local diagonal of (qureg)
        qindex r = n + firstDiagInd;

        // i = global row of nth local diagonal of (str . qureg)
        qindex i = flipBits(r, maskXY);

        // m = local flat index of i
        qindex m = fast_getQuregLocalFlatIndex(i, n, numAmpsPerCol);

        // sign = +-1 induced by Y and Z (excludes Y's imaginary factors)
        int sign = fast_getPlusOrMinusMaskedBitParity(i, maskYZ);
        cpu_qcomp term = sign * amps[m];

        valueRe += real(term);
        valueIm += imag(term);
    }

    // scale by i^numY (because sign above exlcuded i)
    return qcomp(valueRe,valueIm) * util_getPowerOfI(y.size());
}



/*
 * DIAGONAL MATRIX EXPECTATION VALUES
 */


template <bool HasPower, bool UseRealPow> 
qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    assert_quregAndFullStateDiagMatrHaveSameDistrib(qureg, matr);
    assert_exponentMatchesTemplateParam(exponent, HasPower, UseRealPow);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.cpuElems);
    cpu_qcomp expo   = getCpuQcomp(exponent);
    (void) expo; // silence when unused

    // separately reduce real and imag components to make MSVC happy
    qreal valueRe = 0;
    qreal valueIm = 0;

    // every amp, iterated independently, contributes to the expectation value
    qindex numIts = qureg.numAmpsPerNode;

    #pragma omp parallel for reduction(+:valueRe,valueIm) if(qureg.isMultithreaded||matr.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        cpu_qcomp elem = elems[n];

        // compile-time decide if applying power (!= 1) to avoid in-loop branching.
        // beware that pow(qcomp,qcomp) below gives notable error over pow(qreal,qreal) 
        // (by producing an unexpected non-zero imaginary component) when the base is real 
        // and negative, and the exponent is an integer. This is a plausible scenario given 
        // matr is often Hermitian (ergo real) and the exponent may be effecting repeated 
        // application of matr. When UseRealPow is set, we discard the imaginary components
        // of the matrix and exponent (latter gauranteed zero) and use the stable pow().
        if constexpr (HasPower && ! UseRealPow)
            elem = pow(elem, expo); // pow(qcomp)
        if constexpr (HasPower &&   UseRealPow)
            elem = getCpuQcomp(std::pow(real(elem), real(expo)), 0); // pow(qreal)

        cpu_qcomp term = elem * norm(amps[n]);

        valueRe += real(term);
        valueIm += imag(term);
    }

    return qcomp(valueRe, valueIm);
}


template <bool HasPower, bool UseRealPow>
qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    assert_quregAndFullStateDiagMatrHaveSameDistrib(qureg, matr);
    assert_exponentMatchesTemplateParam(exponent, HasPower, UseRealPow);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);
    cpu_qcomp* elems = getCpuQcompPtr(matr.cpuElems);
    cpu_qcomp expo   = getCpuQcomp(exponent);
    (void) expo; // silence when unused

    // separately reduce real and imag components to make MSVC happy
    qreal valueRe = 0;
    qreal valueIm = 0;

    // iterate each column, of which one amp (the diagonal) contributes
    qindex numIts = powerOf2(qureg.logNumColsPerNode);
    qindex numAmpsPerCol = powerOf2(qureg.numQubits);
    qindex firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);

    #pragma omp parallel for reduction(+:valueRe,valueIm) if(qureg.isMultithreaded||matr.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        cpu_qcomp elem = elems[n];

        // compile-time decide if applying power (!= 1) to avoid in-loop branching.
        // Beware that pow(qcomp,qcomp) below gives notable error over pow(qreal,qreal) 
        // (by producing an unexpected non-zero imaginary component) when the base is real 
        // and negative, and the exponent is an integer. This is a plausible scenario given 
        // matr is often Hermitian (ergo real) and the exponent may be effecting repeated 
        // application of matr. So when UseRealPow is set, we discard the imaginary components
        // of the matrix and exponent (latter gauranteed zero) and use the stable pow().
        // An observant reader may realise the erroneous imaginary components introduced by
        // complex-pow do not damage the statevector expectation value; the imaginary component
        // can be discarded entirely at the end without harm to the real component. Alas this
        // is not true of the density matrix calculation here, where the erroneous imaginary
        // components of pow(qcomp,qcomp) will sabotage the final result. Wah!
        if constexpr (HasPower && ! UseRealPow)
            elem = pow(elem, expo);
        if constexpr (HasPower &&   UseRealPow)
            elem = getCpuQcomp(std::pow(real(elem), real(expo)), 0);

        // i = local index of nth local diagonal element
        qindex i = fast_getQuregLocalIndexOfDiagonalAmp(n, firstDiagInd, numAmpsPerCol);
        cpu_qcomp term = elem * amps[i];

        valueRe += real(term);
        valueIm += imag(term);
    }

    return qcomp(valueRe, valueIm);
}


template qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub<true, true >(Qureg, FullStateDiagMatr, qcomp);
template qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub<true, false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub<false,false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub<false,true >(Qureg, FullStateDiagMatr, qcomp); // uncallable

template qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub<true, true >(Qureg, FullStateDiagMatr, qcomp);
template qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub<true, false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub<false,false>(Qureg, FullStateDiagMatr, qcomp);
template qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub<false,true >(Qureg, FullStateDiagMatr, qcomp); // uncallable



/*
 * PROJECTORS
 */


template <int NumQubits>
void cpu_statevec_multiQubitProjector_sub(Qureg qureg, SmallView qubits, SmallView outcomes, qreal prob) {

    // all qubits are in suffix
    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps  = getCpuQcompPtr(qureg.cpuAmps);

    // visit every amp, setting to zero or multiplying it by renorm
    qindex numIts = qureg.numAmpsPerNode;

    // binary value of targeted qubits in basis states which are to be retained
    qindex retainValue = getIntegerFromBits(outcomes.data(), outcomes.size());
    qreal renorm = 1 / std::sqrt(prob);

    // use template param to compile-time unroll loop in getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, qubits.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // val = outcomes corresponding to n-th local amp (all qubits are in suffix)
        qindex val = getValueOfBits(n, qubits.data(), numBits);

        // multiply amp with renorm or zero, if qubit value matches or disagrees
        amps[n] *= renorm * (val == retainValue);
    }
}


template <int NumQubits>
void cpu_densmatr_multiQubitProjector_sub(Qureg qureg, SmallView qubits, SmallView outcomes, qreal prob) {

    // this function is merely an optimisation to avoid calling the above
    // cpu_statevec_multiQubitProjector_sub() twice upon a density matrix;
    // pre- and post-multiply projector versions DO just call above.

    // qubits are unconstrained, and can include prefix qubits
    assert_numTargsMatchesTemplateParam(qubits.size(), NumQubits);

    // use cpu_qcomp arithmetic overloads (avoid qcomp's)
    cpu_qcomp* amps = getCpuQcompPtr(qureg.cpuAmps);

    // visit every amp, setting most to zero and multiplying the remainder by renorm
    qindex numIts = qureg.numAmpsPerNode;

    // binary value of targeted qubits in basis states which are to be retained
    qindex retainValue = getIntegerFromBits(outcomes.data(), outcomes.size());
    qreal renorm = 1 / prob;

    // use template param to compile-time unroll loops in getValueOfBits()
    SET_VAR_AT_COMPILE_TIME(int, numBits, NumQubits, qubits.size());

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = global index of nth local amp
        qindex i = concatenateBits(qureg.rank, n, qureg.logNumAmpsPerNode);

        // r, c = global row and column indices of nth local amp
        qindex r = getBitsRightOfIndex(i, qureg.numQubits);
        qindex c = getBitsLeftOfIndex(i, qureg.numQubits-1);

        qindex v1 = getValueOfBits(r, qubits.data(), numBits);
        qindex v2 = getValueOfBits(c, qubits.data(), numBits);

        // multiply amp with renorm or zero if values disagree with given outcomes
        amps[n] *= renorm * (v1 == v2) * (retainValue == v1);
    }
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, cpu_statevec_multiQubitProjector_sub, (Qureg qureg, SmallView qubits, SmallView outcomes, qreal prob) )
INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS( void, cpu_densmatr_multiQubitProjector_sub, (Qureg qureg, SmallView qubits, SmallView outcomes, qreal prob) )



/*
 * STATE INITIALISATION
 */


void cpu_statevec_initUniformState_sub(Qureg qureg, qcomp amp) {

    // approx-uniformly distribute modified memory pages across threads,
    // in the hope that each std::fill() will touch only memory within 
    // the thread's corresponding NUMA node, for best performance 

    int numAmpsPerPage = cpu_getPageSize() / sizeof(qcomp); // divides evenly

    #pragma omp parallel if(qureg.isMultithreaded)
    {
        const auto [start, end] = util_getBlockMultipleSubRange(
            qureg.numAmpsPerNode, numAmpsPerPage,
            cpu_getOpenmpThreadInd(), cpu_getCurrentNumThreads());

        std::fill(qureg.cpuAmps + start, qureg.cpuAmps + end, amp); // no cpu_qcomp cast needed
    }
}


void cpu_statevec_initDebugState_sub(Qureg qureg) {

    // overwrite all local amps
    qindex numIts = qureg.numAmpsPerNode;

    #pragma omp parallel for if(qureg.isMultithreaded)
    for (qindex n=0; n<numIts; n++) {

        // i = global index of nth local amp
        qindex i = concatenateBits(qureg.rank, n, qureg.logNumAmpsPerNode);
        qureg.cpuAmps[n] = qcomp(2*i/10., (2*i+1)/10.); // no cpu_qcomp cast needed
    }
}


void cpu_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(Qureg qureg) {

    // all amplitudes are re-randomised, one per iteration
    qindex numIts = qureg.numAmpsPerNode;

    // thread seeds uniquely deviate from a random base seed, which may be node-specific
    unsigned seed = rand_getThreadSharedRandomSeed(qureg.isDistributed);

    // create an explicit parallel region to avoid re-initialisation of RNG every iteration
    #pragma omp parallel if(qureg.isMultithreaded)
    {
        int id = cpu_getOpenmpThreadInd(); // zero if OpenMP not compiled

        // prepare uniquely-seeded thread-private generator
        auto gen = rand_getThreadPrivateGenerator(seed, id);
        auto normDist = rand_getThreadPrivateAmpAbsDistribution();
        auto phaseDist = rand_getThreadPrivateAmpPhaseDistribution();

        #pragma omp for
        for (qindex i=0; i<numIts; i++)
            qureg.cpuAmps[i] = rand_getThreadPrivateRandomAmp(gen, normDist, phaseDist); // advances gen, no cpu_qcomp cast needed
    }
}
