/** @file
 * Inlined operations to determine indices of modified amplitudes
 * using efficient bitwise operations, usable in hot-loops.
 * The getNthIndex() routines accept a compile-time template flag 
 * indicating preconditions about the passed control qubits, which 
 * are used to compile-time optimise the invocation.
 * 
 * Note that the getNthIndex() routines are inlined into both
 * OpenMP loops and CUDA kernels, so cannot make use of e.g.
 * vectors, and passed arrays therein might be device pointers.
 * 
 * The INLINE macro is defined in bitwise.hpp
 */

#ifndef INDEXER_HPP
#define INDEXER_HPP

#include "quest/include/types.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/inliner.hpp"
#include "quest/src/core/bitwise.hpp"

#include <tuple>
#include <vector>
#include <algorithm>

using std::tuple;
using std::vector;



/*
 * COMPILE-TIME PRECONDITION FLAGS
 *
 * which are used as template parameters to specify preconditions about the
 * qubits passed to a function, so that a bespoke, optimised version of the
 * function is invoked. The flag determines which logic is used to find
 * the indices of local amplitudes to modify in hot loops, accelerating
 * the loop body and making memory access more predictable so that the compiler
 * can improve caching performance.
 */


namespace index_flags {
    
    enum CtrlFlag { 
        NO_CTRLS, 
        ONE_CTRL, 
        ONE_STATE_CTRL,
        MULTI_CTRLS, 
        MULTI_STATE_CTRLS
    };
}

using namespace index_flags;



/*
 * TEMPLATE INSTANTIATION MACRO
 *
 * which are used by other files to force their CtrlFlag-templated functions
 * to be instantiated and compiled, once for each possible flag value. This
 * is necessary because the invocation of the templated function happens in 
 * a separate compilation unit to where it is defined, and the compiler needs
 * to know which template args will be later invoked and need instantiation. 
 */


#define INSTANTIATE_TEMPLATED_VOID_FUNC_WITH_ALL_CTRL_FLAGS(funcname, signature) \
    template void funcname <NO_CTRLS>          signature; \
    template void funcname <ONE_CTRL>          signature; \
    template void funcname <ONE_STATE_CTRL>    signature; \
    template void funcname <MULTI_CTRLS>       signature; \
    template void funcname <MULTI_STATE_CTRLS> signature;



/*
 * INTERNAL ERROR CHECK
 */


static void indexer_assertValidCtrls(vector<int> ctrls, vector<int> ctrlStates) {

    if (ctrlStates.size() != ctrls.size() && !ctrlStates.empty())
        error_indexerCtrlsInconsistentWithCtrlStates();
}



/*
 * RUNTIME INFERENCE OF QUBIT PRECONDITIONS
 *
 * which is called in advance of hot loops and is only inlined here to 
 * avoid unused function warnings. Note the use of casual 'inline' over the
 * inliner.hpp's forceful INLINE macro; we cannot declare this as INLINE
 * because it accepts vectors incompatible with CUDA kernels.
 */


static inline CtrlFlag indexer_getCtrlFlag(vector<int> ctrls, vector<int> states) {
    indexer_assertValidCtrls(ctrls, states);

    if (ctrls.empty())
        return NO_CTRLS;

    if (ctrls.size() == 1)
        return (states.empty())? ONE_CTRL : ONE_STATE_CTRL;

    return (states.empty())? MULTI_CTRLS : MULTI_STATE_CTRLS;
}



/*
 * INDEX PARAMETERS
 *
 * which contain all the necessary information for hot loops to efficiently
 * compute indices of amplitudes in a manner optimised for the qubit
 * preconditions. Only the fields relevant to the bespoke index logic will be 
 * initialised. Note that these params are used both by OpenMP loops and CUDA 
 * kernels, so must not contain anything incompatible with either (e.g. vectors).
 */


// informs indices where all ctrl qubits are in their active state
typedef struct {

    // when no ctrls are given, params yield all indices

    // populated when a single ctrl is given
    int ctrl;

    // populated when the single ctrl has a given state
    int ctrlState;

    // populated when many ctrls are given, pointing to caller's heap memory
    int numCtrls;
    int* sortedCtrls;
    
    // populated when many ctrl states are given
    qindex ctrlStateMask;

} CtrlIndParams;


// extends CtrlIndParams to constrain the target qubit to be active (i.e. 1)
typedef struct {

    // when no ctrls are given, params yield all indices where target is 1

    // populated when a single ctrl is given, distinguishing ctrl and target
    int lowerQubit;
    int upperQubit;

    // populated when a single ctrl state is given, ordering ctrl=state and target=1 bits
    int lowerBit;
    int upperBit;

    // contains ctrls + target, populated when many ctrls are given, pointing to caller's heap memory
    int numQubits;
    int* sortedQubits;
    
    // populated when many ctrl states are given
    qindex ctrlStateMask;

} CtrlTargIndParams;




/*
 * INDEX PARAM INITIALISERS
 *
 * which are called in advance of hot loops and are only 
 * inlined here to avoid unused function warnings
 */


static inline CtrlTargIndParams getParamsInformingIndsWhereCtrlsAreActiveAndTargIsOne(vector<int> &ctrls, vector<int> states, int targ) {
    indexer_assertValidCtrls(ctrls, states);

    // most fields of params will stay un-initialised
    CtrlTargIndParams params;

    // if there's a single ctrl, determine its relative position to targ
    if (ctrls.size() == 1) {
        params.lowerQubit = (ctrls[0] < targ)? ctrls[0] : targ;
        params.upperQubit = (ctrls[0] > targ)? ctrls[0] : targ;
    }

    // if single ctrl has a given bit state, make sure it's inserted at the ctrl qubit, inserting 1 at target
    if (states.size() == 1) {
        params.lowerBit = (params.lowerQubit == ctrls[0])? states[0] : 1;
        params.upperBit = (params.upperQubit == ctrls[0])? states[0] : 1;
    }

    // if many ctrl states are given, create a mask from them before ctrls gets modified below
    if (states.size() > 1)
        params.ctrlStateMask = getBitMask(ctrls.data(), states.data(), ctrls.size());

    // if many ctrls are given, add targ to them and sort them (modifying passed ctrls vector)
    if (ctrls.size() > 1) {
        ctrls.push_back(targ);
        std::sort(ctrls.begin(), ctrls.end());
        params.sortedQubits = ctrls.data();  // maintain ptr to caller's vector
    }

    return params;
}


static inline CtrlIndParams getParamsInformingIndsWhereCtrlsAreActive(vector<int> &ctrls, vector<int> states) { 
    indexer_assertValidCtrls(ctrls, states);

    // most fields of params will stay un-initialised
    CtrlIndParams params;

    if (ctrls.size() == 1)
        params.ctrl = ctrls[0];

    if (states.size() == 1)
        params.ctrlState = states[0];

    // if many ctrl states are given, create a mask from them before ctrls gets modified below
    if (states.size() > 1)
        params.ctrlStateMask = getBitMask(ctrls.data(), states.data(), ctrls.size());

    // if many ctrls are given, sort them (modifying passed ctrls vector)
    if (ctrls.size() > 1) {
        std::sort(ctrls.begin(), ctrls.end());
        params.sortedCtrls = ctrls.data();  // maintain ptr to caller's vector
    }

    return params;
}



/*
 * INDEX CALCULATIONS
 *
 * which use the compile-time CtrlFlag to choose the fastest strategy for computing
 * amplitude indices, using the information pre-computed in the params structs.
 * These are invoked in hot loops, so must ergo be compatible with / callable by 
 * OpenMP threads and CUDA kernels. Inlining and compile-time type-traiting ensure
 * there is no runtime penalty for function stack preparation nor branching.
 */


template <CtrlFlag flag>
INLINE qindex getNthIndWhereCtrlsAreActive(qindex n, const CtrlIndParams params) {

    // maps |n> to |m, ctrls=active>

    if constexpr (flag == NO_CTRLS)
        return n;

    if constexpr (flag == ONE_CTRL)
        return insertBit(n, params.ctrl, 1);

    if constexpr (flag == ONE_STATE_CTRL)
        return insertBit(n, params.ctrl, params.ctrlState);

    if constexpr (flag == MULTI_CTRLS)
        return insertBits(n, params.sortedCtrls, params.numCtrls, 1);

    if constexpr (flag == MULTI_STATE_CTRLS) {
        qindex k = insertBits(n, params.sortedCtrls, params.numCtrls, 0);
        return activateBits(k, params.ctrlStateMask);
    }
}


template <CtrlFlag flag>
INLINE qindex getNthIndWhereCtrlsAreActiveAndTargIsOne(qindex n, int targ, const CtrlTargIndParams params) {

    // maps |n> to |m, ctrls=active, targ=1>

    if constexpr (flag == NO_CTRLS)
        return insertBit(n, targ, 1);

    if constexpr (flag == ONE_CTRL)
        return insertTwoBits(n, params.upperQubit, 1, params.lowerQubit, 1);

    if constexpr (flag == ONE_STATE_CTRL)
        return insertTwoBits(n, params.upperQubit, params.upperBit, params.lowerQubit, params.lowerBit);

    if constexpr (flag == MULTI_CTRLS)
        return insertBits(n, params.sortedQubits, params.numQubits, 1);

    if constexpr (flag == MULTI_STATE_CTRLS) {
        qindex k0 = insertBits(n, params.sortedQubits, params.numQubits, 0);
        qindex i0 = activateBits(k0, params.ctrlStateMask);
        return flipBit(i0, targ);
    }
}



#endif // INDEXER_HPP