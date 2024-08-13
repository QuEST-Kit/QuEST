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








// only inline to avoid duplication; not essential for performance in any way

INLINE tuple<vector<int>,qindex> getSortedQubitsAndMask(vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, vector<int> targStates) {

    // TODO: assert internal preconditions, i.e. that
    //  - ctrlStates==ctrls
    //  - targs==targStates


    // merge all qubits and states
    vector<int> qubits = ctrls;
    vector<int> states = ctrlStates;
    qubits.insert(qubits.end(), targs.begin(), targs.end());
    states.insert(states.end(), targStates.begin(), targStates.end());

    // create state mask, then sort qubits (strictly in that order)
    qindex mask = getBitMask(qubits.data(), states.data(), states.size());
    std::sort(qubits.begin(), qubits.end());
    
    return {qubits, mask};
}





#endif // INDEXER_HPP