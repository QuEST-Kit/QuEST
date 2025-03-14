/** @file
 * Testing utilities which evaluate measurements upon
 * reference qvector and qmatrix states. These are slow, 
 * serial, un-optimised, defensively-designed routines.
 *
 * @author Tyson Jones
 */

#include "quest/include/quest.h"

#include "qvector.hpp"
#include "qmatrix.hpp"
#include "convert.hpp"
#include "linalg.hpp"
#include "macros.hpp"

#include <algorithm>
#include <complex>
#include <vector>

using std::vector;



/*
 * EXPECTATION VALUES
 */


qcomp getReferenceExpectationValue(qvector state, qmatrix observable) {
    DEMAND( state.size() == observable.size() );

    return getInnerProduct(state, observable * state);
}

qcomp getReferenceExpectationValue(qmatrix state, qmatrix observable) {
    DEMAND( state.size() == observable.size() );

    return getTrace(observable * state);
}


qcomp getRefExpecValInner(auto state, auto paulis) {
    
    int numQubits = getLog2(state.size());
    qmatrix observable = getMatrix(paulis, numQubits);
    return getReferenceExpectationValue(state, observable);
}
qcomp getReferenceExpectationValue(qvector state, PauliStr    str) { return getRefExpecValInner(state, str); }
qcomp getReferenceExpectationValue(qmatrix state, PauliStr    str) { return getRefExpecValInner(state, str); }
qcomp getReferenceExpectationValue(qvector state, PauliStrSum sum) { return getRefExpecValInner(state, sum); }
qcomp getReferenceExpectationValue(qmatrix state, PauliStrSum sum) { return getRefExpecValInner(state, sum); }



/*
 * PROBABILITIES
 */


qreal getRefProbInner(auto& state, vector<int> targets, vector<int> outcomes) {
    DEMAND( getLog2(state.size()) > *std::max_element(targets.begin(), targets.end()) );
    DEMAND( targets.size() == outcomes.size() );

    // <psi| (|o><o| (x) I) |psi> = Tr( (|o><o| (x) I) rho)
    int numQubits = getLog2(state.size());
    qmatrix projector = getProjector(targets, outcomes, numQubits);
    qcomp value = getReferenceExpectationValue(state, projector); // ~0
    return std::real(value);
}
qreal getReferenceProbability(qvector state, vector<int> targets, vector<int> outcomes) { return getRefProbInner(state, targets, outcomes); }
qreal getReferenceProbability(qmatrix state, vector<int> targets, vector<int> outcomes) { return getRefProbInner(state, targets, outcomes); }


qreal getReferenceProbability(qvector state, qindex basisIndex) {
    DEMAND( basisIndex < (qindex) state.size() );

    qcomp elem = state[basisIndex];
    qreal prob = std::norm(elem);
    return prob;
}

qreal getReferenceProbability(qmatrix state, qindex basisIndex) {
    DEMAND( basisIndex < (qindex) state.size() );

    qcomp elem = state[basisIndex][basisIndex];
    qreal prob = std::real(elem);
    return prob;
}


vector<qreal> getAllRefProbsInner(auto& state, vector<int> targets) {

    vector<qreal> out(getPow2(targets.size()));
    vector<int> outcomes(targets.size());

    for (size_t i=0; i<out.size(); i++) {

        for (size_t j=0; j<outcomes.size(); j++)
            outcomes[j] = getBitAt(i, j);

        out[i] = getReferenceProbability(state, targets, outcomes);
    }

    return out;
}
vector<qreal> getAllReferenceProbabilities(qvector state, vector<int> targets) { return getAllRefProbsInner(state, targets); }
vector<qreal> getAllReferenceProbabilities(qmatrix state, vector<int> targets) { return getAllRefProbsInner(state, targets); }


qreal getReferenceProbability(qvector state) {

    qreal out = 0;
    for (auto& elem : state)
        out += std::norm(elem);

    return out;
}

qreal getReferenceProbability(qmatrix state) {

    qreal out = 0;
    for (size_t i=0; i<state.size(); i++)
        out += std::real(state[i][i]);

    return out;
}


qreal getReferencePurity(qmatrix state) {
    
    return std::real(getTrace(state * state));
}
qreal getReferencePurity(qvector state) {

    return getReferencePurity(getOuterProduct(state, state));
}
