/** @file
 * Testing utilities which evolve a reference state
 * (qvector or qmatrix) under the action of a 
 * reference operation. These are slow, serial,
 * un-optimised, defensively-designed routines.
 *
 * @author Tyson Jones
 */

#include "qvector.hpp"
#include "qmatrix.hpp"
#include "macros.hpp"
#include "linalg.hpp"

#include <tuple>
#include <vector>
#include <algorithm>

using std::vector;



/*
 * OPERATOR MATRICES
 */


qmatrix getSwapMatrix(int qb1, int qb2, int numQb) {
    DEMAND( numQb > 1 );
    DEMAND( (qb1 >= 0 && qb1 < numQb) );
    DEMAND( (qb2 >= 0 && qb2 < numQb) );

    if (qb1 == qb2)
        return getIdentityMatrix(getPow2(numQb));
        
    if (qb1 > qb2)
        std::swap(qb1, qb2);

    qmatrix out;
    
    // qubits are either adjacent
    if (qb2 == qb1 + 1) {
        out = qmatrix{{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};
    
    // or distant
    } else {
        int block = getPow2(qb2 - qb1);
        out = getZeroMatrix(block*2);
        qmatrix iden = getIdentityMatrix(block/2);
        
        // Lemma 3.1 of arxiv.org/pdf/1711.09765.pdf
        qmatrix p0{{1,0},{0,0}};
        qmatrix l0{{0,1},{0,0}};
        qmatrix l1{{0,0},{1,0}};
        qmatrix p1{{0,0},{0,1}};
        
        // notating a^(n+1) = identity(getPow2(n)) (otimes) a, we construct the matrix
        // [ p0^(N) l1^N ]
        // [ l0^(N) p1^N ]
        // where N = qb2 - qb1 */
        setSubMatrix(out, getKroneckerProduct(iden, p0), 0, 0);
        setSubMatrix(out, getKroneckerProduct(iden, l0), block, 0);
        setSubMatrix(out, getKroneckerProduct(iden, l1), 0, block);
        setSubMatrix(out, getKroneckerProduct(iden, p1), block, block);
    }
    
    // pad swap with outer identities
    if (qb1 > 0)
        out = getKroneckerProduct(out, getIdentityMatrix(getPow2(qb1)));

    if (qb2 < numQb-1)
        out = getKroneckerProduct(getIdentityMatrix(getPow2(numQb-qb2-1)), out);
        
    return out;
}


auto getSwapAndUnswapMatrices(vector<int> ctrls, vector<int> targs, size_t numQubits) {
    DEMAND( numQubits >= ctrls.size() + targs.size() );

    // matrices which swap targs+ctrls to be contiguous from 0
    qmatrix swaps   = getIdentityMatrix(getPow2(numQubits));
    qmatrix unswaps = getIdentityMatrix(getPow2(numQubits));

    // swap targs to {0, ..., ntargs - 1}
    for (size_t i=0; i<targs.size(); i++) {

        if (i == (size_t) targs[i])
            continue;

        qmatrix m = getSwapMatrix(i, targs[i], numQubits);
        swaps = m * swaps;
        unswaps = unswaps * m;
        
        std::replace(ctrls.begin(), ctrls.end(), (int) i, targs[i]);
        std::replace(targs.begin(), targs.end(), (int) i, targs[i]);
    }

    // swap ctrls to {ntargs, ..., ntargs + nctrls - 1}
    for (size_t i=0; i<ctrls.size(); i++) {

        size_t j = i + targs.size();
        if (j == (size_t) ctrls[i])
            continue;

        qmatrix m = getSwapMatrix(j, ctrls[i], numQubits);
        swaps = m * swaps;
        unswaps = unswaps * m;

        std::replace(ctrls.begin(), ctrls.end(), (int) j, ctrls[i]);
    }

    return std::tuple{swaps, unswaps};
}


qmatrix getFullStateOperator(vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qmatrix matrix, size_t numQubits) {
    DEMAND( numQubits >= ctrls.size() + targs.size() );
    DEMAND( getPow2(targs.size()) == (qindex) matrix.size() );
    DEMAND( (ctrlStates.empty() || ctrlStates.size() == ctrls.size()) );

    // construct controlled-(matrix) upon lowest order qubits
    qmatrix full = getControlledMatrix(matrix, ctrls.size());

    // left-pad 'full' to be numQubits large
    if (numQubits > ctrls.size() + targs.size()) {
        size_t pad = getPow2(numQubits - ctrls.size() - targs.size());
        full = getKroneckerProduct(getIdentityMatrix(pad), full);
    }
    
    // apply swaps to retarget 'full' to given ctrls and targs
    auto [swaps, unswaps] = getSwapAndUnswapMatrices(ctrls, targs, numQubits);
    qmatrix out = unswaps * full * swaps;

    // apply NOT to all zero-controlled qubits (recurses just once)
    qmatrix matrX = {{0,1},{1,0}};
    for (size_t i=0; i<ctrlStates.size(); i++) {
        
        if (ctrlStates[i] == 1)
            continue;

        qmatrix fullX = getFullStateOperator({}, {}, {ctrls[i]}, matrX, numQubits);
        out = fullX * out * fullX;
    }

    return out;
}



/*
 * EVOLUTION
 */


// overloads with no targs (given full operator)

void applyReferenceOperator(qvector& state, qmatrix matrix) {
    DEMAND( state.size() == matrix.size() );

    state = matrix * state;
}
void applyReferenceOperator(qmatrix& state, qmatrix matrix) {
    DEMAND( state.size() == matrix.size() );

    state = matrix * state * getConjugateTranspose(matrix);
}

void multiplyReferenceOperator(qvector& state, qmatrix matrix) {
    DEMAND( state.size() == matrix.size() );

    // for statevectors, multiplying is the same as applying
    applyReferenceOperator(state, matrix);
}

void multiplyReferenceOperator(qmatrix& state, qmatrix matrix) {
    DEMAND( state.size() == matrix.size() );

    // we left-multiply upon density matrices only
    state = matrix * state;
}


// overloads with ctrls, states and targs (given sub-operator)

void applyReferenceOperator(qvector& state, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qmatrix matrix) {

    qmatrix fullOp = getFullStateOperator(ctrls, ctrlStates, targs, matrix, getLog2(state.size()));
    applyReferenceOperator(state, fullOp);
}

void applyReferenceOperator(qmatrix& state, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qmatrix matrix) {
    
    qmatrix fullOp = getFullStateOperator(ctrls, ctrlStates, targs, matrix, getLog2(state.size()));
    applyReferenceOperator(state, fullOp);
}

void multiplyReferenceOperator(qvector& state, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qmatrix matrix) {
    
    applyReferenceOperator(state, ctrls, ctrlStates, targs, matrix);
}

void multiplyReferenceOperator(qmatrix& state, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qmatrix matrix) {
    
    qmatrix left = getFullStateOperator(ctrls, ctrlStates, targs, matrix, getLog2(state.size()));
    multiplyReferenceOperator(state, left);
}


// overloads with only ctrls and targs

void applyReferenceOperator(qvector& state, vector<int> ctrls, vector<int> targs, qmatrix matrix) {

    applyReferenceOperator(state, ctrls, {}, targs, matrix);
}
void applyReferenceOperator(qmatrix& state, vector<int> ctrls, vector<int> targs, qmatrix matrix) {
    
    applyReferenceOperator(state, ctrls, {}, targs, matrix);
}
void multiplyReferenceOperator(qvector& state, vector<int> ctrls, vector<int> targs, qmatrix matrix) {
    
    multiplyReferenceOperator(state, ctrls, {}, targs, matrix);
}
void multiplyReferenceOperator(qmatrix& state, vector<int> ctrls, vector<int> targs, qmatrix matrix) {
    
    multiplyReferenceOperator(state, ctrls, {}, targs, matrix);
}


// overloads with only targs

void applyReferenceOperator(qvector& state, vector<int> targs, qmatrix matrix) {
    
    applyReferenceOperator(state, {}, {}, targs, matrix);
}
void applyReferenceOperator(qmatrix& state, vector<int> targs, qmatrix matrix) {

    applyReferenceOperator(state, {}, {}, targs, matrix);
}
void multiplyReferenceOperator(qvector& state, vector<int> targs, qmatrix matrix) {

    multiplyReferenceOperator(state, {}, {}, targs, matrix);
}
void multiplyReferenceOperator(qmatrix& state, vector<int> targs, qmatrix matrix) {

    multiplyReferenceOperator(state, {}, {}, targs, matrix);
}


// overloads with only targs and kraus operators

void applyReferenceOperator(qmatrix& state, vector<int> targs, vector<qmatrix> matrices) {

    qmatrix in = state;
    qmatrix out = getZeroMatrix(state.size());

    for (auto& matrix : matrices) {
        state = in;
        applyReferenceOperator(state, targs, matrix);
        out += state;
    }

    state = out;
}
