#define _USE_MATH_DEFINES
#include "qvector.hpp"
#include "qmatrix.hpp"
#include "macros.hpp"
#include "linalg.hpp"
#include "lists.hpp"
#include "quest/include/quest.h"

#include <vector>
#include <tuple>
#include <random>
#include <algorithm>

using std::vector;
using std::tuple;



/*
 * RNG
 */


static std::mt19937 RNG;


void setRandomTestStateSeeds() {

    // generate a random seed from hardware rng
    std::random_device cspnrg;
    unsigned seed = cspnrg();
    
    // seed QuEST, using only the root node's seed
    setSeeds(&seed, 1);

    // broadcat root node seed to all nodes
    getSeeds(&seed);

    // seed rand()
    srand(seed);

    // seed RNG
    RNG.seed(seed);
}



/*
 * SCALAR
 */


qreal getRandomReal(qreal min, qreal maxIncl) {
    DEMAND( min <= maxIncl );

    qreal r = rand() / static_cast<qreal>(RAND_MAX);
    return min + r * (maxIncl - min);
}


int getRandomInt(int min, int maxExcl) {
    return (int) round(getRandomReal(min, maxExcl-1));
}


qcomp getRandomComplex() {
    qreal re = getRandomReal(-1,1);
    qreal im = getRandomReal(-1,1);
    return qcomp(re, im);
}



/*
 * LIST
 */


vector<int> getRandomInts(int min, int maxExcl, int len) {
    DEMAND( len >= 0 ); // permit empty

    vector<int> outcomes(len);
    
    for (auto& x : outcomes)
        x = getRandomInt(min, maxExcl);

    return outcomes;
}


vector<int> getRandomSubRange(int start, int endExcl, int numElems) {
    DEMAND( endExcl >= start );
    DEMAND( numElems >= 1 );
    DEMAND( numElems <= endExcl - start );

    // shuffle entire range
    vector<int> range = getRange(start, endExcl);
    std::shuffle(range.begin(), range.end(), RNG);
    
    // return first subrange
    return vector<int>(range.begin(), range.begin() + numElems);
}


vector<qreal> getRandomProbabilities(int numProbs) {
    
    vector<qreal> probs(numProbs, 0);

    // generate random unnormalised scalars
    for (auto& p : probs)
        p = getRandomReal(0, 1);

    // normalise
    qreal total = 0;
    for (auto& p : probs)
        total += p;

    for (auto& p : probs)
        p /= total;
 
    return probs;
}


auto getRandomCtrlsStatesTargs(int numQubits, int minNumTargs, int maxNumTargsIncl) {
    DEMAND( minNumTargs <= maxNumTargsIncl );
    DEMAND( maxNumTargsIncl <= numQubits );

    int numTargs = getRandomInt(minNumTargs, maxNumTargsIncl+1);

    // numCtrls in [0, remainingNumQb]
    int minNumCtrls = 0;
    int maxNumCtrls = numQubits - numTargs;
    int numCtrls = getRandomInt(minNumCtrls, maxNumCtrls+1);

    vector<int> targsCtrls = getRandomSubRange(0, numQubits, numTargs + numCtrls);
    vector<int> targs = getSublist(targsCtrls, 0, numTargs);
    vector<int> ctrls = getSublist(targsCtrls, numTargs, numCtrls);
    vector<int> states = getRandomInts(0, 2, numCtrls);

    return tuple{ctrls,states,targs};
}



/*
 * VECTOR
 */


qvector getRandomVector(size_t dim) { 

    qvector vec = getZeroVector(dim);

    for (auto& elem : vec)
        elem = getRandomComplex();
        
    return vec;
}


vector<qvector> getRandomOrthonormalVectors(size_t dim, int numVecs) {
    DEMAND( dim >= 1 );
    DEMAND( numVecs >= 1);
    
    vector<qvector> vecs(numVecs);

    // produce each vector in-turn
    for (int n=0; n<numVecs; n++) {

        // from a random vector
        vecs[n] = getRandomVector(dim);

        // orthogonalise by substracting projections of existing vectors
        for (int m=0; m<n; m++)
            vecs[n] -= vecs[m] * getInnerProduct(vecs[m], vecs[n]);

        // then re-normalise
        vecs[n] = getNormalised(vecs[n]);
    }

    return vecs;
}



/*
 * MATRIX
 */


qmatrix getRandomMatrix(size_t dim) {
    DEMAND( dim > 1 );
    
    qmatrix out = getZeroMatrix(dim);
    qreal max = RAND_MAX;

    for (auto& row : out) {
        for (auto& elem : row) {
            
            // generate 2 normally-distributed random numbers via Box-Muller
            qreal a = rand()/max;
            qreal b = rand()/max;
            
            // prevent log(0) NaNs
            if (a == 0)
                a = 1/max;

            qreal fa = sqrt(-2 * log(a));
            qreal re = fa * cos(2 * 3.14159265 * b);
            qreal im = fa * sin(2 * 3.14159265 * b);
            elem = qcomp(re, im);
        }
    }
    
    return out;
}



/*
 * STATES
 */


qvector getRandomStateVector(int numQb) {

    return getNormalised(getRandomVector(getPow2(numQb)));
}


vector<qvector> getRandomOrthonormalStateVectors(int numQb, int numStates) {

    return getRandomOrthonormalVectors(getPow2(numQb), numStates);
}


qmatrix getRandomDensityMatrix(int numQb) {
    DEMAND( numQb > 0 );
    
    // generate random probabilities to weight random pure states
    int dim = getPow2(numQb);
    vector<qreal> probs = getRandomProbabilities(dim);
    
    // add random pure states
    qmatrix dens = getZeroMatrix(dim);
    for (int i=0; i<dim; i++) {
        qvector pure = getRandomStateVector(numQb);
        dens += probs[i] * getOuterProduct(pure, pure);
    }
    
    return dens;
}


qmatrix getRandomPureDensityMatrix(int numQb) {

    qvector vec = getRandomStateVector(numQb);
    qmatrix mat = getOuterProduct(vec, vec);
    return mat;
}


void setToRandomState(qvector& state) {
    state = getRandomStateVector(getLog2(state.size()));
}
void setToRandomState(qmatrix& state) {
    state = getRandomDensityMatrix(getLog2(state.size()));
}



/*
 * OPERATORS
 */


qmatrix getRandomUnitary(int numQb) {
    DEMAND( numQb >= 1 );

    // create Z ~ random complex matrix (distribution not too important)
    size_t dim = getPow2(numQb);
    qmatrix matrZ = getRandomMatrix(dim);
    qmatrix matrZT = getTranspose(matrZ);

    // create Z = Q R (via QR decomposition) ...
    qmatrix matrQT = getOrthonormalisedRows(matrZ);
    qmatrix matrQ = getTranspose(matrQT);
    qmatrix matrR = getZeroMatrix(dim);

    // ... where R_rc = (columm c of Z) . (column r of Q) = (row c of ZT) . (row r of QT) = <r|c>
    for (size_t r=0; r<dim; r++)
        for (size_t c=r; c<dim; c++)
            matrR[r][c] = getInnerProduct(matrQT[r], matrZT[c]);

    // create D = normalised diagonal of R
    qmatrix matrD = getZeroMatrix(dim);
    for (size_t i=0; i<dim; i++)
        matrD[i][i] = matrR[i][i] / abs(matrR[i][i]);

    // create U = Q D
    qmatrix matrU = matrQ * matrD;

    DEMAND( isApproxUnitary(matrU) );
    return matrU;
}


qmatrix getRandomDiagonalUnitary(int numQb) {
    DEMAND( numQb >= 1 );

    qmatrix matr = getZeroMatrix(getPow2(numQb));

    for (size_t i=0; i<matr.size(); i++)
        matr[i][i] = getExpI(getRandomReal(0,4*M_PI));

    return matr;
}


vector<qmatrix> getRandomKrausMap(int numQb, int numOps) {
    DEMAND( numOps >= 1 );

    // generate random unitaries
    vector<qmatrix> ops(numOps);
    for (auto& u : ops)
        u = getRandomUnitary(numQb);

    // generate random weights
    vector<qreal> weights(numOps);
    for (auto& w : weights)
        w = getRandomReal(0, 1);
        
    // normalise random weights
    qreal sum = 0;
    for (auto& w : weights)
        sum += w;
    for (auto& w : weights)
        w = sqrt(w/sum);
        
    // normalise unitaries according to weights
    for (int i=0; i<numOps; i++)
        ops[i] *= weights[i];

    DEMAND( isCompletelyPositiveTracePreserving(ops) );
    return ops;
}



/*
 * PAULIS
 */


PauliStr getRandomPauliStr(int numQubits) {

    std::string paulis = "";
    for (int i=0; i<numQubits; i++)
        paulis += "IXYZ"[getRandomInt(0,4)];

    return getPauliStr(paulis);
}


PauliStr getRandomPauliStr(vector<int> targs) {

    std::string paulis = "";
    for (size_t i=0; i<targs.size(); i++)
        paulis += "IXYZ"[getRandomInt(0,4)];

    return getPauliStr(paulis, targs);
}


PauliStr getRandomDiagPauliStr(int numQubits) {

    std::string paulis = "";
    for (int i=0; i<numQubits; i++)
        paulis += "IZ"[getRandomInt(0,2)];

    return getPauliStr(paulis);
}


PauliStrSum createRandomNonHermitianPauliStrSum(int numQubits, int numTerms) {

    vector<PauliStr> strings(numTerms);
    for (auto& str : strings)
        str = getRandomPauliStr(numQubits);

    vector<qcomp> coeffs(numTerms);
    for (auto& c : coeffs)
        c = getRandomComplex();

    return createPauliStrSum(strings, coeffs);
}


PauliStrSum createRandomPauliStrSum(int numQubits, int numTerms) {

    PauliStrSum out = createRandomNonHermitianPauliStrSum(numQubits, numTerms);

    for (qindex i=0; i<numTerms; i++)
        out.coeffs[i] = real(out.coeffs[i]);

    return out;
}
