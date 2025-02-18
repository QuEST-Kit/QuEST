#include "quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/random.hpp"

#include <vector>
using std::vector;

#include <tuple>
using std::tuple;

#define TEST_TAG "[integration]"

using namespace Catch::Matchers;



void testDensityMatrixEvolution(Qureg psi, Qureg rho) {
    DEMAND( psi.numQubits == rho.numQubits );

    initRandomPureState(psi);
    initPureState(rho, psi);

    qreal eps = 1E-5;
    REQUIRE_THAT( calcPurity(rho),        WithinAbs(1, eps) );
    REQUIRE_THAT( calcPurity(psi),        WithinAbs(1, eps) );
    REQUIRE_THAT( calcTotalProb(rho),     WithinAbs(1, eps) );
    REQUIRE_THAT( calcTotalProb(psi),     WithinAbs(1, eps) );
    REQUIRE_THAT( calcFidelity(rho, psi), WithinAbs(1, eps) );

    int numReps = 10 * psi.numQubits;
    vector<int> ctrls;
    vector<int> targs;
    vector<int> states;

    // below, we will apply every operation which invokes a unique backend
    // function, with random parameters which (if they were deterministic
    // and exhaustive) should access every edge-case including all templated
    // optimisations.

    // // apply CompMatr1
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,1);
        CompMatr1 matr = getCompMatr1(getRandomUnitary(1));
        applyMultiStateControlledCompMatr1(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
        applyMultiStateControlledCompMatr1(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
    }

    // apply CompMatr2
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 2,2);
        CompMatr2 matr = getCompMatr2(getRandomUnitary(2));
        applyMultiStateControlledCompMatr2(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
        applyMultiStateControlledCompMatr2(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
    }

    // above causing segfault mpi

    // apply CompMatr upon 1-6 qubits
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,6);
        CompMatr matr = createCompMatr(targs.size());
        setCompMatr(matr, getRandomUnitary(targs.size()));
        applyMultiStateControlledCompMatr(psi, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        applyMultiStateControlledCompMatr(rho, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        destroyCompMatr(matr);
    }
        
    // apply DiagMatr1
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,1);
        DiagMatr1 matr = getDiagMatr1(getDiagonals(getRandomDiagonalUnitary(1)));
        applyMultiStateControlledDiagMatr1(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
        applyMultiStateControlledDiagMatr1(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
    }

    // apply DiagMatr2
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 2,2);
        DiagMatr2 matr = getDiagMatr2(getDiagonals(getRandomDiagonalUnitary(2)));
        applyMultiStateControlledDiagMatr2(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
        applyMultiStateControlledDiagMatr2(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
    }

    // apply DiagMatr upon 1-8 qubits
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,8);
        DiagMatr matr = createDiagMatr(targs.size());
        setDiagMatr(matr, getDiagonals(getRandomDiagonalUnitary(targs.size())));
        applyMultiStateControlledDiagMatr(psi, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        applyMultiStateControlledDiagMatr(rho, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        destroyDiagMatr(matr);
    }

    // apply DiagMatr raised to power (real to retain unitarity) upon 1-8 qubits
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,8);
        DiagMatr matr = createDiagMatr(targs.size());
        qreal expo = getRandomReal(0, 10);
        setDiagMatr(matr, getDiagonals(getRandomDiagonalUnitary(targs.size())));
        applyMultiStateControlledDiagMatrPower(psi, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr, expo);
        applyMultiStateControlledDiagMatrPower(rho, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr, expo);
        destroyDiagMatr(matr);
    }

    // FullStateDiagMatr and FullStateDiagMatrPower omitted for now

    // apply SWAP
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 2,2);
        applyMultiStateControlledSwap(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1]);
        applyMultiStateControlledSwap(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1]);
    }

    // apply PauliStr upon 1-to-all qubits
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,psi.numQubits);
        PauliStr str = getRandomPauliStr(targs);
        applyMultiStateControlledPauliStr(psi, ctrls.data(), states.data(), ctrls.size(), str);
        applyMultiStateControlledPauliStr(rho, ctrls.data(), states.data(), ctrls.size(), str);
    }

    // apply Pauli gadget upon 1-to-all qubits
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,psi.numQubits);
        PauliStr str = getRandomPauliStr(targs);
        qreal phi = getRandomReal(-2*M_PI, 2*M_PI);
        applyMultiStateControlledPauliGadget(psi, ctrls.data(), states.data(), ctrls.size(), str, phi);
        applyMultiStateControlledPauliGadget(rho, ctrls.data(), states.data(), ctrls.size(), str, phi);
    }

    // phase gadet upon 1-to-all qubits
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomCtrlsStatesTargs(psi.numQubits, 1,psi.numQubits);
        qreal phi = getRandomReal(-2*M_PI, 2*M_PI);
        applyMultiStateControlledPhaseGadget(psi, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), phi);
        applyMultiStateControlledPhaseGadget(rho, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), phi);
    }

    // confirm purity and normalisation was maintained
    REQUIRE_THAT( calcPurity(rho),        WithinAbs(1, eps) );
    REQUIRE_THAT( calcPurity(psi),        WithinAbs(1, eps) );
    REQUIRE_THAT( calcTotalProb(rho),     WithinAbs(1, eps) );
    REQUIRE_THAT( calcTotalProb(psi),     WithinAbs(1, eps) );

    // confirm states agree
    REQUIRE_THAT( calcFidelity(rho, psi), WithinAbs(1, eps) );

    // confirm expectation values agree
    for (int r=0; r<numReps; r++) {
        PauliStr str = getRandomPauliStr(rho.numQubits);
        qreal psiExpec = calcExpecPauliStr(psi, str);
        qreal rhoExpec = calcExpecPauliStr(rho, str);
        REQUIRE_THAT( psiExpec, WithinAbs(rhoExpec, eps) );
    }

    // confirm all one-qubit probabilities agree
    for (int q=0; q<psi.numQubits; q++) {
        qreal psiProb = calcProbOfQubitOutcome(psi, q, 0);
        qreal rhoProb = calcProbOfQubitOutcome(rho, q, 0);
        REQUIRE_THAT( psiProb, WithinAbs(rhoProb, eps) );
    }

    // confirm some multi-qubit (1-6) probabilities agree
    for (int r=0; r<numReps; r++) {
        int numTargets = getRandomInt(1, 6+1);
        vector<int> targets = getRandomSubRange(0, psi.numQubits, numTargets);
        vector<int> outcomes = getRandomInts(0, 1+1, targets.size());
        qreal psiProb = calcProbOfMultiQubitOutcome(psi, targets.data(), outcomes.data(), numTargets);
        qreal rhoProb = calcProbOfMultiQubitOutcome(rho, targets.data(), outcomes.data(), numTargets);
        REQUIRE_THAT( psiProb, WithinAbs(rhoProb, eps) );
    }

    // confirm some basis state probabilities agree
    for (int r=0; r<numReps; r++) {
        int ind = getRandomInt(0, psi.numAmps);
        qreal psiProb = calcProbOfBasisState(psi, ind);
        qreal rhoProb = calcProbOfBasisState(rho, ind);
        REQUIRE_THAT( psiProb, WithinAbs(rhoProb, eps) );
    }

    // mix all canonical channels...
    for (int r=0; r<numReps; r++) {
        vector<int> targs = getRandomSubRange(0, rho.numQubits, 2);
        mixDephasing(rho, targs[0], getRandomReal(0,1/2.));
        mixTwoQubitDephasing(rho, targs[0], targs[1], getRandomReal(0,3/4.));
        mixDepolarising(rho, targs[0], getRandomReal(0,3/4.));
        mixTwoQubitDepolarising(rho, targs[0], targs[1], getRandomReal(0,15/16.));
        mixDamping(rho, targs[0], getRandomReal(0,1));
        mixPaulis(rho, targs[0], getRandomReal(0,1/10.), getRandomReal(0,1/20.), getRandomReal(0,1/30.));
    }

    // to confirm total probability unaffected
    REQUIRE_THAT( calcTotalProb(rho), WithinAbs(1, eps) );

    // confirm damping restores outcome probability
    int qubitInd = 0;
    int qubitOutcome = 0;
    qreal dampProb = 1;
    qreal outcomeProb = 1;
    mixDamping(rho, qubitInd, dampProb);
    REQUIRE_THAT( calcProbOfQubitOutcome(rho, qubitInd, qubitOutcome), WithinAbs(outcomeProb, eps) );
}




TEST_CASE( "density evolution", "[integration]" ) {
    
    // TODO:
    // am unsure about how to make these quregs - should it be all combinations for example?!


    int numQubits = 12;
    Qureg psi = createForcedQureg(numQubits);
    Qureg rho = createForcedDensityQureg(numQubits);

    testDensityMatrixEvolution(psi, rho);


    destroyQureg(psi);
    destroyQureg(rho);
}