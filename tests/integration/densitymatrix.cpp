/** @file
 * Integration tests which combine many QuEST API functions in
 * a single test using Quregs which are too large to validate via
 * serial replication (like the unit tests perform), which instead
 * utilise known scale-invariant analytic properties.
 *
 * @author Tyson Jones
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/random.hpp"

#include <algorithm>
#include <string>
#include <vector>
#include <tuple>

using std::vector;
using std::tuple;
using namespace Catch::Matchers;


#define TEST_TAG \
    LABEL_INTEGRATION_TAG



void testDensityMatrixEvolution(Qureg psi, Qureg rho) {
    DEMAND( psi.numQubits == rho.numQubits );

    // set ||rho>> = |psi><psi|
    initRandomPureState(psi);
    initPureState(rho, psi);

    // we will check all alculations produced within 'eps' of expected
    qreal eps = std::max({(qreal) 1E-5, getTestAbsoluteEpsilon()});
    REQUIRE_THAT( calcPurity(rho),        WithinAbs(1, eps) );
    REQUIRE_THAT( calcPurity(psi),        WithinAbs(1, eps) );
    REQUIRE_THAT( calcTotalProb(rho),     WithinAbs(1, eps) );
    REQUIRE_THAT( calcTotalProb(psi),     WithinAbs(1, eps) );
    REQUIRE_THAT( calcFidelity(rho, psi), WithinAbs(1, eps) );

    // maximum size of tested any-target operators
    int maxNumCompMatrTargs = std::min({6, (int) psi.logNumAmpsPerNode, (int) rho.logNumAmpsPerNode});
    int maxNumDiagMatrTargs = std::min({8, psi.numQubits});
    int maxNumPauliStrTargs = psi.numQubits;
    int maxNumPauliGadTargs = psi.numQubits;
    int maxNumPhaseGadTargs = psi.numQubits;

    // below, we will apply every operation which invokes a unique backend
    // function, with random parameters which (if they were deterministic
    // and exhaustive) should access every edge-case including all templated
    // optimisations. Each operation (except a few anomalously demanding ones)
    // will be repeated with re-randomised parameters a total of 'numReps' times
    int numReps = 10 * psi.numQubits;

    vector<int> ctrls;
    vector<int> targs;
    vector<int> states;

    // apply CompMatr1
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,1);
        CompMatr1 matr = getCompMatr1(getRandomUnitary(1));
        applyMultiStateControlledCompMatr1(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
        applyMultiStateControlledCompMatr1(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
    }

    // apply CompMatr2
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 2,2);
        CompMatr2 matr = getCompMatr2(getRandomUnitary(2));
        applyMultiStateControlledCompMatr2(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
        applyMultiStateControlledCompMatr2(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
    }

    // apply CompMatr
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,maxNumCompMatrTargs);
        CompMatr matr = createCompMatr(targs.size());
        setCompMatr(matr, getRandomUnitary(targs.size()));
        applyMultiStateControlledCompMatr(psi, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        applyMultiStateControlledCompMatr(rho, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        destroyCompMatr(matr);
    }
        
    // apply DiagMatr1
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,1);
        DiagMatr1 matr = getDiagMatr1(getDiagonals(getRandomDiagonalUnitary(1)));
        applyMultiStateControlledDiagMatr1(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
        applyMultiStateControlledDiagMatr1(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], matr);
    }

    // apply DiagMatr2
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 2,2);
        DiagMatr2 matr = getDiagMatr2(getDiagonals(getRandomDiagonalUnitary(2)));
        applyMultiStateControlledDiagMatr2(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
        applyMultiStateControlledDiagMatr2(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1], matr);
    }

    // apply DiagMatr
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,maxNumDiagMatrTargs);
        DiagMatr matr = createDiagMatr(targs.size());
        setDiagMatr(matr, getDiagonals(getRandomDiagonalUnitary(targs.size())));
        applyMultiStateControlledDiagMatr(psi, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        applyMultiStateControlledDiagMatr(rho, ctrls.data(), states.data(), ctrls.size(), targs.data(), targs.size(), matr);
        destroyDiagMatr(matr);
    }

    // apply DiagMatr raised to power (real to retain unitarity)
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,maxNumDiagMatrTargs);
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
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 2,2);
        applyMultiStateControlledSwap(psi, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1]);
        applyMultiStateControlledSwap(rho, ctrls.data(), states.data(), ctrls.size(), targs[0], targs[1]);
    }

    // apply PauliStr
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,maxNumPauliStrTargs);
        PauliStr str = getRandomPauliStr(targs);
        applyMultiStateControlledPauliStr(psi, ctrls.data(), states.data(), ctrls.size(), str);
        applyMultiStateControlledPauliStr(rho, ctrls.data(), states.data(), ctrls.size(), str);
    }

    // apply Pauli gadget
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,maxNumPauliGadTargs);
        PauliStr str = getRandomPauliStr(targs);
        qreal phi = getRandomReal(-2 * 3.14, 2 * 3.14);
        applyMultiStateControlledPauliGadget(psi, ctrls.data(), states.data(), ctrls.size(), str, phi);
        applyMultiStateControlledPauliGadget(rho, ctrls.data(), states.data(), ctrls.size(), str, phi);
    }

    // apply phase gadet
    for (int r=0; r<numReps; r++) {
        auto [ctrls,states,targs] = getRandomVariNumCtrlsStatesTargs(psi.numQubits, 1,maxNumPhaseGadTargs);
        qreal phi = getRandomReal(-2 * 3.14, 2 * 3.14);
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
        vector<int> outcomes = getRandomOutcomes(targets.size());
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


/// @ingroup integrationtests
TEST_CASE( "density evolution", TEST_TAG ) {

    auto deployments = getSupportedDeployments();
    
    // try all combination of statevec and density-matrix deploments
    for (auto [rhoDeploy, rhoMPI, rhoGPU, rhoOMP] : deployments) {
        for (auto [psiDeploy, psiMPI, psiGPU, psiOMP] : deployments) {

            // some combinations are illegal
            if (psiMPI && !rhoMPI)
                continue; 

            // Qureg size determined by slowest deployment
            int numQubits = 6;
            if (rhoMPI && rhoMPI) numQubits = 12;
            if (rhoOMP && rhoOMP) numQubits = 12;
            if (rhoGPU && psiGPU) numQubits = 14;

            auto label = (
                "rho = " + rhoDeploy + ", " +
                "psi = " + psiDeploy + ", " +
                "qubits = " + std::to_string(numQubits));

            DYNAMIC_SECTION( label ) {

                Qureg psi = createCustomQureg(numQubits, 0, psiMPI,psiGPU,psiOMP);
                Qureg rho = createCustomQureg(numQubits, 1, rhoMPI,rhoGPU,rhoOMP);

                testDensityMatrixEvolution(psi, rho);

                destroyQureg(psi);
                destroyQureg(rho);
            }
        }
    }
}
