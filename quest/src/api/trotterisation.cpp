/** @file
 * API definitions for functions which involve Trotterising
 * exponential operators, such as PauliStrSum gadgets, and
 * so are inherently approximate.
 * 
 * @author Tyson Jones
 */

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/paulilogic.hpp"
#include "quest/src/core/errors.hpp"

#include <vector>

using std::vector;



/*
 * INTERNAL UTILS
 */

void internal_applyFirstOrderTrotterRepetition(
    Qureg qureg, vector<int>& ketCtrls, vector<int>& braCtrls,
    vector<int>& states, PauliStrSum sum, qcomp angle, bool onlyLeftApply, bool reverse
) {
    // apply each sum term as a gadget, in forward or reverse order
    for (qindex i=0; i<sum.numTerms; i++) {
        int j = reverse? sum.numTerms - i - 1 : i;
        qcomp coeff = sum.coeffs[j];
        PauliStr str = sum.strings[j];

        // effect |psi> -> exp(i angle * coeff * term)|psi>
        qcomp arg = angle * coeff;
        localiser_statevec_anyCtrlPauliGadget(qureg, ketCtrls, states, str, arg);

        // term finished upon statevector 
        if (!qureg.isDensityMatrix)
            continue;

        // Linbladian propagator is only ever pre-multiplied
        if (onlyLeftApply)
            continue;

        // effect rho -> rho exp(i angle * coeff * term)^dagger via linearised
        //    ||rho>> -> conj(exp(i angle * coeff * term)) (x) I ||rho>>
        //             = exp(- i conj(angle coeff) sign term) (x) I ||rho>>
        arg = - std::conj(arg) * paulis_getSignOfPauliStrConj(str);
        str = paulis_getShiftedPauliStr(str, qureg.numQubits);
        localiser_statevec_anyCtrlPauliGadget(qureg, braCtrls, states, str, arg);
    }
}

void internal_applyHigherOrderTrotterRepetition(
    Qureg qureg, vector<int>& ketCtrls, vector<int>& braCtrls,
    vector<int>& states, PauliStrSum sum, qcomp angle, int order, bool onlyLeftApply
) {
    if (order == 1) {
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle, onlyLeftApply, false);
    
    } else if (order == 2) {
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle/2, onlyLeftApply, false);
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle/2, onlyLeftApply, true);
    
    } else {
        qreal p = 1. / (4 - std::pow(4, 1./(order-1)));
        qcomp a = p * angle;
        qcomp b = (1-4*p) * angle;

        int lower = order - 2;
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, onlyLeftApply); // angle -> a
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, onlyLeftApply);
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, b, lower, onlyLeftApply); // angle -> b
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, onlyLeftApply);
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, onlyLeftApply);
    }
}

void internal_applyAllTrotterRepetitions(
    Qureg qureg, int* controls, int* states, int numControls, 
    PauliStrSum sum, qcomp angle, int order, int reps, bool onlyLeftApply
) {
    // exp(i angle sum) = identity when angle=0
    if (angle == qcomp(0,0))
        return;

    // prepare control-qubit lists once for all invoked gadgets below
    auto ketCtrlsVec = util_getVector(controls, numControls);
    auto braCtrlsVec = (qureg.isDensityMatrix)? util_getBraQubits(ketCtrlsVec, qureg) : vector<int>{};
    auto statesVec = util_getVector(states, numControls);

    qcomp arg = angle / reps;

    // perform carefully-ordered sequence of gadgets
    for (int r=0; r<reps; r++)
        internal_applyHigherOrderTrotterRepetition(
            qureg, ketCtrlsVec, braCtrlsVec, statesVec, sum, arg, order, onlyLeftApply);

    /// @todo
    /// the accuracy of Trotterisation is greatly improved by randomisation
    /// or (even sub-optimal) grouping into commuting terms. Should we 
    /// implement these above or into another function?
}

qindex internal_getNumTotalSuperPropagatorTerms(PauliStrSum hamil, PauliStrSum* jumps, int numJumps) {

    // this function returns 0 to indicate an overflow, which will never
    // be confused for the correct non-overflowed output because hamil.numTerms>0
    int OVERFLOW_FLAG = 0;

    if (util_willProdOverflow({2,hamil.numTerms}))
        return OVERFLOW_FLAG;
        
    // I (x) H + conj(H) (x) I
    qindex numTerms = 2 * hamil.numTerms;

    for (int i=0; i<numJumps; i++) {
        qindex n = jumps[i].numTerms;

        if (util_willProdOverflow({n,n,3}))
            return OVERFLOW_FLAG;
        if (util_willSumOverflow({numTerms, 3*n*n}))
            return OVERFLOW_FLAG;

        // conj(J) (x) J has n^2 terms
        numTerms += n * n;

        // I (x) (adj(J) . J ) + conj(...) (x) I is bounded by 2*n^2 terms
        numTerms += 2 * paulis_getNumTermsInPauliStrSumProdOfAdjointWithSelf(jumps[i]);
    }

    // indicate no overflow
    return numTerms;
}



/*
 * PAULI STR SUM GADGETS
 */

extern "C" {

void applyTrotterizedNonUnitaryPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);
    validate_trotterParams(qureg, order, reps, __func__);
    // sum is permitted to be non-Hermitian

    // |psi> -> U |psi>, rho -> U rho U^dagger
    bool onlyLeftApply = false;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, sum, angle, order, reps, onlyLeftApply);
}

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    bool onlyLeftApply = false;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, sum, angle, order, reps, onlyLeftApply);
}

void applyTrotterizedControlledPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlAndPauliStrSumTargets(qureg, control, sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);
    
    bool onlyLeftApply = false;
    internal_applyAllTrotterRepetitions(qureg, &control, nullptr, 1, sum, angle, order, reps, onlyLeftApply);
}

void applyTrotterizedMultiControlledPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlsAndPauliStrSumTargets(qureg, controls, numControls, sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    bool onlyLeftApply = false;
    internal_applyAllTrotterRepetitions(qureg, controls, nullptr, numControls, sum, angle, order, reps, onlyLeftApply);
}

void applyTrotterizedMultiStateControlledPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlsAndPauliStrSumTargets(qureg, controls, numControls, sum, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr
    validate_trotterParams(qureg, order, reps, __func__);

    bool onlyLeftApply = false;
    internal_applyAllTrotterRepetitions(qureg, controls, states, numControls, sum, angle, order, reps, onlyLeftApply);
}

} // end de-mangler

void applyTrotterizedMultiControlledPauliStrSumGadget(Qureg qureg, vector<int> controls, PauliStrSum sum, qreal angle, int order, int reps) {

    applyTrotterizedMultiControlledPauliStrSumGadget(qureg, controls.data(), controls.size(), sum, angle, order, reps);
}

void applyTrotterizedMultiStateControlledPauliStrSumGadget(Qureg qureg, vector<int> controls, vector<int> states, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_controlsMatchStates(controls.size(), states.size(), __func__);

    applyTrotterizedMultiStateControlledPauliStrSumGadget(qureg, controls.data(), states.data(), controls.size(), sum, angle, order, reps);
}



/*
 * CLOSED TIME EVOLUTION
 */

extern "C" {

void applyTrotterizedUnitaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal time, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(hamil, __func__);
    validate_pauliStrSumTargets(hamil, qureg, __func__);
    validate_pauliStrSumIsHermitian(hamil, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    // exp(-i t H) = exp(x i H) | x=-t
    qcomp angle = - time;
    bool onlyLeftApply = false;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, hamil, angle, order, reps, onlyLeftApply);
}

void applyTrotterizedImaginaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal tau, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(hamil, __func__);
    validate_pauliStrSumTargets(hamil, qureg, __func__);
    validate_pauliStrSumIsHermitian(hamil, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    // exp(-tau H) = exp(x i H) | x=tau*i
    qcomp angle = qcomp(0, tau);
    bool onlyLeftApply = false;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, hamil, angle, order, reps, onlyLeftApply);
}

} // end de-mangler



/*
 * OPEN TIME EVOLUTION
 */

extern "C" {

void applyTrotterizedNoisyTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal* damps, PauliStrSum* jumps, int numJumps, qreal time, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_pauliStrSumFields(hamil, __func__);
    validate_pauliStrSumTargets(hamil, qureg, __func__);
    validate_pauliStrSumIsHermitian(hamil, __func__);
    validate_trotterParams(qureg, order, reps, __func__);
    validate_lindbladJumpOps(jumps, numJumps, qureg, __func__);
    validate_lindbladDampingRates(damps, numJumps, __func__);
    
    qindex numSuperTerms = internal_getNumTotalSuperPropagatorTerms(hamil, jumps, numJumps); // 0 indicates overflow
    validate_numLindbladSuperPropagatorTerms(numSuperTerms, __func__);

    // validate memory allocations for all super-propagator terms
    vector<PauliStr> superStrings;
    vector<qcomp> superCoeffs;
    auto callbackString = [&]() { validate_tempAllocSucceeded(false, numSuperTerms, sizeof(PauliStr), __func__); };
    auto callbackCoeff  = [&]() { validate_tempAllocSucceeded(false, numSuperTerms, sizeof(qcomp),    __func__); };
    util_tryAllocVector(superStrings, numSuperTerms, callbackString);
    util_tryAllocVector(superCoeffs,  numSuperTerms, callbackCoeff);

    qindex superTermInd = 0;

    // collect -i[H,rho] terms
    for (qindex n=0; n<hamil.numTerms; n++) {
        PauliStr oldStr = hamil.strings[n];
        qcomp oldCoeff = hamil.coeffs[n];

        // term of -i Id (x) H
        superStrings[superTermInd] = oldStr;
        superCoeffs [superTermInd] = -1_i * oldCoeff;
        superTermInd++;

        // term of i conj(H) (x) I
        superStrings[superTermInd] = paulis_getShiftedPauliStr(oldStr, qureg.numQubits);
        superCoeffs [superTermInd] = 1_i * paulis_getSignOfPauliStrConj(oldStr) * std::conj(oldCoeff);
        superTermInd++;
    }

    // below we bind superStrings/Coeffs to a spoofed PauliStrSum to pass to paulis functions
    PauliStrSum temp;
    int flagForDebugSafety = -1;
    temp.isApproxHermitian = &flagForDebugSafety;

    // collect jump terms
    for (int n=0; n<numJumps; n++) {

        // damp  conj(J) (x) J
        temp.strings = &superStrings[superTermInd];
        temp.coeffs = &superCoeffs[superTermInd];
        temp.numTerms = jumps[n].numTerms * jumps[n].numTerms;
        superTermInd += temp.numTerms;
        paulis_setPauliStrSumToScaledTensorProdOfConjWithSelf(temp, damps[n], jumps[n], qureg.numQubits);

        // -damp/2  I (x) (adj(J) . J)
        temp.strings = &superStrings[superTermInd];
        temp.coeffs = &superCoeffs[superTermInd];
        temp.numTerms = paulis_getNumTermsInPauliStrSumProdOfAdjointWithSelf(jumps[n]);
        superTermInd += temp.numTerms;
        paulis_setPauliStrSumToScaledProdOfAdjointWithSelf(temp, -damps[n]/2, jumps[n]);

        // -damp/2 conj(adj(J) . J) (x) I = conj(above) when damp is real
        PauliStrSum temp2;
        temp2.strings = &superStrings[superTermInd];
        temp2.coeffs = &superCoeffs[superTermInd];
        temp2.numTerms = temp.numTerms;
        superTermInd += temp2.numTerms;
        paulis_setPauliStrSumToShiftedConj(temp2, temp, qureg.numQubits);
    }

    // defensively check we didn't write too few (or too many, though that'd segfault
    // above) Lindblad terms, in case the above code changes when jump ops are generalised
    if (superTermInd != numSuperTerms)
        error_unexpectedNumLindbladSuperpropTerms();

    // pass superpropagator terms as temporary PauliStrSum
    PauliStrSum superSum; 
    superSum.numTerms = numSuperTerms;
    superSum.strings = superStrings.data();
    superSum.coeffs = superCoeffs.data();
    superSum.isApproxHermitian = nullptr; // will not be queried

    // effect exp(t S) = exp(x i S) | x=-i*time, left-multiplying only
    qcomp angle = qcomp(0, -time);
    bool onlyLeftApply = true;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, superSum, angle, order, reps, onlyLeftApply);
}

} // end de-mangler
