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

#include <vector>

using std::vector;



/*
 * INTERNAL UTILS
 */

extern bool paulis_hasOddNumY(PauliStr str);
extern PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift);

void internal_applyFirstOrderTrotterRepetition(
    Qureg qureg, vector<int>& ketCtrls, vector<int>& braCtrls,
    vector<int>& states, PauliStrSum sum, qcomp angle, bool reverse
) {
    // apply each sum term as a gadget, in forward or reverse order
    for (qindex i=0; i<sum.numTerms; i++) {
        int j = reverse? sum.numTerms - i - 1 : i;
        qcomp coeff = sum.coeffs[j];
        PauliStr str = sum.strings[j];

        // effect |psi> -> exp(i angle * sum)|psi>
        qcomp arg = angle * coeff;
        localiser_statevec_anyCtrlPauliGadget(qureg, ketCtrls, states, str, arg);

        if (!qureg.isDensityMatrix)
            continue;

        // effect rho -> rho dagger(i angle * sum)
        arg *= paulis_hasOddNumY(str) ? 1 : -1;
        str = paulis_getShiftedPauliStr(str, qureg.numQubits);
        localiser_statevec_anyCtrlPauliGadget(qureg, braCtrls, states, str, arg);
    }
}

void internal_applyHigherOrderTrotterRepetition(
    Qureg qureg, vector<int>& ketCtrls, vector<int>& braCtrls,
    vector<int>& states, PauliStrSum sum, qcomp angle, int order
) {
    if (order == 1) {
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle, false);
    
    } else if (order == 2) {
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle/2, false);
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle/2, true);
    
    } else {
        qreal p = 1. / (4 - std::pow(4, 1./(order-1)));
        qcomp a = p * angle;
        qcomp b = (1-4*p) * angle;

        int lower = order - 2;
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower);
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower);
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, b, lower);
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower);
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower);
    }
}

void internal_applyAllTrotterRepetitions(
    Qureg qureg, int* controls, int* states, int numControls, 
    PauliStrSum sum, qcomp angle, int order, int reps
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
            qureg, ketCtrlsVec, braCtrlsVec, statesVec, sum, arg, order);

    /// @todo
    /// the accuracy of Trotterisation is greatly improved by randomisation
    /// or (even sub-optimal) grouping into commuting terms. Should we 
    /// implement these above or into another function?
}



/*
 * PAULI STR SUM GADGETS
 */

extern "C" {

void applyNonUnitaryTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qcomp angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);
    validate_trotterParams(qureg, order, reps, __func__);
    // sum is permitted to be non-Hermitian

    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, sum, angle, order, reps);
}

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, sum, angle, order, reps);
}

void applyControlledTrotterizedPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlAndPauliStrSumTargets(qureg, control, sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);
    
    internal_applyAllTrotterRepetitions(qureg, &control, nullptr, 1, sum, angle, order, reps);
}

void applyMultiControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlsAndPauliStrSumTargets(qureg, controls, numControls, sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    internal_applyAllTrotterRepetitions(qureg, controls, nullptr, numControls, sum, angle, order, reps);
}

void applyMultiStateControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlsAndPauliStrSumTargets(qureg, controls, numControls, sum, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr
    validate_trotterParams(qureg, order, reps, __func__);

    internal_applyAllTrotterRepetitions(qureg, controls, states, numControls, sum, angle, order, reps);
}

} // end de-mangler

void applyMultiControlledTrotterizedPauliStrSumGadget(Qureg qureg, vector<int> controls, PauliStrSum sum, qreal angle, int order, int reps) {

    applyMultiControlledTrotterizedPauliStrSumGadget(qureg, controls.data(), controls.size(), sum, angle, order, reps);
}

void applyMultiStateControlledTrotterizedPauliStrSumGadget(Qureg qureg, vector<int> controls, vector<int> states, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_controlsMatchStates(controls.size(), states.size(), __func__);

    applyMultiStateControlledTrotterizedPauliStrSumGadget(qureg, controls.data(), states.data(), controls.size(), sum, angle, order, reps);
}
