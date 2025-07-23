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
    vector<int>& states, PauliStrSum sum, qcomp angle, bool postmultiply, bool reverse
) {
    // apply each sum term as a gadget, in forward or reverse order
    for (qindex i=0; i<sum.numTerms; i++) {
        int j = reverse? sum.numTerms - i - 1 : i;
        qcomp coeff = sum.coeffs[j];
        PauliStr str = sum.strings[j];

        // effect |psi> -> exp(i angle * sum)|psi>
        qcomp arg = angle * coeff;
        localiser_statevec_anyCtrlPauliGadget(qureg, ketCtrls, states, str, arg);

        // term finished upon statevector 
        if (!qureg.isDensityMatrix)
            continue;

        // Linbladian propagator is only ever pre-multiplied
        if (!postmultiply)
            continue;

        // effect rho -> rho dagger(i angle * sum)
        arg *= paulis_hasOddNumY(str) ? 1 : -1;
        str = paulis_getShiftedPauliStr(str, qureg.numQubits);
        localiser_statevec_anyCtrlPauliGadget(qureg, braCtrls, states, str, arg);
    }
}

void internal_applyHigherOrderTrotterRepetition(
    Qureg qureg, vector<int>& ketCtrls, vector<int>& braCtrls,
    vector<int>& states, PauliStrSum sum, qcomp angle, int order, bool postmultiply
) {
    if (order == 1) {
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle, postmultiply, false);
    
    } else if (order == 2) {
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle/2, postmultiply, false);
        internal_applyFirstOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, angle/2, postmultiply, true);
    
    } else {
        qreal p = 1. / (4 - std::pow(4, 1./(order-1)));
        qcomp a = p * angle;
        qcomp b = (1-4*p) * angle;

        int lower = order - 2;
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, postmultiply); // angle -> a
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, postmultiply);
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, b, lower, postmultiply); // angle -> b
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, postmultiply);
        internal_applyHigherOrderTrotterRepetition(qureg, ketCtrls, braCtrls, states, sum, a, lower, postmultiply);
    }
}

void internal_applyAllTrotterRepetitions(
    Qureg qureg, int* controls, int* states, int numControls, 
    PauliStrSum sum, qcomp angle, int order, int reps, bool postmultiply
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
            qureg, ketCtrlsVec, braCtrlsVec, statesVec, sum, arg, order, postmultiply);

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

    bool postmultiply = true;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, sum, angle, order, reps, postmultiply);
}

void applyTrotterizedPauliStrSumGadget(Qureg qureg, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    bool postmultiply = true;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, sum, angle, order, reps, postmultiply);
}

void applyControlledTrotterizedPauliStrSumGadget(Qureg qureg, int control, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlAndPauliStrSumTargets(qureg, control, sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);
    
    bool postmultiply = true;
    internal_applyAllTrotterRepetitions(qureg, &control, nullptr, 1, sum, angle, order, reps, postmultiply);
}

void applyMultiControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int numControls, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlsAndPauliStrSumTargets(qureg, controls, numControls, sum, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    bool postmultiply = true;
    internal_applyAllTrotterRepetitions(qureg, controls, nullptr, numControls, sum, angle, order, reps, postmultiply);
}

void applyMultiStateControlledTrotterizedPauliStrSumGadget(Qureg qureg, int* controls, int* states, int numControls, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);
    validate_controlsAndPauliStrSumTargets(qureg, controls, numControls, sum, __func__);
    validate_controlStates(states, numControls, __func__); // permits states==nullptr
    validate_trotterParams(qureg, order, reps, __func__);

    bool postmultiply = true;
    internal_applyAllTrotterRepetitions(qureg, controls, states, numControls, sum, angle, order, reps, postmultiply);
}

} // end de-mangler

void applyMultiControlledTrotterizedPauliStrSumGadget(Qureg qureg, vector<int> controls, PauliStrSum sum, qreal angle, int order, int reps) {

    applyMultiControlledTrotterizedPauliStrSumGadget(qureg, controls.data(), controls.size(), sum, angle, order, reps);
}

void applyMultiStateControlledTrotterizedPauliStrSumGadget(Qureg qureg, vector<int> controls, vector<int> states, PauliStrSum sum, qreal angle, int order, int reps) {
    validate_controlsMatchStates(controls.size(), states.size(), __func__);

    applyMultiStateControlledTrotterizedPauliStrSumGadget(qureg, controls.data(), states.data(), controls.size(), sum, angle, order, reps);
}



/*
 * CLOSED TIME EVOLUTION
 */

void applyTrotterizedUnitaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal time, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(hamil, __func__);
    validate_pauliStrSumTargets(hamil, qureg, __func__);
    validate_pauliStrSumIsHermitian(hamil, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    // exp(-i t H) = exp(x i H) | x=-t
    qcomp angle = - time;
    bool postmultiply = true;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, hamil, angle, order, reps, postmultiply);
}

void applyTrotterizedImaginaryTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal tau, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(hamil, __func__);
    validate_pauliStrSumTargets(hamil, qureg, __func__);
    validate_pauliStrSumIsHermitian(hamil, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    // exp(-tau H) = exp(x i H) | x=tau*i
    qcomp angle = qcomp(0, tau);
    bool postmultiply = true;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, hamil, angle, order, reps, postmultiply);
}



/*
 * OPEN TIME EVOLUTION
 */

extern PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift);
extern PauliStr paulis_getKetAndBraPauliStr(PauliStr str, Qureg qureg);

extern "C" {

void applyTrotterizedPauliNoisyTimeEvolution(Qureg qureg, PauliStrSum hamil, qreal* damps, PauliStr* jumps, int numJumps, qreal time, int order, int reps) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_pauliStrSumFields(hamil, __func__);
    validate_pauliStrSumTargets(hamil, qureg, __func__);
    validate_pauliStrSumIsHermitian(hamil, __func__);
    validate_trotterParams(qureg, order, reps, __func__);

    // TODO: validate numJumps

    // TODO; validate damps (looped); must be non-negative (0 legally disables one)

    for (int n=0; n<numJumps; n++)
        validate_pauliStrTargets(qureg, jumps[n], __func__);

    // when all jump operators are Paulis, the linblad superop simplifies to:
    // L = -i (Id (x) H - conj(H) (x) I) + sum_k gamma_k conj(L_k) (x) L_k - (sum_k gamma_k) Id
    // where the final term commutes with everything and can just be brought out the front of
    // the Trotter circuit; it becomes just a final scalar factor upon the state.
    //
    // So we need merely prepare a new PauliStrSum which contains
    //   - old hamil times -i
    //   - conj and shifted hamil times -i
    //   - gamma_k conj(L_k) (x) L_k = +- gamma_k L_k (x) L_k depending on L_k Y parity
    // then after effecting that, scale the state

    vector<PauliStr> newStrings;
    vector<qcomp> newCoeffs;

    // premature optimisation
    qindex numNewTerms = 2 * hamil.numTerms + numJumps;
    newStrings.reserve(numNewTerms);
    newCoeffs.reserve(numNewTerms);

    // collect -i[H,rho] terms
    for (qindex n=0; n<hamil.numTerms; n++) {
        PauliStr oldStr = hamil.strings[n];
        qcomp oldCoeff = hamil.coeffs[n];

        // term of -i Id (x) H
        newStrings.push_back(oldStr);
        newCoeffs.push_back(-1_i * oldCoeff);

        // term of i conj(H) (x) I
        newStrings.push_back(paulis_getShiftedPauliStr(oldStr, qureg.numQubits));
        newCoeffs.push_back(1_i * (paulis_hasOddNumY(oldStr) ? -1 : 1) * std::conj(oldCoeff));
    }

    // collect jump terms
    for (int n=0; n<numJumps; n++) {

        // gamma_k conj(L_k) (x) L_k
        newStrings.push_back(paulis_getKetAndBraPauliStr(jumps[n], qureg));
        newCoeffs.push_back(damps[n] * (paulis_hasOddNumY(jumps[n]) ? -1 : 1));
    }

    // spoof a PauliStrSum to avoid superfluous alloc
    PauliStrSum temp; 
    temp.numTerms = numNewTerms;
    temp.strings = newStrings.data();
    temp.coeffs = newCoeffs.data();
    temp.isApproxHermitian = nullptr; // will not be queried

    // effect exp(t S) = exp(x i S) | x=-i*time, premultiplying only
    qcomp angle = qcomp(0, -time);
    bool postmultiply = false;
    internal_applyAllTrotterRepetitions(qureg, nullptr, nullptr, 0, temp, angle, order, reps, postmultiply);

    // scale by exp(- time sum_k gamma_k)
    qreal dampSum = 0;
    for (int n=0; n<numJumps; n++)
        dampSum += damps[n];
    qcomp fac = std::exp(- time * dampSum);
    localiser_statevec_setQuregToSuperposition(fac, qureg, 0, qureg, 0, qureg);
}

}
