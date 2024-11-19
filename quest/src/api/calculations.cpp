/** @file
 * API definitions for calculating properties of quantum states,
 * such as probabilities and expectation values.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/calculations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/bitwise.hpp"

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions



/*
 * C++ ONLY FUNCTIONS
 *
 * which are not directly C-compatible because of limited
 * interoperability of the qcomp type. See calculations.h 
 * for more info. We here define a C++-only signature (with
 * name-mangling), and a C-friendly wrapper which passes by
 * pointer; the C-friendly interface in wrappers.h which itself
 * wrap this.
 */


qcomp calcInnerProduct(Qureg quregA, Qureg quregB) {
    validate_quregFields(quregA, __func__);
    validate_quregFields(quregB, __func__);
    validate_quregsCanBeProducted(quregA, quregB, __func__);

    // <A|B> or Tr(A^dagger B)
    if (quregA.isDensityMatrix == quregB.isDensityMatrix)
        return localiser_statevec_calcInnerProduct(quregA, quregB);

    return (quregA.isDensityMatrix)?
        localiser_densmatr_calcFidelityWithPureState(quregA, quregB, true):  // <B|A^dagger|B>
        localiser_densmatr_calcFidelityWithPureState(quregB, quregA, false); // <A|B|A>
}
extern "C" void _wrap_calcInnerProduct(Qureg bra, Qureg ket, qcomp* out) {

    *out = calcInnerProduct(bra, ket);
}


qcomp calcExpecNonHermitianPauliStrSum(Qureg qureg, PauliStrSum sum) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}
extern "C" void _wrap_calcExpecNonHermitianPauliStrSum(qcomp* out, Qureg qureg, PauliStrSum sum) {

    *out = calcExpecNonHermitianPauliStrSum(qureg, sum);
}


qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr) {

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}
extern "C" void _wrap_calcExpecNonHermitianFullStateDiagMatr(qcomp* out, Qureg qureg, FullStateDiagMatr matr) {

    *out = calcExpecNonHermitianFullStateDiagMatr(qureg, matr);
}



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
extern "C" {



/*
 * EXPECTED VALUES
 */


qreal calcExpecPauliStr(Qureg qureg, PauliStr str) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrTargets(qureg, str, __func__);

    qcomp value = (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecPauliStr(qureg, str):
        localiser_statevec_calcExpecPauliStr(qureg, str);

    validate_expecValIsReal(value, qureg.isDensityMatrix, __func__);
    return std::real(value);
}


qreal calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum) {

    // TODO
    valdidate_pauliStrSumIsHermitian(sum, __func__);

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}


qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr) {

    // TODO
    validate_matrixIsHermitian(matr, __func__);

    // TODO
    error_functionNotImplemented(__func__);
    return -1;
}


/*
 * PROBABILITIES
 */


qreal calcProbOfBasisState(Qureg qureg, qindex index) {
    validate_quregFields(qureg, __func__);
    validate_basisStateIndex(qureg, index, __func__);

    if (qureg.isDensityMatrix)
        index *= 1 + powerOf2(qureg.numQubits);

    qcomp amp = localiser_statevec_getAmp(qureg, index);
    qreal prob = std::norm(amp);
    return prob;
}


qreal calcProbOfQubitOutcome(Qureg qureg, int qubit, int outcome) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, qubit, __func__);
    validate_measurementOutcomeIsValid(outcome, __func__);

    int numQubits = 1;
    return calcProbOfMultiQubitOutcome(qureg, &qubit, &outcome, numQubits);
}


qreal calcProbOfMultiQubitOutcome(Qureg qureg, int* qubits, int* outcomes, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_measurementOutcomesAreValid(outcomes, numQubits, __func__);

    auto qubitVec = util_getVector(qubits, numQubits);
    auto outcomeVec = util_getVector(outcomes, numQubits);

    return (qureg.isDensityMatrix)?
        localiser_densmatr_calcProbOfMultiQubitOutcome(qureg, qubitVec, outcomeVec):
        localiser_statevec_calcProbOfMultiQubitOutcome(qureg, qubitVec, outcomeVec);
}


void calcProbsOfAllMultiQubitOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {
    validate_quregFields(qureg, __func__);
    validate_targets(qureg, qubits, numQubits, __func__);
    validate_measurementOutcomesFitInGpuMem(qureg, numQubits, __func__);

    auto qubitVec = util_getVector(qubits, numQubits);

    (qureg.isDensityMatrix)?
        localiser_densmatr_calcProbsOfAllMultiQubitOutcomes(outcomeProbs, qureg, qubitVec):
        localiser_statevec_calcProbsOfAllMultiQubitOutcomes(outcomeProbs, qureg, qubitVec);
}


qreal calcTotalProb(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    return (qureg.isDensityMatrix)?
        localiser_densmatr_calcTotalProb(qureg): 
        localiser_statevec_calcTotalProb(qureg);
}



/*
 * MEASURES
 */


qreal calcPurity(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // assuming Hermiticity, Tr(rho^2) = sum_ij |rho_ij|^2,
    // and Tr(|x><x|x><x|) = (sum_i |x_i|^2)^2
    qreal prob = localiser_statevec_calcTotalProb(qureg);
    return (qureg.isDensityMatrix)? prob : prob * prob;
}


qreal calcFidelity(Qureg quregA, Qureg quregB) {
    validate_quregFields(quregA, __func__);
    validate_quregFields(quregB, __func__);
    validate_quregsCanBeProducted(quregA, quregB, __func__);

    bool isDensA = quregA.isDensityMatrix;
    bool isDensB = quregB.isDensityMatrix;

    // F(rho,sigma) not yet supported
    if (isDensA && isDensB)
        validate_throwErrorBecauseCalcFidOfDensMatrNotYetImplemented(__func__);

    // |<A|B>|^2
    if (!isDensA && !isDensB) {
        qcomp prod = localiser_statevec_calcInnerProduct(quregA, quregB);
        return std::norm(prod);
    }

    // Re[<B|A|B>] or Re[<A|B|A>]
    qcomp fid = (quregA.isDensityMatrix)?
        localiser_densmatr_calcFidelityWithPureState(quregA, quregB, false): // no conj
        localiser_densmatr_calcFidelityWithPureState(quregB, quregA, false); // no conj

    validate_fidelityIsReal(fid, __func__);
    return std::real(fid);
}


qreal calcDistance(Qureg quregA, Qureg quregB) {
    validate_quregFields(quregA, __func__);
    validate_quregFields(quregB, __func__);
    validate_quregsCanBeProducted(quregA, quregB, __func__);

    bool isDensA = quregA.isDensityMatrix;
    bool isDensB = quregB.isDensityMatrix;

    // Hilbert-Schmidt = sqrt( Tr((A-B)(A-B)^dagger) = sqrt(sum_ij |A_ij - B_ij|^2)
    if (isDensA && isDensB) {
        qreal dif = localiser_densmatr_calcHilbertSchmidtDistance(quregA, quregB);
        return sqrt(dif);
    }

    // Bures = sqrt(2 - 2 |<A|B>|) (even when unnormalised)
    if (!isDensA && !isDensB) {
        qcomp prod = localiser_statevec_calcInnerProduct(quregA, quregB);
        qreal mag = std::abs(prod);

        validate_buresDistanceInnerProdIsNormalised(mag, __func__);
        mag = (mag > 1)? 1 : mag; // forgive eps error to avoid complex
        return std::sqrt(2 - 2 * mag);
    }

    // purified distance = sqrt(1 - <psi|rho|psi>)
    qcomp fid = (quregA.isDensityMatrix)?
        localiser_densmatr_calcFidelityWithPureState(quregA, quregB, false): // no conj
        localiser_densmatr_calcFidelityWithPureState(quregB, quregA, false); // no conj

    validate_purifiedDistanceIsNormalised(fid, __func__); 
    qreal re = std::real(fid);
    re = (re > 1)? 1 : re; // forgive eps error to avoid complex
    return std::sqrt(1 - re);
}



/*
 * PARTIAL TRACE
 */

Qureg calcPartialTrace(Qureg qureg, int* qubits, int numQubits) {

    // BIG TODO
    error_functionNotImplemented(__func__);    
    return createQureg(0);
}


} // end de-name mangler
