/** @file
 * API definitions for calculating properties of quantum states,
 * such as probabilities, expectation values and partial traces.
 * 
 * @author Tyson Jones
 * @author Balint Koczor (prototyped v3 calcHilbertSchmidtDistance and calcDensityInnerProduct)
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/calculations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/errors.hpp"

#include <vector>
using std::vector;



/*
 * INTERNAL FUNCTIONS
 */


extern Qureg validateAndCreateCustomQureg(
    int numQubits, int isDensMatr, int useDistrib, 
    int useGpuAccel, int useMultithread, const char* caller);



/*
 * C++ ONLY FUNCTIONS
 *
 * which are not directly C-compatible because of limited
 * interoperability of the qcomp type. See calculations.h 
 * for more info. We here define a C++-only signature (with
 * name-mangling), and a C-friendly wrapper which passes by
 * pointer; the C-friendly interface in wrappers.h which itself
 * wrap this. This excludes additional convenience C++-only 
 * overloads defined at the bottom of this file.
 */


qcomp calcInnerProduct(Qureg quregA, Qureg quregB) {
    validate_quregFields(quregA, __func__);
    validate_quregFields(quregB, __func__);
    validate_quregsCanBeProducted(quregA, quregB, __func__);

    // <A|B> or Tr(A^dagger B) = <<A|B>>
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
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);

    qcomp value = (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecPauliStrSum(qureg, sum):
        localiser_statevec_calcExpecPauliStrSum(qureg, sum);

    return value;
}
extern "C" void _wrap_calcExpecNonHermitianPauliStrSum(qcomp* out, Qureg qureg, PauliStrSum sum) {

    *out = calcExpecNonHermitianPauliStrSum(qureg, sum);
}


qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, true, __func__);

    return calcExpecNonHermitianFullStateDiagMatrPower(qureg, matrix, 1); // harmlessly re-validates
}
extern "C" void _wrap_calcExpecNonHermitianFullStateDiagMatr(qcomp* out, Qureg qureg, FullStateDiagMatr matr) {

    *out = calcExpecNonHermitianFullStateDiagMatr(qureg, matr);
}


qcomp calcExpecNonHermitianFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, true, __func__);

    // this function never uses the qreal-pow overload (because we make
    // no assumption matrix is Hermitian i.e. real), and instead always
    // uses the relatively numerically inaccurate qcomp overload
    const bool useRealPow = false;

    return (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecFullStateDiagMatr(qureg, matrix, exponent, useRealPow):
        localiser_statevec_calcExpecFullStateDiagMatr(qureg, matrix, exponent, useRealPow);
}
extern "C" void _wrap_calcExpecNonHermitianFullStateDiagMatrPower(qcomp* out, Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    *out = calcExpecNonHermitianFullStateDiagMatrPower(qureg, matr, exponent);
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

    validate_expecPauliStrValueIsReal(value, qureg.isDensityMatrix, __func__);
    return std::real(value);
}


qreal calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum) {
    validate_quregFields(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);
    validate_pauliStrSumIsHermitian(sum, __func__);

    qcomp value = (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecPauliStrSum(qureg, sum):
        localiser_statevec_calcExpecPauliStrSum(qureg, sum);

    validate_expecPauliStrSumValueIsReal(value, qureg.isDensityMatrix, __func__);
    return std::real(value);
}


qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matrix) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, true, __func__);
    validate_matrixIsHermitian(matrix, __func__);

    // unused parameters; see calcExpecFullStateDiagMatrPower() for explanation
    qcomp exponent = qcomp(1,0);
    bool useRealPow = false;

    qcomp value = (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecFullStateDiagMatr(qureg, matrix, exponent, useRealPow):
        localiser_statevec_calcExpecFullStateDiagMatr(qureg, matrix, exponent, useRealPow);

    // the sub-epsilon imaginary components in matrix never damage the real
    // component of the statevector expectation value, so we do not validate
    // imag(value)~0; we can always safely discard it, and only validate for densmatr:
    if (qureg.isDensityMatrix)
        validate_densMatrExpecDiagMatrValueIsReal(value, exponent, __func__);

    return std::real(value);
}


qreal calcExpecFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qreal exponent) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, true, __func__);
    validate_matrixExpIsHermitian(matrix, exponent, __func__);

    // the backend can use either the pow(qcomp,qcomp) or pow(qreal,qreal) overload;
    // the former is significantly less accurate when the base is real & negative and
    // the exponent is real, because complex pow(a,b) = exp(i b Arg(a)) = exp(i b 2 pi),
    // and the numerical error in pi is compounded by the exponent and the exp(). Because
    // this function assumes approx/intended Hermiticity, it will always use the real-pow
    // overload, discarding the imaginary components of 'matrix' during computation - this
    // is a behaviour unique to this function (other functions collect the erroneous
    // imaginary components before a final validation and discarding).
    const bool useRealPow = true;

    qcomp value = (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecFullStateDiagMatr(qureg, matrix, exponent, useRealPow):
        localiser_statevec_calcExpecFullStateDiagMatr(qureg, matrix, exponent, useRealPow);

    // is it impossible for the statevector routine to produce non-zero
    // imaginary components because of our use of real-pow, the result of
    // which is multiplied by abs(amp). Alas, density matrices multiply the
    // result with a complex scalar and can accrue erroneous imaginary
    // components when the density matrix is non-Hermitian, or due to rounding
    // errors. As such, we only post-validate density matrix values.
    if (qureg.isDensityMatrix)
        validate_densMatrExpecDiagMatrValueIsReal(value, exponent, __func__);

    return std::real(value);
}



/*
 * PROBABILITIES
 */


qreal calcProbOfBasisState(Qureg qureg, qindex index) {
    validate_quregFields(qureg, __func__);
    validate_basisStateIndex(qureg, index, __func__);

    // |i><i| = ||(1+2^N)i>>
    if (qureg.isDensityMatrix)
        index *= 1 + powerOf2(qureg.numQubits);

    qcomp amp = localiser_statevec_getAmp(qureg, index);
    qreal prob = (qureg.isDensityMatrix)? 
        std::real(amp): 
        std::norm(amp);

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
        return std::sqrt(dif);
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


Qureg calcPartialTrace(Qureg qureg, int* traceOutQubits, int numTraceQubits) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, traceOutQubits, numTraceQubits, __func__);
    validate_quregCanBeReduced(qureg, numTraceQubits, __func__);

    // attempt to create reduced Qureg with same deployments as input Qureg
    Qureg out = validateAndCreateCustomQureg(
        qureg.numQubits - numTraceQubits, 
        qureg.isDensityMatrix,  qureg.isDistributed, 
        qureg.isGpuAccelerated, qureg.isMultithreaded, __func__);

    // set it to reduced density matrix
    auto targets = util_getVector(traceOutQubits, numTraceQubits);
    localiser_densmatr_partialTrace(qureg, out, targets);

    return out;
}


Qureg calcReducedDensityMatrix(Qureg qureg, int* retainQubits, int numRetainQubits) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_targets(qureg, retainQubits, numRetainQubits, __func__);
    validate_quregCanBeReduced(qureg, qureg.numQubits - numRetainQubits, __func__);

    auto traceQubits = util_getNonTargetedQubits(retainQubits, numRetainQubits, qureg.numQubits);

    // harmlessly re-validates
    return calcPartialTrace(qureg, traceQubits.data(), traceQubits.size());
}


} // end de-name mangler



/*
 * C++ OVERLOADS
 *
 * which are only ever accessible to C++ binaries, and 
 * accept arguments more natural to C++ (e.g. std::vector).
 */


qreal calcProbOfMultiQubitOutcome(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    // C interface discards individual sizes, so we validate
    validate_measurementOutcomesMatchTargets(qubits.size(), outcomes.size(), __func__);

    // but C function performs all other validation
    return calcProbOfMultiQubitOutcome(qureg, qubits.data(), outcomes.data(), qubits.size());
}


vector<qreal> calcProbsOfAllMultiQubitOutcomes(Qureg qureg, vector<int> qubits) {

    // allocate temp vector, and validate successful (since it's exponentially large!)
    vector<qreal> out;
    qindex numOut = powerOf2(qubits.size());
    auto callback = [&]() { validate_tempAllocSucceeded(false, numOut, sizeof(qreal), __func__); };
    util_tryAllocVector(out, numOut, callback);

    calcProbsOfAllMultiQubitOutcomes(out.data(), qureg, qubits.data(), qubits.size());
    return out;
}

Qureg calcPartialTrace(Qureg qureg, vector<int> traceOutQubits) {
    return calcPartialTrace(qureg, traceOutQubits.data(), traceOutQubits.size());
}


Qureg calcReducedDensityMatrix(Qureg qureg, vector<int> retainQubits) {
    return calcReducedDensityMatrix(qureg, retainQubits.data(), retainQubits.size());
}

