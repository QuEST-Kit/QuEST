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



/*
 * INTERNAL FUNCTIONS
 */


extern Qureg validateAndCreateCustomQureg(
    int numQubits, int isDensMatr, int useDistrib, 
    int useGpuAccel, int useMultithread, const char* caller);


vector<int> getNonTargetedQubits(int* targets, int numTargets, int numQubits) {
    
    qindex mask = getBitMask(targets, numTargets);

    vector<int> nonTargets;
    nonTargets.reserve(numQubits - numTargets);

    for (int i=0; i<numQubits; i++)
        if (getBit(mask, i) == 0)
            nonTargets.push_back(i);

    return nonTargets;
}



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

    return (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecFullStateDiagMatr(qureg, matrix, exponent):
        localiser_statevec_calcExpecFullStateDiagMatr(qureg, matrix, exponent);
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

    return calcExpecFullStateDiagMatrPower(qureg, matrix, 1); // harmlessly re-validates
}


qreal calcExpecFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent) {
    validate_quregFields(qureg, __func__);
    validate_matrixFields(matrix, __func__);
    validate_matrixAndQuregAreCompatible(matrix, qureg, true, __func__);
    validate_matrixIsHermitian(matrix, __func__);

    qcomp value = (qureg.isDensityMatrix)?
        localiser_densmatr_calcExpecFullStateDiagMatr(qureg, matrix, exponent):
        localiser_statevec_calcExpecFullStateDiagMatr(qureg, matrix, exponent);

    // demand value is real, despite exponent being complex
    validate_expecFullStateDiagMatrValueIsReal(value, qureg.isDensityMatrix, __func__);
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
        std::real(amp) : 
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

    auto traceQubits = getNonTargetedQubits(retainQubits, numRetainQubits, qureg.numQubits);

    // harmlessly re-validates
    return calcPartialTrace(qureg, traceQubits.data(), traceQubits.size());
}


void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits) {
    validate_quregFields(in, __func__);
    validate_quregFields(out, __func__);
    validate_quregIsDensityMatrix(in, __func__);
    validate_quregIsDensityMatrix(out, __func__);
    validate_targets(in, traceOutQubits, numTraceQubits, __func__);
    validate_quregCanBeSetToReducedDensMatr(out, in, numTraceQubits, __func__);

    auto targets = util_getVector(traceOutQubits, numTraceQubits);
    localiser_densmatr_partialTrace(in, out, targets);
}


void setQuregToReducedDensityMatrix(Qureg out, Qureg in, int* retainQubits, int numRetainQubits) {
    validate_quregFields(in, __func__);
    validate_quregFields(out, __func__);
    validate_quregIsDensityMatrix(in, __func__);
    validate_quregIsDensityMatrix(out, __func__);
    validate_targets(in, retainQubits, numRetainQubits, __func__);
    validate_quregCanBeSetToReducedDensMatr(out, in, in.numQubits - numRetainQubits, __func__);

    auto traceQubits = getNonTargetedQubits(retainQubits, numRetainQubits, in.numQubits);
    localiser_densmatr_partialTrace(in, out, traceQubits);
}



} // end de-name mangler
