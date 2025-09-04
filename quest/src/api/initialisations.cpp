/** @file
 * API definitions for initialising Quregs into 
 * particular states. Note when a Qureg is GPU-
 * accelerated, these functions only update the
 * state in GPU memory; the CPU amps are unchanged.
 * 
 * @author Tyson Jones
 */

#include "quest/include/qureg.h"
#include "quest/include/calculations.h"
#include "quest/include/initialisations.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/localiser.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

using std::vector;



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
extern "C" {


/*
 * INIT
 */


void initBlankState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // |null> = {0, 0, 0, ...}
    qcomp amp = 0;
    localiser_statevec_initUniformState(qureg, amp);
}


void initZeroState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // |0> = |0><0|
    qindex ind = 0;
    localiser_statevec_initClassicalState(qureg, ind);
}


void initPlusState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // |+>    = sum_i 1/sqrt(2^N) |i>  where 2^N = numAmps
    // |+><+| = sum_ij 1/2^N |i><j|    where 2^N = sqrt(numAmps)
    qcomp amp = 1.0 / std::sqrt(qureg.numAmps);
    localiser_statevec_initUniformState(qureg, amp);
}


void initPureState(Qureg qureg, Qureg pure) {
    validate_quregFields(qureg, __func__);
    validate_quregFields(pure, __func__);
    validate_quregCanBeInitialisedToPureState(qureg, pure, __func__);

    (qureg.isDensityMatrix)?
        localiser_densmatr_initPureState(qureg, pure):
        localiser_statevec_setQuregToClone(qureg, pure);
}


void initClassicalState(Qureg qureg, qindex ind) {
    validate_quregFields(qureg, __func__);
    validate_basisStateIndex(qureg, ind, __func__);

    // |ind><ind| = ||ind'>>
    if (qureg.isDensityMatrix)
        ind = util_getGlobalFlatIndex(qureg, ind, ind);

    localiser_statevec_initClassicalState(qureg, ind);
}


void initDebugState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    localiser_statevec_initDebugState(qureg);
}


void initArbitraryPureState(Qureg qureg, qcomp* amps) {
    validate_quregFields(qureg, __func__);

    // set qureg = |amps> or |amps><amps|
    (qureg.isDensityMatrix)?
        localiser_densmatr_initArbitraryPureState(qureg, amps):
        localiser_statevec_initArbitraryPureState(qureg, amps);
}


void initArbitraryMixedState(Qureg qureg, qcomp** amps) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);

    localiser_densmatr_initArbitraryMixedState(qureg, amps);
}


void initRandomPureState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // these invoked localiser functions may harmlessly 
    // re-call the API and re-perform input validation

    if (qureg.isDensityMatrix)
        localiser_densmatr_initUniformlyRandomPureStateAmps(qureg); 
    else {
        localiser_statevec_initUnnormalisedUniformlyRandomPureStateAmps(qureg);
        setQuregToRenormalized(qureg);
    }
}


void initRandomMixedState(Qureg qureg, qindex numPureStates) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_numInitRandomPureStates(numPureStates, __func__);

    localiser_densmatr_initMixtureOfUniformlyRandomPureStates(qureg, numPureStates);
}



/*
 * SET
 */


void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps) {
    validate_quregFields(qureg, __func__);
    validate_quregIsStateVector(qureg, __func__);
    validate_basisStateIndices(qureg, startInd, numAmps, __func__);

    localiser_statevec_setAmps(amps, qureg, startInd, numAmps);
}


void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_basisStateRowCols(qureg, startRow, startCol, numRows, numCols, __func__);

    localiser_densmatr_setAmps(amps, qureg, startRow, startCol, numRows, numCols);
}


void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_basisStateIndices(qureg, startInd, numAmps, __func__); // validation msg correct for density-matrix

    localiser_statevec_setAmps(amps, qureg, startInd, numAmps);
}


void setQuregToClone(Qureg outQureg, Qureg inQureg) {
    validate_quregFields(outQureg, __func__);
    validate_quregFields(inQureg, __func__);
    validate_quregsCanBeCloned(outQureg, inQureg, __func__);

    // we invoke mixing/superposing, which involves superfluous
    // floating-point operators but is not expected to cause an
    // appreciable slowdown since simulation is often memory-bound
    (outQureg.isDensityMatrix)?
        localiser_densmatr_mixQureg(0, outQureg, 1, inQureg):
        localiser_statevec_setQuregToClone(outQureg, inQureg);
}


qreal setQuregToRenormalized(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    qreal prob = calcTotalProb(qureg); // harmlessly re-validates
    validate_quregRenormProbIsNotZero(prob, __func__);

    qreal norm = (qureg.isDensityMatrix)? prob : std::sqrt(prob);
    qreal fac = 1 / norm;
    localiser_statevec_scaleAmps(qureg, fac);

    return fac;
}


void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum) {
    validate_quregFields(qureg, __func__);
    validate_quregIsDensityMatrix(qureg, __func__);
    validate_pauliStrSumFields(sum, __func__);
    validate_pauliStrSumTargets(sum, qureg, __func__);

    // sum is permitted to be non-Hermitian, since Hermiticity
    // is insufficient to ensure qureg would be physical/valid

    localiser_densmatr_setAmpsToPauliStrSum(qureg, sum);
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

    auto traceQubits = util_getNonTargetedQubits(retainQubits, numRetainQubits, in.numQubits);
    localiser_densmatr_partialTrace(in, out, traceQubits);
}


void setQuregToWeightedSum(Qureg out, qcomp* coeffs, Qureg* in, int numIn) {
    validate_quregFields(out, __func__);
    validate_numQuregsInSum(numIn, __func__);
    validate_quregsCanBeSummed(out, in, numIn, __func__); // also validates all init

    auto coeffVec = util_getVector(coeffs, numIn);
    auto inVec = util_getVector(in, numIn);
    localiser_statevec_setQuregToWeightedSum(out, coeffVec, inVec);
}


void setQuregToMixture(Qureg out, qreal* probs, Qureg* in, int numIn) {
    validate_quregFields(out, __func__);
    validate_quregIsDensityMatrix(out, __func__);
    validate_numQuregsInSum(numIn, __func__);
    validate_quregsCanBeMixed(out, in, numIn, __func__); // also validates all init & densmatr
    validate_probabilities(probs, numIn, __func__);

    // convert probs to complex (assume this alloc never fails)
    vector<qcomp> coeffVec(numIn);
    for (int i=0; i<numIn; i++)
        coeffVec[i] = getQcomp(probs[i], 0);

    auto inVec = util_getVector(in, numIn);
    localiser_statevec_setQuregToWeightedSum(out, coeffVec, inVec);
}


} // end de-mangler



/*
 * C++ OVERLOADS
 */


void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, vector<vector<qcomp>> amps) {

    // C++-specific validation
    validate_matrixRowsAllSameSize(amps, __func__);

    // we must pass nested pointers to the C function, so alloc a vector
    // of pointers of amps. We defensively check the temp vector allocates fine
    vector<qcomp*> ptrs;
    size_t len = amps.size();
    auto callback = [&]() { validate_tempAllocSucceeded(false, len, sizeof(qcomp*), __func__); };
    util_tryAllocVector(ptrs, len, callback);

    // then set the pointers
    for (size_t i=0; i<len; i++)
        ptrs[i] = amps[i].data();

    // C function performs main validation
    setDensityQuregAmps(qureg, startRow, startCol, ptrs.data(), len, (len>0)? len : 0); // avoid seg-fault
}

void setQuregAmps(Qureg qureg, qindex startInd, vector<qcomp> amps) {
    setQuregAmps(qureg, startInd, amps.data(), amps.size());
}

void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, vector<qcomp> amps) {
    setDensityQuregFlatAmps(qureg, startInd, amps.data(), amps.size());
}

void setQuregToPartialTrace(Qureg out, Qureg in, vector<int> traceOutQubits) {
    setQuregToPartialTrace(out, in, traceOutQubits.data(), traceOutQubits.size());
}

void setQuregToReducedDensityMatrix(Qureg out, Qureg in, vector<int> retainQubits) {
    setQuregToReducedDensityMatrix(out, in, retainQubits.data(), retainQubits.size());
}

void setQuregToWeightedSum(Qureg out, vector<qcomp> coeffs, vector<Qureg> in) {
    validate_numQuregsMatchesCoeffs(in.size(), coeffs.size(), __func__);

    setQuregToWeightedSum(out, coeffs.data(), in.data(), in.size());
}

void setQuregToMixture(Qureg out, vector<qreal> probs, vector<Qureg> in) {
    validate_numQuregsMatchesProbs(in.size(), probs.size(), __func__);

    setQuregToMixture(out, probs.data(), in.data(), in.size());
}
