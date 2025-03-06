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

#include "quest/src/core/errors.hpp" // only needed for not-implemented functions

#include <algorithm>
#include <cmath>


// enable invocation by both C and C++ binaries
extern "C" {


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
    qcomp amp = 1.0 / sqrt(qureg.numAmps);
    localiser_statevec_initUniformState(qureg, amp);
}


void initPureState(Qureg qureg, Qureg pure) {
    validate_quregFields(qureg, __func__);
    validate_quregFields(pure, __func__);
    validate_quregCanBeInitialisedToPureState(qureg, pure, __func__);

    // when qureg=statevec, we lazily invoke setQuregToSuperposition which
    // invokes superfluous floating-point operations which will be happily
    // occluded by the memory movement costs
    (qureg.isDensityMatrix)?
        localiser_densmatr_initPureState(qureg, pure):
        localiser_statevec_setQuregToSuperposition(0, qureg, 1, pure, 0, pure);
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


void setQuregToClone(Qureg targetQureg, Qureg copyQureg) {
    validate_quregFields(targetQureg, __func__);
    validate_quregFields(copyQureg, __func__);
    validate_quregsCanBeCloned(targetQureg, copyQureg, __func__);

    // we invoke mixing/superposing, which involves superfluous
    // floating-point operators but is not expected to cause an
    // appreciable slowdown since simulation is memory-bound
    (targetQureg.isDensityMatrix)?
        localiser_densmatr_mixQureg(0, targetQureg, 1, copyQureg):
        localiser_statevec_setQuregToSuperposition(0, targetQureg, 1, copyQureg, 0, copyQureg);
}


void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2) {
    validate_quregFields(out, __func__);
    validate_quregFields(qureg1, __func__);
    validate_quregFields(qureg2, __func__);
    validate_quregsCanBeSuperposed(out, qureg1, qureg2, __func__); // asserts statevectors

    localiser_statevec_setQuregToSuperposition(facOut, out, fac1, qureg1, fac2, qureg2);
}


qreal setQuregToRenormalized(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    qreal prob = calcTotalProb(qureg); // harmlessly re-validates
    validate_quregRenormProbIsNotZero(prob, __func__);

    qreal norm = (qureg.isDensityMatrix)? prob : sqrt(prob);
    qreal fac = 1 / norm;
    localiser_statevec_setQuregToSuperposition(fac, qureg, 0, qureg, 0, qureg);

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



} // end de-mangler