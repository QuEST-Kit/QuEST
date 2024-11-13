/** @file
 * API definitions for initialisaing Qureg states.
 *
 * Note many of these can be replaced with bespoke,
 * faster OpenMP and CUDA implementations although
 * that this may prove a premature optimisation
 * since initialisations are expectedly called 
 * relatively very few times in simulation. 
 * Think carefully!
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

    // TODO:
    // this invokes a slow, non-parallel overwrite;
    // bespoke GPU functions (and maybe even OpenMP)
    // are likely much faster

    // fill CPU memory with 0
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(0,0));

    // overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
}


void initZeroState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // TODO:
    // this invokes a slow, non-parallel overwrite;
    // bespoke GPU functions (and maybe even OpenMP)
    // are likely much faster

    // fill CPU memory with 0
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(0,0));

    // set the first global element to 1 (valid for both statevecs and density matrices)
    if (qureg.rank == 0)
        qureg.cpuAmps[0] = 0;

    // overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
}


void initPlusState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    // TODO:
    // this invokes a slow, non-parallel overwrite;
    // bespoke GPU functions (and maybe even OpenMP)
    // are likely much faster

    // |+>    = sum_i 1/sqrt(2^N) |i>  where 2^N = numAmps
    // |+><+| = sum_ij 1/2^N |i><j|    where 2^N = sqrt(numAmps)
    qreal elem = 1.0 / sqrt(qureg.numAmps);
        
    // fill CPU memory with elem
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(elem,0));

    // overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
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

    // TODO:
    // this invokes a slow, non-parallel overwrite;
    // bespoke GPU functions (and maybe even OpenMP)
    // are likely much faster

    // fill CPU memory with 0
    std::fill(qureg.cpuAmps, qureg.cpuAmps + qureg.numAmpsPerNode, qcomp(0,0));

    // |ind><ind| = ||ind'>>
    if (qureg.isDensityMatrix)
        ind = util_getGlobalFlatIndex(qureg, ind, ind);

    // only one node needs to write a 1 elem
    if (qureg.rank == util_getRankContainingIndex(qureg, ind))
        qureg.cpuAmps[ind % qureg.numAmpsPerNode] = 1;

    // all nodes overwrite GPU memory
    if (qureg.isGpuAccelerated)
        gpu_copyCpuToGpu(qureg);
}


void initDebugState(Qureg qureg) {

    // TODO
    error_functionNotImplemented(__func__);
}


void initArbitraryState(Qureg qureg, qcomp* amps) {

    // TODO
    error_functionNotImplemented(__func__); 
}


void initRandomPureState(Qureg qureg) {
    validate_quregFields(qureg, __func__);

    if (qureg.isDensityMatrix)
        localiser_densmatr_setUniformlyRandomPureStateAmps_sub(qureg); // harmlessly recalls API and re-validates
    else {
        localiser_statevec_setUnnormalisedUniformlyRandomPureStateAmps_sub(qureg);
        setQuregToRenormalized(qureg); // harmlessly re-validates
    }
}


void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps) {

    // TODO
    error_functionNotImplemented(__func__);

    // re-use code/logic for FullStateDiagMatr (obviously)
}

void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols) {

    // TODO
    error_functionNotImplemented(__func__);
}


void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum) {

    // TODO
    error_functionNotImplemented(__func__);
}


void setQuregToClone(Qureg targetQureg, Qureg copyQureg) {
    validate_quregFields(targetQureg, __func__);
    validate_quregFields(copyQureg, __func__);
    validate_quregsCanBeCloned(targetQureg, copyQureg, __func__);

    (targetQureg.isDensityMatrix)?
        localiser_densmatr_mixQureg(0, targetQureg, 1, copyQureg):
        localiser_statevec_setQuregToSuperposition(0, targetQureg, 1, copyQureg, 0, copyQureg);
}


void setQuregToSuperposition(qcomp facOut, Qureg out, qcomp fac1, Qureg qureg1, qcomp fac2, Qureg qureg2) {
    validate_quregFields(out, __func__);
    validate_quregFields(qureg1, __func__);
    validate_quregFields(qureg2, __func__);
    validate_quregsCanBeSuperposed(out, qureg1, qureg2, __func__);

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



} // end de-mangler