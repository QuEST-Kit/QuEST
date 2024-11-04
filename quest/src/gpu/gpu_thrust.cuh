/** @file
 * Subroutines which invoke Thrust. This file is only ever included
 * when COMPILE_CUDA=1 so it can safely invoke CUDA signatures without 
 * guards. Further, as it is entirely a header, it can declare templated
 * times without explicitly instantiating them across all parameter values.
 */

#ifndef GPU_THRUST_HPP
#define GPU_THRUST_HPP

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_thrust.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/miscellaneous.hpp"
#include "quest/src/gpu/gpu_types.cuh"

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>



/*
 * QUBIT LISTS
 *
 * are copied to device memory using thrust's device_vector's 
 * copy constructor (devicevec d_vec = hostvec). The pointer 
 * to the data (d_vec.data()) can be cast into a raw pointer
 * and passed directly to CUDA kernels
 */


using devints = thrust::device_vector<int>;

int* getPtr(devints qubits) {

    return thrust::raw_pointer_cast(qubits.data());
}


using devreals = thrust::device_vector<qreal>;

qreal* getPtr(devreals reals) {

    return thrust::raw_pointer_cast(reals.data());
}

void copyFromDeviceVec(devreals reals, qreal* out) {

    thrust::copy(reals.begin(), reals.end(), out);
}



/*
 * AMP POINTERS
 *
 * used to enumerate GPU amps of matrices, quregs and
 * full-state diagonal matrices inside thrust functions
 */


thrust::device_ptr<cu_qcomp> getStartPtr(cu_qcomp* amps) {

    return thrust::device_pointer_cast(amps);
}
auto getStartPtr(qcomp* amps) {

    return getStartPtr(toCuQcomps(amps));
}


auto getStartPtr(Qureg qureg) {

    return getStartPtr(qureg.gpuAmps);
}
auto getEndPtr(Qureg qureg) {

    return getStartPtr(qureg) + qureg.numAmpsPerNode;
}


auto getStartPtr(FullStateDiagMatr matr) {

    return getStartPtr(matr.gpuElems);
}
auto getEndPtr(FullStateDiagMatr matr) {

    return getStartPtr(matr) + matr.numElemsPerNode;
}



/*
 * CUSTOM FUNCTORS
 *
 * used to effect custom transformations upon GPU
 * amps using thrust functions
 */


struct functor_getAmpConj : public thrust::unary_function<cu_qcomp,cu_qcomp> {

    __host__ __device__ cu_qcomp operator()(cu_qcomp amp) {
        return getCompConj(amp);
    }
};

struct functor_getAmpNorm : public thrust::unary_function<cu_qcomp,qreal> {

    __host__ __device__ qreal operator()(cu_qcomp amp) {
        return getCompNorm(amp);
    }
};

struct functor_getAmpReal : public thrust::unary_function<cu_qcomp,qreal> {

    __host__ __device__ qreal operator()(cu_qcomp amp) {
        return getCompReal(amp);
    }
};


struct functor_mixAmps : public thrust::binary_function<cu_qcomp,cu_qcomp,cu_qcomp> {

    // this functor linearly combines the given pair
    // of amplitudes, weighted by the fixed qreals,
    // and is used by mixQureg upon density matrices

    qreal outProb;
    qreal inProb;
    functor_mixAmps(qreal out, qreal in) : outProb(out), inProb(in) {}

    __host__ __device__ cu_qcomp operator()(cu_qcomp outAmp, cu_qcomp inAmp) {
        
        return (outProb * outAmp) + (inProb * inAmp);
    }
};


struct functor_superposeAmps {

    // this functor linearly combines the given trio
    // of amplitudes, weighted by fixed qcomps, and is
    // used by setQuregToSuperposition

    cu_qcomp fac0;
    cu_qcomp fac1;
    cu_qcomp fac2;
    functor_superposeAmps(cu_qcomp f0, cu_qcomp f1, cu_qcomp f2) : fac0(f0), fac1(f1), fac2(f2) {}

    template <typename Tuple> __host__ __device__ void operator()(Tuple t) {
        thrust::get<0>(t) = fac0*thrust::get<0>(t) + fac1*thrust::get<1>(t) + fac2*thrust::get<2>(t);
    }
};


template <bool HasPower>
struct functor_multiplyElemPowerWithAmp : public thrust::binary_function<cu_qcomp,cu_qcomp,cu_qcomp> {

    // this functor multiplies a diagonal matrix element 
    // raised to a power (templated to optimise away the 
    // exponentiation at compile-time when power==1) upon
    // a statevector amp, used to modify the statevector

    cu_qcomp exponent;
    functor_multiplyElemPowerWithAmp(cu_qcomp power) : exponent(power) {}

    __host__ __device__ cu_qcomp operator()(cu_qcomp quregAmp, cu_qcomp matrElem) {

        if constexpr (HasPower)
            matrElem = getCompPower(matrElem, exponent);

        return quregAmp * matrElem;
    }
};


struct functor_getDiagInd : public thrust::unary_function<qindex,qindex> {

    // this functor accepts the index of a statevector 
    // basis-state and produces the index of a density 
    // matrix's corresponding diagonal basis-state

    qindex numAmpsPerCol;
    qindex firstDiagInd;

    functor_getDiagInd(Qureg qureg) {
        firstDiagInd = misc_getLocalIndexOfFirstDiagonalAmp(qureg);
        numAmpsPerCol = powerOf2(qureg.numQubits);
    }

    __host__ __device__ qindex operator()(qindex i) {

        return misc_getLocalIndexOfDiagonalAmp(i, firstDiagInd, numAmpsPerCol);
    }
};


template <int NumBits>
struct functor_insertBits : public thrust::unary_function<qindex,qindex> {

    // this functor inserts bits into a qindex value, and
    // is used to enumerate specific basis-state indices
    // with qubits in the specified bit values

    int* sortedIndsPtr;
    qindex valueMask;
    int numBits;

    functor_insertBits(int* ptr, qindex mask, int nBits) :
        sortedIndsPtr(ptr), valueMask(mask), numBits(nBits)
    {
        assert_numTargsMatchesTemplateParam(nBits, NumBits);
    }

    __host__ __device__ qindex operator()(qindex i) {

        // use the compile-time value if possible, to auto-unroll the insertBits loop
        SET_VAR_AT_COMPILE_TIME(int, nbits, NumBits, numBits);

        // return ith local index where bits have the specified values at the specified indices
        return insertBitsWithMaskedValues(i, sortedIndsPtr, nbits, valueMask);
    }
};


template <int NumTargets>
struct functor_projectStateVec : public thrust::binary_function<qindex,cu_qcomp,cu_qcomp> {

    // this functor multiplies an amp with zero or a 
    // renormalisation codfficient, depending on whether
    // the basis state of the amp has qubits in a particular
    // configuration. This is used to project statevector
    // qubits into a particular measurement outcome

    int* targetsPtr;
    int numTargets;
    int rank;
    qindex logNumAmpsPerNode;
    qindex retainValue;
    qreal renorm;

    functor_projectStateVec(
        int* targetsPtr, int numTargets, int rank, 
        qindex logNumAmpsPerNode, qindex retainValue, qreal renorm
    ) :
        targetsPtr(targetsPtr), numTargets(numTargets), rank(rank), 
        logNumAmpsPerNode(logNumAmpsPerNode), retainValue(retainValue), renorm(renorm)
    { 
        assert_numTargsMatchesTemplateParam(numTargets, NumTargets);
    }

    __host__ __device__ cu_qcomp operator()(qindex n, cu_qcomp amp) {

        // use the compile-time value if possible, to auto-unroll the getValueOfBits() loop below
        SET_VAR_AT_COMPILE_TIME(int, numBits, NumTargets, numTargets);

        // i = global index of nth local amp
        qindex i = concatenateBits(rank, n, logNumAmpsPerNode);

        // return amp scaled by zero or renorm, depending on whether n has projected substate
        qindex val = getValueOfBits(i, targetsPtr, numBits);
        qreal fac = renorm * (val == retainValue);
        return fac * amp;
    }
};


template <int NumTargets>
struct functor_projectDensMatr : public thrust::binary_function<qindex,cu_qcomp,cu_qcomp> {

    // this functor multiplies an amp with zero or a 
    // renormalisation codfficient, depending on whether
    // the basis state of the amp has qubits in a particular
    // configuration. This is used to project density matrix
    // qubits into a particular measurement outcome

    int* targetsPtr;
    int numTargets;
    int rank;
    int numQuregQubits;
    qindex logNumAmpsPerNode;
    qindex retainValue;
    qreal renorm;

    functor_projectDensMatr(
        int* targetsPtr, int numTargets, int rank, int numQuregQubits,
        qindex logNumAmpsPerNode, qindex retainValue, qreal renorm
    ) :
        targetsPtr(targetsPtr), numTargets(numTargets), rank(rank), numQuregQubits(numQuregQubits),
        logNumAmpsPerNode(logNumAmpsPerNode), retainValue(retainValue), renorm(renorm)
    { 
        assert_numTargsMatchesTemplateParam(numTargets, NumTargets);
    }

    __host__ __device__ cu_qcomp operator()(qindex n, cu_qcomp amp) {

        // use the compile-time value if possible, to auto-unroll the getValueOfBits() loop below
        SET_VAR_AT_COMPILE_TIME(int, numBits, NumTargets, numTargets);

        // i = global index of nth local amp
        qindex i = concatenateBits(rank, n, logNumAmpsPerNode);

        // r, c = global row and column indices of nth local amp
        qindex r = getBitsRightOfIndex(i, numQuregQubits);
        qindex c = getBitsLeftOfIndex(i, numQuregQubits-1);

        qindex v1 = getValueOfBits(r, targetsPtr, numBits);
        qindex v2 = getValueOfBits(c, targetsPtr, numBits);

        // multiply amp with renorm or zero if values disagree with given outcomes
        qreal fac = renorm * (v1 == v2) * (retainValue == v1);
        return fac * amp;
    }
};



/*
 * STATE MODIFICATION 
 */


void thrust_setElemsToConjugate(cu_qcomp* matrElemsPtr, qindex matrElemsLen) {

    auto ptr = getStartPtr(matrElemsPtr);
    thrust::transform(ptr, ptr + matrElemsLen, ptr, functor_getAmpConj());
}


void thrust_densmatr_mixQureg_subA(qreal outProb, Qureg outQureg, qreal inProb, Qureg inQureg) {

    thrust::transform(
        getStartPtr(outQureg), getEndPtr(outQureg), 
        getStartPtr(inQureg),  getStartPtr(outQureg), // 4th arg is output pointer
        functor_mixAmps(outProb, inProb));
}


void thrust_statevec_setQuregToSuperposition_sub(cu_qcomp facOut, Qureg outQureg, cu_qcomp fac1, Qureg inQureg1, cu_qcomp fac2, Qureg inQureg2) {

    thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(getStartPtr(outQureg), getStartPtr(inQureg1), getStartPtr(inQureg2))),
        thrust::make_zip_iterator(thrust::make_tuple(getEndPtr(outQureg),   getEndPtr(inQureg1),   getEndPtr(inQureg2))),
        functor_superposeAmps(facOut, fac1, fac2));
}


template <bool HasPower>
void thrust_statevec_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, cu_qcomp exponent) {

    thrust::transform(
        getStartPtr(qureg), getEndPtr(qureg), 
        getStartPtr(matr),  getStartPtr(qureg), // 4th arg is output pointer
        functor_multiplyElemPowerWithAmp<HasPower>(exponent));
}



/*
 * PROBABILITIES
 */


qreal thrust_statevec_calcTotalProb_sub(Qureg qureg) {

    qreal prob = thrust::transform_reduce(
        getStartPtr(qureg), getEndPtr(qureg), 
        functor_getAmpNorm(), 0, thrust::plus<qreal>());

    return prob;
}


qreal thrust_densmatr_calcTotalProb_sub(Qureg qureg) {

    auto rawIter = thrust::make_counting_iterator(0);
    auto indIter = thrust::make_transform_iterator(rawIter, functor_getDiagInd(qureg));
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), indIter);
    auto probIter= thrust::make_transform_iterator(ampIter, functor_getAmpReal());

    qindex numIts = powerOf2(qureg.logNumColsPerNode);
    return thrust::reduce(probIter, probIter + numIts);
}


template <int NumQubits>
qreal thrust_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    devints sortedQubits = util_getSorted(qubits);
    qindex valueMask = util_getBitMask(qubits, outcomes);

    auto indFunctor = functor_insertBits<NumQubits>(getPtr(sortedQubits), valueMask, qubits.size());
    auto probFunctor = functor_getAmpNorm();

    auto rawIter = thrust::make_counting_iterator(0);
    auto indIter = thrust::make_transform_iterator(rawIter, indFunctor);
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), indIter);
    auto probIter = thrust::make_transform_iterator(ampIter, probFunctor);

    qindex numIts = qureg.numAmpsPerNode / powerOf2(qubits.size());
    return thrust::reduce(probIter, probIter + numIts);
}


template <int NumQubits>
qreal thrust_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    // cannot move these into functor_insertBits constructor, since the memory
    // would dangle - and we cannot bind deviceints as an attribute - it's host-only!
    devints sortedQubits = util_getSorted(qubits);
    qindex valueMask = util_getBitMask(qubits, outcomes);

    auto basisIndFunctor = functor_insertBits<NumQubits>(getPtr(sortedQubits), valueMask, qubits.size());
    auto diagIndFunctor = functor_getDiagInd(qureg);
    auto probFunctor = functor_getAmpReal();

    auto rawIter = thrust::make_counting_iterator(0);
    auto indIter = thrust::make_transform_iterator(rawIter, basisIndFunctor);
    auto diagIter = thrust::make_transform_iterator(indIter, diagIndFunctor);
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), diagIter);
    auto probIter = thrust::make_transform_iterator(ampIter, probFunctor);

    qindex numIts = misc_getNumLocalDiagonalsWithBits(qureg, qubits, outcomes);
    return thrust::reduce(probIter, probIter + numIts);
}



/*
 * PROJECTORS
 */


template <int NumQubits>
void thrust_statevec_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal norm) {

    devints devQubits = qubits;
    qindex retainValue = getIntegerFromBits(outcomes.data(), outcomes.size());
    auto projFunctor = functor_projectStateVec<NumQubits>(
        devQubits.data(), qubits.size(), qureg.rank, 
        qureg.logNumAmpsPerNode, retainValue, norm);

    auto indIter = thrust::make_counting_iterator(0);
    auto ampIter = getStartPtr(qureg);

    qindex numIts = qureg.numAmpsPerNode;
    thrust::transform(indIter, indIter + numIts, ampIter, ampIter, projFunctor); // 4th arg gets modified
}


template <int NumQubits>
void thrust_densmatr_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal norm) {

    devints devQubits = qubits;
    qindex retainValue = getIntegerFromBits(outcomes.data(), outcomes.size());
    auto projFunctor = functor_projectDensMatr<NumQubits>(
        devQubits.data(), qubits.size(), qureg.rank, qureg.numQubits,
        qureg.logNumAmpsPerNode, retainValue, norm);

    auto indIter = thrust::make_counting_iterator(0);
    auto ampIter = getStartPtr(qureg);

    qindex numIts = qureg.numAmpsPerNode;
    thrust::transform(indIter, indIter + numIts, ampIter, ampIter, projFunctor); // 4th arg gets modified
}



#endif // GPU_THRUST_HPP