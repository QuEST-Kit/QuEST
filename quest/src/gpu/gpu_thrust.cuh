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

    qindex numRows;
    functor_getDiagInd(Qureg qureg) : numRows(powerOf2(qureg.numQubits)) {}

    __host__ __device__ qindex operator()(qindex i) {

        return i * (numRows + 1);
    }
};


template <int NumBits>
struct functor_insertBits : public thrust::unary_function<qindex,qindex> {

    // this functor inserts bits into a qindex value, and
    // is used to enumerate specific basis-state indices
    // with qubits in the specified bit values 

    vector<int> sortedInds;
    qindex valueMask;
    int numBits;

    functor_insertBits(vector<int> indices, vector<int> bits) {
        assert_numTargsMatchesTemplateParam(indices.size(), NumBits);

        sortedInds = util_getSorted(indices);
        valueMask = util_getBitMask(indices, bits);

        // only consulted when compile-time NumBits==-1
        numBits = indices.size();
    }

    __host__ __device__ qindex operator()(qindex i) {

        // use the compile-time value if possible, to auto-unroll the insertBits loop
        SET_VAR_AT_COMPILE_TIME(int, nbits, NumBits, numBits);

        // return ith local index where bits have the specified values at the specified indices
        return insertBitsWithMaskedValues(i, sortedInds.data(), nbits, valueMask);
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


template <int NumQubits>
qreal thrust_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, bool realOnly) {

    qindex numIters = qureg.numAmpsPerNode / powerOf2(qubits.size());
    auto indFunctor = functor_insertBits<NumQubits>(qubits, outcomes);
    auto valFunctor = (realOnly)? functor_getAmpReal() : functor_getAmpNorm();

    auto rawIter = thrust::make_counting_iterator(0);
    auto indIter = thrust::make_transform_iterator(rawIter, indFunctor);
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), indIter);
    auto probIter= thrust::make_transform_iterator(ampIter, valFunctor);

    qreal prob = thrust::reduce(probIter, probIter + numIters, getQcomp(0,0));
    return prob;
}


qreal thrust_statevec_calcTotalProb_sub(Qureg qureg) {

    qreal prob = thrust::transform_reduce(
        getStartPtr(qureg), getEndPtr(qureg), 
        functor_getAmpNorm(), getQcomp(0,0));
}


qreal thrust_densmatr_calcTotalProb_sub(Qureg qureg) {

    qindex numColsPerNode = powerOf2(qureg.logNumColsPerNode);
    qindex startInd = qureg.rank * numColsPerNode;

    auto rawIter = thrust::make_counting_iterator(startInd);
    auto indIter = thrust::make_transform_iterator(rawIter, functor_getDiagInd(qureg));
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), indIter);
    auto probIter= thrust::make_transform_iterator(ampIter, functor_getAmpReal());

    qreal prob = thrust::reduce(probIter, probIter + numColsPerNode, getQcomp(0,0));
    return prob;
}



#endif // GPU_THRUST_HPP