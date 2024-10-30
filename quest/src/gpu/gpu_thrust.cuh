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
#include "quest/src/gpu/gpu_types.cuh"

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
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



#endif // GPU_THRUST_HPP