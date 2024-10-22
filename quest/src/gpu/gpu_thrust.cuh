/** @file
 * Subroutines which invoke Thrust. This file is only ever included
 * when COMPILE_CUDA=1 so it can safely invoke CUDA
 * signatures without guards.
 */

#ifndef GPU_THRUST_HPP
#define GPU_THRUST_HPP

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_thrust.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"

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


using devicevec = thrust::device_vector<int>;


int* getPtr(devicevec qubits) {

    return thrust::raw_pointer_cast(qubits.data());
}



/*
 * AMP AND QUREG POINTERS
 *
 * used to enumerate GPU amps inside thrust functions
 */


thrust::device_ptr<cu_qcomp> getStartPtr(cu_qcomp* amps) {

    return thrust::device_pointer_cast(amps);
}


thrust::device_ptr<cu_qcomp> getStartPtr(Qureg qureg) {

    return getStartPtr(toCuQcomps(qureg.gpuAmps));
}


thrust::device_ptr<cu_qcomp> getEndPtr(Qureg qureg) {

    return getStartPtr(qureg) + qureg.numAmpsPerNode;
}



/*
 * CUSTOM FUNCTORS
 *
 * used to effect custom transformations upon GPU
 * amps using thrust functions
 */


struct functor_conj : public thrust::unary_function<cu_qcomp,cu_qcomp> {

    __host__ __device__ cu_qcomp operator()(cu_qcomp amp) {
        amp.y *= -1;
        return amp;
    }
};


struct functor_mixDensityQuregs : public thrust::binary_function<cu_qcomp,cu_qcomp,cu_qcomp> {

    qreal inProb;
    qreal outProb;
    functor_mixDensityQuregs(qreal in, qreal out) : inProb(in), outProb(out) {}

    __host__ __device__ cu_qcomp operator()(cu_qcomp inAmp, cu_qcomp outAmp) {
        return (inProb * inAmp) + (outProb * outAmp);
    }
};


/*
 * FUNCTIONS 
 */


void thrust_setElemsToConjugate(cu_qcomp* matrElemsPtr, qindex matrElemsLen) {

    auto ptr = getStartPtr(matrElemsPtr);
    thrust::transform(ptr, ptr + matrElemsLen, ptr, functor_conj());
}


void thrust_densmatr_mixQureg_subA(qreal outProb, Qureg outQureg, qreal inProb, Qureg inQureg) {

    thrust::transform(
        getStartPtr(inQureg),  getEndPtr(inQureg), 
        getStartPtr(outQureg), getStartPtr(outQureg), // 4th arg is output pointer
        functor_mixDensityQuregs(inProb, outProb));
}




#endif // GPU_THRUST_HPP