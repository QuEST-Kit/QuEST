/** @file
 * Definition of gpu_qcomp, an extension of base_qcomp and a
 * compatible alternative to the user-facing qcomp, used
 * exclusively by the GPU backend and which is compatible with
 * both CUDA and HIP.
 * 
 * This header is safe to re-include by multiple files because typedef 
 * redefinition is legal in C++, and all functions herein are inline. 
 * Furthermore, since it is only ever parsed by nvcc, the __host__ symbols 
 * are safely processed by other nvcc-only GPU files, like the cuquantum backend.
 * 
 * @author Tyson Jones
 * @author Oliver Brown (patched former HIP arithmetic overloads)
 * @author Erich Essmann (patched former ROCm build issues)
 */

#ifndef GPU_QCOMP_CUH
#define GPU_QCOMP_CUH

#include "quest/include/config.h"
#include "quest/include/types.h"
#include "quest/include/precision.h"

#include "quest/src/core/inliner.hpp"
#include "quest/src/core/base_qcomp.hpp"

#if ! QUEST_COMPILE_CUDA
    #error "A file being compiled somehow included gpu_qcomp.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#if (QUEST_FLOAT_PRECISION == 4)
    #error "Build bug; precision.h should have prevented non-float non-double qcomp precision on GPU."
#endif

#if defined(__HIP__)
    #include "quest/src/gpu/cuda_to_hip.hpp"
#endif

#include <array>



/*
 * DEFINE GPU_QCOMP
 *
 * which is safe to typdef and define additional overloads
 * below, since never witnessed outside the GPU Backend
 */

typedef base_qcomp gpu_qcomp;



/*
 * CONVERTERS
 *
 * which merely wrap the base_qcomp functions for clarity
 * in the GPU source code, disambiguating from cpu_qcomp
 */

INLINE gpu_qcomp* getGpuQcompPtr(qcomp* list) {
    return getBaseQcompPtr(list);
}

INLINE gpu_qcomp getGpuQcomp(qreal re, qreal im) {
    return getBaseQcomp(re, im);
}

// not INLINE to avoid __device__ because qcomp not supported in CUDA kernels
inline gpu_qcomp getGpuQcomp(const qcomp& a) {
    return getBaseQcomp(a.real(), a.imag());
}

// not INLINE to avoid __device__ because qcomp not supported in CUDA kernels
inline qcomp getQcomp(const gpu_qcomp& a) {
    return qcomp( a.re, a.im );
}

template <int Dim>
__host__ inline std::array<gpu_qcomp,Dim> getGpuQcompArray(qcomp matr[Dim]) {
    static_assert(Dim == 2 || Dim == 4);

    // it's crucial we explicitly copy over the elements,
    // rather than just reinterpret the pointer (like we do
    // for heap-memory), because LLVM-based compilers like HIP
    // use aggressive TBAA on stack memory and break the
    // interoperability, causing segfaults here!

    std::array<gpu_qcomp,Dim> out;
    for (int i=0; i<Dim; i++)
        out[i] = getGpuQcomp(matr[i]);

    return out;
}

template <int Dim>
__host__ inline std::array<gpu_qcomp,Dim*Dim> getFlattenedGpuQcompMatrix(qcomp matr[Dim][Dim]) {
    static_assert(Dim == 2 || Dim == 4);

    std::array<gpu_qcomp,Dim*Dim> out;
    for (int i=0; i<Dim*Dim; i++)
        out[i] = getGpuQcomp(matr[i/Dim][i%Dim]);

    return out;
}



/*
 * GPU-SPECIFIC MATHS
 */

INLINE gpu_qcomp pow(gpu_qcomp base, gpu_qcomp exponent) {

    // using https://mathworld.wolfram.com/ComplexExponentiation.html,
    // and the principal argument of 'base'

    // base = a + b i, exponent = c + d i
    qreal a = base.re;
    qreal b = base.im;
    qreal c = exponent.re;
    qreal d = exponent.im;

    // intermediate quantities (uses CUDA atan2,log,pow,exp,cos,sin)
    qreal arg = atan2(b, a);
    qreal mag = a*a + b*b;
    qreal ln = log(mag);
    qreal fac = pow(mag, c/2) * exp(-d * arg);
    qreal ang = c*arg + d*ln/2;

    // output scalar
    qreal re = fac * cos(ang);
    qreal im = fac * sin(ang);
    return getGpuQcomp(re, im);
}



#endif // GPU_QCOMP_CUH