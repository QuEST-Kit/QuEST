/** @file
 * CUDA-compatible complex types. This file is only ever included
 * when COMPILE_CUDA=1 so it can safely invoke CUDA
 * signatures without guards. This is safe to re-include by
 * multiple files because typedef redefinition is legal in C++,
 * and all functions herein are inline. Furthermore, since it
 * is only ever parsed by nvcc, the __host__ symbols are safely
 * processed by the cuquantum backend.
 */

#ifndef GPU_TYPES_HPP
#define GPU_TYPES_HPP

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/inliner.hpp"

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_types.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#include <array>
#include <cuComplex.h>
#include <thrust/device_vector.h>



/*
 * COPYING VECTORS FROM HOST TO DEVICE
 *
 * is done using thrust's device_vector's copy constructor
 * (devicevec d_vec = hostvec), the pointer of which (.data())
 * can be passed to, and accessed by, a CUDA kernel
 */


using devicevec = thrust::device_vector<int>;



/*
 * CUDA-COMPATIBLE QCOMP ALIAS (cu_qcomp)
 *
 * which we opt to use over a Thrust complex type to gaurantee
 * compatibility with cuQuantum, though this irritatingly 
 * requires explicitly defining operator overloads below
 */


#if (FLOAT_PRECISION == 1)
    typedef cuFloatComplex cu_qcomp;

#elif (FLOAT_PRECISION == 2)
    typedef cuDoubleComplex cu_qcomp;

#else
    #error "Build bug; precision.h should have prevented non-float non-double qcomp precision on GPU."

#endif



/*
 * cu_qcomp ARITHMETIC OVERLOADS
 */


INLINE cu_qcomp operator + (const cu_qcomp& a, const cu_qcomp& b) {
    return (cu_qcomp) {
        .x = a.x + b.x,
        .y = a.y + b.y
    };
}

INLINE cu_qcomp operator - (const cu_qcomp& a, const cu_qcomp& b) {
    return (cu_qcomp) {
        .x = a.x - b.x,
        .y = a.y - b.y
    };
}

INLINE cu_qcomp operator * (const cu_qcomp& a, const cu_qcomp& b) {
    return (cu_qcomp) {
        .x = a.x * b.x - a.y * b.y,
        .y = a.x * b.y + a.y * b.x
    };
}


INLINE cu_qcomp operator + (const cu_qcomp& a, const qreal& b) {
    return (cu_qcomp) {
        .x = a.x + b,
        .y = a.y + b
    };
}

INLINE cu_qcomp operator - (const cu_qcomp& a, const qreal& b) {
    return (cu_qcomp) {
        .x = a.x - b,
        .y = a.y - b
    };
}

INLINE cu_qcomp operator * (const cu_qcomp& a, const qreal& b) {
    return (cu_qcomp) {
        .x = a.x * b,
        .y = a.y * b
    };
}


INLINE void operator += (cu_qcomp& a, const cu_qcomp& b) {
    a = a + b;
}

INLINE void operator -= (cu_qcomp& a, const cu_qcomp& b) {
    a = a - b;
}

INLINE void operator *= (cu_qcomp& a, const cu_qcomp& b) {
    a = a * b;
}



/*
 * CASTS BETWEEN qcomp AND cu_qcomp
 */


__host__ inline cu_qcomp toCuQcomp(qcomp a) {
    return (cu_qcomp) {.x = real(a), .y = imag(a)};
}


__host__ inline cu_qcomp* toCuQcomps(qcomp* a) {
    return reinterpret_cast<cu_qcomp*>(a);
}



/*
 * MATRIX CASTING AND UNPACKING
 */


__host__ inline std::array<cu_qcomp,4> unpackMatrixToCuQcomps(CompMatr1 in) {

    std::array<cu_qcomp,4> arr{};
    for (int i=0; i<4; i++)
        arr[i] = toCuQcomp(in.elems[i/2][i%2]);

    return arr;
}



#endif // GPU_TYPES_HPP