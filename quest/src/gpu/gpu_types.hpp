/** @file
 * CUDA-compatible complex types. This file is only ever included
 * when ENABLE_GPU_ACCELERATION=1 so it can safely invoke CUDA
 * signatures without guards.
 */

#ifndef GPU_TYPES_HPP
#define GPU_TYPES_HPP

#include "quest/include/modes.h"
#include "quest/include/types.h"

#if ! ENABLE_GPU_ACCELERATION
    #error "A file being compiled somehow included gpu_types.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif



/*
 * CUDA-COMPATIBLE QCOMP ALIAS (cu_qcomp)
 */

#include <cuComplex.h>

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

__host__ __device__ inline cu_qcomp operator + (const cu_qcomp& a, const cu_qcomp& b) {
    cu_qcomp res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    return res;
}

__host__ __device__ inline cu_qcomp operator * (const cu_qcomp& a, const cu_qcomp& b) {
    cu_qcomp res;
    res.x = a.x * b.x - a.y * b.y;
    res.y = a.x * b.y + a.y * b.x;
    return res;
}



/*
 * CASTS BETWEEN qcomp AND cu_qcomp
 */

__host__ cu_qcomp toCuQcomp(qcomp x) {
    return reinterpret_cast<cu_qcomp&>(x);
}

__host__ cu_qcomp* toCuQcomps(qcomp* x) {
    return reinterpret_cast<cu_qcomp*>(x);
}



#endif // GPU_TYPES_HPP