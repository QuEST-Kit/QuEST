/** @file
 * CUDA-compatible complex types. This file is only ever included
 * when ENABLE_GPU_ACCELERATION=1 so it can safely invoke CUDA
 * signatures without guards. This is safe to re-include by
 * multiple files because typedef definition is legal in C++,
 * and all functions herein are inline. Furthermore, since it
 * is only ever parsed by nvcc, the __host__ symbols are safely
 * processed by the cuquantum backend.
 */

#ifndef GPU_TYPES_HPP
#define GPU_TYPES_HPP

#include "quest/include/modes.h"
#include "quest/include/types.h"

#if ! ENABLE_GPU_ACCELERATION
    #error "A file being compiled somehow included gpu_types.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#include <cuComplex.h>
#include <vector>



/*
 * CUDA-COMPATIBLE QCOMP ALIAS (cu_qcomp)
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

__host__ inline cu_qcomp toCuQcomp(qcomp x) {
    return reinterpret_cast<cu_qcomp&>(x);
}

__host__ inline cu_qcomp* toCuQcomps(qcomp* x) {
    return reinterpret_cast<cu_qcomp*>(x);
}



/*
 * MATRIX CASTING AND UNPACKING
 */

__host__ inline void unpackMatrixToCuQcomps(CompMatr1 in, cu_qcomp &m00, cu_qcomp &m01, cu_qcomp &m10, cu_qcomp &m11) {

    m00 = toCuQcomp(in.elems[0][0]);
    m01 = toCuQcomp(in.elems[0][1]);
    m10 = toCuQcomp(in.elems[1][0]);
    m11 = toCuQcomp(in.elems[1][1]);
}

__host__ inline std::vector<cu_qcomp> unpackMatrixToCuQcomps(CompMatr1 in) {

    std::vector<cu_qcomp> vec(4);
    vec[0] = toCuQcomp(in.elems[0][0]);
    vec[1] = toCuQcomp(in.elems[0][1]);
    vec[2] = toCuQcomp(in.elems[1][0]);
    vec[3] = toCuQcomp(in.elems[1][1]);
    return vec;
}


#endif // GPU_TYPES_HPP