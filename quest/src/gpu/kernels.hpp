/** @file
 * Custom CUDA kernels invoked by gpu.cpp, usually only necessary when
 * there is no equivalent utility in Thrust (or cuQuantum, when it is
 * targeted).
 */

#ifndef KERNELS_HPP
#define KERNELS_HPP

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include <cuComplex.h>



/*
 * CREATE CUDA-COMPATIBLE QCOMP ALIAS
 */

#if (FLOAT_PRECISION == 1)
    typedef cuFloatComplex cu_qcomp;

#elif (FLOAT_PRECISION == 2)
    typedef cuDoubleComplex cu_qcomp;

#else
    #error "Build bug; precision.h should have prevented non-float non-double qcomp precision on GPU."

#endif



/*
 * CREATE CUDA-COMPATIBLE QCOMP OVERLOADS
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



#endif // KERNELS_HPP