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
#include "quest/include/precision.h"

#include "quest/src/core/inliner.hpp"

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_types.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#if defined(__NVCC__)
    #include <cuComplex.h>
#elif defined(__HIPCC__)
    #include "quest/src/gpu/cuda_to_hip.hpp"
#endif

#include <array>
#include <vector>



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
 * CASTS BETWEEN qcomp AND cu_qcomp
 */


INLINE cu_qcomp getCuQcomp(qreal re, qreal im) {

#if (FLOAT_PRECISION == 1)
    return make_cuFloatComplex(re, im);
#else
    return make_cuDoubleComplex(re, im);
#endif
}


__host__ inline cu_qcomp toCuQcomp(qcomp a) {
    return getCuQcomp(real(a), imag(a));
}
__host__ inline qcomp toQcomp(cu_qcomp a) {
    return getQcomp(a.x, a.y);
}


__host__ inline cu_qcomp* toCuQcomps(qcomp* a) {
    return reinterpret_cast<cu_qcomp*>(a);
}



/*
 * cu_qcomp ARITHMETIC OVERLOADS
 *
 * which are only needed by NVCC because
 * HIP defines them for us. This good deed
 * goes punished; a HIP bug disables our
 * use of *= and += overloads, so kernels.cuh
 * has disgusting (x = x * y) statements. Bah!
 */


// TODO:
// - clean this up (with templates?)
// - use getCuQcomp() rather than struct creation,
//   to make the algebra implementation-agnostic


#if defined(__NVCC__)

INLINE cu_qcomp operator + (const cu_qcomp& a, const cu_qcomp& b) {
    cu_qcomp out = {
        .x = a.x + b.x,
        .y = a.y + b.y
    };
    return out;
}

INLINE cu_qcomp operator - (const cu_qcomp& a, const cu_qcomp& b) {
    cu_qcomp out = {
        .x = a.x - b.x,
        .y = a.y - b.y
    };
    return out;
}

INLINE cu_qcomp operator * (const cu_qcomp& a, const cu_qcomp& b) {
    cu_qcomp out = {
        .x = a.x * b.x - a.y * b.y,
        .y = a.x * b.y + a.y * b.x
    };
    return out;
}


INLINE cu_qcomp operator + (const cu_qcomp& a, const qreal& b) {
    cu_qcomp out = {
        .x = a.x + b,
        .y = a.y + b
    };
    return out;
}
INLINE cu_qcomp operator + (const qreal& b, const cu_qcomp& a) {
    cu_qcomp out = {
        .x = a.x + b,
        .y = a.y + b
    };
    return out;
}

INLINE cu_qcomp operator - (const cu_qcomp& a, const qreal& b) {
    cu_qcomp out = {
        .x = a.x - b,
        .y = a.y - b
    };
    return out;
}
INLINE cu_qcomp operator - (const qreal& b, const cu_qcomp& a) {
    cu_qcomp out = {
        .x = a.x - b,
        .y = a.y - b
    };
    return out;
}

INLINE cu_qcomp operator * (const cu_qcomp& a, const qreal& b) {
    cu_qcomp out = {
        .x = a.x * b,
        .y = a.y * b
    };
    return out;
}
INLINE cu_qcomp operator * (const qreal& b, const cu_qcomp& a) {
    cu_qcomp out = {
        .x = a.x * b,
        .y = a.y * b
    };
    return out;
}

#endif



/*
 * cu_qcomp UNARY FUNCTIONS
 */


INLINE qreal getCompReal(cu_qcomp num) {
    return num.x;
}

INLINE cu_qcomp getCompConj(cu_qcomp num) {
    num.y *= -1;
    return num;
}

INLINE qreal getCompNorm(cu_qcomp num) {
    return (num.x * num.x) + (num.y * num.y);
}

INLINE cu_qcomp getCompPower(cu_qcomp base, cu_qcomp exponent) {

    // using https://mathworld.wolfram.com/ComplexExponentiation.html,
    // and the principal argument of 'base'

    // base = a + b i, exponent = c + d i
    qreal a = base.x;
    qreal b = base.y;
    qreal c = exponent.x;
    qreal d = exponent.y;

    // intermediate quantities
    qreal arg = atan2(b, a);
    qreal mag = a*a + b*b;
    qreal ln = log(mag);
    qreal fac = pow(mag, c/2) * exp(-d * arg);
    qreal ang = c*arg + d*ln/2;

    // output scalar
    qreal re = fac * cos(ang);
    qreal im = fac * sin(ang);
    return getCuQcomp(re, im);
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


__host__ inline std::array<cu_qcomp,16> unpackMatrixToCuQcomps(CompMatr2 in) {

    std::array<cu_qcomp,16> arr{};
    for (int i=0; i<16; i++)
        arr[i] = toCuQcomp(in.elems[i/4][i%4]);

    return arr;
}



#endif // GPU_TYPES_HPP