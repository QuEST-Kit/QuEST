/** @file
 * CUDA and HIP-compatible complex types. This file is only ever included
 * when COMPILE_CUDA=1 so it can safely invoke CUDA signatures without guards. 
 * 
 * This header is safe to re-include by multiple files because typedef 
 * redefinition is legal in C++, and all functions herein are inline. 
 * Furthermore, since it is only ever parsed by nvcc, the __host__ symbols 
 * are safely processed by other nvcc-only GPU files, like the cuquantum backend.
 * 
 * @author Tyson Jones
 * @author Oliver Brown (patched HIP arithmetic overloads)
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
#elif defined(__HIP__)
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
 * TRANSFORMING qcomp AND cu_qcomp
 */


INLINE cu_qcomp getCuQcomp(qreal re, qreal im) {

#if (FLOAT_PRECISION == 1)
    return make_cuFloatComplex(re, im);
#else
    return make_cuDoubleComplex(re, im);
#endif
}


__host__ inline cu_qcomp toCuQcomp(qcomp a) {
    return getCuQcomp(std::real(a), std::imag(a));
}
__host__ inline qcomp toQcomp(cu_qcomp a) {
    return getQcomp(a.x, a.y);
}


__host__ inline cu_qcomp* toCuQcomps(qcomp* a) {

    // reinterpret a qcomp ptr as a cu_qcomp ptr,
    // which is ONLY SAFE when comp and cu_qcomp 
    // have identical memory layouts. Be very
    // careful; HIP stack arrays (e.g. qcomp[])
    // seg-fault when passed here, so this funciton
    // should only ever be used on malloc'd data!
    // Stack objects should use the below unpacks.

    return reinterpret_cast<cu_qcomp*>(a);
}


__host__ inline std::array<cu_qcomp,2> unpackMatrixToCuQcomps(DiagMatr1 in) {

    // it's crucial we explicitly copy over the elements,
    // rather than just reinterpret the pointer, to avoid
    // segmentation faults when memory misaligns (like on HIP)

    return {toCuQcomp(in.elems[0]), toCuQcomp(in.elems[1])};
}


__host__ inline std::array<cu_qcomp,4> unpackMatrixToCuQcomps(DiagMatr2 in) {

    return {
        toCuQcomp(in.elems[0]), toCuQcomp(in.elems[1]),
        toCuQcomp(in.elems[2]), toCuQcomp(in.elems[3])};
}


__host__ inline std::array<cu_qcomp,4> unpackMatrixToCuQcomps(CompMatr1 in) {

    std::array<cu_qcomp,4> out{};
    for (int i=0; i<4; i++)
        out[i] = toCuQcomp(in.elems[i/2][i%2]);

    return out;
}


__host__ inline std::array<cu_qcomp,16> unpackMatrixToCuQcomps(CompMatr2 in) {

    std::array<cu_qcomp,16> out{};
    for (int i=0; i<16; i++)
        out[i] = toCuQcomp(in.elems[i/4][i%4]);

    return out;
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


/// @todo
/// - clean this up (with templates?)
/// - use getCuQcomp() rather than struct creation,
///   to make the algebra implementation-agnostic


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

    // intermediate quantities (uses CUDA atan2,log,pow,exp,cos,sin)
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



#endif // GPU_TYPES_HPP