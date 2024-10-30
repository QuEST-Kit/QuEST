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
#include <vector>
#include <cuComplex.h>



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


// TODO:
// clean this up with templates!


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
INLINE cu_qcomp operator + (const qreal& b, const cu_qcomp& a) {
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
INLINE cu_qcomp operator - (const qreal& b, const cu_qcomp& a) {
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
INLINE cu_qcomp operator * (const qreal& b, const cu_qcomp& a) {
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


INLINE void operator += (cu_qcomp& a, const qreal& b) {
    a = a + b;
}

INLINE void operator -= (cu_qcomp& a, const qreal& b) {
    a = a - b;
}

INLINE void operator *= (cu_qcomp& a, const qreal& b) {
    a = a * b;
}



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
    cu_qcomp out = {.x = re, .y = im};
    return out;
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


__host__ inline std::array<cu_qcomp,16> unpackMatrixToCuQcomps(CompMatr2 in) {

    std::array<cu_qcomp,16> arr{};
    for (int i=0; i<16; i++)
        arr[i] = toCuQcomp(in.elems[i/4][i%4]);

    return arr;
}



#endif // GPU_TYPES_HPP