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

#include "quest/include/config.h"
#include "quest/include/types.h"
#include "quest/include/precision.h"

#include "quest/src/core/inliner.hpp"
#include "quest/src/core/basetypes.hpp"

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_types.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#if (FLOAT_PRECISION == 4)
    #error "Build bug; precision.h should have prevented non-float non-double qcomp precision on GPU."
#endif

#if defined(__HIP__)
    #include "quest/src/gpu/cuda_to_hip.hpp"
#endif

#include <array>


typedef base_qcomp gpu_qcomp;


INLINE gpu_qcomp* getGpuQcompPtr(qcomp* list) {
    return getBaseQcompPtr(list);
}
INLINE gpu_qcomp getGpuQcomp(qreal re, qreal im) {
    return getBaseQcomp(re, im);
}
INLINE gpu_qcomp getGpuQcomp(const qcomp& a) {
    return getBaseQcomp(a);
}


// backend specific maths functions
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


// check the memory layout of gpu_qcomp agrees with qcomp, since
// it is not formally gauranteed, unlike _Complex and std::complex
static_assert(sizeof (gpu_qcomp) == sizeof (qcomp));
static_assert(alignof(gpu_qcomp) == alignof(qcomp));
static_assert(std::is_standard_layout_v   <gpu_qcomp>);
static_assert(std::is_trivially_copyable_v<gpu_qcomp>);


// TODO:
// the above checks are potentially inadequate to identify an
// insidious incompatibility between qcomp and gpu_qcomp - perhaps
// we should perform a compile-time duck-check, casting a small
// array between them and checking no data is corrupted? Perhaps
// a runtime check in initQuESTEnv() is also necessary, checking the
// casting is safe for all circumstances (e.g. heap mem, static lists)





__host__ inline std::array<gpu_qcomp,2> unpackMatrixToGpuQcomps(DiagMatr1 in) {

    // it's crucial we explicitly copy over the elements,
    // rather than just reinterpret the pointer, to avoid
    // segmentation faults when memory misaligns (like on HIP)

    return {getGpuQcomp(in.elems[0]), getGpuQcomp(in.elems[1])};
}


__host__ inline std::array<gpu_qcomp,4> unpackMatrixToGpuQcomps(DiagMatr2 in) {

    return {
        getGpuQcomp(in.elems[0]), getGpuQcomp(in.elems[1]),
        getGpuQcomp(in.elems[2]), getGpuQcomp(in.elems[3])};
}


__host__ inline std::array<gpu_qcomp,4> unpackMatrixToGpuQcomps(CompMatr1 in) {

    std::array<gpu_qcomp,4> out{};
    for (int i=0; i<4; i++)
        out[i] = getGpuQcomp(in.elems[i/2][i%2]);

    return out;
}


__host__ inline std::array<gpu_qcomp,16> unpackMatrixToGpuQcomps(CompMatr2 in) {

    std::array<gpu_qcomp,16> out{};
    for (int i=0; i<16; i++)
        out[i] = getGpuQcomp(in.elems[i/4][i%4]);

    return out;
}


#endif // GPU_TYPES_HPP