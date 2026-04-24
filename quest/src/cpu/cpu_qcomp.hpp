/** @file
 * Definition of cpu_qcomp, an extension of base_qcomp and a
 * compatible alternative to the user-facing qcomp, used
 * exclusively by the CPU backend. This custom complex type
 * avoids performance pitfalls of qcomp (i.e. std::complex),
 * e.g. due to NaN checks, without resorting to platform
 * and compiler-specific flags.
 * 
 * @author Tyson Jones
 */

#ifndef CPU_QCOMP_HPP
#define CPU_QCOMP_HPP

#include "quest/include/types.h"

#include "quest/src/core/inliner.hpp"
#include "quest/src/core/base_qcomp.hpp"

#include <array>



/*
 * DEFINE CPU_QCOMP
 *
 * which is safe to typdef and define additional overloads
 * below, since never witnessed outside the CPU Backend
 */

typedef base_qcomp cpu_qcomp;



/*
 * CONVERTERS
 *
 * which merely wrap the base_qcomp functions for clarity
 * in the CPU source code, disambiguating from gpu_qcomp
 */

INLINE cpu_qcomp* getCpuQcompPtr(qcomp* ptr) {
    return getBaseQcompPtr(ptr);
}

INLINE cpu_qcomp getCpuQcomp(qreal re, qreal im) {
    return getBaseQcomp(re, im);
}

INLINE cpu_qcomp getCpuQcomp(const qcomp& a) {
    return getBaseQcomp(a.real(), a.imag());
}

INLINE qcomp getQcomp(const cpu_qcomp& a) {
    return qcomp( a.re, a.im );
}

template <int Dim>
std::array<std::array<cpu_qcomp,Dim>,Dim> getCpuQcompsMatrix(qcomp matr[Dim][Dim]) {

    // Creator for fixed-size dense matrices CompMatr1 and CompMatr2,
    // which are respectively 2x2 and 4x4 - deliberately not inlined!
    // We create new cpu_qcomp in lieu of reinterpreting a 2D pointer
    // in fear of alignment and static array nightmares
    static_assert(Dim == 2 || Dim == 4);

    std::array<std::array<cpu_qcomp,Dim>,Dim> out;

    for (int i=0; i<Dim; i++)
        for (int j=0; j<Dim; j++)
            out[i][j] = getCpuQcomp(matr[i][j]);

    return out;
}



/*
 * CPU-SPECIFIC MATHS
 */

INLINE cpu_qcomp pow(cpu_qcomp base, cpu_qcomp expo) noexcept {

    // Here, we re-use std::pow(std::complex) to avoid a custom definition,
    // and so accept NaN-check performance penalties. Notice too we also
    // create new qcomp(), rather than just reinterpreting the given cpu_qcomp,
    // just to avoid any insiduous issues alignment/aliasing issues (since the
    // creation time iss occluded by std::pow time).

    qcomp base_ = getQcomp(base);
    qcomp expo_ = getQcomp(expo);
    qcomp out_ = std::pow(base_, expo_);
    return getCpuQcomp(out_);
}



#endif // CPU_QCOMP_HPP