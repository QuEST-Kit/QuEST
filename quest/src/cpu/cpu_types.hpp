/** @file
 * Custom types used exclusively by the CPU backend.
 * 
 * @author Tyson Jones
 */

#ifndef CPU_TYPES_HPP
#define CPU_TYPES_HPP

#include "quest/include/types.h"

#include "quest/src/core/inliner.hpp"
#include "quest/src/core/basetypes.hpp"

#include <array>


typedef base_qcomp cpu_qcomp;


INLINE cpu_qcomp* getCpuQcompPtr(qcomp* list) {
    return getBaseQcompPtr(list);
}
INLINE cpu_qcomp getCpuQcomp(qreal re, qreal im) {
    return getBaseQcomp(re, im);
}
INLINE cpu_qcomp getCpuQcomp(const qcomp& a) {
    return getBaseQcomp(a);
}


// backend specific maths functions
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


// creator for fixed-size dense matrices (CompMatr1 and CompMatr2) ((not inlined!))
template <int dim>
std::array<std::array<cpu_qcomp,dim>,dim> getCpuQcomps(qcomp matr[dim][dim]) {

    // detect brain-dead compiler inferencing (looking at you MSVC...)
    static_assert(dim == 2 || dim == 4, "getCpuQcomps called with unexpected dim");

    std::array<std::array<cpu_qcomp,dim>,dim> out;

    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
            out[i][j] = getCpuQcomp(matr[i][j]);

    return out;
}


// check the memory layout of cpu_qcomp agrees with qcomp, since
// it is not formally gauranteed, unlike _Complex and std::complex
static_assert(sizeof (cpu_qcomp) == sizeof (qcomp));
static_assert(alignof(cpu_qcomp) == alignof(qcomp));
static_assert(std::is_standard_layout_v   <cpu_qcomp>);
static_assert(std::is_trivially_copyable_v<cpu_qcomp>);


#endif // CPU_TYPES_HPP