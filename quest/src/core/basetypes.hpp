/** @file
 * Custom types used exclusively by the hot loops in the
 * hardware accelerated backends
 * 
 * @author Tyson Jones
 */

#ifndef BASETYPES_HPP
#define BASETYPES_HPP

#include "quest/include/types.h"

#include "quest/src/core/inliner.hpp"


struct base_qcomp {

    // memory layout
    qreal re;
    qreal im;

    // in-place complex arithmetic overloads
    INLINE base_qcomp& operator += (const base_qcomp& a) noexcept {
        re += a.re;
        im += a.im;
        return *this;
    }
    INLINE base_qcomp& operator -= (const base_qcomp& a) noexcept {
        re -= a.re;
        im -= a.im;
        return *this;
    }
    INLINE base_qcomp& operator *= (const base_qcomp& a) noexcept {
        qreal re_ = re;
        qreal im_ = im;
        re = (re_ * a.re) - (im_ * a.im);
        im = (re_ * a.im) + (im_ * a.re);
        return *this;
    }

    // in-place mixed-type arithmetic overloads
    INLINE base_qcomp& operator *= (const int& a) noexcept {
        re *= a;
        im *= a;
        return *this;
    }
    INLINE base_qcomp& operator *= (const qreal& a) noexcept {
        re *= a;
        im *= a;
        return *this;
    }
};


// out-of-place complex arithmetic overloads (optimised)
INLINE base_qcomp operator + (base_qcomp a, const base_qcomp& b) noexcept {
    a += b;
    return a;
}
INLINE base_qcomp operator - (base_qcomp a, const base_qcomp& b) noexcept {
    a -= b;
    return a;
}
INLINE base_qcomp operator * (base_qcomp a, const base_qcomp& b) noexcept {
    a *= b;
    return a;
}


// out-of-place mixed-type arithmetic overloads
INLINE base_qcomp operator * (base_qcomp a, const int& b) noexcept {
    a *= b;
    return a;
}
INLINE base_qcomp operator * (base_qcomp a, const qreal& b) noexcept {
    a *= b;
    return a;
}


// reverse order of out-of-place mixed-type arithmetic (via commutation)
INLINE base_qcomp operator * (const int& a, const base_qcomp& b) noexcept {
    return b * a;
}
INLINE base_qcomp operator * (const qreal& a, const base_qcomp& b) noexcept {
    return b * a;
}


// backend agnostic maths
INLINE qreal real(const base_qcomp& a) {
    return a.re;
}
INLINE qreal imag(const base_qcomp& a) {
    return a.im;
}
INLINE base_qcomp conj(const base_qcomp& a) {
    return {a.re, - a.im};
}
INLINE qreal norm(const base_qcomp& a) noexcept {
    return (a.re * a.re) + (a.im * a.im);
}


// backend specific maths must be defined elsewhere

    // INLINE base_qcomp pow(base_qcomp base, base_qcomp expo) noexcept {

    //     // Here, we re-use std::pow(std::complex) to avoid a custom definition,
    //     // and so accept NaN-check performance penalties. Notice too we also
    //     // create new qcomp(), rather than just reinterpreting the given base_qcomp,
    //     // just to avoid any insiduous issues alignment/aliasing issues (since the
    //     // creation time iss occluded by std::pow time).
    //     qcomp base_ = getQcomp(base);
    //     qcomp expo_ = getQcomp(expo);
    //     qcomp out_ = std::pow(base_, expo_);
    //     return getCpuQcomp(out_);
    // }


// (base) creators
INLINE base_qcomp* getBaseQcompPtr(qcomp* list) {
    return reinterpret_cast<base_qcomp*>(list);
}
INLINE base_qcomp getBaseQcomp(qreal re, qreal im) {
    return { re, im };
}
INLINE base_qcomp getBaseQcomp(const qcomp& a) {
    return { a.real(), a.imag() };
}
INLINE qcomp getQcomp(const base_qcomp& a) {
    return qcomp( a.re, a.im );
}


// creator for fixed-size dense matrices (CompMatr1 and CompMatr2) ((not inlined!))
    // template <int dim>
    // std::array<std::array<base_qcomp,dim>,dim> getCpuQcomps(qcomp matr[dim][dim]) {

    //     // detect brain-dead compiler inferencing (looking at you MSVC...)
    //     static_assert(dim == 2 || dim == 4, "getCpuQcomps called with unexpected dim");

    //     std::array<std::array<base_qcomp,dim>,dim> out;

    //     for (int i=0; i<dim; i++)
    //         for (int j=0; j<dim; j++)
    //             out[i][j] = getCpuQcomp(matr[i][j]);

    //     return out;
    // }








// check the memory layout of base_qcomp agrees with qcomp, since
// it is not formally gauranteed, unlike _Complex and std::complex
static_assert(sizeof (base_qcomp) == sizeof (qcomp));
static_assert(alignof(base_qcomp) == alignof(qcomp));
static_assert(std::is_standard_layout_v   <base_qcomp>);
static_assert(std::is_trivially_copyable_v<base_qcomp>);


// TODO:
// the above checks are potentially inadequate to identify an
// insidious incompatibility between qcomp and base_qcomp - perhaps
// we should perform a compile-time duck-check, casting a small
// array between them and checking no data is corrupted? Perhaps
// a runtime check in initQuESTEnv() is also necessary, checking the
// casting is safe for all circumstances (e.g. heap mem, static lists)


#endif // BASETYPES_HPP