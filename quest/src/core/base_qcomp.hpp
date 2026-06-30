/** @file
 * Definition of base_qcomp, which is extended by the CPU and GPU
 * backends (into cpu_qcomp and gpu_qcomp) and used in hot loops
 * and kernels.
 * 
 * The user-facing qcomp (which in the QuEST middle-end, resolves to 
 * std::complex) is not used by the CPU backend, since it creates
 * performance pitfalls (e.g. expensive NaN checks within arithmetic
 * operators) in some compilers, and is furthermore illegal in the
 * GPU backend (i.e. within CUDA kernels). So the backends instead
 * use custom complex types with identical memory layouts/alignment
 * to qcomp. Those types extend base_qcomp defined in this file,
 * since they otherwise share all the same arithmetic boilerplate.
 *
 * Beware that this file is parsed by both the CPU and GPU compiler,
 * for which the meaning of INLINE is different, and so all INLINE
 * functions must be both OpenMP and CUDA compatible. Non-inline
 * functions are not permitted since this header is included by
 * multiple src files.
 * 
 * @author Tyson Jones
 */

#ifndef BASE_QCOMP_HPP
#define BASE_QCOMP_HPP

#include "quest/include/types.h"

#include "quest/src/core/inliner.hpp"



/*
 * BASE DEFINITION
 *
 * which must remain POD (a simple {re,im}) and with an identical
 * memory layout and alignment to qcomp (i.e. std::complex). Only
 * the in-place arithmetic overloads are defined below which are
 * reused by the subsequent out-of-place overloads, to avoid
 * code duplication.
 */

struct alignas(qcomp) base_qcomp {

    qreal re;
    qreal im;


    /*
     * IN-PLACE COMPLEX ARITHMETIC
     */
    
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


    /*
     * IN-PLACE MIXED-TYPE ARITHMETIC
     */

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

    INLINE base_qcomp& operator *= (const size_t& a) noexcept {
        re *= a;
        im *= a;
        return *this;
    }

}; // base_qcomp



/*
 * OUT-OF-PLACE COMPLEX ARITHMETIC
 * 
 * which avoid code duplication by re-using the
 * in-place arithmetic operator overloads above
 */

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



/*
 * OUT-OF-PLACE MIXED-TYPE ARITHMETIC
 * 
 * which avoid code duplication by re-using the
 * in-place arithmetic operator overloads above
 */


// base_qcomp * other

INLINE base_qcomp operator * (base_qcomp a, const int& b) noexcept {
    a *= b;
    return a;
}

INLINE base_qcomp operator * (base_qcomp a, const qreal& b) noexcept {
    a *= b;
    return a;
}

INLINE base_qcomp operator * (base_qcomp a, const size_t& b) noexcept {
    a *= b;
    return a;
}


// other * base_qcomp (via commutation)

INLINE base_qcomp operator * (const int& a, const base_qcomp& b) noexcept {
    return b * a;
}

INLINE base_qcomp operator * (const qreal& a, const base_qcomp& b) noexcept {
    return b * a;
}

INLINE base_qcomp operator * (const size_t& a, const base_qcomp& b) noexcept {
    return b * a;
}



/*
 * BACKEND-AGNOSTIC MATHS
 */

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



/*
 * CONVERTERS
 */

INLINE base_qcomp* getBaseQcompPtr(qcomp* list) {
    return reinterpret_cast<base_qcomp*>(list);
}

INLINE base_qcomp getBaseQcomp(qreal re, qreal im) {
    return { re, im };
}



/*
 * CHECK COMPATIBILITY WITH QCOMP
 */


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



#endif // BASE_QCOMP_HPP