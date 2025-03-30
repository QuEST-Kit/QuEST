/** @file
 * C-compatible functions which are alternatives to C++-only API
 * functions, ultimately providing an identical interface. This is 
 * necessary because these functions otherwise return qcomps by-value
 * which is prohibited between C and C++ compiled binaries (because
 * complex numbers are not agreed upon in their ABI, despite having 
 * identical memory layouts in the C and C++ standard libraries).
 * Ergo this file defines no new API functions as far as the user/
 * documentation is aware, but secretly ensures the backend C++ 
 * binaries return qcomps to the user's C code only by pointer.
 * 
 * (( _passing_ qcomps by value to a function seems to be okay,
 *    although I am not entirely sure why ))
 * 
 * Note that matrix getters and setters (like getCompMatr1()) are
 * excluded, and instead defined directly in matrices.h/.cpp .
 * 
 * An unimportant by-product of this method of achieving interoperability
 * is that the internal wrapped functions are exposed to the C user; so
 * we prefix them with "_wrap_" to imply privacy. Note the "extern"
 * declarations are superfluous (the behaviour is default), but used
 * to explicitly distinguish the intendedly-private internal functions
 * from the C API functions herein defined.
 * 
 * @author Tyson Jones
 * 
 * (no doxygen doc)
 */

#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

/// @cond EXCLUDE_FROM_DOXYGEN

// these definitions are only exposed to C, 
// since they duplicate existing C++ functions
#ifndef __cplusplus



extern void _wrap_calcInnerProduct(Qureg bra, Qureg ket, qcomp* out);

qcomp calcInnerProduct(Qureg bra, Qureg ket) {

    qcomp out;
    _wrap_calcInnerProduct(bra, ket, &out);
    return out;
}


extern void _wrap_calcExpecNonHermitianPauliStrSum(qcomp*, Qureg, PauliStrSum);

qcomp calcExpecNonHermitianPauliStrSum(Qureg qureg, PauliStrSum sum) {

    qcomp out;
    _wrap_calcExpecNonHermitianPauliStrSum(&out, qureg, sum);
    return out;
}


extern void _wrap_calcExpecNonHermitianFullStateDiagMatr(qcomp*, Qureg, FullStateDiagMatr);

qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr) {

    qcomp out;
    _wrap_calcExpecNonHermitianFullStateDiagMatr(&out, qureg, matr);
    return out;
}


extern void _wrap_calcExpecNonHermitianFullStateDiagMatrPower(qcomp*, Qureg, FullStateDiagMatr, qcomp);

qcomp calcExpecNonHermitianFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matr, qcomp expo) {

    qcomp out;
    _wrap_calcExpecNonHermitianFullStateDiagMatrPower(&out, qureg, matr, expo);
    return out;
}


extern void _wrap_getQuregAmp(qcomp* out, Qureg qureg, qindex index);

qcomp getQuregAmp(Qureg qureg, qindex index) {

    qcomp out;
    _wrap_getQuregAmp(&out, qureg, index);
    return out;
}


extern void _wrap_getDensityQuregAmp(qcomp* out, Qureg qureg, qindex row, qindex column);

qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column) {

    qcomp out;
    _wrap_getDensityQuregAmp(&out, qureg, row, column);
    return out;
}



#endif // !__cplusplus

/// @endcond // EXCLUDE_FROM_DOXYGEN

#endif // WRAPPERS_H