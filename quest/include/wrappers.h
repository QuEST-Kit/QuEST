/** @file
 * C-compatible functions which are alternatives to C++-only API
 * functions, ultimately providing an identical interface. This is 
 * necessary because these functions otherwise pass qcomps by-value
 * which is prohibited between C and C++ compiled binaries (because
 * complex numbers are not agreed upon in their ABI, despite having 
 * identical memory layouts in the C and C++ standard libraries).
 * Ergo this file defines no new API functions as far as the user/
 * documentation is aware, but secretly ensures the backend C++ 
 * binaries are send qcomps from the user's C code only by pointer.
 * 
 * Note that CompMatr getters and setters (like getCompMatr1()) are
 * missing from this file, even though they contain qcomp[2][2] (which
 * are NOT pointers) and cannot be directly passed between binaries.
 * Those functions are instead defined in structures.h/.cpp because
 * those structs are declared 'const'; they must be intiailised inline, 
 * and can never be modified by pointer. We shouldn't even address them!
 */

#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "quest/include/types.h"
#include "quest/include/structures.h"

// these definitions are only exposed to C, 
// since they duplicate existing C++ functions
#ifndef __cplusplus



// TODO:
//      this file will contain wrappers of functions like C++'s 
//          qcomp getAmp(Qureg, i);
//      which will (for C users) be secretly invoking something like:
//
// void wrap_getAmp(Qureg qureg, qindex i, qcomp* amp); // defined by C++ backend
//
// void getAmp(Qureg qureg, qindex i) {
//     qcomp amp;
//     wrap_getAmp(qureg, i, &amp);
//     return amp;
// }



#endif // !__cplusplus

#endif // WRAPPERS_H