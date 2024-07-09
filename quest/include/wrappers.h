/** @file
 * C-compatible functions which are alternatives to C++-only API
 * functions, ultimately providing an identical interface. This is 
 * necessary because these functions otherwise pass qcomps by-value
 * which is prohibited between C and C++ compiled binaries (because
 * complex numbers are not specified in the ABI, despite having 
 * identical memory layouts in the C and C++ standard libraries).
 * Ergo this file defines no new API functions as far as the user/
 * documentation is aware, but secretly ensures C binaries receive
 * qcomps from the C++ backend only by-reference. It must be used
 * by C compilers which otherwise lack the C++-only API signatures.
 */

#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "types.h"
#include "structures.h"

// these definitions are only exposed to C, 
// since they duplicate existing C++ functions
#ifndef __cplusplus



void wrap_getCompMatr1FromArr(CompMatr1* out, qcomp in[2][2]);

CompMatr1 getCompMatr1FromArr(qcomp in[2][2]) {

    CompMatr1 out;
    wrap_getCompMatr1FromArr(&out, in);
    return out;
}


void wrap_getCompMatr1FromPtr(CompMatr1* out, qcomp** in);

CompMatr1 getCompMatr1FromPtr(qcomp** in) {

    CompMatr1 out;
    wrap_getCompMatr1FromPtr(&out, in);
    return out;
}


void wrap_getCompMatr2FromArr(CompMatr2* out, qcomp in[4][4]);

CompMatr2 getCompMatr2FromArr(qcomp in[4][4]) {

    CompMatr2 out;
    wrap_getCompMatr2FromArr(&out, in);
    return out;
}


void wrap_getCompMatr2FromPtr(CompMatr2* out, qcomp** in);

CompMatr2 getCompMatr2FromPtr(qcomp** in) {

    CompMatr2 out;
    wrap_getCompMatr2FromPtr(&out, in);
    return out;
}



#endif // !__cplusplus

#endif // WRAPPERS_H