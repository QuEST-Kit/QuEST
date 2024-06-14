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

#include "quest/include/types.h"
#include "quest/include/structures.h"

// these definitions are only exposed to C, 
// since they duplicate existing C++ functions
#ifndef __cplusplus



void wrap_getCompMatr1(CompMatr1* out, qcomp in[2][2]);

CompMatr1 getCompMatr1(qcomp in[2][2]) {

    CompMatr1 out;
    wrap_getCompMatr1(&out, in);
    return out;
}


void wrap_getCompMatr2(CompMatr2* out, qcomp in[4][4]);

CompMatr2 getCompMatr2(qcomp in[4][4]) {

    CompMatr2 out;
    wrap_getCompMatr2(&out, in);
    return out;
}



#endif // !__cplusplus

#endif // WRAPPERS_H