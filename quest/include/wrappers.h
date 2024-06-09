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



#endif // WRAPPERS_H