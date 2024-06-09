/** @file
 * API signatures for debugging QuEST behaviour, or controlling
 * input validation.
 */

#ifndef DEBUG_H
#define DEBUG_H

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



void invalidQuESTInputError(const char* msg, const char* func);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // DEBUG_H