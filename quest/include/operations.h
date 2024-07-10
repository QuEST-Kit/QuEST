/** @file
 * API signatures for effecting operators (such as gates and unitaries) 
 * upon Quregs which are instantiated as either statevectors or 
 * density matrices.
 */

#ifndef OPERATIONS_H
#define OPERATIONS_H

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/*
 * GATES
 */

void unitary(Qureg qureg, int target, CompMatr1 matrix);



// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // OPERATIONS_H