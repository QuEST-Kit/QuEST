/** @file
 * API definitions for effecting operators (such as gates and unitaries) 
 * upon Quregs which are instantiated as either statevectors or 
 * density matrices.
 */

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/localiser.hpp"


// enable invocation by both C and C++ binaries
extern "C" {


void unitary(Qureg qureg, int target, CompMatr1 matrix) {
    validate_quregFields(qureg, __func__);
    validate_target(qureg, target, __func__);
    validate_matrixIsUnitary(matrix, __func__);

    statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, target, matrix);
    if (qureg.isDensityMatrix)
        statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, util_getShifted(target, qureg), util_getConj(matrix));
}




// TODO: prefix with localiser, c'mon man


// DEBUG - just for benchmarking!

void noCtrlGate(Qureg qureg, int target, CompMatr1 matrix) {

    statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, target, matrix);
}

void oneCtrlGate(Qureg qureg, int control, int target, CompMatr1 matrix) {

    statevec_anyCtrlOneTargDenseMatr(qureg, {control}, {}, target, matrix);
} 
void oneStateCtrlGate(Qureg qureg, int control, int controlState, int target, CompMatr1 matrix) {

    statevec_anyCtrlOneTargDenseMatr(qureg, {control}, {controlState}, target, matrix);
}






} // end de-mangler