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

    localiser_statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, target, matrix);
    if (qureg.isDensityMatrix)
        localiser_statevec_anyCtrlOneTargDenseMatr(qureg, {}, {}, util_getBraQubit(target, qureg), util_getConj(matrix));
}


} // end de-mangler