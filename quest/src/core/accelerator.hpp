/** @file
 * Internal functions for choosing which accelerator backend
 * (CPU or GPU) to call, and which qubit preconditions are
 * satisfied in order to accelerate simulation.
 */

#ifndef ACCELERATOR_HPP
#define ACCELERATOR_HPP

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include <vector>

using std::vector;



qindex accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates);

void accel_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
void accel_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr);
void accel_statevec_anyCtrlOneTargDenseMatrix_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);



#endif // ACCELERATOR_HPP