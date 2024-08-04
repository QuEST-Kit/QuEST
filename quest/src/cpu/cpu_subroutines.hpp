/** @file
 * CPU signatures of the subroutines called by accelerator.cpp. 
 */

#ifndef CPU_SUBROUTINES_HPP
#define CPU_SUBROUTINES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include <vector>

using std::vector;
using namespace index_flags;


template <CtrlFlag flag> qindex cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates);

template <CtrlFlag flag> void cpu_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
template <CtrlFlag flag> void cpu_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr);
template <CtrlFlag flag> void cpu_statevec_anyCtrlOneTargDenseMatrix_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);



#endif // CPU_SUBROUTINES_HPP