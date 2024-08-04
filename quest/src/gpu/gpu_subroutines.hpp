/** @file
 * CUDA GPU-accelerated signatures of the subroutines called by accelerator.cpp.
 */

#ifndef GPU_SUBROUTINES_HPP
#define GPU_SUBROUTINES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/indexer.hpp"

#include <vector>

using std::vector;
using namespace index_flags;



template <CtrlFlag flag> qindex gpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates);

template <CtrlFlag flag> void gpu_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
template <CtrlFlag flag> void gpu_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr);
template <CtrlFlag flag> void gpu_statevec_anyCtrlOneTargDenseMatrix_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);



#endif // GPU_SUBROUTINES_HPP