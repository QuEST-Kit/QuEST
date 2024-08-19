/** @file
 * CUDA GPU-accelerated signatures of the subroutines called by accelerator.cpp.
 */

#ifndef GPU_SUBROUTINES_HPP
#define GPU_SUBROUTINES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include <vector>

using std::vector;


/*
 * COMMUNICATION BUFFER PACKING
 */

template <int NumQubits> void gpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubits, vector<int> qubitStates);


/*
 * SWAPS
 */

template <int NumCtrls> void gpu_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2);
template <int NumCtrls> void gpu_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates);
template <int NumCtrls> void gpu_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState);


/*
 * DENSE MATRIX
 */

template <int NumCtrls> void gpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
template <int NumCtrls> void gpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);

template <int NumCtrls, int NumTargs> void gpu_statevec_anyCtrlAnyTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr);


/*
 * DIAGONAL MATRIX
 */

template <int NumCtrls, int NumTargs> void gpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr);


#endif // GPU_SUBROUTINES_HPP