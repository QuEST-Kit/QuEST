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


/*
 * COMMUNICATION BUFFER PACKING
 */

template <int NumQubits> void cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubitInds, vector<int> qubitStates);


/*
 * SWAPS
 */

template <int NumCtrls> void cpu_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2);
template <int NumCtrls> void cpu_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates);
template <int NumCtrls> void cpu_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState);


/*
 * DENSE MATRIX
 */

template <int NumCtrls> void cpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
template <int NumCtrls> void cpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);

template <int NumCtrls, int NumTargs> void cpu_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr);


/*
 * DIAGONAL MATRIX
 */

template <int NumCtrls, int NumTargs> void cpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr);


/*
 * PAULI TENSOR AND GADGET
 */

template <int NumCtrls, int NumTargs> void cpu_statevector_anyCtrlPauliTensorOrGadget_subA(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> suffixTargsXY, 
    qindex suffixMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
);

template <int NumCtrls> void cpu_statevector_anyCtrlPauliTensorOrGadget_subB(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates,
    qindex suffixMaskXY, qindex bufferMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
);

template <int NumCtrls> void cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, 
    qcomp fac0, qcomp fac1
);


/*
 * DECOHERENCE
 */

void cpu_densmatr_oneQubitDephasing_subA(Qureg qureg, int qubit, qreal prob);
void cpu_densmatr_oneQubitDephasing_subB(Qureg qureg, int qubit, qreal prob);

void cpu_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubitA, int qubitB, qreal prob);
void cpu_densmatr_twoQubitDephasing_subB(Qureg qureg, int qubitA, int qubitB, qreal prob);

void cpu_densmatr_oneQubitDepolarising_subA(Qureg qureg, int ketQubit, qreal prob);
void cpu_densmatr_oneQubitDepolarising_subB(Qureg qureg, int ketQubit, qreal prob);


#endif // CPU_SUBROUTINES_HPP