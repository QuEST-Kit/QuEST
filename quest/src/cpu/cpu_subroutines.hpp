/** @file
 * CPU OpenMP-accelerated signatures of the main backend simulation routines,
 * as mirrored by gpu_subroutines.cpp, and called by accelerator.cpp. 
 * 
 * @author Tyson Jones
 */

#ifndef CPU_SUBROUTINES_HPP
#define CPU_SUBROUTINES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/utilities.hpp"

#include <vector>

using std::vector;


/*
 * GETTERS
 */

qcomp cpu_statevec_getAmp_sub(Qureg qureg, qindex ind);


/*
 * SETTERS
 */

void cpu_densmatr_setAmpsToPauliStrSum_sub(Qureg qureg, PauliStrSum sum);

void cpu_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in);

 // CPU only
void cpu_fullstatediagmatr_setElemsFromMultiDimLists(FullStateDiagMatr out, void* lists, int* numQubitsPerDim, int numDims);
void cpu_fullstatediagmatr_setElemsFromMultiVarFunc(FullStateDiagMatr out, qcomp (*callbackFunc)(qindex*), int* numQubitsPerVar, int numVars, int areSigned);


/*
 * COMMUNICATION BUFFER PACKING
 */

template <int NumQubits> qindex cpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubitInds, vector<int> qubitStates);

qindex cpu_statevec_packPairSummedAmpsIntoBuffer(Qureg qureg, int qubit1, int qubit2, int qubit3, int bit2);


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

template <int NumCtrls> void cpu_statevec_anyCtrlTwoTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, CompMatr2 matr);

template <int NumCtrls, int NumTargs, bool ApplyConj> void cpu_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr);


/*
 * DIAGONAL MATRIX
 */

template <int NumCtrls> void cpu_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr);

template <int NumCtrls> void cpu_statevec_anyCtrlTwoTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, DiagMatr2 matr);

template <int NumCtrls, int NumTargs, bool ApplyConj, bool HasPower> void cpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, qcomp exponent);

template <bool HasPower> void cpu_statevec_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);

template <bool HasPower, bool MultiplyOnly> void cpu_densmatr_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);


/*
 * PAULI TENSOR AND GADGET
 */

template <int NumCtrls, int NumTargs> void cpu_statevector_anyCtrlPauliTensorOrGadget_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> x, vector<int> y, vector<int> z, qcomp ampFac, qcomp pairAmpFac);

template <int NumCtrls> void cpu_statevector_anyCtrlPauliTensorOrGadget_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> x, vector<int> y, vector<int> z, qcomp ampFac, qcomp pairAmpFac, qindex bufferMaskXY);

template <int NumCtrls> void cpu_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qcomp fac0, qcomp fac1);


/*
 * QUREG COMBINATION
 */

void cpu_statevec_setQuregToSuperposition_sub(qcomp facOut, Qureg outQureg, qcomp fac1, Qureg inQureg1, qcomp fac2, Qureg inQureg2);

void cpu_densmatr_mixQureg_subA(qreal outProb, Qureg outQureg, qreal inProb, Qureg inDensMatr);
void cpu_densmatr_mixQureg_subB(qreal outProb, Qureg outQureg, qreal inProb, Qureg inStateVec);
void cpu_densmatr_mixQureg_subC(qreal outProb, Qureg outQureg, qreal inProb);


/*
 * DECOHERENCE
 */

void cpu_densmatr_oneQubitDephasing_subA(Qureg qureg, int qubit, qreal prob);
void cpu_densmatr_oneQubitDephasing_subB(Qureg qureg, int qubit, qreal prob);

void cpu_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubit1, int qubit2, qreal prob);
void cpu_densmatr_twoQubitDephasing_subB(Qureg qureg, int qubit1, int qubit2, qreal prob);

void cpu_densmatr_oneQubitDepolarising_subA(Qureg qureg, int qubit, qreal prob);
void cpu_densmatr_oneQubitDepolarising_subB(Qureg qureg, int qubit, qreal prob);

void cpu_densmatr_twoQubitDepolarising_subA(Qureg qureg, int qubit1, int qubit2, qreal prob);
void cpu_densmatr_twoQubitDepolarising_subB(Qureg qureg, int qubit1, int qubit2, qreal prob);
void cpu_densmatr_twoQubitDepolarising_subC(Qureg qureg, int qubit1, int qubit2, qreal prob);
void cpu_densmatr_twoQubitDepolarising_subD(Qureg qureg, int qubit1, int qubit2, qreal prob);
void cpu_densmatr_twoQubitDepolarising_subE(Qureg qureg, int qubit1, int qubit2, qreal prob);
void cpu_densmatr_twoQubitDepolarising_subF(Qureg qureg, int qubit1, int qubit2, qreal prob);

void cpu_densmatr_oneQubitPauliChannel_subA(Qureg qureg, int qubit, qreal pI, qreal pX, qreal pY, qreal pZ);
void cpu_densmatr_oneQubitPauliChannel_subB(Qureg qureg, int qubit, qreal pI, qreal pX, qreal pY, qreal pZ);

void cpu_densmatr_oneQubitDamping_subA(Qureg qureg, int qubit, qreal prob);
void cpu_densmatr_oneQubitDamping_subB(Qureg qureg, int qubit, qreal prob);
void cpu_densmatr_oneQubitDamping_subC(Qureg qureg, int qubit, qreal prob);
void cpu_densmatr_oneQubitDamping_subD(Qureg qureg, int qubit, qreal prob);


/*
 * PARTIAL TRACE
 */

template <int NumTargs> void cpu_densmatr_partialTrace_sub(Qureg inQureg, Qureg outQureg, vector<int> targs, vector<int> pairTargs);


/*
 * PROBABILITIES
 */

qreal cpu_statevec_calcTotalProb_sub(Qureg qureg);
qreal cpu_densmatr_calcTotalProb_sub(Qureg qureg);

template <int NumQubits> qreal cpu_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes);
template <int NumQubits> qreal cpu_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes);

template <int NumQubits> void cpu_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits);
template <int NumQubits> void cpu_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits);


/*
 * INNER PRODUCTS
 */

qcomp cpu_statevec_calcInnerProduct_sub(Qureg quregA, Qureg quregB);

qreal cpu_densmatr_calcHilbertSchmidtDistance_sub(Qureg quregA, Qureg quregB);

template <bool Conj> qcomp cpu_densmatr_calcFidelityWithPureState_sub(Qureg rho, Qureg psi);


/*
 * EXPECTATION VALUES
 */

qreal cpu_statevec_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs);
qcomp cpu_densmatr_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs);

qcomp cpu_statevec_calcExpecPauliStr_subA(Qureg qureg, vector<int> x, vector<int> y, vector<int> z);
qcomp cpu_statevec_calcExpecPauliStr_subB(Qureg qureg, vector<int> x, vector<int> y, vector<int> z);
qcomp cpu_densmatr_calcExpecPauliStr_sub (Qureg qureg, vector<int> x, vector<int> y, vector<int> z);

template <bool HasPower, bool UseRealPow> qcomp cpu_statevec_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);
template <bool HasPower, bool UseRealPow> qcomp cpu_densmatr_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);


/*
 * PROJECTORS
 */

template <int NumQubits> void cpu_statevec_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob);
template <int NumQubits> void cpu_densmatr_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob);


/*
 * STATE INITIALISATION
 */

void cpu_statevec_initUniformState_sub(Qureg qureg, qcomp amp);

void cpu_statevec_initDebugState_sub(Qureg qureg);

void cpu_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(Qureg qureg);


#endif // CPU_SUBROUTINES_HPP