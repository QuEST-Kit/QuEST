/** @file
 * Internal functions which localize the data needed for simulation.
 * That is, they determine whether performing a simulation requires
 * Querg amplitudes from other distributed nodes and if so, invoke
 * the necessary communication, before finally calling the 
 * embarrassingly parallel subroutines in accelerator.cpp. This is
 * done agnostically of whether amplitudes of the Qureg are being
 * stored in RAM (CPU) or VRAM (GPU).
 */

#ifndef LOCALISER_HPP
#define LOCALISER_HPP

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include <vector>

using std::vector;



/*
 * SWAP
 */


void localiser_statevec_anyCtrlSwap(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2);


/*
 * MATRICES
 */

void localiser_statevec_anyCtrlOneTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);

void localiser_statevec_anyCtrlAnyTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr);

void localiser_statevec_anyCtrlAnyTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr);


/*
 * PAULI TENSORS AND GADGETS
 */

void localiser_statevec_anyCtrlPauliTensor(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str);

void localiser_statevec_anyCtrlPauliGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str, qreal angle);

void localiser_statevec_anyCtrlAnyTargZ(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs);

void localiser_statevec_anyCtrlPhaseGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qreal angle);


#endif // LOCALISER_HPP