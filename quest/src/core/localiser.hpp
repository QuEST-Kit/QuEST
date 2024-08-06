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
#include "quest/include/matrices.h"

#include <vector>

using std::vector;



/*
 * ANY-CTRL ONE-TARG MATRIX
 */

void statevec_anyCtrlOneTargMatrix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
void statevec_anyCtrlOneTargMatrix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr);





// NEW STUFF

void NEW_statevec_anyCtrlOneTargMatrix(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);




#endif // LOCALISER_HPP