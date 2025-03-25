/** @file
 * Internal functions which localize the data needed for simulation.
 * That is, they determine whether performing a simulation requires
 * Querg amplitudes from other distributed nodes and if so, invoke
 * the necessary communication, before finally calling the 
 * embarrassingly parallel subroutines in accelerator.cpp. This is
 * done agnostically of whether amplitudes of the Qureg are being
 * stored in RAM (CPU) or VRAM (GPU).
 * 
 * @author Tyson Jones
 */

#ifndef LOCALISER_HPP
#define LOCALISER_HPP

#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"

#include <vector>

using std::vector;


/*
 * GETTERS
 */

qcomp localiser_statevec_getAmp(Qureg qureg, qindex globalInd);

void localiser_statevec_getAmps(qcomp* outAmps, Qureg qureg, qindex globalStartInd, qindex globalNumAmps);
void localiser_densmatr_getAmps(qcomp** outAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);

void localiser_fullstatediagmatr_getElems(qcomp* outElems, FullStateDiagMatr matr, qindex globalStartInd, qindex globalNumElems);


/*
 * SETTERS
 */

void localiser_statevec_setAmps(qcomp* inAmps, Qureg qureg, qindex globalStartInd, qindex globalNumAmps);
void localiser_densmatr_setAmps(qcomp** inAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);

void localiser_densmatr_setAmpsToPauliStrSum(Qureg qureg, PauliStrSum sum);

void localiser_fullstatediagmatr_setElems(FullStateDiagMatr matr, qindex startInd, qcomp* in, qindex numElems);

void localiser_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in);


/*
 * STATE INITIALISATION
 */

void localiser_statevec_initUniformState(Qureg qureg, qcomp amp);

void localiser_statevec_initDebugState(Qureg qureg);

void localiser_statevec_initClassicalState(Qureg qureg, qindex globalInd);

void localiser_densmatr_initPureState(Qureg qureg, Qureg pure);

void localiser_statevec_initArbitraryPureState(Qureg qureg,  qcomp* amps);
void localiser_densmatr_initArbitraryPureState(Qureg qureg,  qcomp* amps);
void localiser_densmatr_initArbitraryMixedState(Qureg qureg, qcomp** amps);

void localiser_statevec_initUnnormalisedUniformlyRandomPureStateAmps(Qureg qureg);

void localiser_densmatr_initUniformlyRandomPureStateAmps(Qureg qureg);

void localiser_densmatr_initMixtureOfUniformlyRandomPureStates(Qureg qureg, qindex numPureStates);


/*
 * SWAP
 */

void localiser_statevec_anyCtrlSwap(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2);


/*
 * DENSE MATRICES
 */

void localiser_statevec_anyCtrlOneTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr, bool conj);

void localiser_statevec_anyCtrlTwoTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, CompMatr2 matr, bool conj);

void localiser_statevec_anyCtrlAnyTargDenseMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr, bool conj);


/*
 * DIAGONAL MATRICES
 */

void localiser_statevec_anyCtrlOneTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr, bool conj);

void localiser_statevec_anyCtrlTwoTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, DiagMatr2 matr, bool conj);

void localiser_statevec_anyCtrlAnyTargDiagMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, qcomp exponent, bool conj);

void localiser_statevec_allTargDiagMatr(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);
void localiser_densmatr_allTargDiagMatr(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool multiplyOnly);


/*
 * CONVENIENCE GATEWAY TO MATRICES
 */

template <class T>
void localiser_statevec_anyCtrlAnyTargAnyMatr(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, T matr, bool conj);


/*
 * PAULI TENSORS AND GADGETS
 */

void localiser_statevec_anyCtrlPauliTensor(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str, qcomp globalFactor=1);

void localiser_statevec_anyCtrlPauliGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, PauliStr str, qreal phase);

void localiser_statevec_anyCtrlPhaseGadget(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, qreal phase);


/*
 * QUREG COMBINATION
 */

void localiser_statevec_setQuregToSuperposition(qcomp facOut, Qureg outQureg, qcomp fac1, Qureg inQureg1, qcomp fac2, Qureg inQureg2);

void localiser_densmatr_mixQureg(qreal outProb, Qureg out, qreal inProb, Qureg in);


/*
 * DECOHERENCE
 */

void localiser_densmatr_oneQubitDephasing(Qureg qureg, int qubit, qreal prob);
void localiser_densmatr_twoQubitDephasing(Qureg qureg, int qubitA, int qubitB, qreal prob);

void localiser_densmatr_oneQubitDepolarising(Qureg qureg, int qubit, qreal prob);
void localiser_densmatr_twoQubitDepolarising(Qureg qureg, int qubitA, int qubitB, qreal prob);

void localiser_densmatr_oneQubitPauliChannel(Qureg qureg, int qubit, qreal pX, qreal pY, qreal pZ);

void localiser_densmatr_oneQubitDamping(Qureg qureg, int qubit, qreal prob);

void localiser_densmatr_superoperator(Qureg qureg, SuperOp op, vector<int> ketTargs);

void localiser_densmatr_krausMap(Qureg qureg, KrausMap map, vector<int> qubits);


/*
 * PARTIAL TRACE
 */

void localiser_densmatr_partialTrace(Qureg inQureg, Qureg outQureg, vector<int> targs);


/*
 * PROBABILITIES
 */

qreal localiser_statevec_calcTotalProb(Qureg qureg);
qreal localiser_densmatr_calcTotalProb(Qureg qureg);

qreal localiser_statevec_calcProbOfMultiQubitOutcome(Qureg qureg, vector<int> qubits, vector<int> outcomes);
qreal localiser_densmatr_calcProbOfMultiQubitOutcome(Qureg qureg, vector<int> qubits, vector<int> outcomes);

void localiser_statevec_calcProbsOfAllMultiQubitOutcomes(qreal* outProbs, Qureg qureg, vector<int> qubits);
void localiser_densmatr_calcProbsOfAllMultiQubitOutcomes(qreal* outProbs, Qureg qureg, vector<int> qubits);


/*
 * INNER PRODUCTS
 */

qcomp localiser_statevec_calcInnerProduct(Qureg quregA, Qureg quregB);

qcomp localiser_densmatr_calcFidelityWithPureState(Qureg rho, Qureg psi, bool conj);

qreal localiser_densmatr_calcHilbertSchmidtDistance(Qureg quregA, Qureg quregB);


/*
 * EXPECTATION VALUES
 */

qcomp localiser_statevec_calcExpecPauliStr(Qureg qureg, PauliStr str);
qcomp localiser_densmatr_calcExpecPauliStr(Qureg qureg, PauliStr str);

qcomp localiser_statevec_calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum);
qcomp localiser_densmatr_calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum);

qcomp localiser_statevec_calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool useRealPow);
qcomp localiser_densmatr_calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool useRealPow);


/*
 * PROJECTORS 
 */

void localiser_statevec_multiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob);
void localiser_densmatr_multiQubitProjector(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob);


#endif // LOCALISER_HPP