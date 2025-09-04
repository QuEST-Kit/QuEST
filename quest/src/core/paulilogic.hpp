/** @file
 * Internal signatures which process Pauli strings
 * and their weighted sums
 * 
 * @author Tyson Jones
 */

#ifndef PAULILOGIC_HPP
#define PAULILOGIC_HPP

#include "quest/include/precision.h"
#include "quest/include/paulis.h"
#include "quest/include/qureg.h"

#include <utility>
#include <vector>
#include <array>

using std::vector;


/*
 * CONSTANTS
 */

static const int MAX_NUM_PAULIS_PER_MASK = sizeof(PAULI_MASK_TYPE) * 8 / 2;
static const int MAX_NUM_PAULIS_PER_STR  = MAX_NUM_PAULIS_PER_MASK * 2;


/*
 * PauliStr
 */

bool paulis_isIdentity(PauliStr str);

bool paulis_containsXOrY(PauliStr str);

int paulis_getPauliAt(PauliStr str, int ind);

int paulis_getIndOfLefmostNonIdentityPauli(PauliStr str);
int paulis_getIndOfLefmostNonIdentityPauli(PauliStr* strings, qindex numStrings);

int paulis_getSignOfPauliStrConj(PauliStr str);

int paulis_getPrefixZSign(Qureg qureg, vector<int> prefixZ);

qcomp paulis_getPrefixPaulisElem(Qureg qureg, vector<int> prefixY, vector<int> prefixZ);

vector<int> paulis_getTargetInds(PauliStr str);

std::array<vector<int>,3> paulis_getSeparateInds(PauliStr str);

qindex paulis_getTargetBitMask(PauliStr str);

PauliStr paulis_getShiftedPauliStr(PauliStr str, int pauliShift);

PauliStr paulis_getKetAndBraPauliStr(PauliStr str, Qureg qureg);

PAULI_MASK_TYPE paulis_getKeyOfSameMixedAmpsGroup(PauliStr str);


// below are not currently used outside of paulilogic.cpp but are natural methods

PauliStr paulis_getTensorProdOfPauliStr(PauliStr left, PauliStr right, int numQubits);

std::pair<qcomp,PauliStr> paulis_getPauliStrProd(PauliStr strA, PauliStr strB);


/*
 * PauliStrSum
 */

bool paulis_containsXOrY(PauliStrSum sum);

int paulis_getIndOfLefmostNonIdentityPauli(PauliStrSum sum);

qindex paulis_getTargetBitMask(PauliStrSum sum);


// below are used exclusively by Trotterisation

qindex paulis_getNumTermsInPauliStrSumProdOfAdjointWithSelf(PauliStrSum in);

void paulis_setPauliStrSumToScaledTensorProdOfConjWithSelf(PauliStrSum out, qreal factor, PauliStrSum in, int numQubits);

void paulis_setPauliStrSumToScaledProdOfAdjointWithSelf(PauliStrSum out, qreal factor, PauliStrSum in);

void paulis_setPauliStrSumToShiftedConj(PauliStrSum out, PauliStrSum in, int numQubits);


#endif // PAULILOGIC_HPP