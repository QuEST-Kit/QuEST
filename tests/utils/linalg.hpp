/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilslinalg Linalg
 * @ingroup testutils
 * @brief 
 * Testing utilities which perform linear algebra
 * routines upon reference qvector and qmatrix. 
 * These are slow, serial, un-optimised, defensively-
 * designed routines.
 * @{
 */

#ifndef LINALG_HPP
#define LINALG_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;

int getLog2(qindex);
int getBitAt(qindex num, int ind);
int getNumPermutations(int n, int k);
qindex getPow2(int);
qcomp getExpI(qreal);

qvector getNormalised(qvector);
qvector getDisceteFourierTransform(qvector);

qcomp   getInnerProduct(qvector bra, qvector ket);
qmatrix getOuterProduct(qvector ket, qvector bra);

qvector operator * (const qmatrix&, const qvector&);

bool isDiagonal(qmatrix);
bool isApproxUnitary(qmatrix);

qcomp getTrace(qmatrix m);
qmatrix getTranspose(qmatrix);
qmatrix getConjugateTranspose(qmatrix);
qmatrix getPowerOfDiagonalMatrix(qmatrix diag, qcomp power);
qmatrix getExponentialOfDiagonalMatrix(qmatrix diag);
qmatrix getExponentialOfPauliMatrix(qreal arg, qmatrix pauli);
qmatrix getExponentialOfNormalisedPauliVector(qreal arg, qreal x, qreal y, qreal z);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getKroneckerProduct(qmatrix, qmatrix);
qmatrix getKroneckerProduct(qmatrix, int);
qmatrix getKroneckerProduct(vector<qmatrix> matrices);
qmatrix getProjector(int outcome);
qmatrix getProjector(vector<int> targets, vector<int> outcomes, int numQubits);
qmatrix getPartialTrace(qmatrix matrix, vector<int> targets);
qmatrix getControlledMatrix(qmatrix matrix, int numCtrls);

bool isApproxCPTP(vector<qmatrix>);


#endif // LINALG_HPP

/** @} (end defgroup) */
