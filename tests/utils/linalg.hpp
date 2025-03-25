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

int getNumPermutations(int n, int k);
int getLog2(qindex);
int getBitAt(qindex num, int ind);
vector<int> getBits(qindex num, int numBits);
qindex getBitsAt(qindex num, vector<int> inds);
qindex setBitAt(qindex num, int ind, int bit);
qindex setBitsAt(qindex num, vector<int> inds, qindex bits);
qindex getPow2(int);

qreal getSum(vector<qreal> vec);
qcomp getSum(qvector);
qvector getNormalised(qvector);
qvector getDisceteFourierTransform(qvector);
qvector getDisceteFourierTransform(qvector in, vector<int> targs);

qcomp   getInnerProduct(qvector bra, qvector ket);
qmatrix getOuterProduct(qvector ket, qvector bra);

qvector operator * (const qmatrix&, const qvector&);

bool isDiagonal(qmatrix);
bool isApproxUnitary(qmatrix);

qcomp getTrace(qmatrix);
qmatrix getTranspose(qmatrix);
qmatrix getConjugate(qmatrix);
qmatrix getConjugateTranspose(qmatrix);
qmatrix getPowerOfDiagonalMatrix(qmatrix diag, qcomp power);
qmatrix getExponentialOfDiagonalMatrix(qmatrix);
qmatrix getExponentialOfPauliMatrix(qreal arg, qmatrix pauli);
qmatrix getExponentialOfNormalisedPauliVector(qreal arg, qreal x, qreal y, qreal z);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getOrthonormalisedRows(qmatrix);
qmatrix getKroneckerProduct(qmatrix, qmatrix);
qmatrix getKroneckerProduct(qmatrix, int count);
qmatrix getKroneckerProduct(vector<qmatrix>);
qmatrix getProjector(int outcome);
qmatrix getProjector(vector<int> targets, vector<int> outcomes, int numQubits);
qmatrix getPartialTrace(qmatrix matrix, vector<int> targets);
qmatrix getControlledMatrix(qmatrix matrix, int numCtrls);
qmatrix getMixture(vector<qvector> statevecs, vector<qreal> probs);
qmatrix getMixture(vector<qmatrix> densmatrs, vector<qreal> probs);
qmatrix getSuperOperator(vector<qmatrix>);

bool isApproxCPTP(vector<qmatrix>);


#endif // LINALG_HPP

/** @} (end defgroup) */
