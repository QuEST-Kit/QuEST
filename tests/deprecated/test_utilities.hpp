/** @file
 * Utility signatures used by the ported v3 tests of QuEST's deprecated v3 API.
 *
 * @author Tyson Jones
 * @author Oliver Thomson Brown (ported to Catch2 v3)
 * @author Ali Rezaei (tested porting to QuEST v4)
 */
 
#ifndef QUEST_TEST_UTILS_H
#define QUEST_TEST_UTILS_H

#include "quest/include/quest.h"

#include <catch2/generators/catch_generators.hpp>

#include <string>
#include <vector>

using std::string;
using std::vector;

/** The default number of qubits in the registers created for unit testing 
 * (both statevectors and density matrices). Creation of non-NUM_QUBITS sized 
 * Quregs should be justified in a comment. 
 * Note that the smaller this number is, the fewer nodes can be employed in 
 * distribution testing, since each node must contain at least one amplitude.
 * Furthermore, the larger this number is, the greater the deviation of correct 
 * results from their expected value, due to numerical error; this is especially 
 * apparent for density matrices.
 */
#define NUM_QUBITS 5


// replace REAL_EPS macro with constant
#undef REAL_EPS
#if FLOAT_PRECISION == 1
    constexpr qreal REAL_EPS = 1E-1;
#elif FLOAT_PRECISION == 2
    constexpr qreal REAL_EPS = 1E-8;
#elif FLOAT_PRECISION == 4
    constexpr qreal REAL_EPS = 1E-10;
#endif


#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

/** A complex square matrix. 
 * Should be initialised with getZeroMatrix().
 * These have all the natural linear-algebra operator overloads, including 
 * left-multiplication onto a vector.
 *
 * This data-structure is not partitioned between nodes in distributed mode.
 * That is, every node has a complete copy, allowing for safe comparisons.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
typedef vector<vector<qcomp>> QMatrix;

/** A complex vector, which can be zero-initialised with QVector(numAmps).
 * These have all the natural linear-algebra operator overloads.
 *
 * This data-structure is not partitioned between nodes in distributed mode.
 * That is, every node has a complete copy, allowing for safe comparisons.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
typedef vector<qcomp> QVector;

/** Seed the C and C++ RNGs using hardware CSPRNG
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void setRandomTestStateSeeds();

#ifdef absReal
#undef absReal
#endif
qreal absReal(qreal x);
qreal absComp(qcomp x);

/** Asserts the given statevector qureg and reference agree, and are properly initialised in the debug state.
 * Assertion uses the DEMAND() macro, calling Catch2's FAIL() if unsatisfied, so does not contribute
 * toward unit test statistics. This should be called within every PREPARE_TEST macro, to ensure that
 * the test states themselves are initially correct, and do not accidentally agree by (e.g.) being all-zero.
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void assertQuregAndRefInDebugState(Qureg qureg, QVector ref);

/** Asserts the given density qureg and reference agree, and are properly initialised in the debug state.
 * Assertion uses the DEMAND() macro, calling Catch2's FAIL() if unsatisfied, so does not contribute
 * toward unit test statistics. This should be called within every PREPARE_TEST macro, to ensure that
 * the test states themselves are initially correct, and do not accidentally agree by (e.g.) being all-zero.
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void assertQuregAndRefInDebugState(Qureg qureg, QMatrix ref);

/* (Excluded from Doxygen doc)
 *
 * Define QVector and QMatrix operator overloads.
 * Note that QMatrix overloads don't simply use QVector 
 * overloads, since the complex vector dot product involves 
 * conjugation, which doesn't occur in complex matrix multiplication.
 * Note too we also avoid defining operators in terms of other operators
 * (e.g. minus is plus(negative times)) since compiler optimisations 
 * may change the order of operations and confuse the overloads invoked.
 * Definition of division using multiplication can furthermore 
 * heighten numerical errors.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QVector operator + (const QVector& v1, const QVector& v2);
QVector operator - (const QVector& v1, const QVector& v2);
QVector operator * (const qcomp& a, const QVector& v);
QVector operator * (const QVector& v, const qcomp& a);
QVector operator / (const QVector& v, const qcomp& a);
qcomp operator * (const QVector &v1, const QVector& v2);
void operator += (QVector& v1, const QVector& v2);
void operator -= (QVector& v1, const QVector& v2);
void operator *= (QVector& v1, const qcomp& a);
void operator /= (QVector& v1, const qcomp& a);
QMatrix operator + (const QMatrix& m1, const QMatrix& m2);
QMatrix operator - (const QMatrix& m1, const QMatrix& m2);
QMatrix operator * (const qcomp& a, const QMatrix& m);
QMatrix operator * (const QMatrix& m, const qcomp& a);
QMatrix operator / (const QMatrix& m, const qcomp& a);
QMatrix operator * (const QMatrix& m1, const QMatrix& m2);
void operator += (QMatrix& m1, const QMatrix& m2);
void operator -= (QMatrix& m1, const QMatrix& m2);
void operator *= (QMatrix& m1, const qreal& a);
void operator /= (QMatrix& m1, const qreal& a);
void operator *= (QMatrix& m1, const QMatrix& m2);
QVector operator * (const QMatrix& m, const QVector& v);

/** Returns an equal-size copy of the given state-vector \p qureg.
 * In GPU mode, this function involves a copy of \p qureg from GPU memory to RAM.
 * In distributed mode, this involves an all-to-all broadcast of \p qureg.
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QVector toQVector(Qureg qureg);

/** Returns a vector with the given diagonal's elements.
 * In distributed mode, this involves an all-to-all broadcast of \p op.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QVector toQVector(DiagMatr op);

/** Returns a vector with the given diagonal's elements.
 * In distributed mode, this involves an all-to-all broadcast of \p op.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QVector toQVector(FullStateDiagMatr op);

/** Returns an equal-size copy of the given density matrix \p qureg.
 * In GPU mode, this function involves a copy of \p qureg from GPU memory to RAM.
 * In distributed mode, this involves an all-to-all broadcast of \p qureg.
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(Qureg qureg);

/** Returns a copy of the given 2-by-2 matrix.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(CompMatr1 src);

/** Returns a copy of the given 4-by-4 matrix.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(CompMatr2 src);

/** Returns a copy of the given matrix
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(CompMatr src);

/** Returns a 2^\p N-by-2^\p N Hermitian matrix form of the specified 
 * weighted sum of Pauli products
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(qreal* coeffs, pauliOpType* paulis, int numQubits, int numTerms);

/** Returns a 2^\p N-by-2^\p N Hermitian matrix form of the PauliHamil
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(PauliHamil hamil);

/** Returns a 2^\p N-by-2^\p N Hermitian Z-basis
 * matrix of the given complex-weighted sum of Pauli 
 * strings, where N is the number of non-Identity
 * operators.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(PauliStrSum sum);

/** Returns a dense matrix equivalent to the given diagonal
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(DiagMatr matr);

/** Returns a dense matrix equivalent to the given diagonal
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix toQMatrix(FullStateDiagMatr matr);

/** Returns a \p ComplexMatrix2 copy of QMatix \p qm.
 * Demands that \p qm is a 2-by-2 matrix.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
ComplexMatrix2 toComplexMatrix2(QMatrix qm);

/** Returns a \p ComplexMatrix4 copy of QMatix \p qm.
 * Demands that \p qm is a 4-by-4 matrix.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
ComplexMatrix4 toComplexMatrix4(QMatrix qm);

/** Populates the ComplexMatrixN with the contents of a QMatrix. In
 * GPU-mode, this will then sync the elements ot the matrix's
 * persistent GPU memory
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void toComplexMatrixN(QMatrix qm, ComplexMatrixN cm);

/** Returns a diagonal complex matrix formed by the given vector 
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones 
 */
QMatrix toDiagonalQMatrix(QVector vec);

/** Initialises the state-vector \p qureg to have the same amplitudes as \p vec.
 * Demands \p qureg is a state-vector of an equal size to \p vec.
 * In GPU mode, this function involves a copy from RAM to GPU memory.
 * This function has no communication cost in distributed mode.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void toQureg(Qureg qureg, QVector vec);

/** Initialises the density matrix \p qureg to have the same amplitudes as \p mat.
 * Demands \p qureg is a density matrix of equal dimensions to \p mat.
 * In GPU mode, this function involves a copy from RAM to GPU memory.
 * This function has no communication cost in distributed mode.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void toQureg(Qureg qureg, QMatrix mat);

/** Returns b (otimes) a. If b and a are state-vectors, the resulting kronecker 
 * product is the seperable state formed by joining the qubits in the state-vectors, 
 * producing |b>|a> (a is least significant)
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QVector getKroneckerProduct(QVector b, QVector a);

/** Returns a dim-by-dim square complex matrix, initialised to all zeroes.
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix getZeroMatrix(size_t dim);

/** Returns a dim-by-dim identity matrix
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getIdentityMatrix(size_t dim);

/** Returns the matrix exponential of a diagonal, square, complex matrix.
 * This method explicitly checks that the passed matrix \p a is diagonal.
 * 
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getExponentialOfDiagonalMatrix(QMatrix a);

/** Returns the matrix exponential of a kronecker product of pauli matrices 
 * (or of any involutory matrices), with exponent factor (-i \p angle / 2).
 * This method will not explicitly check that the passed matrix \p a is 
 * kronecker product of involutory matrices, but will otherwise return an 
 * incorrect exponential.
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix getExponentialOfPauliMatrix(qreal angle, QMatrix a);

/** Returns the kronecker product of \p a and \p b, where \p a and \p b are 
 * square but possibly differently-sized complex matrices.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getKroneckerProduct(QMatrix a, QMatrix b);

/** Returns the 2^\p numQb-by-2^\p numQb unitary matrix which swaps qubits 
 * \p qb1 and \p qb2; the SWAP gate of not-necessarily-adjacent qubits.
 * If \p qb1 == \p qb2, returns the identity matrix.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getSwapMatrix(int qb1, int qb2, int numQb);

/** Takes a 2^\p numTargs-by-2^\p numTargs matrix \p op and a returns a 
 * 2^\p numQubits-by-2^\p numQubits matrix where \p op is controlled on the given 
 * \p ctrls qubits. The union of {\p ctrls} and {\p targs} must be unique (though
 * this is not explicitly checked), and every element must be >= 0 (not checked).
 * The passed {\p ctrls} and {\p targs} arrays are unmodified.
 *
 * This funciton works by first swapping {\p targs} and {\p ctrls} (via swap unitaries) 
 * to be strictly increasing {0,1,...}, building controlled(\p op), tensoring it to 
 * the full Hilbert space, and then 'unswapping'. The returned matrix has form:
 * swap1 ... swapN . controlled(\p op) . swapN ... swap1
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
QMatrix getFullOperatorMatrix(int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op, int numQubits);

/** Returns the matrix |\p ket><\p bra|, with ith-jth element \p ket(i) conj(\p bra(j)), since
 * |\p ket><\p bra| = sum_i a_i|i> sum_j b_j* <j| = sum_{ij} a_i b_j* |i><j|.
 * The dimensions of bra and ket must agree, and the returned square complex matrix 
 * has dimensions size(bra) x size(bra).
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getKetBra(QVector ket, QVector bra);

/** Returns the conjugate transpose of the complex square matrix \p a
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getConjugateTranspose(QMatrix a);

/** Returns a random integer between \p min (inclusive) and \p max (exclusive),
 * from the uniform distribution.
 * Demands that \p max > \p min.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
int getRandomInt(int min, int max);

/** Returns a random real between \p min (inclusive) and \p max (exclusive),
 * from the uniform distribution.
 * Demands that \p max > \p min.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
qreal getRandomReal(qreal min, qreal max);

/** Returns a random complex number within the square closing (-1-i) and (1+i),
 * from a distribution uniformly randomising the individual real and imaginary 
 * components in their domains.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
qcomp getRandomComplex();

/** Returns a \p dim-length vector with random complex amplitudes in the 
 * square joining {-1-i, 1+i}, of an undisclosed distribution. The resulting 
 * vector is NOT L2-normalised.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QVector getRandomQVector(int dim);

/** Returns a \p dim-by-\p dim complex matrix, where the real and imaginary value of 
 * each element are independently random, under the standard normal distribution 
 * (mean 0, standard deviation 1).
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getRandomQMatrix(int dim);

/** Returns a uniformly random (under Haar) 2^\p numQb-by-2^\p numQb unitary matrix.
 * This function works by first generating a complex matrix where 
 * each element is independently random; the real and imaginary component thereof are 
 * independent standard normally-distributed (mean 0, standard-dev 1).
 * Then, the matrix is orthonormalised via the Gram Schmidt algorithm. 
 * The resulting unitary matrix MAY be uniformly distributed under the Haar 
 * measure, but we make no assurance. 
 * This routine may return an identity matrix if it was unable to sufficiently 
 * precisely produce a unitary of the given size.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getRandomUnitary(int numQb);

/** Returns a random \p numQb-length L2-normalised state-vector from an 
 * undisclosed distribution. This function works by randomly generating each 
 * complex amplitude, then L2-normalising.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QVector getRandomStateVector(int numQb);

/** Returns a random \p numQb-by-\p numQb density matrix, from an undisclosed 
 * distribution, in a very mixed state. This function works by generating 
 * 2^\p numQb random pure states, and mixing them with random probabilities.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getRandomDensityMatrix(int numQb);

/** Returns a random \p numQb-by-\p numQb density matrix, from an undisclosed 
 * distribution, which is pure (corresponds to a random state-vector)
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QMatrix getRandomPureDensityMatrix(int numQb);

/** Returns a density matrix initialised into the given pure state 
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones 
 */
QMatrix getPureDensityMatrix(QVector state);

/** Returns the diagonal vector of the given matrix
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QVector getMatrixDiagonal(QMatrix matr);

/** Returns a random Kraus map of #`numOps` 2^\p numQb-by-2^\p numQb operators, 
 * from an undisclosed distribution.
 * Note this method is very simple and cannot generate all possible Kraus maps. 
 * It works by generating \p numOps random unitary matrices, and randomly 
 * re-normalising them, such that the sum of ops[j]^dagger ops[j] = 1
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
vector<QMatrix> getRandomKrausMap(int numQb, int numOps);

/** Returns a list of random real scalars, each in [0, 1], which sum to unity. 
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
vector<qreal> getRandomProbabilities(int numProbs);

/** Returns a list of random orthonormal complex vectors, from an undisclosed 
 * distribution. 
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
vector<QVector> getRandomOrthonormalVectors(int numQb, int numStates);

/** Returns a mixed density matrix formed from mixing the given pure states, 
 * which are assumed normalised, but not necessarily orthogonal.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones 
 */
QMatrix getMixedDensityMatrix(vector<qreal> probs, vector<QVector> states);

/** Returns an L2-normalised copy of \p vec, using Kahan summation for improved accuracy.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
QVector getNormalised(QVector vec);

/** Returns the discrete fourier transform of vector in 
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones 
 */
QVector getDFT(QVector in);

/** Returns the discrete fourier transform of a sub-partition of the vector in.
 * 
 * @ingroup deprecatedutils
 * @author Tyson Jones 
 */
QVector getDFT(QVector in, int* targs, int numTargs);

/** Returns the integer value of the targeted sub-register for the given 
 * full state index \p ind. 
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones 
 */
long long int getValueOfTargets(long long int ind, int* targs, int numTargs);

/** Modifies \p dest by overwriting its submatrix (from top-left corner 
 * (\p r, \p c) to bottom-right corner (\p r + \p dest.size(), \p c + \p dest.size()) 
 * with the complete elements of sub.
 * This demands that dest.size() >= sub.size() + max(r,c).
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void setSubMatrix(QMatrix &dest, QMatrix sub, size_t r, size_t c);

 /** Modifies the density matrix \p state to be the result of applying the multi-target operator 
  * matrix \p op, with the specified control and target qubits (in \p ctrls and \p targs 
  * respectively). This updates \p state under
  * \f[
            \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
  * \f]
  * even if \p op is not unitary (which is useful for applying Kraus operators).
  * 
  * \p op must be a 2^\p numTargs-by-2^\p numTargs matrix. Furthermore, every element of \p targs 
  * must not appear in \p ctrls (and vice-versa), though this is not explicitly checked.
  * Elements of \p targs and \p ctrls should be unique.
  *
  * This function works by computing getFullOperatorMatrix() from the given 
  * arguments, and left-multipling it to \p state, then right-multiplying its 
  * conjugate transpose onto the result.
  *
  * @ingroup deprecatedutils 
  * @author Tyson Jones
  */
void applyReferenceOp(QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);

/** Modifies the density matrix \p state to be the result of applying the two-target operator 
 * matrix \p op, with the specified control qubits (in \p ctrls). 
 * This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
 * \f]
 * even if \p op is not unitary (which is useful for applying Kraus operators).
 * 
 * \p op must be a 4-by-4 matrix. Both \p targ1 and \p targ2 must not appear in \p ctrls, 
 * though this is not explicitly checked. Elements of \p ctrls, and \p targ1 and \p targ2,
 * should be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multipling it to \p state, then right-multiplying its 
 * conjugate transpose onto the result.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QMatrix &state, int* ctrls, int numCtrls, int targ1, int targ2, QMatrix op);

/** Modifies the density matrix \p state to be the result of applying the single-target 
 * operator matrix \p op, with the specified control qubits (in \p ctrls). 
 * This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
 * \f]
 * even if \p op is not unitary (which is useful for applying Kraus operators).
 * 
 * \p op must be a 2-by-2 matrix. \p target must not appear in \p ctrls, 
 * though this is not explicitly checked.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multipling it to \p state, then right-multiplying its 
 * conjugate transpose onto the result.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QMatrix &state, int* ctrls, int numCtrls, int target, QMatrix op);

/** Modifies the density matrix \p state to be the result of applying the multi-target 
 * operator matrix \p op, with no control qubits.
 * This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
 * \f]
 * even if \p op is not unitary (which is useful for applying Kraus operators).
 * 
 * \p op must be a 2^\p numTargs-by-2^\p numTargs matrix. 
 * Every element in \p targs should be unique, though this is not explicitly checked.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multipling it to \p state, then right-multiplying its 
 * conjugate transpose onto the result.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QMatrix &state, int *targs, int numTargs, QMatrix op);

/** Modifies the density matrix \p state to be the result of applying the single-control 
 * single-target operator matrix \p op.
 * This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
 * \f]
 * even if \p op is not unitary (which is useful for applying Kraus operators).
 * 
 * \p op must be a 2-by-2 matrix, and \p ctrl and \p targ should be different.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multipling it to \p state, then right-multiplying its 
 * conjugate transpose onto the result.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QMatrix &state, int ctrl, int targ, QMatrix op);

/** Modifies the density matrix \p state to be the result of applying the multi-target 
 * operator matrix \p op, with a single control qubit \p ctrl.
 * This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
 * \f]
 * even if \p op is not unitary (which is useful for applying Kraus operators).
 * 
 * \p op must be a 2^\p numTargs-by-2^\p numTargs matrix, and \p ctrl must not 
 * appear in \p targs (though this is not explicitly checked).
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multipling it to \p state, then right-multiplying its 
 * conjugate transpose onto the result.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QMatrix &state, int ctrl, int* targs, int numTargs, QMatrix op);

/** Modifies the density matrix \p state to be the result of applying the two-target 
 * operator matrix \p op, with a single control qubit \p ctrl.
 * This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
 * \f]
 * even if \p op is not unitary (which is useful for applying Kraus operators).
 * 
 * \p op must be a 4-by-4 matrix, and \p ctrl, \p targ1 and \p targ2 must be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multipling it to \p state, then right-multiplying its 
 * conjugate transpose onto the result.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QMatrix &state, int ctrl, int targ1, int targ2, QMatrix op);

/** Modifies the density matrix \p state to be the result of applying the single-target 
 * operator matrix \p op, with no control qubit.
 * This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state} \, \text{op}^\dagger
 * \f]
 * even if \p op is not unitary (which is useful for applying Kraus operators).
 * 
 * \p op must be a 2-by-2 matrix.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multipling it to \p state, then right-multiplying its 
 * conjugate transpose onto the result.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QMatrix &state, int targ, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the multi-target operator 
 * matrix \p op, with the specified control and target qubits (in \p ctrls and \p targs 
 * respectively). This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 2^\p numTargs-by-2^\p numTargs matrix. Furthermore, every element of \p targs 
 * must not appear in \p ctrls (and vice-versa), though this is not explicitly checked.
 * Elements of \p targs and \p ctrls should be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the two-target operator 
 * matrix \p op, with the specified control qubits (in \p ctrls). This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 4-by-4 matrix. Furthermore, \p ctrls, \p targ1 and \p targ2 should 
 * all be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int* ctrls, int numCtrls, int targ1, int targ2, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the single-target operator 
 * matrix \p op, with the specified control qubits (in \p ctrls). This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 2-by-2 matrix. Furthermore, elements in \p ctrls and \p target should 
 * all be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int* ctrls, int numCtrls, int target, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the multi-target operator 
 * matrix \p op, with no contorl qubits. This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 2^\p numTargs-by-2^\p numTargs matrix. Furthermore, elements in \p targs should be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int *targs, int numTargs, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the single-target operator 
 * matrix \p op, with a single control qubit (\p ctrl). This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 2-by-2 matrix. Furthermore, \p ctrl and \p targ must be different.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int ctrl, int targ, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the multi-target operator 
 * matrix \p op, with a single control qubit (\p ctrl) This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 2^\p numTargs-by-2^\p numTargs matrix. Furthermore, elements 
 * in \p targs and \p ctrl should be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int ctrl, int* targs, int numTargs, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the two-target operator 
 * matrix \p op, with a single control qubit (\p ctrl). This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 4-by-4 matrix. Furthermore, \p ctrl, \p targ1 and \p targ2 should 
 * all be unique.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int ctrl, int targ1, int targ2, QMatrix op);

/** Modifies the state-vector \p state to be the result of applying the single-target operator 
 * matrix \p op, with no contorl qubits. This updates \p state under
 * \f[
           \text{state} \to \text{op} \, \text{state}
 * \f]
 * even if \p op is not unitary.
 * 
 * \p op must be a 2-by-2 matrix.
 *
 * This function works by computing getFullOperatorMatrix() from the given 
 * arguments, and left-multiplying it onto \p state.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceOp(QVector &state, int targ, QMatrix op);

/** Modifies the state-vector \p state to be the result of left-multiplying the multi-target operator 
 * matrix \p op, with the specified control and target qubits (in \p ctrls and \p targs 
 * respectively). This is an alias of applyReferenceOp(), since operators are always 
 * left-multiplied as matrices onto state-vectors.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceMatrix(QVector &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);

/** Modifies the state-vector \p state to be the result of left-multiplying the multi-target operator 
 * matrix \p op, with the specified target qubits (assuming no control qubits). T
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceMatrix(QVector &state, int *targs, int numTargs, QMatrix op);

/** Modifies the density matrix \p state to be the result of left-multiplying the multi-target operator 
 * matrix \p op, with the specified control and target qubits (in \p ctrls and \p targs 
 * respectively). Here, \p op is treated like a simple matrix and is hence left-multiplied 
 * onto the state once.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceMatrix(QMatrix &state, int* ctrls, int numCtrls, int *targs, int numTargs, QMatrix op);

/** Modifies the density matrix \p state to be the result of left-multiplying the multi-target operator 
 * matrix \p op, with the target qubits (assuming no control qubits). 
 * Here, \p op is treated like a simple matrix and is hence left-multiplied onto the state once.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void applyReferenceMatrix(QMatrix &state, int *targs, int numTargs, QMatrix op); 

/** Performs a hardware-agnostic comparison of the given quregs, checking 
 * whether the difference between the real and imaginary components of every amplitude
 * is smaller than the QuEST_PREC-specific REAL_EPS (defined in QuEST_precision) precision.
 * This function demands that \p qureg1 and \p qureg2 are of the same type 
 * (i.e. both state-vectors or both density matrices), and of an equal number 
 * of qubits.
 *
 * In GPU mode, this function involves a GPU to CPU memory copy overhead.
 * In distributed mode, it involves a all-to-all single-int broadcast.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(Qureg qureg1, Qureg qureg2);

/** Performs a hardware-agnostic comparison of state-vector \p qureg to \p vec, checking 
 * whether the difference between the real and imaginary components of every amplitude
 * is smaller than the QuEST_PREC-specific REAL_EPS (defined in QuEST_precision) precision.
 * This function demands \p qureg is a state-vector, and that \p qureg and 
 * \p vec have the same number of amplitudes.
 *
 * In GPU mode, this function involves a GPU to CPU memory copy overhead.
 * In distributed mode, it involves a all-to-all single-int broadcast.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(Qureg qureg, QVector vec);

/** Performs a hardware-agnostic comparison of density-matrix \p qureg to \p matr, checking 
 * whether the difference between the real and imaginary components of every amplitude
 * is smaller than the QuEST_PREC-specific REAL_EPS (defined in QuEST_precision) precision.
 * This function demands \p qureg is a density matrix, and that \p qureg and \p matr have
 * equal dimensions.
 *
 * In GPU mode, this function involves a GPU to CPU memory copy overhead.
 * In distributed mode, it involves a all-to-all single-int broadcast.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(Qureg qureg, QMatrix matr);

/** Performs a hardware-agnostic comparison of the given quregs, checking 
 * whether the difference between the real and imaginary components of every amplitude
 * is smaller than \p precision.
 * This function demands that \p qureg1 and \p qureg2 are of the same type 
 * (i.e. both state-vectors or both density matrices), and of an equal number 
 * of qubits.
 *
 * In GPU mode, this function involves a GPU to CPU memory copy overhead.
 * In distributed mode, it involves a all-to-all single-int broadcast.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(Qureg qureg1, Qureg qureg2, qreal precision);

/** Performs a hardware-agnostic comparison of state-vector \p qureg to \p vec, checking 
 * whether the difference between the real and imaginary components of every amplitude
 * is smaller than \p precision.
 * This function demands \p qureg is a state-vector, and that \p qureg and 
 * \p vec have the same number of amplitudes.
 *
 * In GPU mode, this function involves a GPU to CPU memory copy overhead.
 * In distributed mode, it involves a all-to-all single-int broadcast.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(Qureg qureg, QVector vec, qreal precision);

/** Performs a hardware-agnostic comparison of density-matrix \p qureg to \p matr, checking 
 * whether the difference between the real and imaginary components of every amplitude
 * is smaller than \p precision.
 * This function demands \p qureg is a density matrix, and that \p qureg and \p matr have
 * equal dimensions.
 *
 * In GPU mode, this function involves a GPU to CPU memory copy overhead.
 * In distributed mode, it involves a all-to-all single-int broadcast.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(Qureg qureg, QMatrix matr, qreal precision);

/** Returns true if the absolute value of the difference between every amplitude in 
 * vectors \p a and \p b is less than \p REAL_EPS.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(QVector a, QVector b);

/** Returns true if the absolute value of the difference between every amplitude in 
 * matrices \p a and \p b is less than \p REAL_EPS.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(QMatrix a, QMatrix b);

/** Returns true if the absolute value of the difference between every element in 
 * \p vec and those implied by \p reals and \p imags, is less than \p REAL_EPS.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(QVector vec, qreal* reals, qreal* imags);

/** Returns true if the absolute value of the difference between every element in 
 * \p vec (which must be strictly real) and those implied by \p reals, is less 
 * than \p REAL_EPS.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
bool areEqual(QVector vec, qreal* reals);

/** Returns the unit-norm complex number exp(i*\p phase). This function uses the 
 * Euler formula, and avoids problems with calling exp(__complex__) in a platform 
 * agnostic way 
 */
qcomp expI(qreal phase);

/** Returns log2 of numbers which must be gauranteed to be 2^n 
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
unsigned int calcLog2(long unsigned int res);

/** Populates the \p coeffs array with random qreals in (-5, 5), and 
 * populates \p codes with random Pauli codes
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void setRandomPauliSum(qreal* coeffs, pauliOpType* codes, int numQubits, int numTerms);

/** Populates \p hamil with random coefficients and pauli codes
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void setRandomPauliSum(PauliHamil hamil, int numQubits);

/** Populates \p hamil with random coefficients and a random amount number of 
 * PAULI_I and PAULI_Z operators.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void setRandomDiagPauliHamil(PauliHamil hamil, int numQubits);

/** Populates \p targs with a random selection of \p numTargs elements from [0,\p numQb-1].
 * List \p targs does not need to be initialised and its elements are overwritten.
 * 
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void setRandomTargets(int* targs, int numTargs, int numQb);

/** Populates \p targs with a random selection of elements from [0,\p numQb-1].
 * List \p targs does not need to be initialised and its elements are overwritten.
 * 
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
void setRandomTargets(vector<int> &targs, int numQb);

/** Returns the two's complement signed encoding of the unsigned number decimal, 
 * which must be a number between 0 and 2^numBits (exclusive). The returned number 
 * lies in [-2^(numBits-1), 2^(numBits-1)-1]
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
long long int getTwosComplement(long long int decimal, int numBits);

/** Return the unsigned value of a number, made of `#numBits` bits, which under 
 * two's complement, encodes the signed number twosComp. The returned number 
 * lies in [0, 2^(numBits)-1]
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
long long int getUnsigned(long long int twosComp, int numBits);

// /** Modifies the given diagonal matrix such that the diagonal elements which 
//  * correspond to the coordinates in overrideInds are replaced with exp(i phase), as
//  * prescribed by overridePhases. This function assumes that the given registers 
//  * are contiguous, are in order of increasing significance, and that the matrix 
//  * is proportionately sized and structured to act on the space of all registers 
//  * combined. Overrides can be repeated, and only the first encountered for a given 
//  * index will be effected (much like applyMultiVarPhaseFuncOverrides()).
//  *
//  * @ingroup deprecatedutils
//  * @author Tyson Jones
//  */
// void setDiagMatrixOverrides(QMatrix &matr, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, long long int* overrideInds, qreal* overridePhases, int numOverrides);

/** Modifies outFn to be a filename of format prefix_NUM.txt where NUM 
 * is a new unique integer so far. This is useful for getting unique filenames for 
 * independent test cases of functions requiring reading/writing to file, to 
 * avoid IO locks (especially common in distributed mode).
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void setUniqueFilename(char* outFn, int maxlen, char* prefix);

/** Writes contents to the file with filename fn, which is created and/or overwritten.
 * In distributed mode, the master node writes while the other nodes wait until complete. 
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void writeToFileSynch(char* fn, const string& contents);

/** Deletes all files with filename starting with prefix. In distributed mode, the 
 * master node deletes while the other nodes wait until complete.
 *
 * @ingroup deprecatedutils
 * @author Tyson Jones
 */
void deleteFilesWithPrefixSynch(char* prefix);

// makes below signatures more concise
template<class T> using CatchGen = Catch::Generators::GeneratorWrapper<T>;

 /** Returns a Catch2 generator of every length-\p sublen sublist of length-\p len 
  * \p list, in increasing lexographic order. This generates every fixed-length 
  * combination of the given list and every permutation of each. 
  & If the sublist length is the full list length, this generator produces every 
  * permutation correctly. Note that the sublist must not be modified, else further 
  & generation may break (QuEST's internal functions will indeed modify but restore 
  * the qubit index lists given to them, which is ok). 
  * Assumes \p list contains no duplicates, otherwise the generated sublists may be 
  * duplicated. 
  *
  * This function can be used like
  *
  *     int list[4] = {1,2,3,4};
  *     int sublen = 2;
  *     int* sublist = GENERATE_COPY( sublists(list, 4, sublen) );
  *
  * to generate {1,2}, {1,3}, {1,4}, {2,1}, {2,3}, {2,4}, {3,1}, {3,2}, {3, 4},
  * {4,1}, {4,2}, {4, 3}.
  *
  * @ingroup deprecatedutils 
  * @author Tyson Jones
  */
CatchGen<int*> sublists(int* list, int len, int sublen);

/** Returns a Catch2 generator of every length-\p sublen sublist of the elements
 * generated by \p gen, which exclude all elements in \p exclude, in increasing lexographic order. 
 * This generates every fixed-length combination of \p gen's elements the nominated exclusions,
 * and every permutation of each. 
 *
 * There is on need for the elements of \p exclude to actually appear in those of \p gen.
 * \p sublen must less than or equal to the number of elements in \p gen, after 
 * the nominated exclusions.
 *
 * Note that the sublist must not be modified, else further generation may break (QuEST's 
 * internal functions will indeed modify but restore the qubit index lists given 
 * to them, which is ok). Assumes \p list contains no duplicates, 
 * otherwise the generated sublists may be duplicated. 
 *
 * This function can be used like
 *
 *     int sublen = 2;
 *     int exclude[2] = {3,4};
 *     int* sublist = GENERATE_COPY( sublists(range(1,6), sublen, exclude, 2) );
 *
 * to generate {1,2}, {1,5}, {2,1}, {2,5}, {5,1}, {5,2}
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
CatchGen<int*> sublists(CatchGen<int>&& gen, int numSamps, const int* exclude, int numExclude);

/** Returns a Catch2 generator of every length-\p sublen sublist of the elements
 * generated by \p gen which exclude element \p excluded, in increasing lexographic order. 
 * This generates every fixed-length combination of \p gen's elements the nominated exclusions,
 * and every permutation of each. 
 *
 * \p sublen must less than or equal to the number of elements in \p gen, after 
 * the nominated exclusion. There is no need for \p excluded to actually appear 
 * in the elements of \p gen.
 *
 * Note that the sublist must not be modified, else further generation may break (QuEST's 
 * internal functions will indeed modify but restore the qubit index lists given 
 * to them, which is ok). Assumes \p list contains no duplicates, 
 * otherwise the generated sublists may be duplicated. 
 *
 * This function can be used like
 *
 *     int sublen = 2;
 *     int excluded = 1;
 *     int* sublist = GENERATE_COPY( sublists(range(1,4), sublen, excluded) );
 *
 * to generate {2,3}, {3,2}.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
CatchGen<int*> sublists(CatchGen<int>&& gen, int numSamps, int excluded);

/** Returns a Catch2 generator of every length-\p sublen sublist of the elements 
 * generated by \p gen, in increasing lexographic order. This generates every fixed-length 
 * combination of \p gen's elements, and every permutation of each. 
 * Note that the produced sublist must not be modified, else further 
 * generation may break (QuEST's internal functions will indeed modify but restore 
 * the qubit index lists given to them, which is ok). 
 * Assumes \p list contains no duplicates, otherwise the generated sublists may be 
 * duplicated. 
 *
 * This function can be used like
 *
 *     int sublen = 2;
 *     int* sublist = GENERATE_COPY( sublists(list, 4, sublen) );
 *
 * to generate {1,2}, {1,3}, {1,4}, {2,1}, {2,3}, {2,4}, {3,1}, {3,2}, {3, 4},
 * {4,1}, {4,2}, {4, 3}.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
CatchGen<int*> sublists(CatchGen<int>&& gen, int sublen);

/** Returns a Catch2 generator of every \p numBits-length bit-set,
 * in increasing lexographic order,
 * where left-most (zero index) bit is treated as LEAST significant (opposite 
 * typical convention). Note that the produced bitset must not be modified during generation.
 *
 * This function can be used like
 *
 *     int* bits = GENERATE( bitsets(3) );
 *
 * to produce {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
CatchGen<int*> bitsets(int numBits);

/** Returns a Catch2 generator of every \p numDigits-length sequence in the given 
 * \p base, in increasing lexographic order,
 * where left-most (zero index) bit is treated as LEAST significant (opposite 
 * typical convention). Note that the sequence must not be modified during generation.
 *
 * This function can be used like
 *
 *     int base = 3;
 *     int numDigits = 2;
 *     int* seq = GENERATE_COPY( sequences(base, numDigits) );
 *
 * to produce {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2}.
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
CatchGen<int*> sequences(int base, int numDigits);

/** Returns a Catch2 generator of every \p numPaulis-length set of Pauli-matrix 
 * types (or base-4 integers).
 * Generates in increasing lexographic order, where the left-most (zero index) 
 * pauli is treated as LEAST significant. 
 * Note that the sequence must not be modified during generation.
 *
 * This function can be used like
 *
 *     pauliOpType* set = GENERATE( pauliseqs(2) );
 *
 * to produce {I,I}, {X,I}, {Y,I}, {Z,I}, {I,X}, {X,X}, {Y,X}, {Z,X}, {I,Y},
 * {X,Y}, {Y,Y}, {Z,Y}, {I,Z}, {X,Z}, {Y,Z}, {Z,Z}/
 *
 * @ingroup deprecatedutils 
 * @author Tyson Jones
 */
CatchGen<pauliOpType*> pauliseqs(int numPaulis);

#endif // QUEST_TEST_UTILS_H
