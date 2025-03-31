/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsrandom Random
 * @ingroup testutils
 * @brief
 * Testing utilities which generate random objects
 * independently of QuEST's internal generators. 
 * @{
 */

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"
#include "macros.hpp"
#include "linalg.hpp"
#include "lists.hpp"
#include "quest/include/quest.h"

#include <vector>
using std::vector;

#include <tuple>
using std::tuple;


void setRandomTestStateSeeds();

int getRandomInt(int min, int maxExcl);
qreal getRandomReal(qreal min, qreal maxExcl);
qreal getRandomPhase();
qcomp getRandomComplex();
qcomp getRandomUnitComplex();

listpair getRandomFixedNumCtrlsTargs(int numQubits, int numCtrls, int numTargs);
listtrio getRandomVariNumCtrlsStatesTargs(int numQubits, int minNumTargs, int maxNumTargsIncl);

qvector getRandomVector(size_t dim);
qmatrix getRandomMatrix(size_t dim);
qmatrix getRandomDiagonalMatrix(size_t dim);
qmatrix getRandomNonSquareMatrix(size_t numRows, size_t numCols);

qvector getRandomStateVector(int numQb);
qmatrix getRandomDensityMatrix(int numQb);
qmatrix getRandomPureDensityMatrix(int numQb);

void setToRandomState(qvector& state);
void setToRandomState(qmatrix& state);

qmatrix getRandomUnitary(int numQb);
qmatrix getRandomDiagonalUnitary(int numQb);
qmatrix getRandomDiagonalHermitian(int numQb);
vector<qmatrix> getRandomKrausMap(int numQb, int numOps);

PauliStr getRandomPauliStr(int numQubits);
PauliStr getRandomPauliStr(vector<int> targs);
PauliStr getRandomDiagPauliStr(int numQubits);

vector<int> getRandomInts(int min, int maxExcl, int len);
vector<int> getRandomOutcomes(int len);
vector<int> getRandomSubRange(int start, int endExcl, int numElems);
vector<qreal> getRandomProbabilities(int numProbs);

vector<qvector> getRandomOrthonormalVectors(size_t dim, int numVecs);
vector<qvector> getRandomOrthonormalStateVectors(int numQb, int numStates);

PauliStrSum createRandomPauliStrSum(int numQubits, int numTerms);
PauliStrSum createRandomNonHermitianPauliStrSum(int numQubits, int numTerms);


#endif // RANDOM_HPP

/** @} (end defgroup) */
