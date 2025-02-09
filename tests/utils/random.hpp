#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"
#include "macros.hpp"
#include "linalg.hpp"
#include "quest.h"

#include <vector>
using std::vector;


void setRandomTestStateSeeds();

int getRandomInt(int min, int maxExcl);
qreal getRandomReal(qreal min, qreal max);
qcomp getRandomComplex();

qvector getRandomVector(size_t dim);
qmatrix getRandomMatrix(size_t dim);

qvector getRandomStateVector(int numQb);
qmatrix getRandomDensityMatrix(int numQb);
qmatrix getRandomPureDensityMatrix(int numQb);

qmatrix getRandomUnitary(int numQb);
qmatrix getRandomDiagonalUnitary(int numQb);
vector<qmatrix> getRandomKrausMap(int numQb, int numOps);

PauliStr getRandomPauliStr(int numQubits);
PauliStr getRandomDiagPauliStr(int numQubits);

vector<int> getRandomInts(int min, int maxExcl, int len);
vector<int> getRandomSubRange(int start, int endExcl, int numElems);
vector<qreal> getRandomProbabilities(int numProbs);

vector<qvector> getRandomOrthonormalVectors(size_t dim, int numVecs);
vector<qvector> getRandomOrthonormalStateVectors(int numQb, int numStates);


#endif // RANDOM_HPP