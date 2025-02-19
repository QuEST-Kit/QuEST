#include "quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/evolve.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/lists.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"
#include "tests/utils/cache.hpp"

#include <vector>
using std::vector;


#define TEST_TAG "[calculations]"



qreal calcExpecPauliStrSum(Qureg qureg, PauliStrSum sum);

qreal calcExpecFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);

qreal calcExpecFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);




qreal calcTotalProb(Qureg qureg);

qreal calcProbOfBasisState(Qureg qureg, qindex index);

qreal calcProbOfQubitOutcome(Qureg qureg, int qubit, int outcome);

qreal calcProbOfMultiQubitOutcome(Qureg qureg, int* qubits, int* outcomes, int numQubits);

void  calcProbsOfAllMultiQubitOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits);


qreal calcPurity(Qureg qureg);

qreal calcFidelity(Qureg qureg, Qureg other);

qreal calcDistance(Qureg qureg1, Qureg qureg2);


Qureg calcPartialTrace(Qureg qureg, int* traceOutQubits, int numTraceQubits);

Qureg calcReducedDensityMatrix(Qureg qureg, int* retainQubits, int numRetainQubits);

void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits);

void setQuregToReducedDensityMatrix(Qureg out, Qureg in, int* retainQubits, int numRetainQubits);


qcomp calcInnerProduct(Qureg qureg1, Qureg qureg2);

qcomp calcExpecNonHermitianPauliStrSum(Qureg qureg, PauliStrSum sum); 

qcomp calcExpecNonHermitianFullStateDiagMatr(Qureg qureg, FullStateDiagMatr matr);

qcomp calcExpecNonHermitianFullStateDiagMatrPower(Qureg qureg, FullStateDiagMatr matrix, qcomp exponent);
