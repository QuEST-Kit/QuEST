/** @file
 * Unit tests of the decoherence module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitdeco Decoherence unit tests
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/evolve.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/lists.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/random.hpp"



/*
 * TODO:
 * UNTESTED FUNCTIONS
 */

void mixDephasing(Qureg qureg, int qubit, qreal prob);

void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob);

void mixDepolarising(Qureg qureg, int qubit, qreal prob);

void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob);

void mixDamping(Qureg qureg, int qubit, qreal prob);

void mixPaulis(Qureg qureg, int qubit, qreal probX, qreal probY, qreal probZ);

void mixQureg(Qureg qureg, Qureg other, qreal prob);

void mixKrausMap(Qureg qureg, int* qubits, int numQubits, KrausMap map);

