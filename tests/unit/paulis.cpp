#include "quest.h"

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
 * UNTESTED FUNCTIONS BELOW
 */

PauliStr getPauliStr(const char* paulis, int* indices, int numPaulis);
PauliStr getPauliStr(int* paulis, int* indices, int numPaulis);
PauliStr getPauliStr(std::string paulis, int* indices, int numPaulis);
PauliStr getPauliStr(std::string paulis, std::vector<int> indices);
PauliStr getPauliStr(std::string paulis);


// (macro) getInlinePauliStr


PauliStrSum createPauliStrSum(PauliStr* strings, qcomp* coeffs, qindex numTerms);
PauliStrSum createPauliStrSum(std::vector<PauliStr> strings, std::vector<qcomp> coeffs);


PauliStrSum createInlinePauliStrSum(const char* str);
PauliStrSum createInlinePauliStrSum(std::string str);

PauliStrSum createPauliStrSumFromFile(const char* fn);
PauliStrSum createPauliStrSumFromFile(std::string fn);

PauliStrSum createPauliStrSumFromReversedFile(const char* fn);
PauliStrSum createPauliStrSumFromReversedFile(std::string fn);


void destroyPauliStrSum(PauliStrSum sum);

void reportPauliStr(PauliStr str);

void reportPauliStrSum(PauliStrSum str);