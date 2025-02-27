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

static inline CompMatr1 getCompMatr1(qcomp** in);
CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in);

static inline CompMatr2 getCompMatr2(qcomp** in);
CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in);

static inline DiagMatr1 getDiagMatr1(qcomp* in);
DiagMatr1 getDiagMatr1(std::vector<qcomp> in);

static inline DiagMatr2 getDiagMatr2(qcomp* in);
DiagMatr2 getDiagMatr2(std::vector<qcomp> in);

// (macros)
// getInlineCompMatr1
// getInlineCompMatr2
// getInlineDiagMatr1
// getInlineDiagMatr2

CompMatr createCompMatr(int numQubits);

DiagMatr createDiagMatr(int numQubits);

FullStateDiagMatr createFullStateDiagMatr(int numQubits);

FullStateDiagMatr createCustomFullStateDiagMatr(int numQubits, int useDistrib, int useGpuAccel);


void destroyCompMatr(CompMatr matrix);

void destroyDiagMatr(DiagMatr matrix);

void destroyFullStateDiagMatr(FullStateDiagMatr matrix);


void syncCompMatr(CompMatr matr);

void syncDiagMatr(DiagMatr matr);

void syncFullStateDiagMatr(FullStateDiagMatr matr);


void setCompMatr(CompMatr matr, qcomp** vals);
void setCompMatr(CompMatr out, std::vector<std::vector<qcomp>> in);

void setDiagMatr(DiagMatr out, qcomp* in);
void setDiagMatr(DiagMatr out, std::vector<qcomp> in);

void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems);
void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, std::vector<qcomp> in);


void setInlineCompMatr(CompMatr matr, int numQb, std::vector<std::vector<qcomp>> in);

void setInlineDiagMatr(DiagMatr matr, int numQb, std::vector<qcomp> in);

void setInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems, std::vector<qcomp> in);


CompMatr createInlineCompMatr(int numQb, std::vector<std::vector<qcomp>> elems);

DiagMatr createInlineDiagMatr(int numQb, std::vector<qcomp> elems);


void setDiagMatrFromMultiVarFunc(DiagMatr out, qcomp (*func)(qindex*), int* numQubitsPerVar, int numVars, int areSigned);

void setDiagMatrFromMultiDimLists(DiagMatr out, void* lists, int* numQubitsPerDim, int numDims);


FullStateDiagMatr createFullStateDiagMatrFromPauliStrSum(PauliStrSum in);

void setFullStateDiagMatrFromPauliStrSum(FullStateDiagMatr out, PauliStrSum in);

void setFullStateDiagMatrFromMultiVarFunc(FullStateDiagMatr out, qcomp (*func)(qindex*), int* numQubitsPerVar, int numVars, int areSigned);

void setFullStateDiagMatrFromMultiDimLists(FullStateDiagMatr out, void* lists, int* numQubitsPerDim, int numDims);


void reportCompMatr1(CompMatr1 matrix);

void reportCompMatr2(CompMatr2 matrix);

void reportCompMatr(CompMatr matrix);


void reportDiagMatr1(DiagMatr1 matrix);

void reportDiagMatr2(DiagMatr2 matrix);

void reportDiagMatr(DiagMatr matrix);


void reportFullStateDiagMatr(FullStateDiagMatr matr);