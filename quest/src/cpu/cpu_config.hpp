/** @file
 * Utility signatures for querying CPU hardware.
 */

#ifndef CPU_CONFIG_HPP
#define CPU_CONFIG_HPP

#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include <vector>

using std::vector;



/*
 * ENABLE OPENMP REDUCTION OF qcomp
 * 
 * which is incompatible with MSVC; we have not yet decided
 * how to remedy this - force Windows users who seek multithreading
 * to compile using another compiler?
 */

#if defined(COMPILE_OPENMP) && !defined(_MSC_VER)
    #pragma omp declare reduction(+ : qcomp : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#endif



/*
 * OPENMP CONFIG
 */

bool cpu_isOpenmpCompiled();

int cpu_getCurrentNumThreads();

int cpu_getNumOpenmpProcessors();



/*
 * OPENMP SUBROUTINES
 */

int cpu_getOpenmpThreadInd();



/*
 * MEMORY ALLOCATION
 */

qcomp* cpu_allocArray(qindex length);
void cpu_deallocArray(qcomp* arr);

qcomp** cpu_allocMatrix(qindex dim);
void cpu_deallocMatrix(qcomp** matrix, qindex dim);

qcomp*** cpu_allocMatrixList(int numMatrices, qindex numRows);
void cpu_deallocMatrixList(qcomp*** matrices, int numMatrices, qindex numRows);

int* cpu_allocHeapFlag();
void cpu_deallocHeapFlag(int* ptr);

PauliStr* cpu_allocPauliStrings(qindex numStrings);
void cpu_deallocPauliStrings(PauliStr* strings);



/*
 * MEMORY MOVEMENT
 */

void cpu_copyArray(qcomp* out, qcomp* in, qindex dim);

void cpu_copyMatrix(qcomp** out, qcomp** in, qindex dim);
void cpu_copyMatrix(qcomp** out, vector<vector<qcomp>> in, qindex dim);

void cpu_copyPauliStrSum(PauliStrSum out, PauliStr* strings, qcomp* coeffs); 



#endif // CPU_CONFIG_HPP