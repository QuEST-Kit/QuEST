/** @file
 * Utility signatures for querying the CPU multithreadng
 * configuration, and allocating and copying RAM data.
 * 
 * @author Tyson Jones
 */

#ifndef CPU_CONFIG_HPP
#define CPU_CONFIG_HPP

#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include <vector>

using std::vector;



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

qcomp** cpu_allocAndInitMatrixWrapper(qcomp* arr, qindex dim);
void cpu_deallocMatrixWrapper(qcomp** wrapper);

qcomp** cpu_allocMatrix(qindex dim);
void cpu_deallocMatrix(qcomp** matrix, qindex dim);

qcomp*** cpu_allocMatrixList( qindex numRows, int numMatrices);
void cpu_deallocMatrixList(qcomp*** matrices, qindex numRows, int numMatrices);

int* cpu_allocHeapFlag();
void cpu_deallocHeapFlag(int* ptr);

PauliStr* cpu_allocPauliStrings(qindex numStrings);
void cpu_deallocPauliStrings(PauliStr* strings);



/*
 * MEMORY MOVEMENT
 */

void cpu_copyArray(qcomp* dest, qcomp* src, qindex dim);

void cpu_copyMatrix(qcomp** dest, qcomp** src, qindex dim);
void cpu_copyMatrix(qcomp** dest, vector<vector<qcomp>> src, qindex dim);

void cpu_copyPauliStrSum(PauliStrSum out, PauliStr* strings, qcomp* coeffs); 



#endif // CPU_CONFIG_HPP