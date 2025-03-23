/** @file
 * Utility functions for querying the CPU multithreadng
 * configuration, and allocating and copying RAM data.
 * 
 * @author Tyson Jones
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include "quest/src/core/errors.hpp"

#include <vector>
#include <cstring>
#include <cstdlib>

using std::vector;


// when COMPILE_OPENMP=1, the compiler expects arguments like -fopenmp
// which cause _OPENMP to be defined, which we check to ensure that
// COMPILE_OPENMP has been set correctly. Note that HIP compilers do
// not define _OPENMP even when parsing OpenMP, and it's possible that
// the user is compiling all the source code (including this file) with
// HIP; we tolerate _OPENMP being undefined in that instance

#if COMPILE_OPENMP && !defined(_OPENMP) && !defined(__HIP__)
    #error "Attempted to compile in multithreaded mode without enabling OpenMP in the compiler flags."
#endif


#if COMPILE_OPENMP
    #include <omp.h>
#endif



/*
 * OPENMP CONFIG
 */


bool cpu_isOpenmpCompiled() {
    return (bool) COMPILE_OPENMP;
}


int cpu_getCurrentNumThreads() {
#if COMPILE_OPENMP
    int n = -1;

    #pragma omp parallel shared(n)
    n = omp_get_num_threads();

    return n;
#else
    error_cpuThreadsQueriedButEnvNotMultithreaded();
    return -1;
#endif
}


int cpu_getNumOpenmpProcessors() {
#if COMPILE_OPENMP
    return omp_get_num_procs();
#else
    error_cpuThreadsQueriedButEnvNotMultithreaded();
    return -1;
#endif
}



/*
 * OPENMP SUBROUTINES
 *
 * which must be queried within OpenMP parallel
 * regions to get reliable results, but which are
 * safely invoked when OpenMP is not compiled
 */


int cpu_getOpenmpThreadInd() {
#if COMPILE_OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}



/*
 * MEMORY ALLOCATION
 */


qcomp* cpu_allocArray(qindex length) {

    /// @todo
    /// here, we calloc the entire array in a serial setting, rather than one malloc 
    /// followed by threads subsequently memset'ing their own partitions. The latter
    /// approach would distribute the array pages across NUMA nodes, accelerating 
    /// their subsequent access by the same threads (via NUMA's first-touch policy).
    /// We have so far foregone this optimisation since a thread's memory-access pattern
    /// in many of the QuEST functions is non-trivial, and likely to be inconsistent 
    /// with the memset pattern. As such, I expect the benefit is totally occluded
    /// and only introduces potential new bugs - but this should be tested and confirmed!

    // we call calloc over malloc in order to fail immediately if mem isn't available;
    // caller must handle nullptr result

    return (qcomp*) calloc(length, sizeof(qcomp));
}


void cpu_deallocArray(qcomp* arr) {

    // arr can safely be nullptr
    free(arr);
}


qcomp** cpu_allocAndInitMatrixWrapper(qcomp* arr, qindex dim) {

    // do not allocate if arr alloc failed (caller will handle)
    if (arr == nullptr)
        return nullptr;

    // allocate only the outer memory (i.e. one row's worth)
    qcomp** out = (qcomp**) malloc(dim * sizeof *out);

    // caller will handle malloc failure
    if (out == nullptr)
        return out;

    // populate out with offsets of arr
    for (qindex i=0; i<dim; i++)
        out[i] = &arr[i*dim];

    return out; // may be nullptr
}


void cpu_deallocMatrixWrapper(qcomp** wrapper) {

    // only the outer pointer is freed; the
    // inner pointers are offsets to another
    // malloc which is separately freed. 
    // Safe to call even when nullptr
    free(wrapper);
}


qcomp** cpu_allocMatrix(qindex dim) {

    // NOTE:
    // this function creates a matrix where rows are not necessarily
    // contiguous in memory, which can incur gratuitous caching penalties
    // when accessed in hot loops. As such, we do not use this function
    // to allocate memory for CompMatr (instead, cpu_allocAndInitMatrixWrapper()),
    // but instead use it for the individual Kraus matrices of a KrausMap,
    // which are each quadratically smaller than the important superoperator.

    // allocate outer array
    qcomp** rows = (qcomp**) malloc(dim * sizeof *rows); // nullptr if failed

    // if that did not fail, allocate each inner array
    if (rows != nullptr)
        for (qindex r=0; r<dim; r++)
            rows[r] = cpu_allocArray(dim); // nullptr if failed

    // caller will validate whether mallocs were successful
    return rows;
}


void cpu_deallocMatrix(qcomp** matrix, qindex dim) {

    // we attempt to deallocate every row (assuming the outer array was
    // successfully allocated), regardless of whether they are actually
    // allocated; it is legal to call free() on nullptr

    if (matrix != nullptr)
        for (qindex r=0; r<dim; r++)
            cpu_deallocArray(matrix[r]);

    free(matrix);
}


qcomp*** cpu_allocMatrixList(qindex numRows, int numMatrices) {

    // attempt to allocate the outer list
    qcomp*** matrices = (qcomp***) malloc(numMatrices * sizeof *matrices); // nullptr if failed

    // attempt to allocate each matrix
    if (matrices != nullptr)
        for (int n=0; n<numMatrices; n++)
            matrices[n] = cpu_allocMatrix(numRows); // nullptr if failed

    return matrices; // may be or contain nullptrs, user will handle
}


void cpu_deallocMatrixList(qcomp*** matrices, qindex numRows, int numMatrices) {

    // free everything that allocated (but permit anything to have failed)
    if (matrices != nullptr)
        for (int n=0; n<numMatrices; n++)
            cpu_deallocMatrix(matrices[n], numRows);

    // legal to free nullptr
    free(matrices);
}


int* cpu_allocHeapFlag() {

    // we use int over bool for the flag, because often we use
    // value -1 as a third value
    return (int*) malloc(sizeof(int)); // may be nullptr, caller will handle
}


void cpu_deallocHeapFlag(int* ptr) {

    // safe to free if nullptr
    free(ptr);
}


PauliStr* cpu_allocPauliStrings(qindex numStrings) {

    return (PauliStr*) malloc(numStrings * sizeof(PauliStr)); // may be nullptr, caller will handle
}


void cpu_deallocPauliStrings(PauliStr* strings) {

    // safe to free if nullptr
    free(strings);
}



/*
 * MEMORY COPYING
 */


void cpu_copyArray(qcomp* dest, qcomp* src, qindex dim) {

    memcpy(dest, src, dim * sizeof(qcomp));
}


void cpu_copyMatrix(qcomp** dest, qcomp** src, qindex dim) {

    /// @todo
    /// there may be a faster, asynchronous way to perform
    /// these memcpys then do a final wait

    // note that we cannot call a single memcpy to copy all rows at once,
    // because dest/src may not be contiguous stack arrays; instead, each
    // row is likely a unique, discontiguous span of heap memory. So we
    // memcpy each row in-turn
    for (qindex r=0; r<dim; r++)
        cpu_copyArray(dest[r], src[r], dim);
}

void cpu_copyMatrix(qcomp** dest, vector<vector<qcomp>> src, qindex dim) {

    /// @todo
    /// there may be a faster, asynchronous way to perform
    /// these memcpys then do a final wait

    for (qindex r=0; r<dim; r++)
        cpu_copyArray(dest[r], src[r].data(), dim);
}


void cpu_copyPauliStrSum(PauliStrSum out, PauliStr* strings, qcomp* coeffs) {

    // serially copy data over to new heap memory
    for (int i=0; i<out.numTerms; i++) {
        out.strings[i] = strings[i];
        out.coeffs[i] = coeffs[i];
    }
}
