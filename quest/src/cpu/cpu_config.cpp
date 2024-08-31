/** @file
 * Utility definitions for querying CPU hardware.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include "quest/src/core/errors.hpp"

#include <vector>
#include <cstring>
#include <cstdlib>

using std::vector;


#if COMPILE_OPENMP && !defined(_OPENMP)
    #error "Attempted to compile in multithreaded mode without enabling OpenMP in the compiler flags."
#endif


#if COMPILE_OPENMP
    #include <omp.h>
#endif



/*
 * ENABLE OPENMP REDUCTION OF qcomp (except on MSVC compilers)
 */

#if defined(COMPILE_OPENMP) && !defined(_MSC_VER)
     #pragma omp declare reduction(+ : qcomp : omp_out += omp_in ) initializer( omp_priv = omp_orig )
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
 * MEMORY ALLOCATION
 */


qcomp* cpu_allocArray(qindex length) {

    // we call calloc over malloc in order to fail immediately if mem isn't available;
    // caller must handle NULL result
    return (qcomp*) calloc(length, sizeof(qcomp));
}


void cpu_deallocArray(qcomp* arr) {

    // arr can safely be NULL
    free(arr);
}


qcomp** cpu_allocMatrix(qindex dim) {

    // TODO:
    // the design of storing the CPU matrix elements as a 2D structure will impede
    // performance for many qubits; the allocated heap memories for each row
    // have no gaurantee to reside near other, so that their access/iteration in
    // hot loops may incur unnecessary caching penalties. Consider storing the
    // elements as a flat array, like we do for the GPU memory. This makes manual
    // modification by the user trivially harder (changing [r][c] to [r*n+c]),
    // but should improve caching, and significantly simplify allocation and its
    // validation; no more enumerating nested pointers! Benchmark this scenario.

    // allocate outer array
    qcomp** rows = (qcomp**) malloc(dim * sizeof *rows); // NULL if failed

    // if that did not fail, allocate each inner array
    if (rows != NULL)
        for (qindex r=0; r<dim; r++)
            rows[r] = cpu_allocArray(dim); // NULL if failed

    // caller will validate whether mallocs were successful
    return rows;
}


void cpu_deallocMatrix(qcomp** matrix, qindex dim) {

    // we attempt to deallocate every row (assuming the outer array was
    // successfully allocated), regardless of whether they are actually
    // allocated; it is legal to call free() on NULL

    if (matrix != NULL)
        for (qindex r=0; r<dim; r++)
            cpu_deallocArray(matrix[r]);

    free(matrix);
}


int* cpu_allocHeapFlag() {

    // we use int over bool for the flag, because often we use
    // value -1 as a third value
    return (int*) malloc(sizeof(int)); // may be NULL, caller will handle
}


void cpu_deallocHeapFlag(int* ptr) {

    // safe to free if NULL
    free(ptr);
}


PauliStr* cpu_allocPauliStrings(qindex numStrings) {

    return (PauliStr*) malloc(numStrings * sizeof(PauliStr)); // may be NULL, caller will handle
}


void cpu_deallocPauliStrings(PauliStr* strings) {

    // safe to free if NULL
    free(strings);
}


/*
 * MEMORY COPYING
 */


void cpu_copyArray(qcomp* out, qcomp* in, qindex dim) {

    memcpy(out, in, dim * sizeof(qcomp));
}


void cpu_copyMatrix(qcomp** out, qcomp** in, qindex dim) {

    // note that we cannot call a single memcpy to copy all rows at once,
    // because in/out may not be contiguous stack arrays; instead, each
    // row is likely a unique, discontiguous span of heap memory. So we
    // memcpy each row in-turn
    for (qindex r=0; r<dim; r++)
        cpu_copyArray(out[r], in[r], dim);
}

void cpu_copyMatrix(qcomp** out, vector<vector<qcomp>> in, qindex dim) {

    for (qindex r=0; r<dim; r++)
        cpu_copyArray(out[r], in[r].data(), dim);
}


void cpu_copyPauliStrSum(PauliStrSum out, PauliStr* strings, qcomp* coeffs) {

    // serially copy data over to new heap memory
    for (int i=0; i<out.numTerms; i++) {
        out.strings[i] = strings[i];
        out.coeffs[i] = coeffs[i];
    }
}