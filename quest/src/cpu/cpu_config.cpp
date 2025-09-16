/** @file
 * Utility functions for querying the CPU multithreadng
 * configuration, and allocating and copying RAM data.
 * 
 * @author Tyson Jones
 * @author Luc Jaulmes (NUMA awareness)
 */

#include "quest/include/config.h"
#include "quest/include/types.h"
#include "quest/include/paulis.h"

#include "quest/src/core/memory.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"

#include <vector>
#include <cstring>
#include <cstdlib>
#include <cstdint>

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


/// @todo
/// Windows provides a NUMA API we could access in theory, although we 
/// forego the hassle for now - who is running QuEST on big multi-core 
/// Windows? This validation protects against enabling NUMA awareness
/// on Windows but silently recieving no benefit due to no NUMA API calls

#if NUMA_AWARE && defined(_WIN32)
    #error "NUMA awareness is not currently supported on non-POSIX systems like Windows."
#endif


#if COMPILE_OPENMP
    #include <omp.h>
#endif

#if NUMA_AWARE && ! defined(_WIN32)
    #include <sys/mman.h>
    #include <numaif.h>
    #include <numa.h>
#endif

#if defined(_WIN32)
    #define NOMINMAX
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#else
    #include <unistd.h>
#endif



/*
 * OPENMP CONFIG
 */


bool cpu_isOpenmpCompiled() {
    return (bool) COMPILE_OPENMP;
}


int cpu_getAvailableNumThreads() {
#if COMPILE_OPENMP
    int n = -1;

    #pragma omp parallel shared(n)
    #pragma omp single
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


int cpu_getCurrentNumThreads() {
#if COMPILE_OPENMP
    return omp_get_num_threads();
#else
    return 1;
#endif
}



/*
 * MEMORY ALLOCATION
 */


qindex getNumPagesToContainArray(long pageLen, qindex arrLen) {

    // round up to the nearest page
    return static_cast<qindex>(std::ceil(arrLen / (qreal) pageLen));
}


long cpu_getPageSize() {

    // avoid repeated queries to this fixed value
    static long pageSize = 0;
    if (pageSize > 0)
        return pageSize;

    // obtain pageSize for the first time
#if defined(_WIN32)
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);
    pageSize = sysInfo.dwPageSize;
#else
    pageSize = sysconf(_SC_PAGESIZE);
#endif

    // rigorously check the found pagesize is valid
    // and consistent with preconditions assumed by
    // callers, to avoid extremely funky bugs on
    // esoteric future systems

    if (pageSize <= 0)
        error_gettingPageSizeFailed();

    if (!isPowerOf2(pageSize))
        error_pageSizeNotAPowerOf2();

    if (pageSize % sizeof(qcomp) != 0)
        error_pageSizeNotAMultipleOfQcomp();

    return pageSize;
}


qcomp* cpu_allocArray(qindex length) {
    return (qcomp*) calloc(length, sizeof(qcomp));
}


qcomp* cpu_allocNumaArray(qindex length) {
#if ! NUMA_AWARE
    return cpu_allocArray(length);

#elif defined(_WIN32)
    error_numaAllocOrDeallocAttemptedOnWindows();

#else
    // we will divide array's memory into pages
    long pageSize = cpu_getPageSize();
    qindex arraySize = length * sizeof(qcomp); // gauranteed no overflow

    // if entire array fits within a single page, alloc like normal
    if (arraySize <= pageSize)
        return cpu_allocArray(length);

    // otherwise we will bind pages across NUMA nodes
    static int numNodes = numa_num_configured_nodes();
    if (numNodes < 1)
        error_gettingNumNumaNodesFailed();

    qindex numPages = getNumPagesToContainArray(pageSize, arraySize);
    qindex numBytes = numPages * pageSize; // prior validation gaurantees no overflow
    
    // allocate memory, potentially more than arraySize (depending on page divisibility)
    void *rawAddr = mmap(NULL, numBytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    
    // indicate memory alloc failure to caller (no NUMA-specific validation error message)
    if (rawAddr == MAP_FAILED)
        return nullptr;

    // if there is only a single NUMA node, then all memory access will occur within it
    qcomp* outAddr = reinterpret_cast<qcomp*>(rawAddr);
    if (numNodes == 1)
        return outAddr;

    // otherwise, we bind continguous pages to NUMA nodes, distributing the pages 
    // attemptedly uniformly and spreading remaining pages maximally apart
    qindex baseNumPagesPerNode = numPages / numNodes; // floors
    qindex remainingNumPagesTotal = numPages % numNodes;

    // use integer type for safe address arithmetic below
    uintptr_t offsetAddr = reinterpret_cast<uintptr_t>(rawAddr);

    for (int node=0, shift=numNodes; node < numNodes; ++node) {

        // decide number of pages to bind to NUMA node
        shift -= remainingNumPagesTotal;
        qindex numPagesInNode = baseNumPagesPerNode + (shift <= 0);
        qindex numBytesInNode = numPagesInNode * pageSize; // validation prevents overflow

        // bind those pages from the offset address to the node (identified by mask)
        unsigned long nodeMask = 1UL << node;
        unsigned long numBitsInMask = 8 * sizeof(nodeMask);
        void* nodeAddr = reinterpret_cast<void*>(offsetAddr);
        long success = mbind(nodeAddr, numBytesInNode, MPOL_BIND, &nodeMask, numBitsInMask, 0);

        // treat bind failure as internal error (even though it can result from insufficient kernel mem),
        // rather than permitting silent fallback to non-NUMA awareness which might be astonishingly slow
        if (success == -1)
            error_numaBindingFailed();

        // prepare next node's address
        offsetAddr += numPagesInNode * pageSize;
        if (shift <= 0)
            shift += numNodes;
    }

    return outAddr;
#endif
}


void cpu_deallocArray(qcomp* arr) {

    // arr can safely be nullptr
    free(arr);
}


void cpu_deallocNumaArray(qcomp* arr, qindex length) {

    // musn't pass nullptr to munmap() below
    if (arr == nullptr)
        return;

#if ! NUMA_AWARE
    cpu_deallocArray(arr);

#elif defined(_WIN32)
    error_numaAllocOrDeallocAttemptedOnWindows();

#else
    qindex arrSize = length * sizeof(qcomp);
    long pageSize = cpu_getPageSize();

    // sub-page arrays were allocated with calloc()
    if (arrSize <= pageSize)
        return cpu_deallocArray(arr);

    qindex numPages = getNumPagesToContainArray(pageSize, arrSize);
    qindex numBytes = numPages * pageSize; // gauranteed no overflow
    int success = munmap(arr, numBytes);

    if (success == -1)
        error_numaUnmappingFailed();
#endif
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
