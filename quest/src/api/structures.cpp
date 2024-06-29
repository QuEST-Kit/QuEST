/** @file
 * Definitions of API data structures like gate matrices. 
 * Note QuESTEnv and Qureg structs have their own definitions 
 * in environment.cpp and qureg.cpp respectively.
 */

#include "quest/include/structures.h"
#include "quest/include/environment.h"
#include "quest/include/types.h"

#include "quest/src/comm/comm_config.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/core/formatter.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <cstdlib>
#include <iostream>
#include <string>


// for concise reporters
using namespace form_substrings;



/*
 * C++-ONLY MATRIX 1&2 CONSTRUCTORS
 * 
 * and their corresponding C-compatible alternatives which are wrapped by wrappers.h.
 * These are necessary because of the non-interchangeable C and C++ qcomp types.
 * See structures.h for an explanation.
 */


CompMatr1 getCompMatr1(qcomp in[2][2]) {

    CompMatr1 out = {
        .numQubits = 1,
        .numRows = 2,
        .elems = {
            in[0][0], in[0][1],
            in[1][0], in[1][1]
        }
    };

    return out;
}
extern "C" void wrap_getCompMatr1(CompMatr1* out, qcomp in[2][2]) {

    *out = getCompMatr1(in);
}


CompMatr2 getCompMatr2(qcomp in[4][4]) {

    CompMatr2 out = {
        .numQubits = 2,
        .numRows = 4,
        .elems = {
            in[0][0], in[0][1], in[0][2], in[0][3],
            in[1][0], in[1][1], in[1][2], in[1][3],
            in[2][0], in[2][1], in[2][2], in[2][3],
            in[3][0], in[3][1], in[3][2], in[3][3]
        }
    };
    
    return out;
}
extern "C" void wrap_getCompMatr2(CompMatr2* out, qcomp in[4][4]) {

    *out = getCompMatr2(in);
}



/*
 * C & C++ MATRIX N CONSTRUCTORS
 *
 * all of which are de-mangled for C-compatibility
 */


extern "C" CompMatrN createCompMatrN(int numQubits) {
    validate_newMatrixNumQubits(numQubits, __func__);

    // allocate GPU memory if the env is GPU-accelerated, regardless of (unknown) Qureg allocations
    bool isGpuAccel = getQuESTEnv().isGpuAccelerated;

    // validation ensures these (and below mem sizes) never overflow
    qindex numRows = powerOf2(numQubits);
    qindex numElems = numRows * numRows;

    CompMatrN out = {
        .numQubits = numQubits,
        .numRows = numRows,

        // 2D CPU memory (NULL if failed)
        .elems = (qcomp**) malloc(numRows * sizeof *out.elems),

        // 1D GPU memory (NULL if failed or not needed)
        .gpuElems = (isGpuAccel)? gpu_allocAmps(numElems) : NULL // first amp will be un-sync'd flag
    };

    // only if outer CPU allocation succeeded, attempt to allocate each row array
    if (out.elems != NULL)
        for (qindex r=0; r < numRows; r++)
            out.elems[r] = (qcomp*) calloc(numRows, sizeof **out.elems); // NULL if failed

    // check all CPU & GPU malloc and calloc's succeeded (no attempted freeing if not)
    bool isNewMatr = true;
    validate_newOrExistingMatrixAllocs(out, isNewMatr, __func__);

    return out;
}


extern "C" void destroyCompMatrN(CompMatrN matrix) {
    validate_matrixInit(matrix, __func__);

    // free each CPU row array
    for (qindex r=0; r < matrix.numRows; r++)
        free(matrix.elems[r]);

    // free CPU array of rows
    free(matrix.elems);

    // free flat GPU array if it exists
    if (matrix.gpuElems != NULL)
        gpu_deallocAmps(matrix.gpuElems);
}



/*
 * MATRIX N INITIALISER
 */

extern "C" void syncCompMatrN(CompMatrN matr) {
    validate_matrixInit(matr, __func__);
    validate_matrixElemsDontContainUnsyncFlag(matr.elems, __func__);

    gpu_copyCpuToGpu(matr);
}

extern "C" void setCompMatrN(CompMatrN matr, qcomp** vals) {
    validate_matrixInit(matr, __func__);
    validate_matrixElemsDontContainUnsyncFlag(vals, __func__);

    // serially copy values to CPU memory
    populateCompMatrElems(matr.elems, vals, matr.numRows);

    // overwrite GPU elements
    if (matr.gpuElems != NULL)
        gpu_copyCpuToGpu(matr);
}



/*
 * C & C++ MATRIX REPORTERS
 *
 * and private (permittedly name-mangled) inner functions
 */


void printMatrixHeader(int numQubits) {

    // find memory used by matrix; equal to that of a non-distributed density matrix
    bool isMatr = true;
    int numNodes = 1;
    size_t mem = mem_getLocalMemoryRequired(numQubits, isMatr, numNodes);
    std::string memStr = form_str(mem) + by;

    // prepare dim substring
    qindex dim = powerOf2(numQubits);
    std::string dimStr = form_str(dim) + mu + form_str(dim);

    // e.g. CompMatr2 (4 x 4, 256 bytes):
    std::cout << "CompMatr" << numQubits << " (" << dimStr << ", " << memStr << "):" << std::endl;
}


template<class T> 
void rootPrintMatrix(T matrix) {
    
    if (comm_getRank() != 0)
        return;

    printMatrixHeader(matrix.numQubits);
    form_printMatrix(matrix);
}


extern "C" void reportCompMatr1(CompMatr1 matr) {

    rootPrintMatrix(matr);
}

extern "C" void reportCompMatr2(CompMatr2 matr) {

    rootPrintMatrix(matr);
}

extern "C" void reportCompMatrN(CompMatrN matr) {

    rootPrintMatrix(matr);
}

