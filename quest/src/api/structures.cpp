/** @file
 * Definitions of API data structures like gate matrices. 
 * Note QuESTEnv and Qureg structs have their own definitions 
 * in environment.cpp and qureg.cpp respectively.
 */

#include "quest/include/structures.h"
#include "quest/include/types.h"

#include "quest/src/comm/comm_config.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/core/formatter.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"

#include <cstdlib>
#include <iostream>
#include <string>


// for concise reporters
using namespace form_substrings;



/*
 * C++-ONLY MATRIX CONSTRUCTORS
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
 * C & C++ MATRIX CONSTRUCTORS
 *
 * all of which are de-mangled for C-compatibility
 */


extern "C" CompMatrN createCompMatrN(int numQubits) {
    validate_newMatrixNumQubits(numQubits, __func__);

    // validation ensures this and mem sizes below don't overflow
    qindex numRows = powerOf2(numQubits);

    // allocate array for rows
    CompMatrN out = {
        .numQubits = numQubits,
        .numRows = numRows,
        .elems = (qcomp**) malloc(numRows * sizeof *out.elems) // NULL if failed
    };

    // only if that allocation succeeded, allocate each row array
    if (out.elems != NULL)
        for (qindex r=0; r < numRows; r++)
            out.elems[r] = (qcomp*) calloc(numRows, sizeof **out.elems);

    // check malloc and calloc's succeeded
    bool isNewMatr = true;
    validate_newOrExistingMatrixAllocs(out, isNewMatr, __func__);

    return out;
}


extern "C" void destroyCompMatrN(CompMatrN matrix) {
    validate_matrixInit(matrix, __func__);

    for (qindex r=0; r < matrix.numRows; r++)
        free(matrix.elems[r]);

    free(matrix.elems);
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

