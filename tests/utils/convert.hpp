/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsconvert Convert
 * @ingroup testutils
 * @brief
 * Testing utilities for converting QuEST API structures
 * (like Qureg, CompMatr, PauliStr) to/from testing types 
 * (like qvector and qmatrix).
 * @{
 */

#ifndef CONVERT_HPP
#define CONVERT_HPP

#include "quest/include/quest.h"
#include "qvector.hpp"
#include "qmatrix.hpp"


void setQuregToReference(Qureg, qvector);
void setQuregToReference(Qureg, qmatrix);

qvector getVector(Qureg);
qmatrix getMatrix(Qureg);

template <typename T> qmatrix getMatrix(T);

qmatrix getMatrix(PauliStr str,     vector<int> targs);
qmatrix getMatrix(PauliStr str,    int numQubits);
qmatrix getMatrix(PauliStrSum sum, int numQubits);


#endif // CONVERT_HPP

/** @} (end defgroup) */
