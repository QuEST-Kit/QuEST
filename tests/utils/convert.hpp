#ifndef CONVERT_HPP
#define CONVERT_HPP

#include "quest/include/quest.h"
#include "qvector.hpp"
#include "qmatrix.hpp"


void setQureg(Qureg, qvector);
void setQureg(Qureg, qmatrix);

qvector getVector(Qureg);
qmatrix getMatrix(Qureg);

template <typename T> qmatrix getMatrix(T);


#endif // CONVERT_HPP