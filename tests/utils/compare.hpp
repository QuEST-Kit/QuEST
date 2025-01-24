#ifndef COMPARE_HPP
#define COMPARE_HPP

#include "quest.h"
#include "qvector.hpp"
#include "qmatrix.hpp"


bool operator == (const qvector&, const qvector&);
bool operator == (const qmatrix&, const qmatrix&);
bool operator == (const qvector&, const Qureg&);
bool operator == (const Qureg&,   const qvector&);
bool operator == (const qmatrix&, const Qureg&);
bool operator == (const Qureg&,   const qmatrix&);
bool operator == (const Qureg&,   const Qureg&);


#endif // COMPARE_HPP