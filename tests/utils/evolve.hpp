#ifndef EVOLVE_HPP
#define EVOLVE_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;


void applyReferenceOperator(   qvector& state, vector<int> ctrls, vector<int> targs, qmatrix matrix);
void applyReferenceOperator(   qmatrix& state, vector<int> ctrls, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qmatrix& state, vector<int> ctrls, vector<int> targs, qmatrix matrix);

void applyReferenceOperator(   qvector& state, vector<int> targs, qmatrix matrix);
void applyReferenceOperator(   qmatrix& state, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qmatrix& state, vector<int> targs, qmatrix matrix);


#endif // EVOLVE_HPP