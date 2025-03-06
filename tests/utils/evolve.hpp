/** @file
 * Testing utilities which evolve a reference state
 * (qvector or qmatrix) under the action of a 
 * reference operation. These are slow, serial,
 * un-optimised, defensively-designed routines.
 *
 * @author Tyson Jones
 */

#ifndef EVOLVE_HPP
#define EVOLVE_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;

qmatrix getControlledMatrix(qmatrix matrix, int numCtrls);

void applyReferenceOperator(   qvector& state, vector<int> ctrls, vector<int> states, vector<int> targs, qmatrix matrix);
void applyReferenceOperator(   qmatrix& state, vector<int> ctrls, vector<int> states, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qvector& state, vector<int> ctrls, vector<int> states, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qmatrix& state, vector<int> ctrls, vector<int> states, vector<int> targs, qmatrix matrix);

void applyReferenceOperator(   qvector& state, vector<int> ctrls, vector<int> targs, qmatrix matrix);
void applyReferenceOperator(   qmatrix& state, vector<int> ctrls, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qvector& state, vector<int> ctrls, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qmatrix& state, vector<int> ctrls, vector<int> targs, qmatrix matrix);

void applyReferenceOperator(   qvector& state, vector<int> targs, qmatrix matrix);
void applyReferenceOperator(   qmatrix& state, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qvector& state, vector<int> targs, qmatrix matrix);
void multiplyReferenceOperator(qmatrix& state, vector<int> targs, qmatrix matrix);


#endif // EVOLVE_HPP