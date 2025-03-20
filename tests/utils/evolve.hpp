/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsevolve Evolve
 * @ingroup testutils
 * @brief
 * Testing utilities which evolve a reference state
 * (qvector or qmatrix) under the action of a 
 * reference operation. These are slow, serial,
 * un-optimised, defensively-designed routines.
 * @{
 */

#ifndef EVOLVE_HPP
#define EVOLVE_HPP

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;


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

void applyReferenceOperator(   qvector& state, qmatrix matrix);
void applyReferenceOperator(   qmatrix& state, qmatrix matrix);
void multiplyReferenceOperator(qvector& state, qmatrix matrix);
void multiplyReferenceOperator(qmatrix& state, qmatrix matrix);

void applyReferenceOperator(qmatrix& state, vector<int> targs, vector<qmatrix> matrices);


#endif // EVOLVE_HPP

/** @} (end defgroup) */
