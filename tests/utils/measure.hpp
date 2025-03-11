/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsmeasure Measure
 * @ingroup testutils
 * @brief
 * Testing utilities which evaluate measurements upon
 * reference qvector and qmatrix states. These are slow, 
 * serial, un-optimised, defensively-designed routines.
 * @{
 */

#ifndef MEASURE_HPP
#define MEASURE_HPP

#include "quest/include/quest.h"

#include "qvector.hpp"
#include "qmatrix.hpp"


qcomp getReferenceExpectationValue(qvector state, qmatrix observable);
qcomp getReferenceExpectationValue(qmatrix state, qmatrix observable);

qcomp getReferenceExpectationValue(qvector state, PauliStr str);
qcomp getReferenceExpectationValue(qmatrix state, PauliStr str);

qcomp getReferenceExpectationValue(qvector state, PauliStrSum sum);
qcomp getReferenceExpectationValue(qmatrix state, PauliStrSum sum);

qreal getReferenceProbability(qvector state);
qreal getReferenceProbability(qmatrix state);

qreal getReferenceProbability(qvector state, qindex basisIndex);
qreal getReferenceProbability(qmatrix state, qindex basisIndex);

qreal getReferenceProbability(qvector state, vector<int> targets, vector<int> outcomes);
qreal getReferenceProbability(qmatrix state, vector<int> targets, vector<int> outcomes);

vector<qreal> getAllReferenceProbabilities(qvector state, vector<int> targets);
vector<qreal> getAllReferenceProbabilities(qmatrix state, vector<int> targets);

qreal getReferencePurity(qvector state);
qreal getReferencePurity(qmatrix state);


#endif // MEASURE_HPP

/** @} (end defgroup) */
