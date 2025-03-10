/** @file
 * Testing utilities which compare scalars produced by the
 * QuEST API to those produced by other test utilities, and
 * Quregs modified by the API to qvector qmatrix references.
 *
 * @author Tyson Jones
 * 
 * @defgroup compare Compare
 * @ingroup utils
 * @{
 */

#ifndef COMPARE_HPP
#define COMPARE_HPP

#include "quest/include/quest.h"
#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;


bool doScalarsAgree(qcomp a, qcomp b);


void REQUIRE_AGREE( Qureg qureg, qvector reference );
void REQUIRE_AGREE( Qureg qureg, qmatrix reference );

void REQUIRE_AGREE( qvector reference, Qureg qureg );
void REQUIRE_AGREE( qmatrix reference, Qureg qureg );

void REQUIRE_AGREE( qreal apiScalar, qreal refScalar );
void REQUIRE_AGREE( qcomp apiScalar, qcomp refScalar );

void REQUIRE_AGREE( vector<qreal> apiList, vector<qreal> refList );


#endif // COMPARE_HPP

/** @} (end defgroup) */
