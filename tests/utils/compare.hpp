/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilscompare Compare
 * @ingroup testutils
 * @brief
 * Testing utilities which compare scalars produced by the
 * QuEST API to those produced by other test utilities, and
 * Quregs modified by the API to qvector qmatrix references.
 * @{
 */

#ifndef COMPARE_HPP
#define COMPARE_HPP

#include "quest/include/quest.h"
#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;


qreal getTestAbsoluteEpsilon();
qreal getTestRelativeEpsilon();

bool doScalarsAgree(qcomp a, qcomp b);
bool doMatricesAgree(qmatrix a, qmatrix b);

void REQUIRE_AGREE( Qureg qureg, Qureg other );

void REQUIRE_AGREE( Qureg qureg, qvector reference );
void REQUIRE_AGREE( Qureg qureg, qmatrix reference );

void REQUIRE_AGREE( qreal scalar, qreal reference );
void REQUIRE_AGREE( qcomp scalar, qcomp reference );

void REQUIRE_AGREE( vector<qreal> list, vector<qreal> reference );
void REQUIRE_AGREE( vector<qcomp> list, vector<qcomp> reference );

void REQUIRE_AGREE( qmatrix   matrix, qmatrix reference );
void REQUIRE_AGREE( CompMatr1 matrix, qmatrix reference );
void REQUIRE_AGREE( CompMatr2 matrix, qmatrix reference );
void REQUIRE_AGREE( CompMatr  matrix, qmatrix reference );
void REQUIRE_AGREE( DiagMatr1 matrix, qmatrix reference );
void REQUIRE_AGREE( DiagMatr2 matrix, qmatrix reference );
void REQUIRE_AGREE( DiagMatr  matrix, qmatrix reference );
void REQUIRE_AGREE( SuperOp   matrix, qmatrix reference );


#endif // COMPARE_HPP

/** @} (end defgroup) */
