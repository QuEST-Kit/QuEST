#ifndef COMPARE_HPP
#define COMPARE_HPP

#include "quest/include/quest.h"
#include "qvector.hpp"
#include "qmatrix.hpp"

#include <vector>
using std::vector;


void REQUIRE_AGREE( Qureg qureg, qvector reference );
void REQUIRE_AGREE( Qureg qureg, qmatrix reference );

void REQUIRE_AGREE( qvector reference, Qureg qureg );
void REQUIRE_AGREE( qmatrix reference, Qureg qureg );

void REQUIRE_AGREE( qreal apiScalar, qreal refScalar );
void REQUIRE_AGREE( qcomp apiScalar, qcomp refScalar );

void REQUIRE_AGREE( vector<qreal> apiList, vector<qreal> refList );


#endif // COMPARE_HPP