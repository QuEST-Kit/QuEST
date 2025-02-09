#ifndef COMPARE_HPP
#define COMPARE_HPP

#include "quest.h"
#include "qvector.hpp"
#include "qmatrix.hpp"


void REQUIRE_AGREE( Qureg qureg, qvector reference );
void REQUIRE_AGREE( Qureg qureg, qmatrix reference );

void REQUIRE_AGREE( qvector reference, Qureg qureg );
void REQUIRE_AGREE( qmatrix reference, Qureg qureg );


#endif // COMPARE_HPP