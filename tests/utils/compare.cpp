#include "qvector.hpp"
#include "qmatrix.hpp"
#include "convert.hpp"
#include "linalg.hpp"
#include "macros.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <complex>
#include <vector>

using std::abs;
using std::vector;
using namespace Catch::Matchers;



/*
 * We compare a Qureg and its reference qvector/qmatrix
 * in a manner which avoids polluting the Catch2 console 
 * output, and which reports a single disgreeing amplitude
 * when the comparison fails.
 */


void REPORT_AND_FAIL( size_t index, qcomp amplitude, qcomp reference ) {
    CAPTURE( index, amplitude, reference );
    FAIL( );
}


void REQUIRE_AGREE( Qureg q, qvector v1 ) {
    DEMAND( !q.isDensityMatrix );
    DEMAND( q.numAmps == (qindex) v1.size() );

    qvector v2 = getVector(q);   

    for (size_t i=0; i<v1.size(); i++)
        if (abs(v1[i] - v2[i]) > TEST_EPSILON)
            REPORT_AND_FAIL(i, v2[i], v1[i]);

    SUCCEED( );
}

void REQUIRE_AGREE( qvector v1, Qureg q ) {
    REQUIRE_AGREE( q, v1 );
}


void REQUIRE_AGREE( Qureg q, qmatrix m1 ) {
    DEMAND( q.isDensityMatrix );
    DEMAND( getPow2(q.numQubits) == (qindex) m1.size() );

    qmatrix m2 = getMatrix(q);

    for (size_t i=0; i<m1.size(); i++)
        for (size_t j=0; j<m1.size(); j++)
            if (abs(m1[i][j] - m2[i][j]) > TEST_EPSILON)
                REPORT_AND_FAIL(j*m1.size()+i, m2[i][j], m1[i][j]);

    SUCCEED( );
}

void REQUIRE_AGREE( qmatrix m1, Qureg q ) {
    REQUIRE_AGREE( q, m1 );
}



/*
 * REAL AND COMPLEX SCALARS
 */


void REQUIRE_AGREE( qreal apiScalar, qreal refScalar ) {

    REQUIRE_THAT( apiScalar, WithinAbs(refScalar, TEST_EPSILON) );
}


void REQUIRE_AGREE( qcomp apiScalar, qcomp refScalar ) {

    REQUIRE_THAT( real(apiScalar), WithinAbs(real(refScalar), TEST_EPSILON) );
    REQUIRE_THAT( imag(apiScalar), WithinAbs(imag(refScalar), TEST_EPSILON) );
}



/*
 * LISTS
 */


void REQUIRE_AGREE( vector<qreal> apiList, vector<qreal> refList ) {
    DEMAND( apiList.size() == refList.size() );

    for (size_t i=0; i<apiList.size(); i++)
        REQUIRE_THAT( apiList[i], WithinAbs(refList[i], TEST_EPSILON) );
}
