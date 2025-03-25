/** @file
 * Testing utilities which compare scalars produced by the
 * QuEST API to those produced by other test utilities, and
 * Quregs modified by the API to qvector qmatrix references.
 *
 * @author Tyson Jones
 */

#include "qvector.hpp"
#include "qmatrix.hpp"
#include "convert.hpp"
#include "linalg.hpp"
#include "macros.hpp"

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <complex>
#include <vector>
#include <algorithm>

using std::abs;
using std::vector;
using namespace Catch::Matchers;



/*
 * Maximum tolerated difference between API and reference
 * scalar in order to be considered equivalent.
 */


#if FLOAT_PRECISION == 1
    const qreal ABSOLUTE_EPSILON = 1E-2;
    const qreal RELATIVE_EPSILON = 1E-2;
#elif FLOAT_PRECISION == 2
    const qreal ABSOLUTE_EPSILON = 1E-8;
    const qreal RELATIVE_EPSILON = 1E-8;
#elif FLOAT_PRECISION == 4
    const qreal ABSOLUTE_EPSILON = 1E-10;
    const qreal RELATIVE_EPSILON = 1E-10;
#endif


qreal getTestAbsoluteEpsilon() {

    return ABSOLUTE_EPSILON;
}

qreal getTestRelativeEpsilon() {

    return RELATIVE_EPSILON;
}


qreal getAbsDif(qcomp a, qcomp b) {

    return std::abs(a - b);
}

qreal getRelDif(qcomp a, qcomp b) {

    qreal denom = std::min({std::abs(a), std::abs(b)});
    return getAbsDif(a,b) / denom;
}


bool doScalarsAgree(qcomp a, qcomp b) {

    // permit absolute OR relative agreement

    if (getAbsDif(a, b) <= ABSOLUTE_EPSILON)
        return true;

   return (getRelDif(a, b) <= RELATIVE_EPSILON);
}

bool doMatricesAgree(qmatrix a, qmatrix b) {
    DEMAND( a.size() == b.size() );

    // assumed square and equal-size
    size_t dim = a.size();

    for (size_t i=0; i<dim; i++)
        for (size_t j=0; j<dim; j++)
            if (!doScalarsAgree(a[i][j], b[i][j]))
                return false;

    return true;
}



/*
 * We compare a Qureg and its reference qvector/qmatrix
 * in a manner which avoids polluting the Catch2 console 
 * output, and which reports a single disgreeing amplitude
 * when the comparison fails.
 */


void REPORT_AMP_AND_FAIL( size_t index, qcomp amplitude, qcomp reference ) {
    qreal absolute_difference = getAbsDif(amplitude, reference);
    qreal relative_difference = getRelDif(amplitude, reference);
    CAPTURE( 
        index, amplitude, reference, 
        absolute_difference, ABSOLUTE_EPSILON,
        relative_difference, RELATIVE_EPSILON
    );
    FAIL( );
}


void REQUIRE_AGREE( Qureg q, qvector v1 ) {
    DEMAND( !q.isDensityMatrix );
    DEMAND( q.numAmps == (qindex) v1.size() );

    qvector v2 = getVector(q);   

    for (size_t i=0; i<v1.size(); i++)
        if (!doScalarsAgree(v1[i], v2[i]))
            REPORT_AMP_AND_FAIL(i, v2[i], v1[i]);

    SUCCEED( );
}


void REQUIRE_AGREE( Qureg q, qmatrix m1 ) {
    DEMAND( q.isDensityMatrix );
    DEMAND( getPow2(q.numQubits) == (qindex) m1.size() );

    qmatrix m2 = getMatrix(q);

    for (size_t i=0; i<m1.size(); i++)
        for (size_t j=0; j<m1.size(); j++)
            if (!doScalarsAgree(m1[i][j], m2[i][j]))
                REPORT_AMP_AND_FAIL(j*m1.size()+i, m2[i][j], m1[i][j]);

    SUCCEED( );
}



/*
 * REAL AND COMPLEX SCALARS
 */


void REPORT_SCALAR_AND_FAIL( qcomp scalar, qcomp reference ) {
    qreal absolute_difference = getAbsDif(scalar, reference);
    qreal relative_difference = getRelDif(scalar, reference);
    CAPTURE( 
        scalar, reference, 
        absolute_difference, ABSOLUTE_EPSILON,
        relative_difference, RELATIVE_EPSILON
    );
    FAIL( );
}

void REPORT_SCALAR_AND_FAIL( qreal scalar, qreal reference ) {

    // like above but does not display redundant imag-components
    qreal absolute_difference = getAbsDif(qcomp(scalar,0), qcomp(reference,0));
    qreal relative_difference = getRelDif(qcomp(scalar,0), qcomp(reference,0));
    CAPTURE( 
        scalar, reference, 
        absolute_difference, ABSOLUTE_EPSILON,
        relative_difference, RELATIVE_EPSILON
    );
    FAIL( );
}


void REQUIRE_AGREE( qcomp scalar, qcomp reference ) {

    if (!doScalarsAgree(scalar, reference))
        REPORT_SCALAR_AND_FAIL(scalar, reference);

    SUCCEED( );
}

void REQUIRE_AGREE( qreal scalar, qreal reference ) {

    if (!doScalarsAgree(qcomp(scalar,0), qcomp(reference,0)))
        REPORT_SCALAR_AND_FAIL(scalar, reference);
    
    SUCCEED( );
}



/*
 * LISTS
 */


void REQUIRE_AGREE( vector<qreal> apiList, vector<qreal> refList ) {
    DEMAND( apiList.size() == refList.size() );

    for (size_t i=0; i<apiList.size(); i++)
        if (!doScalarsAgree(apiList[i], refList[i]))
            REPORT_SCALAR_AND_FAIL(apiList[i], refList[i]);

    SUCCEED( );
}


void REQUIRE_AGREE( vector<qcomp> apiList, vector<qcomp> refList ) {
    DEMAND( apiList.size() == refList.size() );

    for (size_t i=0; i<apiList.size(); i++)
        if (!doScalarsAgree(apiList[i], refList[i]))
            REPORT_SCALAR_AND_FAIL(apiList[i], refList[i]);

    SUCCEED( );
}



/*
 * MATRICES
 */


void REPORT_ELEM_AND_FAIL( size_t row, size_t col, qcomp elem, qcomp reference ) {
    qreal absolute_difference = getAbsDif(elem, reference);
    qreal relative_difference = getRelDif(elem, reference);
    CAPTURE( 
        row, col, elem, reference, 
        absolute_difference, ABSOLUTE_EPSILON,
        relative_difference, RELATIVE_EPSILON
    );
    FAIL( );
}


void REQUIRE_AGREE( qmatrix matrix, qmatrix reference ) {
    DEMAND( matrix.size() == reference.size() );

    size_t dim = matrix.size();

    for (size_t r=0; r<dim; r++)
        for (size_t c=0; c<dim; c++)
            if (!doScalarsAgree(matrix[r][c], reference[r][c]))
                REPORT_ELEM_AND_FAIL(r, c, matrix[r][c], reference[r][c]);

    SUCCEED( );
}


void REQUIRE_AGREE( CompMatr1 matrix, qmatrix reference ) { REQUIRE_AGREE( getMatrix(matrix), reference ); }
void REQUIRE_AGREE( CompMatr2 matrix, qmatrix reference ) { REQUIRE_AGREE( getMatrix(matrix), reference ); }
void REQUIRE_AGREE( CompMatr  matrix, qmatrix reference ) { REQUIRE_AGREE( getMatrix(matrix), reference ); }
void REQUIRE_AGREE( DiagMatr1 matrix, qmatrix reference ) { REQUIRE_AGREE( getMatrix(matrix), reference ); }
void REQUIRE_AGREE( DiagMatr2 matrix, qmatrix reference ) { REQUIRE_AGREE( getMatrix(matrix), reference ); }
void REQUIRE_AGREE( DiagMatr  matrix, qmatrix reference ) { REQUIRE_AGREE( getMatrix(matrix), reference ); }
void REQUIRE_AGREE( SuperOp   matrix, qmatrix reference ) { REQUIRE_AGREE( getMatrix(matrix), reference ); }



/*
 * we lazily compare quregs by converting one into qvector
 * or qmatrix, rather than an optimised direct comparison,
 * so that we can compare the states of distinct deployments
 */


void REQUIRE_AGREE( Qureg qureg, Qureg other ) {
    DEMAND( qureg.numQubits == other.numQubits );
    DEMAND( qureg.isDensityMatrix == other.isDensityMatrix );

    return (qureg.isDensityMatrix)?
        REQUIRE_AGREE( qureg, getMatrix(other) ):
        REQUIRE_AGREE( qureg, getVector(other) );
}
