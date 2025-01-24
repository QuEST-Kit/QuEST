#include "qvector.hpp"
#include "qmatrix.hpp"
#include "convert.hpp"
#include "linalg.hpp"
#include "macros.hpp"

#include <complex>
using std::abs;


/*
 * VECTOR & VECTOR
 */

bool operator == (const qvector& v1, const qvector& v2) {
    DEMAND( v1.size() == v2.size() );

    for (size_t i=0; i<v1.size(); i++)
        if (abs(v1[i] - v2[i]) > TEST_EPSILON)
            return false;

    return true;
}


/*
 * MATRIX & MATRIX
 */

bool operator == (const qmatrix& m1, const qmatrix& m2) {
    DEMAND( m1.size() == m2.size() );

    for (size_t i=0; i<m1.size(); i++)
        for (size_t j=0; j<m1.size(); j++)
            if (abs(m1[i][j] - m2[i][j]) > TEST_EPSILON)
                return false;

    return true;
}


/*
 * VECTOR & QUREG
 */

bool operator == (const qvector& v, const Qureg& q) {
    DEMAND( !q.isDensityMatrix );
    DEMAND( q.numAmps == v.size() );

    return v == getVector(q);
}

bool operator == (const Qureg& q, const qvector& v) {
    return v == q;
}


/*
 * MATRIX & QUREG
 */

bool operator == (const qmatrix& m, const Qureg& q) {
    DEMAND( q.isDensityMatrix );
    DEMAND( pow2(q.numQubits) == m.size() );

    return m == getMatrix(q);
}

bool operator == (const Qureg& q, const qmatrix& m) {
    return m == q;
}


/*
 * QUREG & QUREG
 */

bool operator == (const Qureg& q1, const Qureg& q2) {
    DEMAND( q1.isDensityMatrix == q2.isDensityMatrix );
    DEMAND( q1.numQubits == q2.numQubits );

    return q1.isDensityMatrix?
        getMatrix(q1) == getMatrix(q2) :
        getVector(q1) == getVector(q2);
}
