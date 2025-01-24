#include "qvector.hpp"
#include "qmatrix.hpp"
#include "linalg.hpp"
#include "macros.hpp"
#include "quest.h"

#include <type_traits>
using std::is_same_v;



/*
 * TO QUREG
 */


void setQureg(Qureg qureg, qvector vector) {
    DEMAND( !qureg.isDensityMatrix );
    DEMAND( qureg.numAmps == vector.size() );

    setQuregAmps(qureg, 0, vector.data(), vector.size());
}


void setQureg(Qureg qureg, qmatrix matrix) {
    DEMAND( qureg.isDensityMatrix );
    DEMAND( pow2(qureg.numQubits) == matrix.size() );

    qindex numRows = 1;
    qindex numCols = pow2(qureg.numQubits);

    for (size_t r=0; r<matrix.size(); r++) {
        qcomp* arr[] = {matrix[r].data()};
        setDensityQuregAmps(qureg, r, 0, arr, numRows, numCols);
    }
}



/*
 * FROM QUREG
 */


qvector getVector(Qureg qureg) {
    DEMAND( !qureg.isDensityMatrix );

    qvector out = getZeroVector(qureg.numAmps);
    getQuregAmps(out.data(), qureg, 0, qureg.numAmps);
    return out;
}


qmatrix getMatrix(Qureg qureg) {
    DEMAND( qureg.isDensityMatrix );

    qindex numRows = 1;
    qindex numCols = pow2(qureg.numQubits);
    qmatrix out = getZeroMatrix(numCols);

    for (size_t r=0; r<out.size(); r++) {
        qcomp* arr[] = {out[r].data()};
        getDensityQuregAmps(arr, qureg, r, 0, numRows, numCols);

    }

    return out;
}



/*
 * TO MATRIX
 */


template <typename T> 
qcomp getElem(T m, size_t r, size_t c) {

    if constexpr (is_same_v<T, CompMatr>)
        return m.cpuElems[r][c];
    
    if constexpr (is_same_v<T, CompMatr1> || is_same_v<T, CompMatr2>)
        return m.elems[r][c];

    if constexpr (is_same_v<T, DiagMatr>)
        return (r==c)? m.cpuElems[r] : 0;

    if constexpr (is_same_v<T, DiagMatr1> || is_same_v<T, DiagMatr2>)
        return (r==c)? m.elems[r] : 0;
}


template <typename T> 
qmatrix getMatrix(T m) {

    qmatrix out = getZeroMatrix(pow2(m.numQubits));

    for (size_t r=0; r<out.size(); r++)
        for (size_t c=0; c<out.size(); c++)
            out[r][c] = getElem<T>(m, r, c);

    return out;
}

template qmatrix getMatrix(CompMatr1);
template qmatrix getMatrix(CompMatr2);
template qmatrix getMatrix(CompMatr );
template qmatrix getMatrix(DiagMatr1);
template qmatrix getMatrix(DiagMatr2);
template qmatrix getMatrix(DiagMatr );

// FullStateDiagMatr excluded
