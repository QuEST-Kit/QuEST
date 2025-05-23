/** @file
 * Examples of initialising KrausMaps in C11. Note that
 * MSVC's C11 (which is already weird) doesn't support
 * assigning any non-complex literal to complex variables
 * nor any complex arithmetic operators, so it doesn't
 * get to play with the other children.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <stdlib.h>

#if !defined(_MSC_VER)



void demo_createInlineKrausMap() {

    // inline literal without C99 compound-literal syntax (non-MSVC only)
    KrausMap a = createInlineKrausMap(1, 3, {
        {{1,2},{3,4}},
        {{5,5},{6,6}},
        {{1i,2i},{-3i,-4i}}
    });
    reportKrausMap(a);
    destroyKrausMap(a);

    // unspecified elements/rows/matrices will be defaulted to all 0 (non-MSVC C only)
    KrausMap b = createInlineKrausMap(2, 5, {
        {
            {1,2,3,4},
            {1i,2i,3i},
            {4i,5i},
            {6}
        }, {
            {6,7,8,9},
            {9,8,7,6},
        }, {
            {1}
        }
    });
    reportKrausMap(b);
    destroyKrausMap(b);
}


void demo_setInlineKrausMap() {

    // inline literal without C99 compound-literal syntax (non-MSVC only)
    KrausMap a = createKrausMap(1, 3);
    setInlineKrausMap(a, 1, 3, {
        {{1,2},{3,4}},
        {{5,5},{6,6}},
        {{1i,2i},{-3i,-4i}}
    });
    reportKrausMap(a);
    destroyKrausMap(a);

    // unspecified elements/rows/matrices will be defaulted to all 0 (non-MSVC C only)
    KrausMap b = createKrausMap(2, 5);
    setInlineKrausMap(b, 2, 5, {
        {
            {1,2,3,4},
            {1i,2i,3i},
            {4i,5i},
            {6}
        }, {
            {6,7,8,9},
            {9,8,7,6},
        }, {
            {1}
        }
    });
    reportKrausMap(b);
    destroyKrausMap(b);
}


void demo_setKrausMap() {

    // 3D compile-time array passed to VLA arg (non-MSVC C only)
    qcomp arr[2][4][4] = {
        {
            {1,2,3,4},
            {5,6,7,8},
            {9,8,7,6},
            {5,4,3,2},
        }, {
            {1i,2i,3i},
            {5i}
        }
    };
    KrausMap a = createKrausMap(2, 2);
    setKrausMap(a, arr);
    reportKrausMap(a);
    destroyKrausMap(a);

    // 3D VLA (non-MSVC C only)
    int nQb = 2;
    int nOps = 3;
    int dim = 1 << nQb;
    qcomp elems[nOps][dim][dim];
    for (int n=0; n<nOps; n++)
        for (int r=0; r<dim; r++)
            for (int c=0; c<dim; c++)
                elems[n][r][c] = (1+n+r+c) * 1i;
    KrausMap b = createKrausMap(nQb, nOps);
    setKrausMap(b, elems);
    reportKrausMap(b);
    destroyKrausMap(b);

    // 3D nested pointers
    nQb = 2;
    nOps = 4;
    dim = 1 << nQb;
    qcomp*** ptrs = malloc(nOps * sizeof *ptrs);
    for (int n=0; n<nOps; n++) {
        ptrs[n] = malloc(dim * sizeof **ptrs);
        for (int i=0; i<dim; i++) {
            ptrs[n][i] = malloc(dim * sizeof ***ptrs);
            for (int j=0; j<dim; j++)
                ptrs[n][i][j] = (n==j) * (n+1) * -10i;
        }
    }
    KrausMap c = createKrausMap(nQb, nOps);
    setKrausMap(c, ptrs);
    reportKrausMap(c);
    destroyKrausMap(c);

    // 1D array of 2D pointers
    qcomp** ptrArr[4];
    for (int n=0; n<nOps; n++)
        ptrArr[n] = ptrs[n];
    KrausMap d = createKrausMap(nQb, nOps);
    setKrausMap(d, ptrArr);
    reportKrausMap(d);
    destroyKrausMap(d);

    // inline C99 temporary array -> VLA (non-MSVC C only)
    KrausMap e = createKrausMap(1, 3);
    setKrausMap(e, (qcomp[3][2][2]) {
        {{1,2},{3,4}},
        {{5,6},{7,8}},
        {{9}}
    });
    reportKrausMap(e);
    destroyKrausMap(e);

    // cleanup
    for (int n=0; n<nOps; n++) {
        for (int i=0; i<dim; i++)
            free(ptrs[n][i]);
        free(ptrs[n]);
    }
    free(ptrs);
}


void demo_syncKrausMap() {

    KrausMap a = createKrausMap(2, 4);

    // manually modify elements
    for (int n=0; n<a.numMatrices; n++)
        for (int r=0; r<a.numRows; r++)
            for (int c=0; c<a.numRows; c++)
                a.matrices[n][r][c] = (n+r+c+1)*3i;

    syncKrausMap(a);
    reportKrausMap(a);
    destroyKrausMap(a);
}


int main() {
    
    initQuESTEnv();

    demo_createInlineKrausMap();
    demo_setInlineKrausMap();
    demo_setKrausMap();
    demo_syncKrausMap();

    finalizeQuESTEnv();
    return 0;
}



// MSVC's naughty corner
#else
int main() { return 0; }
#endif
