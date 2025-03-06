/** @file
 * Examples of initialising KrausMaps in C++14. Note that
 * when compiling QuEST in C++ with 'double' or 'long double'
 * precision, you can freely use <complex>'s double-precision 
 * literals, e.g.
 * 
 *      qcomp x = 1.3i
 * 
 * However, such literals are forbidden when compiling in 'single' 
 * precision; the compiler will not allow the double-precision 
 * literal to be autocast to a float. If you want to use 
 * precision-agnostic complex literals, you can use "_i", e.g.
 * 
 *      qcomp x = 1.3_i
 * 
 * This file will use the latter literal.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <stdlib.h>
#include <vector>

using std::vector;



void demo_createInlineKrausMap() {

    // inline literal
    KrausMap a = createInlineKrausMap(1, 3, {
        {{1,2},{3,4}},
        {{5,5},{6,6}},
        {{1_i,2_i},{-3_i,-4_i}}
    });
    reportKrausMap(a);
    destroyKrausMap(a);

    // existing vector (C++ only)
    vector<vector<vector<qcomp>>> elems = {
        {{1,2},{3,4}},
        {{5,5},{6,6}},
        {{1_i,2_i},{-3_i,-4_i}}
    };
    KrausMap b = createInlineKrausMap(1, 3, elems);
    reportKrausMap(b);
    destroyKrausMap(b);

    // must specify every element (unlike in C) otherwise runtime validation is triggered
}


void demo_setInlineKrausMap() {

    // inline literal; equivalent to to setKrausMap() but needed for consistencty with C API
    KrausMap a = createKrausMap(1, 3);
    setInlineKrausMap(a, 1, 3, {
        {{1,2_i},{3,4}},
        {{5,5},{6_i,6}},
        {{1_i,2_i},{-3_i,-4_i}}
    });
    reportKrausMap(a);
    destroyKrausMap(a);

    // must specify every element (unlike in C) otherwise runtime validation is triggered
}


void demo_setKrausMap() {

    // inline literal (C++ only)
    KrausMap a = createKrausMap(1, 3);
    setKrausMap(a, {
        {{1,2},{3,4}},
        {{5,5},{6,6}},
        {{1_i,2_i},{-3_i,-4_i}}
    });
    reportKrausMap(a);
    destroyKrausMap(a);

    // 3D vector (C++ only)
    vector<vector<vector<qcomp>>> vec {
        {
            {-9,-8, -8, -9},
            {7, 7, 6, 6},
            {0, -1, -2, -3},
            {-4_i, -5_i, 0, 0}
        }, {
            {1,2,3,4},
            {5,6,7,8},
            {1,1,1,1},
            {2,2,2,2}
        }
    };
    KrausMap b = createKrausMap(2, 2);
    setKrausMap(b, vec);
    reportKrausMap(b);
    destroyKrausMap(b);

    // 3D nested pointers
    int nQb = 2;
    int nOps = 4;
    int dim = 1 << nQb;
    qcomp*** ptrs = (qcomp***) malloc(nOps * sizeof *ptrs);
    for (int n=0; n<nOps; n++) {
        ptrs[n] = (qcomp**) malloc(dim * sizeof **ptrs);
        for (int i=0; i<dim; i++) {
            ptrs[n][i] = (qcomp*) malloc(dim * sizeof ***ptrs);
            for (int j=0; j<dim; j++)
                ptrs[n][i][j] = (n==j) * (n+1) * -10_i;
        }
    }
    KrausMap c = createKrausMap(nQb, nOps);
    setKrausMap(c, ptrs);
    reportKrausMap(c);
    destroyKrausMap(c);

    // 1D array of 2D pointers (decays, so C++ supported)
    qcomp** ptrArr[4];
    for (int n=0; n<nOps; n++)
        ptrArr[n] = ptrs[n];
    KrausMap d = createKrausMap(nQb, nOps);
    setKrausMap(d, ptrArr);
    reportKrausMap(d);
    destroyKrausMap(d);

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
                a.matrices[n][r][c] = (n+r+c+1)*3_i;

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
