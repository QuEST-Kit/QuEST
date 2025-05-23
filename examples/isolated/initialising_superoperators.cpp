/** @file
 * Examples of initialising a SuperOp in C++14. 
 * Note that when compiling QuEST in C++ with 'double' 
 * or 'long double' precision, you can freely use <complex>'s 
 * double-precision literals, e.g.
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



void demo_createInlineSuperOp() {

    // inline literal
    SuperOp a = createInlineSuperOp(1, {
        {1,2,3,4},
        {5,0.0000000006-(10E-11) * 3.14_i,7,8},
        {9,10,11,12},
        {13,14,15,16+1.23_i}
    });
    reportSuperOp(a);
    destroySuperOp(a);

    // existing vector (C++ only)
    vector<vector<qcomp>> elems(16, vector<qcomp>(16,5_i));
    SuperOp b = createInlineSuperOp(2, elems);
    reportSuperOp(b);
    destroySuperOp(b);
}


void demo_setInlineSuperOp() {

    // inline literal; achieves the same as setSuperOp(), and exists for consistencty with C API
    SuperOp a = createSuperOp(1);
    setInlineSuperOp(a, 1, {
        {1,2,3,4},
        {5,0.0000000006-(10E-11) * 3.14_i,7,8},
        {9,10,11,12},
        {13,14,15,16+1.23_i}
    });
    reportSuperOp(a);
    destroySuperOp(a); 
}


void demo_setSuperOp() {

    // inline literal (C++ only)
    SuperOp a = createSuperOp(1);
    setSuperOp(a, {
        {1,2,3,4},
        {5,0.0000000006-(10E-11) * 3.14_i,7,8},
        {9,10,11,12},
        {13,14,15,16+1.23_i}});
    reportSuperOp(a);
    destroySuperOp(a); 

    // 2D vector (C++ only)
    vector<vector<qcomp>> vec {
        {-9,-8, -8, -9},
        {7, 7, 6, 6},
        {0, -1, -2, -3},
        {-4_i, -5_i, 0, 0}
    };
    SuperOp b = createSuperOp(1);
    setSuperOp(b, vec);
    reportSuperOp(b);
    destroySuperOp(b);

    // nested pointers
    int n = 3;
    int d = 1 << (2*n);
    qcomp** ptrs = (qcomp**) malloc(d * sizeof *ptrs);
    for (int i=0; i<d; i++) {
        ptrs[i] = (qcomp*) malloc(d * sizeof **ptrs);
        for (int j=0; j<d; j++)
            ptrs[i][j] = (j==5)*(i + j*1_i);
    }
    SuperOp c = createSuperOp(n);
    setSuperOp(c, ptrs);
    reportSuperOp(c);
    destroySuperOp(c);

    // array of pointers (decays, so C++ supported too)
    qcomp* ptrArr[64];
    for (int i=0; i<d; i++)
        ptrArr[i] = ptrs[i];
    SuperOp e = createSuperOp(n);
    setSuperOp(e, ptrArr);
    reportSuperOp(e);
    destroySuperOp(e);

    // cleanup
    for (int i=0; i<d; i++)
        free(ptrs[i]);
    free(ptrs);
}


void demo_syncSuperOp() {

    // manually modify the superoperator elements
    SuperOp a = createSuperOp(1);
    a.cpuElems[0][0] = 1;
    a.cpuElems[1][1] = 2_i;
    a.cpuElems[2][2] = 3_i;
    syncSuperOp(a);
    reportSuperOp(a);
    destroySuperOp(a);
}


int main() {
    
    initQuESTEnv();

    demo_createInlineSuperOp();
    demo_setInlineSuperOp();
    demo_setSuperOp();
    demo_syncSuperOp();

    finalizeQuESTEnv();
    return 0;
}
