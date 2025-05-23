/** @file
 * Examples of initialising a SuperOp in C11. Note that
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



void demo_createInlineSuperOp() {

    // inline literal without gross C99 compound-literal syntax (non-MSVC only)
    SuperOp a = createInlineSuperOp(1, {
        {1,2,3,4},
        {5,0.0000000006-(10E-11) * 3.14i,7,8},
        {9,10,11,12},
        {13,14,15,16+1.23i}
    });
    reportSuperOp(a);
    destroySuperOp(a);

    // unspecified elements default to 0 (C only)
    SuperOp b = createInlineSuperOp(3, {
        {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16},
        {5},
        {7,8,9}
    });
    reportSuperOp(b);
    destroySuperOp(b);
}


void demo_setInlineSuperOp() {

    // inline literal without gross C99 compound-literal syntax (non-MSVC only)
    SuperOp a = createSuperOp(1);
    setInlineSuperOp(a, 1, {
        {1,2,3,4},
        {5,0.0000000006-(10E-11) * 3.14i,7,8},
        {9,10,11,12},
        {13,14,15,16+1.23i}
    });
    reportSuperOp(a);
    destroySuperOp(a);

    // unspecified elements default to 0 (non-MSVC C only)
    SuperOp b = createSuperOp(3);
    setInlineSuperOp(b, 3, {
        {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16},
        {5},
        {7,8,9}
    });
    reportSuperOp(b);
    destroySuperOp(b);
}


void demo_setSuperOp() {

    // 2D compile-time array passed to VLA arg (non-MSVC C only)
    qcomp arr[4][4] = {
        {1,2,3,4},
        {5,6,7,8},
        {9,8,7,6},
        {5,4,3,2}
    };
    SuperOp a = createSuperOp(1);
    setSuperOp(a, arr);
    reportSuperOp(a);
    destroySuperOp(a);

    // 2D VLA (non-MSVC C only)
    int n = 2;
    int d = 1 << (2*n);
    qcomp elems[d][d];
    for (int i=0; i<d; i++)
        for (int j=0; j<d; j++)
            elems[i][j] = (i==j) * ((i+1)*1i);
    SuperOp b = createSuperOp(n);
    setSuperOp(b, elems);
    reportSuperOp(b);
    destroySuperOp(b);

    // nested pointers
    n = 3;
    d = 1 << (2*n);
    qcomp** ptrs = malloc(d * sizeof *ptrs);
    for (int i=0; i<d; i++) {
        ptrs[i] = malloc(d * sizeof **ptrs);
        for (int j=0; j<d; j++)
            ptrs[i][j] = (j==5)*(i + j*1i);
    }
    SuperOp c = createSuperOp(n);
    setSuperOp(c, ptrs);
    reportSuperOp(c);
    destroySuperOp(c);

    // array of pointers
    qcomp* ptrArr[1 << (2*2)];
    d = 1 << (2*2);
    for (int i=0; i<d; i++)
        ptrArr[i] = ptrs[i];
    SuperOp e = createSuperOp(n);
    setSuperOp(e, ptrArr);
    reportSuperOp(e);
    destroySuperOp(e);

    // inline C99 temporary array -> VLA (non-MSVC only)
    SuperOp f = createSuperOp(1);
    setSuperOp(f, (qcomp[4][4]) {
        {1,2,3,4},
        {1,2,3},
        {1,2},
        {1}
    });
    reportSuperOp(f);
    destroySuperOp(f);

    // cleanup
    for (int i=0; i<d; i++)
        free(ptrs[i]);
    free(ptrs);
}


void demo_syncSuperOp() {

    // manually modify the superoperator elements
    SuperOp a = createSuperOp(1);
    a.cpuElems[0][0] = 1;
    a.cpuElems[1][1] = 2i;
    a.cpuElems[2][2] = 3i;
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



// MSVC's naughty corner
#else
int main() { return 0; }
#endif
