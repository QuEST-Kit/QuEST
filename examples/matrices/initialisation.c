#include "quest.h"
#include <stdlib.h>


void demo_getInlineCompMatr() {

    // inline literal without gross C compound-literal syntax
    CompMatr1 a = getInlineCompMatr1( {{1,2i},{3i+.1,-4}} );
    reportCompMatr1(a);

    // unspecified elements default to 0 (C only)
    CompMatr2 b = getInlineCompMatr2({
        {1, 2, 3, 4},
        {5, 6},
        {10}
    });
    reportCompMatr2(b);
}


void demo_getCompMatr() {

    // 2D compile-time array
    qcomp arr[2][2] = {{5, 4},{3, 2}};
    CompMatr1 a = getCompMatr1(arr);
    reportCompMatr1(a);

    // 2D VLA (C only)
    int len = 2;
    qcomp elems[len][len];
    elems[0][0] = .1;
    elems[0][1] = 2i;
    elems[1][0] = -.7i;
    elems[1][1] = 1E-5;
    CompMatr1 b = getCompMatr1(elems);
    reportCompMatr1(b);

    // nested pointers
    int dim = 4;
    qcomp** ptrs = malloc(dim * sizeof(qcomp));
    for (int i=0; i<dim; i++) {
        ptrs[i] = malloc(dim * sizeof(qcomp));
        for (int j=0; j<dim; j++)
            ptrs[i][j] = i + j*1i;
    }
    CompMatr2 c = getCompMatr2(ptrs);
    reportCompMatr2(c);

    // array of pointers
    qcomp* ptrArr[dim];
    for (int i=0; i<dim; i++)
        ptrArr[i] = ptrs[i];
    CompMatr2 d = getCompMatr2(ptrArr);
    reportCompMatr2(d);

    // temporary array as compound literal (C only)
    CompMatr2 e = getCompMatr2( (qcomp[4][4]) {{42}} );
    reportCompMatr2(e);

    // cleanup
    for (int i=0; i<dim; i++)
        free(ptrs[i]);
    free(ptrs);
}


void demo_setInlineCompMatr() {

    // inline literal without gross C compound-literal syntax
    CompMatrN a = createCompMatrN(1);
    setInlineCompMatrN(a, 1, {{.3,.4},{.6,.7}});
    reportCompMatrN(a);
    destroyCompMatrN(a);

    // unspecified elements default to 0 (C only)
    CompMatrN b = createCompMatrN(3);
    setInlineCompMatrN(b, 3, {
        {1,2,3,4,5,6,7,8},
        {8i, 7i, 6i, 5i},
        {9,9,9},
        {10}      
    });
    reportCompMatrN(b);
    destroyCompMatrN(b);
}


void demo_setCompMatr() {

    // 2D compile-time array passed to VLA arg (C only) 
    qcomp arr[2][2] = {{5, 4},{3, 2}};
    CompMatrN a = createCompMatrN(1);
    setCompMatrN(a, arr);
    reportCompMatrN(a);
    destroyCompMatrN(a);

    // 2D VLA (C only)
    int len = 2;
    qcomp elems[len][len];
    elems[0][0] = .1;
    elems[0][1] = 2i;
    elems[1][0] = -.7i;
    elems[1][1] = 1E-5;
    CompMatrN b = createCompMatrN(1);
    setCompMatrN(b, elems);
    reportCompMatrN(b);
    destroyCompMatrN(b);

    // nested pointers
    int dim = 8;
    qcomp** ptrs = malloc(dim * sizeof(qcomp));
    for (int i=0; i<dim; i++) {
        ptrs[i] = malloc(dim * sizeof(qcomp));
        for (int j=0; j<dim; j++)
            ptrs[i][j] = i + j*1i;
    }
    CompMatrN c = createCompMatrN(3);
    setCompMatrN(c, ptrs);
    reportCompMatrN(c);
    destroyCompMatrN(c);

    // array of pointers
    qcomp* ptrArr[dim];
    for (int i=0; i<dim; i++)
        ptrArr[i] = ptrs[i];
    CompMatrN d = createCompMatrN(3);
    setCompMatrN(d, ptrArr);
    reportCompMatrN(d);
    destroyCompMatrN(d);

    // temporary array as compound literal (C only)
    CompMatrN e = createCompMatrN(1);
    setCompMatrN(e, (qcomp[2][2]) {{3,2},{1,0}});
    reportCompMatrN(e);
    destroyCompMatrN(e);

    // cleanup
    for (int i=0; i<dim; i++)
        free(ptrs[i]);
    free(ptrs);
}


int main() {
    
    initQuESTEnv();

    demo_getInlineCompMatr();
    demo_getCompMatr();
    demo_setInlineCompMatr();
    demo_setCompMatr();

    finalizeQuESTEnv();
    return 0;
}
