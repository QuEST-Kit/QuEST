/** @file
 * Examples of initialising matrices (CompMatr1,
 * CompMatr1, CompMatr, DiagMatr1, DiagMatr2, 
 * DiagMatr, FullStateDiagMatr) in C11. Note that
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



/*
 * CompMatr1, CompMatr2
 */


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
    qcomp** ptrs = malloc(dim * sizeof *ptrs);
    for (int i=0; i<dim; i++) {
        ptrs[i] = malloc(dim * sizeof **ptrs);
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



/*
 * CompMatr
 */


void demo_createInlineCompMatr() {

    // inline literal without gross C compound literal syntax (non-MSVC only)
    CompMatr a = createInlineCompMatr(2, {
        {1,2,3i,4},
        {4,5,6,7},
        {9,8,7,6},
        {1i,2i,0,0}
    });
    reportCompMatr(a);
    destroyCompMatr(a);

    // unspecified elements default to 0 (non-MSVC C only)
    CompMatr b = createInlineCompMatr(3, {
        {1,2,3,4,5,6,7,8},
        {8i, 7i, 6i, 5i},
        {9,9,9},
        {10}      
    });
    reportCompMatr(b);
    destroyCompMatr(b);
}


void demo_setInlineCompMatr() {

    // inline literal without gross C compound-literal syntax (non-MSVC only)
    CompMatr a = createCompMatr(1);
    setInlineCompMatr(a, 1, {{.3,.4},{.6,.7}});
    reportCompMatr(a);
    destroyCompMatr(a);

    // unspecified elements default to 0 (non-MSVC C only)
    CompMatr b = createCompMatr(3);
    setInlineCompMatr(b, 3, {
        {1,2,3,4,5,6,7,8},
        {8i, 7i, 6i, 5i},
        {9,9,9},
        {10}      
    });
    reportCompMatr(b);
    destroyCompMatr(b);
}


void demo_setCompMatr() {

    // 2D compile-time array passed to VLA arg (non-MSVC C only)
    qcomp arr[2][2] = {{5, 4},{3, 2}};
    CompMatr a = createCompMatr(1);
    setCompMatr(a, arr);
    reportCompMatr(a);
    destroyCompMatr(a);

    // 2D VLA (non-MSVC C only)
    int len = 2;
    qcomp elems[len][len];
    elems[0][0] = .1;
    elems[0][1] = 2i;
    elems[1][0] = -.7i;
    elems[1][1] = 1E-5;
    CompMatr b = createCompMatr(1);
    setCompMatr(b, elems);
    reportCompMatr(b);
    destroyCompMatr(b);

    // nested pointers
    int dim = 8;
    qcomp** ptrs = malloc(dim * sizeof *ptrs);
    for (int i=0; i<dim; i++) {
        ptrs[i] = malloc(dim * sizeof **ptrs);
        for (int j=0; j<dim; j++)
            ptrs[i][j] = i + j*1i;
    }
    CompMatr c = createCompMatr(3);
    setCompMatr(c, ptrs);
    reportCompMatr(c);
    destroyCompMatr(c);

    // array of pointers
    qcomp* ptrArr[dim];
    for (int i=0; i<dim; i++)
        ptrArr[i] = ptrs[i];
    CompMatr d = createCompMatr(3);
    setCompMatr(d, ptrArr);
    reportCompMatr(d);
    destroyCompMatr(d);

    // temporary array as compound literal (C only)
    CompMatr e = createCompMatr(1);
    setCompMatr(e, (qcomp[2][2]) {{3,2},{1,0}});
    reportCompMatr(e);
    destroyCompMatr(e);

    // cleanup
    for (int i=0; i<dim; i++)
        free(ptrs[i]);
    free(ptrs);
}


void demo_syncCompMatr() {

    CompMatr a = createCompMatr(2);

    // manually modify the elems
    a.cpuElems[0][0] = 1;
    a.cpuElems[1][1] = 2i;
    a.cpuElems[2][2] = -3i;
    a.cpuElems[3][3] = -2+4i;
    
    syncCompMatr(a);
    reportCompMatr(a);
    destroyCompMatr(a);
}



/*
 * DiagMatr1, DiagMatr2
 */


void demo_getInlineDiagMatr() {

    // inline literal without gross C compound-literal syntax
    DiagMatr1 a = getInlineDiagMatr1({.1, .2});
    reportDiagMatr1(a);

    // unspecified elements default to 0 (C only)
    DiagMatr2 b = getInlineDiagMatr2({1i});
    reportDiagMatr2(b);
}


void demo_getDiagMatr() {

    // compile-time array
    qcomp arr[2] = {5,6};
    DiagMatr1 a = getDiagMatr1(arr);
    reportDiagMatr1(a);

    // VLA (C only)
    int len = 4;
    qcomp elems[len];
    elems[0] = .1;
    elems[1] = 2i;
    elems[2] = -.7i;
    elems[3] = 1E-5;
    DiagMatr1 b = getDiagMatr1(elems);
    reportDiagMatr1(b);

    // pointer
    int dim = 4;
    qcomp* ptr = malloc(dim * sizeof *ptr);
    for (int i=0; i<dim; i++)
        ptr[i] = i + 1i*i;
    DiagMatr2 c = getDiagMatr2(ptr);
    reportDiagMatr2(c);

    // temporary array as compound literal (C only)
    DiagMatr2 d = getDiagMatr2( (qcomp[4]) {5, 6, 8}); // defaults 0
    reportDiagMatr2(d);

    // cleanup
    free(ptr);
}



/*
 * DiagMatr
 */


void demo_createInlineDiagMatr() {

    // inline literal without gross C compound-literal syntax (non-MSVC only)
    DiagMatr a = createInlineDiagMatr(1, {3i,5i});
    reportDiagMatr(a);
    destroyDiagMatr(a);

    // unspecified elemenrts default to 0 (non-MSVC C only)
    DiagMatr b = createInlineDiagMatr(4, {1, 2, 3});
    reportDiagMatr(b);
    destroyDiagMatr(b);
}


void demo_setInlineDiagMatr() {

    // inline literal without gross C compound-literal syntax
    DiagMatr a = createDiagMatr(2);
    setInlineDiagMatr(a, 2, {.7, .8, .9, .9});
    reportDiagMatr(a);
    destroyDiagMatr(a);

    // unspecified elements default to 0 (C only)
    DiagMatr b = createDiagMatr(3);
    setInlineDiagMatr(b, 3, {1, 2, 3, 4, -1});
    reportDiagMatr(b);
    destroyDiagMatr(b);
}


void demo_setDiagMatr() {

    // compile-time array passed to VLA arg (C only) 
    qcomp arr[2] = {5, 4};
    DiagMatr a = createDiagMatr(1);
    setDiagMatr(a, arr);
    reportDiagMatr(a);
    destroyDiagMatr(a);

    // VLA (C only)
    int len = 2;
    qcomp elems[len];
    elems[0] = .1;
    elems[1] = 2i;
    DiagMatr b = createDiagMatr(1);
    setDiagMatr(b, elems);
    reportDiagMatr(b);
    destroyDiagMatr(b);

    // heap pointer
    int dim = 1 << 3;
    qcomp* ptr = malloc(dim * sizeof *ptr);
    for (int i=0; i<dim; i++)
        ptr[i] = i + 1i*i;
    DiagMatr c = createDiagMatr(3);
    setDiagMatr(c, ptr);
    reportDiagMatr(c);
    destroyDiagMatr(c);

    // temporary array as compound literal (C only)
    DiagMatr d = createDiagMatr(1);
    setDiagMatr(d, (qcomp[2]) {9,8});
    reportDiagMatr(d);
    destroyDiagMatr(d);

    // cleanup
    free(ptr);
}


void demo_syncDiagMatr() {

    DiagMatr a = createDiagMatr(2);

    // manually modify the elems
    a.cpuElems[0] = 1;
    a.cpuElems[1] = 2i;
    a.cpuElems[2] = -3i;
    a.cpuElems[3] = -2+4i;

    syncDiagMatr(a);
    reportDiagMatr(a);
    destroyDiagMatr(a);
}



/*
 * FullStateDiagMatr
 */


void demo_setInlineFullStateDiagMatr() {

    // force use of available hardware acceleration
    QuESTEnv env = getQuESTEnv();
    FullStateDiagMatr matr = createCustomFullStateDiagMatr(5, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded);

    // inline literal; identical to setFullStateDiagMatr() for consistencty with C API
    setInlineFullStateDiagMatr(matr, 3, 10, {3,4,5,6,7,8,9,1,12,13});
    setInlineFullStateDiagMatr(matr, 15, 6, {6,5,4,3,2,1});
    
    reportFullStateDiagMatr(matr);
    destroyFullStateDiagMatr(matr);
}


void demo_setFullStateDiagMatr() {

    // force use of available hardware acceleration
    QuESTEnv env = getQuESTEnv();
    FullStateDiagMatr matr = createCustomFullStateDiagMatr(5, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded);

    // VLA (C only)
    int n = 10;
    qcomp vla[n];
    for (int i=0; i<n; i++)
        vla[i] =.123*i*1i;
    setFullStateDiagMatr(matr, 2, vla, n);

    // compile-time array
    qcomp arr[4] = {7i, 8i, 8, 7};
    setFullStateDiagMatr(matr, 15, arr, 4);

    // heap pointer
    int dim = 6;
    qcomp* ptr = (qcomp*) malloc(dim * sizeof *ptr);
    for (int i=0; i<dim; i++)
        ptr[i] = -i*.5i;
    setFullStateDiagMatr(matr, 25, ptr, 6);

    reportFullStateDiagMatr(matr);
    destroyFullStateDiagMatr(matr);
}


void demo_syncFullStateDiagMatr() {

    // force use of available hardware acceleration
    QuESTEnv env = getQuESTEnv();
    FullStateDiagMatr a = createCustomFullStateDiagMatr(5, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded);

    // every node modifies its first local element
    a.cpuElems[0] = -10i * (1+getQuESTEnv().rank);

    syncFullStateDiagMatr(a);
    reportFullStateDiagMatr(a);
    destroyFullStateDiagMatr(a);
}



/*
 * main
 */


int main() {
    
    initQuESTEnv();

    demo_getInlineCompMatr();
    demo_getCompMatr();

    demo_createInlineCompMatr();
    demo_setInlineCompMatr();
    demo_setCompMatr();
    demo_syncCompMatr();

    demo_getInlineDiagMatr();
    demo_getDiagMatr();
    
    demo_createInlineDiagMatr();
    demo_setInlineDiagMatr();
    demo_setDiagMatr();
    demo_syncDiagMatr();

    demo_setInlineFullStateDiagMatr();
    demo_setFullStateDiagMatr();
    demo_syncFullStateDiagMatr();

    finalizeQuESTEnv();
    return 0;
}



// MSVC's naughty corner
#else
int main() { return 0; }
#endif
