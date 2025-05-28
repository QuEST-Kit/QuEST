/** @file
 * Examples of initialising matrices (CompMatr1, CompMatr1, 
 * CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr) 
 * in C++14. Note that when compiling QuEST in C++ with 'double' 
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



/*
 * CompMatr1, CompMatr2
 */


void demo_getInlineCompMatr() {

    // inline literal; identical to getCompMatr1() for consistencty with C API
    CompMatr1 a = getInlineCompMatr1( {{1,2_i},{3_i+.1,-4}} );
    reportCompMatr1(a);

    // we must specify all elements (only necessary in C++), else runtime validation triggered
    CompMatr2 b = getInlineCompMatr2({
        {1, 2, 3, 4},
        {5, 6, 8, 8},
        {5, 7, 6, 1},
        {4, 4, 3, 3}
    });
    reportCompMatr2(b);
}


void demo_getCompMatr() {

    // inline literal (C++ only)
    CompMatr1 a = getCompMatr1( {{1,2_i},{3_i+.1,-4}} );
    reportCompMatr1(a);

    // 2D vector (C++ only)
    vector<vector<qcomp>> vec {{-9,-8},{-7,-6}};
    CompMatr1 b = getCompMatr1(vec);
    reportCompMatr1(b);

    // 2D compile-time array
    qcomp arr[2][2] = {{5, 4},{3, 2}};
    CompMatr1 c = getCompMatr1(arr);
    reportCompMatr1(c);

    // nested pointers
    int dim = 4;
    qcomp** ptrs = (qcomp**) malloc(dim * sizeof *ptrs);
    for (int i=0; i<dim; i++) {
        ptrs[i] =(qcomp*) malloc(dim * sizeof **ptrs);
        for (int j=0; j<dim; j++)
            ptrs[i][j] = i + j*1i;
    }
    CompMatr2 d = getCompMatr2(ptrs);
    reportCompMatr2(d);

    // array of pointers
    qcomp* ptrArr[4];
    for (int i=0; i<dim; i++)
        ptrArr[i] = ptrs[i];
    CompMatr2 e = getCompMatr2(ptrArr);
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

    // inline literal
    CompMatr a = createInlineCompMatr(2, {
        {1,2,3_i,4},
        {4,5,6,7},
        {9,8,7,6},
        {1_i,2_i,0,0}
    });
    reportCompMatr(a);
    destroyCompMatr(a);

    // existing vector (C++ only)
    vector<vector<qcomp>> elems = {
        {1,2},
        {3,4}
    };
    CompMatr b = createInlineCompMatr(1, elems);
    reportCompMatr(b);
    destroyCompMatr(b);

    // must specify every element (unlike in C) otherwise runtime validation is triggered
}


void demo_setInlineCompMatr() {

    // inline literal; identical to setCompMatr() for consistencty with C API
    CompMatr a = createCompMatr(1);
    setInlineCompMatr(a, 1, {{.3,.4},{.6,.7}});
    reportCompMatr(a);
    destroyCompMatr(a);

    // we must specify all elements (only necessary in C++)
    CompMatr b = createCompMatr(3);
    setInlineCompMatr(b, 3, {
        {1,2,3,4,5,6,7,8},
        {8_i,7_i,6_i,5_i,0,0,0,0},
        {0,0,0,0,9,9,9,9},
        {8,7,6,5,4,3,2,1},
        {3,3,3,3,3,3,3,3},
        {0,0,0,0,1_i,1_i,1_i,1_i},
        {9,8,7,6,5,4,3,2},
        {-1,-2,-3,-4,4,3,2,1}    
    });
    reportCompMatr(b);
    destroyCompMatr(b);
}


void demo_setCompMatr() {

    // inline literal (C++ only)
    CompMatr a = createCompMatr(1);
    setCompMatr(a, {{1,2_i},{3_i+.1,-4}});
    reportCompMatr(a);
    destroyCompMatr(a);

    // 2D vector (C++ only)
    vector<vector<qcomp>> vec {
        {-9,-8, -8, -9},
        {7, 7, 6, 6},
        {0, -1, -2, -3},
        {-4_i, -5_i, 0, 0}
    };
    CompMatr b = createCompMatr(2);
    setCompMatr(b, vec);
    reportCompMatr(b);
    destroyCompMatr(b);

    // nested pointers
    int dim = 8;
    qcomp** ptrs = (qcomp**) malloc(dim * sizeof *ptrs);
    for (int i=0; i<dim; i++) {
        ptrs[i] = (qcomp*) malloc(dim * sizeof **ptrs);
        for (int j=0; j<dim; j++)
            ptrs[i][j] = i + j*1i;
    }
    CompMatr c = createCompMatr(3);
    setCompMatr(c, ptrs);
    reportCompMatr(c);
    destroyCompMatr(c);

    // array of pointers (decays, so C++ supported)
    qcomp* ptrArr[8];
    for (int i=0; i<dim; i++)
        ptrArr[i] = ptrs[i];
    CompMatr d = createCompMatr(3);
    setCompMatr(d, ptrArr);
    reportCompMatr(d);
    destroyCompMatr(d);

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
    a.cpuElems[3][3] = -2+4_i;
    
    syncCompMatr(a);
    reportCompMatr(a);
    destroyCompMatr(a);
}



/*
 * DiagMatr1, DiagMatr2
 */

void demo_getInlineDiagMatr() {

    // inline literal; identical to getDiagMatr1() for consistencty with C API
    DiagMatr1 a = getInlineDiagMatr1( {9,8} );
    reportDiagMatr1(a);

    // we must specify all elements (only necessary in C++)
    DiagMatr2 b = getInlineDiagMatr2({1_i, 2_i, 9_i, 8_i});
    reportDiagMatr2(b);
}


void demo_getDiagMatr() {

    // inline literal (C++ only)
    DiagMatr1 a = getDiagMatr1({-5_i, -10_i});
    reportDiagMatr1(a);

    // vector (C++ only)
    vector<qcomp> vec {20_i, 30_i};
    DiagMatr1 b = getDiagMatr1(vec);
    reportDiagMatr1(b);

    // compile-time array
    qcomp arr[4] = {7_i, 8_i, 8, 7};
    DiagMatr1 c = getDiagMatr1(arr);
    reportDiagMatr1(c);

    // heap pointer
    int dim = 4;
    qcomp* ptr = (qcomp*) malloc(dim * sizeof *ptr);
    for (int i=0; i<dim; i++)
        ptr[i] = 10_i*i;
    DiagMatr2 d = getDiagMatr2(ptr);
    reportDiagMatr2(d);

    // cleanup
    free(ptr);
}


void demo_setInlineDiagMatr() {

    // inline literal; identical to setDiagMatr() for consistencty with C API
    DiagMatr a = createDiagMatr(1);
    setInlineDiagMatr(a, 1, {333, 444});
    reportDiagMatr(a);
    destroyDiagMatr(a);

    // we must specify all elements (only necessary in C++)
    DiagMatr b = createDiagMatr(3);
    setInlineDiagMatr(b, 3, {4,5,4,5,6,7,6,7_i});
    reportDiagMatr(b);
    destroyDiagMatr(b);
}



/*
 * DiagMatr
 */


void demo_createInlineDiagMatr() {

    // inline literal
    DiagMatr a = createInlineDiagMatr(2, {1,2,3_i,4});
    reportDiagMatr(a);
    destroyDiagMatr(a);

    // existing vector (C++ only)
    vector<qcomp> elems = {1,2,3,4,5,6,7,8};
    DiagMatr b = createInlineDiagMatr(3, elems);
    reportDiagMatr(b);
    destroyDiagMatr(b);

    // must specify every element (unlike in C) otherwise runtime validation is triggered
}


void demo_setDiagMatr() {

    // inline literal (C++ only)
    DiagMatr a = createDiagMatr(1);
    setDiagMatr(a, {6_i,5_i});
    reportDiagMatr(a);
    destroyDiagMatr(a);

    // vector (C++ only)
    vector<qcomp> vec {11, 22, 33, 44};
    DiagMatr b = createDiagMatr(2);
    setDiagMatr(b, vec);
    reportDiagMatr(b);
    destroyDiagMatr(b);

    // compile-time array
    qcomp arr[4] = {7_i, 8_i, 8, 7};
    DiagMatr c = createDiagMatr(2);
    setDiagMatr(c, arr);
    reportDiagMatr(c);
    destroyDiagMatr(c);

    // heap pointer
    int dim = 8;
    qcomp* ptr = (qcomp*) malloc(dim * sizeof *ptr);
    for (int i=0; i<dim; i++)
        ptr[i] = -i*.5_i;
    DiagMatr d = createDiagMatr(3);
    setDiagMatr(d, ptr);
    reportDiagMatr(d);
    destroyDiagMatr(d);

    // cleanup
    free(ptr);
}


void demo_syncDiagMatr() {

    DiagMatr a = createDiagMatr(2);

    // manually modify the elems
    a.cpuElems[0] = 1;
    a.cpuElems[1] = 2_i;
    a.cpuElems[2] = -3_i;
    a.cpuElems[3] = -2+4_i;

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

    // inline literal (C++ only)
    setFullStateDiagMatr(matr, 0, {1,2,3,4});
    setFullStateDiagMatr(matr, 4, {5,6,7,8});

    // vector (C++ only)
    vector<qcomp> vec {11, 22, 33, 44};
    setFullStateDiagMatr(matr, 7, vec);

    // compile-time array must pass size
    qcomp arr[4] = {7_i, 8_i, 8, 7};
    setFullStateDiagMatr(matr, 15, arr, 4);

    // heap pointer
    int dim = 6;
    qcomp* ptr = (qcomp*) malloc(dim * sizeof *ptr);
    for (int i=0; i<dim; i++)
        ptr[i] = -i*.5_i;
    setFullStateDiagMatr(matr, 25, ptr, 6);

    reportFullStateDiagMatr(matr);
    destroyFullStateDiagMatr(matr);
}


void demo_syncFullStateDiagMatr() {

    // force use of available hardware acceleration
    QuESTEnv env = getQuESTEnv();
    FullStateDiagMatr a = createCustomFullStateDiagMatr(5, env.isDistributed, env.isGpuAccelerated, env.isMultithreaded);

    // every node modifies its first local element
    a.cpuElems[0] = -10_i * (1+getQuESTEnv().rank);

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
