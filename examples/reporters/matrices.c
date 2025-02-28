#include "quest.h"

#include <stdio.h>



// MSVC's C11 (which is already weird) doesn't support
// assigning any non-complex literal to complex variables
// nor any complex arithmetic operators, so it doesn't
// get to play with the other children.
#if !defined(_MSC_VER)



/*
 * This demo is most interesting when run distributed with up to 64 nodes
 */



void rootPrint(qindex num) {

    if (getQuESTEnv().rank != 0)
        return;

    printf("\n\n>> set number of reported items to %lld", num);

    if (num == 0)
        printf(" (which means report all)");
         
    printf("\n\n");
}



void demo_CompMatr() {

    CompMatr matr = createCompMatr(4);
    for (qindex r=0; r<matr.numRows; r++)
        for (qindex c=0; c<matr.numRows; c++)
            matr.cpuElems[r][c] = r*matr.numRows + c + 1;
    syncCompMatr(matr);

    int len = 4;
    qindex numReportElems[] = {0, 12, 5, 1};

    for (int i=0; i<len; i++) {
        qindex num = numReportElems[i];
        rootPrint(num);
        setMaxNumReportedItems(num, num);
        reportCompMatr(matr);
    }
}



void demo_DiagMatr() {

    DiagMatr matr = createDiagMatr(5);
    for (qindex i=0; i<matr.numElems; i++)
        matr.cpuElems[i] = i + 1;
    syncDiagMatr(matr);

    const int len = 2;
    qindex numReportElems[] = {0, 10};

    for (int i=0; i<len; i++) {
        qindex num = numReportElems[i];
        rootPrint(num);
        setMaxNumReportedItems(num, num);
        reportDiagMatr(matr);
    }
}



void demo_FullStateDiagMatr() {

    // force use of available hardware acceleration
    QuESTEnv env = getQuESTEnv();
    FullStateDiagMatr matr = createCustomFullStateDiagMatr(6, env.isDistributed, env.isGpuAccelerated);

    // lazily set each element individually
    for (qindex i=0; i<matr.numElems; i++)
        setInlineFullStateDiagMatr(matr, i, 1, {i+1});

    qindex numReportElems[] = {0, 50, 30, 10, 5, 1};

    for (int i=0; i<6; i++) {
        qindex num = numReportElems[i];
        rootPrint(num);
        setMaxNumReportedItems(num, num);
        reportFullStateDiagMatr(matr);
    }
}



int main() {

    initQuESTEnv();

    demo_CompMatr();
    demo_DiagMatr();
    demo_FullStateDiagMatr();

    finalizeQuESTEnv();
    return 0;
}



// MSVC's naughty corner
#else
int main() { return 0; }
#endif
