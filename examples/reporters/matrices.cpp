#include "quest.h"

#include <iostream>



/*
 * This demo is most interesting when run distributed with up to 64 nodes
 */



void rootPrint(qindex num) {

    if (getQuESTEnv().rank != 0)
        return;

    std::cout << "\n\n>> set number of reported items to " << num;

    if (num == 0)
        std::cout << " (which means report all)";
         
    std::cout << "\n" << std::endl; // flush
}



void demo_CompMatr() {

    CompMatr matr = createCompMatr(4);
    for (qindex r=0; r<matr.numRows; r++)
        for (qindex c=0; c<matr.numRows; c++)
            matr.cpuElems[r][c] = r*matr.numRows + c + 1;
    syncCompMatr(matr);

    for (int num : {0, 12, 5, 1}) {
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

    for (int num : {0, 10}) {
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
        setInlineFullStateDiagMatr(matr, i, 1, {(qreal) i+1});

    for (int num : {0, 50, 30, 10, 5, 1}) {
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
