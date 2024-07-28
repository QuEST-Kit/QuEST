#include "quest.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>



/*
 * distributed printing
 */


void rootPrint(qindex num) {

    if (getQuESTEnv().rank != 0)
        return;

    printf("\n\n>> set number of reported items to %lld", num);

    if (num == 0)
        printf(" (which means report all)");
         
    printf("\n\n");
}



/*
 * PauliStr
 */


void demo_PauliStr() {

    // leftmost identities are not printed
    reportPauliStr(
        getInlinePauliStr("XYZ", {2,1,0})
    );
    reportPauliStr(
        getInlinePauliStr("XYZ", {63,62,61})
    );
}



/*
 * PauliStrSum
 */


PauliStrSum prepareRandomPauliStrSum(int numQubits, int numTerms) {

    char paulis[numQubits];
    int qubits[numQubits];
    for (int i=0; i<numQubits; i++)
        qubits[i] = i;

    qcomp coeffs[numTerms];
    PauliStr strings[numTerms];
    
    for (int i=0; i<numTerms; i++) {
        coeffs[i] = qcomp(
            rand() / (qreal) RAND_MAX, 
            rand() / (qreal) RAND_MAX);

        for (int j=0; j<numQubits; j++)
            paulis[j] = "IXYZ"[rand() % 4];

        strings[i] = getPauliStr(paulis, qubits, numQubits);
    }

    return createPauliStrSum(strings, coeffs, numTerms);
}


void demo_PauliStrSum() {

    int numQubits = 64; // max 64
    int numTerms = 100;
    PauliStrSum sum = prepareRandomPauliStrSum(numQubits, numTerms);

    int len = 3;
    int numReportElems[] = {0, 10, 2};

    for (int i=0; i<len; i++) {
        qindex num = numReportElems[i];
        rootPrint(num);

        setNumReportedItems(num);
        reportPauliStrSum(sum);
    }

    destroyPauliStrSum(sum);
}



/*
 * main
 */


int main() {

    initQuESTEnv();

    // seed our RNG (not QuEST's)
    srand(time(NULL));

    demo_PauliStr();
    demo_PauliStrSum();

    finalizeQuESTEnv();
    return 0;
}
