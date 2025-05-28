/** @file
 * Examples of using reportPauliStr and 
 * reportPauliStrSum in C11.
 * 
 * @author Tyson Jones
*/

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

    char paulis[64];
    int qubits[64];
    for (int i=0; i<numQubits; i++)
        qubits[i] = i;

    qcomp coeffs[100];
    PauliStr strings[100];
    
    for (int i=0; i<numTerms; i++) {
        coeffs[i] = getQcomp(
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
    int numTerms = 100; // max 100, only imposed by hardcoded arrays
    PauliStrSum sum = prepareRandomPauliStrSum(numQubits, numTerms);

    int len = 3;
    int numReportElems[] = {0, 10, 2};

    for (int i=0; i<len; i++) {
        qindex num = numReportElems[i];
        rootPrint(num);

        setMaxNumReportedItems(num, num);
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
