#include "quest.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <time.h>



/*
 * distributed printing
 */


void rootPrint(qindex num) {

    if (getQuESTEnv().rank != 0)
        return;

    std::cout << std::endl << std::endl << ">> set number of reported items to " << num;

    if (num == 0)
        std::cout << " (which means report all)";
         
    std::cout << std::endl << std::endl;
}



/*
 * PauliStr
 */


void demo_PauliStr() {

    // leftmost identities are not printed
    reportPauliStr(
        getPauliStr("XYZ", {2,1,0})
    );
    reportPauliStr(
        getPauliStr("XYZ", {63,62,61})
    );
}



/*
 * PauliStrSum
 */


PauliStrSum prepareRandomPauliStrSum(int numQubits, int numTerms) {

    std::vector<char> paulis(numQubits);
    std::vector<int> qubits(numQubits);
    for (int i=0; i<numQubits; i++)
        qubits[i] = i;

    std::vector<qcomp> coeffs(numTerms);
    std::vector<PauliStr> strings(numTerms);
    
    for (int i=0; i<numTerms; i++) {
        coeffs[i] = qcomp(
            rand() / (qreal) RAND_MAX, 
            rand() / (qreal) RAND_MAX);

        for (int j=0; j<numQubits; j++)
            paulis[j] = "IXYZ"[rand() % 4];

        std::string seq = std::string(paulis.begin(), paulis.end());
        strings[i] = getPauliStr(seq, qubits);
    }

    return createPauliStrSum(strings, coeffs);
}


void demo_PauliStrSum() {

    int numQubits = 15; // max 64
    int numTerms = 100;
    PauliStrSum sum = prepareRandomPauliStrSum(numQubits, numTerms);

    std::vector<int> numReportElems {0, 10, 2};

    for (int num : numReportElems) {
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
