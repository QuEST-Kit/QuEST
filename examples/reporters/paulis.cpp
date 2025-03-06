/** @file
 * Examples of using reportPauliStr and 
 * reportPauliStrSum in C++14.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <time.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;



/*
 * distributed printing
 */


void rootPrint(qindex num) {

    if (getQuESTEnv().rank != 0)
        return;

    cout << endl << endl << ">> set number of reported items to " << num;

    if (num == 0)
        cout << " (which means report all)";
         
    cout << endl << endl;
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

    vector<char> paulis(numQubits);
    vector<int> qubits(numQubits);
    for (int i=0; i<numQubits; i++)
        qubits[i] = i;

    vector<qcomp> coeffs(numTerms);
    vector<PauliStr> strings(numTerms);
    qreal randMax = static_cast<qreal>(RAND_MAX);
    
    for (int i=0; i<numTerms; i++) {
        coeffs[i] = getQcomp(
            rand() / randMax,
            rand() / randMax);

        for (int j=0; j<numQubits; j++)
            paulis[j] = "IXYZ"[rand() % 4];

        string seq = string(paulis.begin(), paulis.end());
        strings[i] = getPauliStr(seq, qubits);
    }

    return createPauliStrSum(strings, coeffs);
}


void demo_PauliStrSum() {

    int numQubits = 15; // max 64
    int numTerms = 100;
    PauliStrSum sum = prepareRandomPauliStrSum(numQubits, numTerms);

    vector<int> numReportElems {0, 10, 2};

    for (int num : numReportElems) {
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
