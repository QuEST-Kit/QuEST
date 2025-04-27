/** @file
 * Examples of initialising PauliStr and PauliStrSum
 * from strings, types or files, in C++14.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

using std::vector;



/*
 * PauliStr
 */


void demo_getInlinePauliStr() {

    reportStr("[demo_getInlinePauliStr]");

    // can specify indices as {...} without temporary-array syntax 
    PauliStr a = getInlinePauliStr("XYZII", {4,3,2,1,0});
    reportPauliStr(a);

    // additionally accept 0-3 and lowercase
    reportPauliStr(
        getInlinePauliStr("0123ixyzIXYZ", {11,10,9,8, 7,6,5,4, 3,2,1,0})
    );

    // can specify arbitrary qubit indices of Paulis
    reportPauliStr(
        getInlinePauliStr("XXYYZZ", {5,50, 10,60, 30,40})
    );
}


void demo_getPauliStr() {

    reportStr("[demo_getPauliStr]");

    // C++ can pass no indices to set Paulis upon rightmost qubits
    PauliStr a = getPauliStr("XYZ");
    reportPauliStr(a);

    // C++ can pass vector initialisers
    reportPauliStr(
        getPauliStr("XYZ", {5,6,7})
    );

    // string literal
    int num = 3;
    std::string paulis = "xxx";
    int  qubits[] = {2,5,8};
    reportPauliStr(
        getPauliStr(paulis, qubits, num)
    );

    // heap arrays
    int n = 64;
    char* p = (char*) malloc(n * sizeof *p); // will have no terminal char
    int*  q = (int*)  malloc(n * sizeof *q);
    for (int i=0; i<n; i++) {
        p[i] = (char) ("IXYZ"[i%4]);
        q[i] = i;
    }
    reportPauliStr(
        getPauliStr(p, q, n)
    );
    free(p);
    free(q);

    // can also pass integer Pauli codes, alas using the C API, because
    // integer initializer literals resemble both vector<int> and 
    // std::string types, leading to ambiguous overload invocation :(
    int codes[] = {2,0,3};
    int targs[] = {5,6,7};
    reportPauliStr( getPauliStr(codes, targs, 3) );
}



/*
 * PauliStrSum
 */


void demo_createInlinePauliStrSum() {

    reportStr("[demo_createInlinePauliStrSum]");

    // coeffs can be real, imag, or complex (via C++ raw multilines)
    PauliStrSum a = createInlinePauliStrSum(R"(
        0.123 XXIXX
        1.23i XYZXZ
        -1-6i IIIII
    )");
    reportPauliStrSum(a);
    destroyPauliStrSum(a);

    // and use scientific notation, with any amount of whitespace
    PauliStrSum b = createInlinePauliStrSum(R"(
        + 5E2-1i     XYZ 
        - 1E-50i     IXY 
        + 1 - 6E-5i  IIX 
          0          III 
          5.         XXX 
          .5         ZZZ 
    )");
    reportPauliStrSum(b);
    destroyPauliStrSum(b);

    // paulis can be lowercase and 0-3, with any whitespace
    PauliStrSum c = createInlinePauliStrSum(R"(
        2.5     0123 ixyz IXYZ 
        3.5     0iII 1xXX 2yYY 
                                                               
        -1.5E-15   -   5.123E-30i    0 0 0 0 0 0 0 0 1 1 1 1 
        1.5        +   5i            3 3 3 3 2 2 2 2 I X Y Z 
    )");
    reportPauliStrSum(c);
    destroyPauliStrSum(c);
}


void demo_createPauliStrSum() {

    reportStr("[demo_createPauliStrSum]");

    // inline using C++ vector initialisers
    PauliStrSum a = createPauliStrSum(
        {getPauliStr("XYZ", {1,2,3}),
         getPauliStr("XYZ", {1,2,3}),  // duplication allowed
         getPauliStr("XYZ", {1,2,5}),
         getPauliStr("xxx", {4,2,3}),
         getPauliStr("yyy", {1,2,3}),
         getPauliStr("zzz", {9,2,3}),
         getPauliStr("000", {1,2,3}),
         getPauliStr("111", {1,8,3}),
         getPauliStr("222", {8,2,3}),
         getPauliStr("333", {1,2,8}),
         getPauliStr("1xX", {1,6,3})}, 
        {0.5, 0.2_i, 0.9 + 2.131_i, 1, 2, 3, 4, 5, 6, 7, 8} // use 1_i over 1i for precision agnosticism
    );
    reportPauliStrSum(a);
    destroyPauliStrSum(a);

    // pre-allocated
    int n = 64;
    char* p = (char*) malloc(n * sizeof *p); // will have no terminal char
    int*  q = (int*)  malloc(n * sizeof *q);
    for (int i=0; i<n; i++) {
        p[i] = "IXYZ"[i%4];
        q[i] = i;
    }
    vector<PauliStr> strings = {getPauliStr(p, q, n)};
    vector<qcomp> coeffs = {-1E-40 + .45E2_i};
    PauliStrSum b = createPauliStrSum(strings, coeffs);
    reportPauliStrSum(b);
    destroyPauliStrSum(b);
}


void demo_createPauliStrSumFromFile() {

    reportStr("[demo_createPauliStrSumFromFile]");

    std::string fn = "test.txt";

    // file contents can be identical to createInlinePauliStrSum input
    if (getQuESTEnv().rank == 0) {
        std::ofstream file;
        file.open(fn);
        file << R"(
            + 5E2-1i     XYZ 
            - 1E-50i     IXY 
            + 1 - 6E-5i  IIX 
            0            III 
            5.           IXX
            .5           ZYX 
        )";
        file.close();
    }
    syncQuESTEnv();

    PauliStrSum a = createPauliStrSumFromFile(fn);
    reportPauliStrSum(a);
    destroyPauliStrSum(a);
}


void demo_createPauliStrSumFromReversedFile() {

    reportStr("[demo_createPauliStrSumFromReversedFile]");

    PauliStrSum a = createPauliStrSumFromReversedFile("test.txt");
    reportPauliStrSum(a);
    destroyPauliStrSum(a);
}



/*
 * main
 */


int main() {

    initQuESTEnv();

    demo_getInlinePauliStr();
    demo_getPauliStr();
    demo_createInlinePauliStrSum();
    demo_createPauliStrSum();
    demo_createPauliStrSumFromFile();
    demo_createPauliStrSumFromReversedFile();

    finalizeQuESTEnv();
    return 0;
}
