/** @file
 * Examples of initialising PauliStr and PauliStrSum
 * from strings, types or files, in C11. Note that
 * MSVC's C11 (which is already weird) doesn't support
 * assigning any non-complex literal to complex variables
 * nor any complex arithmetic operators, so it doesn't
 * get to play with the other children.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <stdio.h>
#include <stdlib.h>

#if !defined(_MSC_VER)



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

    // in-line via temporary array
    PauliStr a = getPauliStr("XYZ", (int[]) {5,6,7}, 3);
    reportPauliStr(a);

    // string literal
    int num = 3;
    char paulis[] = "xxx";
    int  qubits[] = {2,5,8};
    reportPauliStr(
        getPauliStr(paulis, qubits, num)
    );

    // heap arrays
    int n = 64;
    char* p = malloc(n * sizeof *p); // will have no terminal char
    int*  q = malloc(n * sizeof *q);
    for (int i=0; i<n; i++) {
        p[i] = "IXYZ"[i%4];
        q[i] = i;
    }
    reportPauliStr(
        getPauliStr(p, q, n)
    );

    // as above, but with pauli-code integers; they must not be an inline temporary array
    int codes[] = {3,2,1};
    reportPauliStr(
        getPauliStr(codes, (int[]) {5,6,7}, 3)
    );

    int* c = malloc(n * sizeof *c);
    for (int i=0; i<n; i++)
        c[i] = i % 4;
    reportPauliStr(
        getPauliStr(c, q, n)
    );

    free(p);
    free(q);
    free(c);
}



/*
 * PauliStrSum
 */


void demo_createInlinePauliStrSum() {

    reportStr("[demo_createInlinePauliStrSum]");

    // coeffs can be real, imag, or complex
    PauliStrSum a = createInlinePauliStrSum(
        "0.123 XXIXX   \n"
        "1.23i XYZXZ   \n"
        "-1-6i IIIII   \n" // no multiline strings in C, grr
    );
    reportPauliStrSum(a);
    destroyPauliStrSum(a);

    // and use scientific notation, with any amount of whitespace
    PauliStrSum b = createInlinePauliStrSum("\
        + 5E2-1i     XYZ    \n\
        - 1E-50i     IXY    \n\
        + 1 - 6E-5i  IIX    \n\
          0          III    \n\
          5.         XXX    \n\
           .5        ZZZ    \n\
    ");
    reportPauliStrSum(b);
    destroyPauliStrSum(b);

    // paulis can be lowercase and 0-3, with any whitespace
    PauliStrSum c = createInlinePauliStrSum("\
        2.5     0123 ixyz IXYZ                                 \n\
        3.5     0iII 1xXX 2yYY                                 \n\
                                                               \n\
        -1.5E-15   -   5.123E-30i    0 0 0 0 0 0 0 0 1 1 1 1   \n\
        1.5        +   5i            3 3 3 3 2 2 2 2 I X Y Z   \n\
    ");
    reportPauliStrSum(c);
    destroyPauliStrSum(c);
}


void demo_createPauliStrSum() {

    reportStr("[demo_createPauliStrSum]");

    // inline using temporary arrays
    PauliStrSum a = createPauliStrSum(
        (PauliStr[]) {
            getInlinePauliStr("XYZ", {1,2,3}),
            getInlinePauliStr("XYZ", {1,2,3}),  // duplication allowed
            getInlinePauliStr("XYZ", {1,2,5}),
            getInlinePauliStr("xxx", {4,2,3}),
            getInlinePauliStr("yyy", {1,2,3}),
            getInlinePauliStr("zzz", {9,2,3}),
            getInlinePauliStr("000", {1,2,3}),
            getInlinePauliStr("111", {1,8,3}),
            getInlinePauliStr("222", {8,2,3}),
            getInlinePauliStr("333", {1,2,8}),
            getInlinePauliStr("1xX", {1,6,3}) }, 
        (qcomp[]) {
            0.5, 0.2i, 0.9 + 2.131i, 1, 2, 3, 4, 5, 6, 7, 8 },
        11
    );
    reportPauliStrSum(a);
    destroyPauliStrSum(a);

    // pre-allocated
    int n = 64;
    char* p = malloc(n * sizeof *p); // will have no terminal char
    int*  q = malloc(n * sizeof *q);
    for (int i=0; i<n; i++) {
        p[i] = "IXYZ"[i%4];
        q[i] = i;
    }
    int numTerms = 1;
    PauliStr strings[] = {getPauliStr(p, q, n)};
    qcomp coeffs[] = {-1E-40 + .45E2i};
    PauliStrSum b = createPauliStrSum(strings, coeffs, numTerms);
    reportPauliStrSum(b);
    destroyPauliStrSum(b);
}


void demo_createPauliStrSumFromFile() {

    reportStr("[demo_createPauliStrSumFromFile]");

    char *fn = "test.txt";

    // file contents can be identical to createInlinePauliStrSum input
    if (getQuESTEnv().rank == 0) {
        FILE *f = fopen(fn, "w");
        fprintf(f, "\
            + 5E2-1i     XYZ    \n\
            - 1E-50i     IXY    \n\
            + 1 - 6E-5i  IIX    \n\
            0            III    \n\
            5.           IXX    \n\
            .5           ZYX    \n\
        ");
        fclose(f);
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



// MSVC's naughty corner
#else
int main() { return 0; }
#endif
