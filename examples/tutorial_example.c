/** @file 
 * A demo of QuEST
 *
 * @author Ania Brown
 * @author Tyson Jones
 */

#include <stdio.h>
#include "QuEST.h"

int main (int narg, char *varg[]) {



    /*
     * PREPARE QuEST environment
     * (Required only once per program)
     */

    QuESTEnv env = createQuESTEnv();

    printf("-------------------------------------------------------\n");
    printf("Running QuEST tutorial:\n\t Basic circuit involving a system of 3 qubits.\n");
    printf("-------------------------------------------------------\n");



    /*
     * PREPARE QUBIT SYSTEM
     */

    Qureg qubits = createQureg(3, env);
    initZeroState(qubits);



    /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    printf("\nThis is our environment:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);



    /*
     * APPLY CIRCUIT
     */

    hadamard(qubits, 0);
    controlledNot(qubits, 0, 1);
    rotateY(qubits, 2, .1);

    int targs[] = {0,1,2};
    multiControlledPhaseFlip(qubits, targs, 3);

    ComplexMatrix2 u = {
        .real={{.5,.5},{.5,.5}},
        .imag={{.5,-.5},{-.5,.5}}
    };
    unitary(qubits, 0, u);

    Complex a, b;
    a.real = .5; a.imag =  .5;
    b.real = .5; b.imag = -.5;
    compactUnitary(qubits, 1, a, b);

    Vector v;
    v.x = 1; v.y = 0; v.z = 0;
    rotateAroundAxis(qubits, 2, 3.14/2, v);

    controlledCompactUnitary(qubits, 0, 1, a, b);

    int ctrls[] = {0,1};
    multiControlledUnitary(qubits, ctrls, 2, 2, u);
    
    ComplexMatrixN toff = createComplexMatrixN(3);
    toff.real[6][7] = 1;
    toff.real[7][6] = 1;
    for (int i=0; i<6; i++)
        toff.real[i][i] = 1;
    multiQubitUnitary(qubits, targs, 3, toff);
    
    
    
    /*
     * STUDY QUANTUM STATE
     */

    printf("\nCircuit output:\n");

    qreal prob;
    prob = getProbAmp(qubits, 7);
    printf("Probability amplitude of |111>: " REAL_STRING_FORMAT "\n", prob);

    prob = calcProbOfOutcome(qubits, 2, 1);
    printf("Probability of qubit 2 being in state 1: " REAL_STRING_FORMAT "\n", prob);

    int outcome = measure(qubits, 0);
    printf("Qubit 0 was measured in state %d\n", outcome);

    outcome = measureWithStats(qubits, 2, &prob);
    printf("Qubit 2 collapsed to %d with probability " REAL_STRING_FORMAT "\n", outcome, prob);



    /*
     * FREE MEMORY
     */

    destroyQureg(qubits, env); 
    destroyComplexMatrixN(toff);



    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return 0;
}
