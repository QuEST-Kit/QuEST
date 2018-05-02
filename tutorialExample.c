#include <stdio.h>
#include "QuEST.h"

int main (int narg, char *varg[]) {

    /*
     * PREPARE QuEST environment
     * (Required only once per program)
     */

    QuESTEnv env;
    initQuESTEnv(&env);

    printf("-------------------------------------------------------\n");
    printf("Running QuEST tutorial:\n\t Basic circuit involving a system of 3 qubits.\n");
    printf("-------------------------------------------------------\n");

    /*
     * PREPARE QUBIT SYSTEM
     */

    MultiQubit qubits; 
    createMultiQubit(&qubits, 3, env);
    initStateZero(&qubits);


    /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    printf("\nThis is our environment:\n");
    reportMultiQubitParams(qubits);
    reportQuESTEnv(env);

    /*
     * APPLY CIRCUIT
     */

    hadamard(qubits, 0);
    controlledNot(qubits, 0, 1);
    rotateY(qubits, 2, .1);

    multiControlledPhaseGate(qubits, (int []){0, 1, 2}, 3);

    ComplexMatrix2 u;
    u.r0c0 = (Complex) {.real=.5, .imag= .5};
    u.r0c1 = (Complex) {.real=.5, .imag=-.5}; 
    u.r1c0 = (Complex) {.real=.5, .imag=-.5};
    u.r1c1 = (Complex) {.real=.5, .imag= .5};
    unitary(qubits, 0, u);

    Complex a, b;
    a.real = .5; a.imag =  .5;
    b.real = .5; b.imag = -.5;
    compactUnitary(qubits, 1, a, b);

    Vector v;
    v.x = 1; v.y = 0; v.z = 0;
    rotateAroundAxis(qubits, 2, 3.14/2, v);

    controlledCompactUnitary(qubits, 0, 1, a, b);

    multiControlledUnitary(qubits, (int []){0, 1}, 2, 2, u);

    printf("\nCircuit output:\n");

    REAL prob = getProbEl(qubits, 7);
    printf("Probability amplitude of |111>: %f\n", prob);

    prob = findProbabilityOfOutcome(qubits, 2, 1);
    printf("Probability of qubit 2 being in state 1: %f\n", prob);

    int outcome = measure(qubits, 0);
    printf("Qubit 0 was measured in state %d\n", outcome);

    outcome = measureWithStats(qubits, 2, &prob);
    printf("Qubit 2 collapsed to %d with probability %f\n", outcome, prob);


    /*
     * FREE MEMORY
     */

    destroyMultiQubit(qubits, env); 


    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    closeQuESTEnv(env);
    return 0;
}
