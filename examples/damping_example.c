#include <stdio.h>
#include "QuEST.h"

int main(int narg, char *varg[])
{

    /*
     * PREPARE QuEST environment
     * (Required only once per program)
     */

    QuESTEnv env = createQuESTEnv();

    printf("-------------------------------------------------------\n");
    printf("Running QuEST damping example:\n\t Basic circuit involving damping of a qubit.\n");
    printf("-------------------------------------------------------\n");

    /*
     * PREPARE QUBIT SYSTEM in density matrix form so that we can apply decoherence/errors
     */

    Qureg qubits = createDensityQureg(1, env);
    initPlusState(qubits);

    /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    printf("\n Reporting the qubit stat to screen:\n");
    reportStateToScreen(qubits, env, 0);

    /*
     * APPLY Damping 10 times with probability 0.1
     */
    printf("\n Applying damping 10 times with probability 0.1 \n");

    int counter;
    for (counter = 0; counter < 10; counter++)
    {
        applyOneQubitDampingError(qubits, 0, 0.1);
        printf("\n Qubit state after applying damping %d times:\n", counter+1);
        reportStateToScreen(qubits, env, 0);
    }

    /*
     * FREE MEMORY
     */

    destroyQureg(qubits, env);

    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    destroyQuESTEnv(env);
    return 0;
}
