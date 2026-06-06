/** @file
 * Minimal end-to-end example of a standard QuEST simulation.
 *
 * This example demonstrates the typical lifecycle of a QuEST program:
 *   1. initialise the QuEST environment,
 *   2. create a quantum register,
 *   3. initialise its state,
 *   4. apply operations,
 *   5. report or calculate results,
 *   6. destroy allocated objects,
 *   7. finalise the environment.
 *
 * @author Nathan Delcid
 */

#include "quest.h"

int main(void) {

    // Initialise the QuEST environment. This prepares available hardware
    // backends such as MPI, GPU acceleration and multithreading if enabled.
    initQuESTEnv();

    // Create a 2-qubit statevector Qureg.
    //
    // Its amplitudes correspond to:
    //   alpha_0 |00>
    //   alpha_1 |01>
    //   alpha_2 |10>
    //   alpha_3 |11>
    //
    // Qubit 0 is the rightmost, least-significant qubit.
    Qureg qureg = createQureg(2);

    // Initialise the register to |00>.
    initZeroState(qureg);

    // Apply some example operations.
    //
    // Hadamard on qubit 0 prepares:
    //   (|00> + |01>) / sqrt(2)
    //
    // Pauli-X on qubit 1 then flips the left qubit, giving:
    //   (|10> + |11>) / sqrt(2)
    applyHadamard(qureg, 0);
    applyPauliX(qureg, 1);

    // Report the resulting statevector amplitudes.
    reportQureg(qureg);

    // Free memory allocated for the Qureg.
    destroyQureg(qureg);

    // Finalise the QuEST environment.
    finalizeQuESTEnv();

    return 0;
}
