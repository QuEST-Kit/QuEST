/** @file
 * 
 * An example of using QuEST's experimental
 * setQuESTNumGpuThreadsPerBlock() function
 * to change the parallelisation granularity
 * of GPU simulation
 * 
 * @author Tyson Jones
 */

#include "quest.h"
#include <stdio.h>
#include <time.h>


const int NUM_REPS = 10;
const int NUM_QUBITS = 25;  // 512 MiB (at double precision)


void simulation(Qureg qureg)
{
    // put your favourite QuEST simulation here
    initRandomPureState(qureg);
    applyFullQuantumFourierTransform(qureg, /*inverse=*/false);
    calcTotalProb(qureg);
}


void benchmark(Qureg qureg, int numThreadsPerBlock)
{
    printf("Using %d threads per block... ", numThreadsPerBlock);
    fflush(stdout);

    setQuESTNumGpuThreadsPerBlock(numThreadsPerBlock);

    // warmup
    for (int r=0; r<NUM_REPS; r++)
        simulation(qureg);
    syncQuESTEnv();

    double start = (double) clock();

    for (int r = 0; r < NUM_REPS; r++)
        simulation(qureg);
    syncQuESTEnv();

    double end = (double) clock();
    double dur = (end - start) / CLOCKS_PER_SEC;
    double av = dur / NUM_REPS;

    printf("took %fs\n", av);
}


int main(void)
{
    initQuESTEnv();

    // This example is pointless without a GPU!
    if (!getQuESTEnv().isGpuAccelerated)
    {
        printf(
            "GPU acceleration is not enabled, and so changing the number "
            "of threads per block has no effect. Exiting...\n");
        finalizeQuESTEnv();
        return 0;
    }

    int initNumTPB = getQuESTNumGpuThreadsPerBlock();
    printf("Initial numThreadsPerBlock: %d\n\n", initNumTPB);

    // Create a statevector parallelised only by the GPU
    Qureg qureg = createCustomQureg(NUM_QUBITS, 0, 0, 1, 0);
    reportQuregParams(qureg);

    // Benchmark sensible parameters
    int goodTPB[] = {64, 128, 256, 512, 1024};
    for (int i = 0; i < 5; i++)
        benchmark(qureg, goodTPB[i]);

    // Try silly parameters
    setQuESTValidationOff();
    int badTPB[] = {31, 15, 5, 1};
    for (int i = 0; i < 4; i++)
        benchmark(qureg, badTPB[i]);

    destroyQureg(qureg);
    finalizeQuESTEnv();

    return 0;
}
