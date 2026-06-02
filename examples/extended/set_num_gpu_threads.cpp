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
#include <iostream>
#include <chrono>


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
    std::cout << "Using " << numThreadsPerBlock << " threads per block... " << std::flush;

    setQuESTNumGpuThreadsPerBlock(numThreadsPerBlock);

    // warmup
    for (int r=0; r<NUM_REPS; r++)
        simulation(qureg);
    syncQuESTEnv();

    using clock = std::chrono::steady_clock;
    auto start = clock::now();

    for (int r=0; r<NUM_REPS; r++)
        simulation(qureg);
    syncQuESTEnv();

    auto end = clock::now();
    auto dur = std::chrono::duration<double>(end - start).count();
    auto av  = dur / NUM_REPS;

    std::cout << " took " << av << "s" << std::endl;
}


int main()
{
    initQuESTEnv();

    // This example is pointless without a GPU!
    if (!getQuESTEnv().isGpuAccelerated) {
        std::cout 
            << "GPU acceleration is not enabled, and so changing the number "
            << "of threads per block has no effect. Exiting..."
            << std::endl;
        finalizeQuESTEnv();
        return 0;
    }

    // The initial number of threads per block is informed by the optional environment
    // variable QUEST_DEFAULT_NUM_GPU_THREADS_PER_BLOCK. If not specified, QuEST will
    // use the value of the CMake option of the same name passed during compilation,
    // which itself will has a default of 128
    auto initNumTPB = getQuESTNumGpuThreadsPerBlock();
    std::cout << "Initial numThreadsPerBlock: " << initNumTPB << "\n\n";

    // Create a statevector parallelised only by the GPU
    Qureg qureg = createCustomQureg(NUM_QUBITS, 0, 0, 1, 0);
    reportQuregParams(qureg);

    // Benchmark QuEST with sensible numbers of threads per block (multiples of warp size)
    for (auto numTPB : {64, 128, 256, 512, 1024})
        benchmark(qureg, numTPB);

    // Try silly parameters ¯\_(ツ)_/¯
    setQuESTValidationOff();
    for (auto numTPB : {31, 15, 5, 1})
        benchmark(qureg, numTPB);
    
    finalizeQuESTEnv();
    return 0;
}
