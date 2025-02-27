#include "quest.h"

#include <stdio.h>


/*
 * This demo is most interesting when using simultaneous
 * multithreaded, GPU-accelerated and distributed deployments
 */



/*
 * distributed printing
 */


void rootPrint(const char* mode) {

    if (getQuESTEnv().rank == 0)
        printf("\n\n>> [%s]\n\n", mode);
}



/*
 * Qureg deployments
 */


void demo_serial() {

    rootPrint("serial");

    int numQubits = 10;
    int isDensMatr = 0;

    int useDistrib = 0;
    int useGpuAccel = 0;
    int useMultithread = 0;

    Qureg qureg = createCustomQureg(
        numQubits, isDensMatr, 
        useDistrib, useGpuAccel, useMultithread);

    reportQuregParams(qureg);
    reportQureg(qureg);
    destroyQureg(qureg);
}


void demo_multithreaded() {

    rootPrint("multithreaded");

    int numQubits = 10;
    int isDensMatr = 0;

    int useDistrib = 0;
    int useGpuAccel = 0;
    int useMultithread = getQuESTEnv().isMultithreaded;

    Qureg qureg = createCustomQureg(
        numQubits, isDensMatr, 
        useDistrib, useGpuAccel, useMultithread);

    reportQuregParams(qureg);
    reportQureg(qureg);
    destroyQureg(qureg);
}


void demo_gpuAccelerated() {

    rootPrint("GPU-accelerated");

    int numQubits = 10;
    int isDensMatr = 0;

    int useDistrib = 0;
    int useGpuAccel = getQuESTEnv().isGpuAccelerated;
    int useMultithread = 0;

    Qureg qureg = createCustomQureg(
        numQubits, isDensMatr, 
        useDistrib, useGpuAccel, useMultithread);

    reportQuregParams(qureg);
    reportQureg(qureg);
    destroyQureg(qureg);
}


void demo_distributed() {

    rootPrint("distributed");

    int numQubits = 10;
    int isDensMatr = 0;

    int useDistrib = getQuESTEnv().isDistributed;
    int useGpuAccel = 0;
    int useMultithread = 0;

    Qureg qureg = createCustomQureg(
        numQubits, isDensMatr, 
        useDistrib, useGpuAccel, useMultithread);

    reportQuregParams(qureg);
    reportQureg(qureg);
    destroyQureg(qureg);
}


void demo_auto() {

    rootPrint("auto (statevector)");
    
    Qureg pure = createQureg(10);
    reportQuregParams(pure);
    reportQureg(pure);
    destroyQureg(pure);

    rootPrint("auto (density matrix)");

    Qureg mixed = createDensityQureg(10);
    reportQuregParams(mixed);
    reportQureg(mixed);
    destroyQureg(mixed);
}


void demo_all() {

    rootPrint("all");

    int numQubits = 10;
    int isDensMatr = 0;

    int useDistrib = getQuESTEnv().isDistributed;
    int useGpuAccel = getQuESTEnv().isGpuAccelerated;
    int useMultithread = getQuESTEnv().isMultithreaded;

    // throws if both GPU and multithreaded
    Qureg qureg = createCustomQureg(
        numQubits, isDensMatr, 
        useDistrib, useGpuAccel, useMultithread);

    reportQuregParams(qureg);
    reportQureg(qureg);
    destroyQureg(qureg);
}



/*
 * main
 */


int main() {

    initQuESTEnv();
    QuESTEnv env = getQuESTEnv();

    demo_serial();

    if (env.isMultithreaded)
        demo_multithreaded();

    if (env.isGpuAccelerated)
        demo_gpuAccelerated();

    if (env.isDistributed)
        demo_distributed();
        
    demo_auto();
    demo_all(); // may throw

    finalizeQuESTEnv();
    return 0;
}
