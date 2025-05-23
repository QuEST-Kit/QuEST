/** @file
 * Examples of using reportQureg in C11. This
 * example is most interesting when using 
 * simultaneously multithreaded, GPU-accelerated 
 * and distributed deployments
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <stdio.h>



/*
 * Qureg deployments
 */


void demo_serial() {

    reportStr("\n[serial]");

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

    reportStr("\n[multithreaded]");

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

    reportStr("\n[GPU-accelerated]");

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

    reportStr("\n[distributed]");

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

    reportStr("\n[auto (statevector)]");
    
    Qureg pure = createQureg(10);
    reportQuregParams(pure);
    reportQureg(pure);
    destroyQureg(pure);

    reportStr("\n[auto (density matrix)]");

    Qureg mixed = createDensityQureg(10);
    reportQuregParams(mixed);
    reportQureg(mixed);
    destroyQureg(mixed);
}


void demo_all() {

    reportStr("\n[all]");

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
