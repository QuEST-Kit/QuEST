#include "quest.h"

#include <stdio.h>
#include <stdlib.h>


/*
 * This demo is most interesting when compiling to simultaneously
 * enable multithreaded, GPU-accelerated and distributed deployments.
 * When launching, pass an integer 1-6 as a command-line argument to
 * choose a deployment.
 */



/*
 * distributed printing
 */


void rootPrint(const char* mode) {

    // cannot call this before initialising QuEST
    if (getQuESTEnv().rank == 0)
        printf("[%s]\n\n", mode);
}



/*
 * environment deployments
 */


void demo_serial() {

    int useDistrib = 0;
    int useGpuAccel = 0;
    int useMultithread = 0;

    initCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread);

    rootPrint("serial");
    reportQuESTEnv();
}


void demo_multithreaded() {

    int useDistrib = 0;
    int useGpuAccel = 0;
    int useMultithread = 1;

    initCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread);

    rootPrint("multithreaded");
    reportQuESTEnv();
}


void demo_gpuAccelerated() {

    int useDistrib = 0;
    int useGpuAccel = 1;
    int useMultithread = 0;

    initCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread);

    rootPrint("GPU-accelerated");
    reportQuESTEnv();
}


void demo_distributed() {

    int useDistrib = 1;
    int useGpuAccel = 0;
    int useMultithread = 0;

    initCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread);

    rootPrint("distributed");
    reportQuESTEnv();
}


void demo_all() {

    int useDistrib = 1;
    int useGpuAccel = 1;
    int useMultithread = 1;

    initCustomQuESTEnv(useDistrib, useGpuAccel, useMultithread);

    rootPrint("all");
    reportQuESTEnv();
}


void demo_auto() {

    initQuESTEnv();

    rootPrint("auto");
    reportQuESTEnv();
}



/*
 * main
 */


int printCmdArgInfo() {

    // we can only initialise the QuEST environment once, so
    // we accept command-line arg 1-6 to choose the deployment

    // all nodes must print because we've not yet setup MPI
    printf(
        "Must pass single cmd-line argument:\n"
        "  1 = serial\n"
        "  2 = multithreaded\n"
        "  3 = GPU-accelerated\n"
        "  4 = distributed\n"
        "  5 = all\n"
        "  6 = auto\n");

    // process return code, indicate error
    return 1;
}


int main(int argc, char* argv[]) {

    if (argc != 2)
        return printCmdArgInfo();

    int opt = atoi(argv[1]);
    if (opt < 1 || opt > 6)
        return printCmdArgInfo();

    if (opt == 1) demo_serial();
    if (opt == 2) demo_multithreaded();
    if (opt == 3) demo_gpuAccelerated();
    if (opt == 4) demo_distributed();
    if (opt == 5) demo_all();
    if (opt == 6) demo_auto();

    finalizeQuESTEnv();
    return 0;
}
