/** @file 
 * Implements Grover's algorithm for unstructured search,
 * using only X, H and multi-controlled Z gates
 *
 * Compile and run from within the build folder, using:
cmake .. -DUSER_SOURCE=../examples/grovers_search.c \
        -DOUTPUT_EXE=grovers
make
./grovers
 *
 *
 * @author Tyson Jones
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "QuEST.h"



/* effect |solElem> -> -|solElem> via a 
 * multi-controlled phase flip gate 
 */
void applyOracle(Qureg qureg, int numQubits, int solElem) {
    
    // apply X to transform |111> into |solElem>
    for (int q=0; q<numQubits; q++)
        if (((solElem >> q) & 1) == 0)
            pauliX(qureg, q);
        
    // effect |111> -> -|111>    
    int ctrls[numQubits];
    for (int q=0; q<numQubits; q++)
        ctrls[q] = q;
    multiControlledPhaseFlip(qureg, ctrls, numQubits);
    
    // apply X to transform |solElem> into |111>
    for (int q=0; q<numQubits; q++)
        if (((solElem >> q) & 1) == 0)
            pauliX(qureg, q);
}



/* apply 2|+><+|-I by transforming into the Hadamard basis 
 * and effecting 2|0><0|-I. We do this, by observing that 
 *   c..cZ = diag{1,..,1,-1} 
 *         = I - 2|1..1><1..1|
 * and hence 
 *   X..X c..cZ X..X = I - 2|0..0><0..0|
 * which differs from the desired 2|0><0|-I state only by 
 * the irrelevant global phase pi
 */
void applyDiffuser(Qureg qureg, int numQubits) {
    
    // apply H to transform |+> into |0>
    for (int q=0; q<numQubits; q++)
        hadamard(qureg, q);

    // apply X to transform |11..1> into |00..0>
    for (int q=0; q<numQubits; q++)
        pauliX(qureg, q);
    
    // effect |11..1> -> -|11..1>
    int ctrls[numQubits];
    for (int q=0; q<numQubits; q++)
        ctrls[q] = q;
    multiControlledPhaseFlip(qureg, ctrls, numQubits);
    
    // apply X to transform |00..0> into |11..1>
    for (int q=0; q<numQubits; q++)
        pauliX(qureg, q);
    
    // apply H to transform |0> into |+>
    for (int q=0; q<numQubits; q++)
        hadamard(qureg, q);
}



int main() {
    
    // prepare the hardware-agnostic QuEST environment
    QuESTEnv env = createQuESTEnv();
    
    // choose the system size
    int numQubits = 15;
    int numElems = (int) pow(2, numQubits);
    int numReps = ceil(M_PI/4 * sqrt(numElems));
    
    printf("numQubits: %d, numElems: %d, numReps: %d\n", 
        numQubits, numElems, numReps);
    
    // randomly choose the element for which to search
    srand(time(NULL));
    int solElem = rand() % numElems;
    
    // prepare |+>
    Qureg qureg = createQureg(numQubits, env);
    initPlusState(qureg);
    
    // apply Grover's algorithm
    for (int r=0; r<numReps; r++) {
        applyOracle(qureg, numQubits, solElem);
        applyDiffuser(qureg, numQubits);
        
        // monitor the probability of the solution state
        printf("prob of solution |%d> = %g\n", 
            solElem, getProbAmp(qureg, solElem));
    }
    
    // free memory 
    destroyQureg(qureg, env);
    destroyQuESTEnv(env);
    return 0;
}