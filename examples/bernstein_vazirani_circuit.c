/** @file 
 * Implements the Bernstien--Vazirani circuit
 *
 * @author Tyson Jones
 */

 # include <stdio.h>
 # include <math.h>
 # include <time.h>
 # include <stdlib.h>

 # include "QuEST.h" 



 void applyOracle(Qureg qureg, int numQubits, int secret) {

     int bits = secret;

     for (int q=1; q<numQubits; q++) {

         // extract the (q-1)-th bit of secret
         int bit = bits % 2;
         bits /= 2;
         
         // NOT the ancilla, controlling on the q-th qubit
         if (bit)
             controlledNot(qureg, q, 0);
     }
 }



 void measureResult(Qureg qureg, int secret) {

     // |ind> = |s>|1>
     int ind = 2*secret + 1;

     qreal prob = getProbAmp(qureg, ind);

     printf("success probability: " REAL_QASM_FORMAT " \n", prob);
 }



 void applyBernsteinVazirani(Qureg qureg, int numQubits, int secret) {

     // start in |0>
     initZeroState(qureg);

     // NOT the ancilla
     pauliX(qureg, 0);

     // H all qubits, including the ancilla
     for (int q=0; q<numQubits; q++)
         hadamard(qureg, q);

     applyOracle(qureg, numQubits, secret);

     for (int q=0; q<numQubits; q++)
         hadamard(qureg, q);

     // infer the output basis state
     measureResult(qureg, secret);
 }



 int main() {
     
     // prepare the hardware-agnostic QuEST environment
     QuESTEnv env = createQuESTEnv();

     // choose the register size
     int numQubits = 15;

     // randomly choose the secret parameter
     srand(time(NULL));
     int secret = rand() % (int) pow(2, numQubits - 1);

     // prepare our register in the |0> state
     Qureg qureg = createQureg(numQubits, env);

     // search for s using BV's algorithm
     applyBernsteinVazirani(qureg, numQubits, secret);

     // tidy up
     destroyQureg(qureg, env);
     destroyQuESTEnv(env);
     return 0;
 }
