/** @file
 * An example of using QuEST to perform closed
 * dynamical simulation via Trotterisation of the
 * unitary-time evolution operator.
 * 
 * @author Tyson Jones
*/

#include "quest.h"
#include <vector>
#include <string>

using std::vector;
using std::string;



/*
 * Prepare a Hamiltonian H under which dyanmical
 * evolution will be simulated via Trotterisation
 * of unitary-time evolution operator e^(-itH).
 * If the Hamiltonian was fixed/known in advance,
 * we could instead use createInlinePauliStrSum()
 */

PauliStrSum createMyHamiltonian(int numQubits) {

    // we prepare a Heisenberg XYZ spin-ring Hamiltonian,
    // i.e. H = -1/2 sum( Jx XX + Jy YY + Jz ZZ + h Z )
    // upon all nearest neighbour qubits, with periodicity.
    // The coefficients must be real for H to be Hermitian
    // and ergo its time-evolution operator to be unitary,
    // although they must be represented with a qcomp type.
    vector<string> operators = {"XX", "YY", "ZZ", "Z"};
    vector<qcomp> coefficients = {.1, .2, .3, .4}; // Jx,Jy,Jz,h

    // we will populate the below vectors with 4*numQubits
    // elements which we could pre-allocate with .reserve,
    // but we might incur Donald Knuth's justified wrath.
    vector<PauliStr> allStrings;
    vector<qcomp> allCoeffs;

    // prepare all XX + YY + ZZ
    for (int p=0; p<3; p++) {
        for (int i=0; i<numQubits; i++) {

            // A_i, A_i+1
            vector<int> targs = {i, (i+1)%numQubits};
            PauliStr str = getPauliStr(operators[p], targs);

            allStrings.push_back(str);
            allCoeffs.push_back(coefficients[p]);
        }
    }

    // prepare Z
    for (int i=0; i<numQubits; i++) {
        allStrings.push_back(getPauliStr(operators[3], {i}));
        allCoeffs.push_back(coefficients[3]);
    }

    // must be freed by caller
    return createPauliStrSum(allStrings, allCoeffs);
}



/*
 * Prepare the observable operator O under which the
 * evolved state (under H above) will be measured.
 * If this were one term (a single ensor product of
 * Pauli operators), we could return instead a PauliStr
 * but we here return an arbitrary weighted sum thereof.
 */

PauliStrSum createMyObservable(int numQubits) {

    // we prepare a weighted sum of alternating Paulis
    // upon each qubit, i.e. 1 X0 + 2 Y1 + 3 Z2 + 1 X3 + ...
    // where the coefficients are real such that the
    // output observable is Hermitian.

    vector<PauliStr> strings(numQubits);
    vector<qcomp> coeffs(numQubits);

    for (int i=0; i<numQubits; i++) {
        strings[i] = getPauliStr({"XYZ"[i%3]}, {i});
        coeffs[i] = getQcomp(i%4 + 1, 0);
    }

    // must be freed by caller
    return createPauliStrSum(strings, coeffs);
}



/*
 * Preview the pre-evolved system and operators thereupon.
 */

void reportMyStructs(Qureg qureg, PauliStrSum hamil, PauliStrSum observ) {

    setMaxNumReportedSigFigs(6);   // sig-figs in scalars
    setNumReportedNewlines(2);     // spacing between reports
    setReportedPauliChars(".XYZ"); // print I as .
    setReportedPauliStrStyle(0);   // print XYZ (0) or Z3 Y2 X1 (1)
    setMaxNumReportedItems(8, 8);  // show max 8 qureg amplitudes

    reportStr("[Initial state]");
    reportQureg(qureg);

    setMaxNumReportedItems(0, 0); // show 0=all Pauli operators

    reportStr("[Hamiltonian]");
    reportPauliStrSum(hamil);

    reportStr("[Observable]");
    reportPauliStrSum(observ);
}



/*
 * Simulate evolution of |+> to |psi(t)> = exp(-itH)|+> 
 * via Trotterisation, calculating the expected value
 * <O> = <psi(t)|O|psi(t)> where Hamiltonian H and
 * observable O are as prepared above.
 */

int main() {
    initQuESTEnv();

    // prepare qureg=|0>, H and O
    int numQubits = 20;
    Qureg qureg = createQureg(numQubits);
    PauliStrSum hamil = createMyHamiltonian(numQubits);
    PauliStrSum observ = createMyObservable(numQubits);

    // init qureg=|+> and report qureg, H, O
    initPlusState(qureg);
    reportMyStructs(qureg, hamil, observ);

    // tidy reporting of below expectation values
    setMaxNumReportedSigFigs(3);
    setNumReportedNewlines(1);

    // evolve by repeatedly (each is a "step") Trotterising
    // exp(-i dt H) with the specified order and repetitions.
    qreal dt  = 0.1;
    int order = 4;
    int reps  = 5;
    int steps = 20;

    for (int i=0; i<steps; i++) {

        // evolve qureg under (approx) exp(-i dt H)
        applyTrotterizedUnitaryTimeEvolution(qureg, hamil, dt, order, reps);

        // calculate and report <O>
        qreal time = dt * (i+1);
        qreal expec = calcExpecPauliStrSum(qureg, observ);
        reportScalar("<O(t=" + std::to_string(time) + ")>", expec);
    }

    reportStr("");
    
    // preview the final state...
    setNumReportedNewlines(2);
    setMaxNumReportedItems(25, 25);
    reportStr("[Final state]");
    reportQureg(qureg);

    // and some of its properties...
    reportScalar("Normalisation", calcTotalProb(qureg));
    reportScalar("Probability of first qubit", calcProbOfQubitOutcome(qureg, 0, 0));

    // including the overlap with |+>
    for (int i=0; i<numQubits; i++) 
        applyHadamard(qureg, i);
    reportScalar("Probability of initial state", calcProbOfBasisState(qureg, 0));

    // verify results by uninterrupted higher-order simulation to target time
    initPlusState(qureg);
    applyTrotterizedUnitaryTimeEvolution(qureg, hamil, dt*steps, order+2, reps*steps);
    reportScalar("final <O>", calcExpecPauliStrSum(qureg, observ));

    // clean up
    destroyQureg(qureg);
    destroyPauliStrSum(hamil);
    destroyPauliStrSum(observ);
    finalizeQuESTEnv();
    return 0;
}
