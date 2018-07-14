#include <QuEST.h>

int main (int narg, char *varg[]) {

    QuESTEnv env;
	MultiQubit qubits; 
    initQuESTEnv(&env);
    createMultiQubit(&qubits, 3, env);
    initStateZero(qubits);

    hadamard(qubits, 0);
    controlledNot(qubits, 0, 1);
    rotateY(qubits, 2, .1);

    destroyMultiQubit(qubits, env); 
    closeQuESTEnv(env);

    return 0;
}
