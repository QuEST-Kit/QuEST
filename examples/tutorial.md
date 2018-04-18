QuEST Tutorial
======

**Table of Contents**
- [Coding](#coding)
- [Compiling](#compiling)
- [Running QuEST through a job submission system](#running-quest-through-a-job-submission-system)


# Coding

Independent of which platform you'll run your simulation on (multicore machine, distributed machines, a GPU), your QuEST code will look the same. View the API [here](https://aniabrown.github.io/QuEST/qubits_8h.html)

> Currently, a limited set of functions are available on the GPU, listed [here](https://aniabrown.github.io/QuEST_GPU/qubits_8h.html)

A simulation of a very simple circuit looks like
```C
#include "qubits.h"

int main(int narg, char *varg[]) {

  // load QuEST
  QuESTEnv env;
  initQuESTEnv(&env);
  
  // create 2 qubits in the zero state
  MultiQubit qubits; 
  createMultiQubit(&qubits, 2, env);
  initStateZero(&qubits);
	
  // apply circuit
  hadamard(qubits, 0);
  controlledNot(qubits, 0, 1);
  measure(qubits, 1);
	
  // unload QuEST
  destroyMultiQubit(qubits, env); 
  closeQuESTEnv(env);
  return 0;
}
```
though of course this doesn't output anything!


----------------------

Let's walk through a more sophisticated circuit.

We first construct a quest environment, which abstracts away all the MPI and GPU preparation.
```C
QuESTEnv env;
initQuESTEnv(&env);
```

We then create a quantum register, in this case of 3 qubits.
```C
MultiQubit qubits; 
createMultiQubit(&qubits, 3, env);
```
and set it to be in the zero state.
```C
initStateZero(&qubits);
```

We're now ready to apply some gates to our qubits, which have indices 0, 1 and 2.
When applying a gate, we pass along the quantum register.
```C
hadamard(qubits, 0);
controlledNot(qubits, 0, 1);
rotateY(qubits, 2, .1);
```

Some gates allow us to specify a general number of control qubits
```C
multiControlledPhaseGate(qubits, (int []){0, 1, 2}, 3);
```

We can specify general single-qubit unitary operations as 2x2 matrices
```C
ComplexMatrix2 u;
u.r0c0 = (Complex) {.real=.5, .imag= .5};
u.r0c1 = (Complex) {.real=.5, .imag=-.5}; 
u.r1c0 = (Complex) {.real=.5, .imag=-.5};
u.r1c1 = (Complex) {.real=.5, .imag= .5};
unitary(qubits, 0, u);
```
or more compactly, foregoing the global phase factor,
```C
Complex a, b;
a.real = .5; a.imag =  .5;
b.real = .5; b.imag = -.5;
compactUnitary(qubits, 1, a, b);
```
or even more compactly, as a rotation around an arbitrary axis on the Bloch-sphere
```C
Vector v;
v.x = 1; v.y = 0; v.z = 0;
rotateAroundAxis(qubits, 2, 3.14/2, v);
```

We can controlled-apply general unitaries
```C
controlledCompactUnitary(qubits, 0, 1, a, b);
```
even with multiple control qubits
```C
multiControlledUnitary(qubits, (int []){0, 1}, 2, 2, u);
```

What has this done to the probability of the basis state |111> = |7>?
```C
REAL prob = getProbEl(qubits, 7);
printf("Probability amplitude of |111>: %f\n", prob);
```

How probable is measuring our final qubit (2) in outcome `1`?
```C
prob = findProbabilityOfOutcome(qubits, 2, 1);
printf("Probability of qubit 2 being in state 1: %f\n", prob);
```

Let's measure the first qubit, collapsing it to the classical 0, 1 basis
```C
int outcome = measure(qubits, 0);
printf("Qubit 0 was measured in state %d\n", outcome);
```
and now measure our final qubit, while also learning of the probability of its outcome.
```
outcome = measureWithStats(qubits, 2, &prob);
printf("Qubit 2 collapsed to %d with probability %f\n", outcome, prob);
```

At the conclusion of our circuit, we should free up our state-vector.
```C
destroyMultiQubit(qubits, env); 
closeQuESTEnv(env);
return 0;
```

Executing all the [code above](tutorialExample.c) simulates the below circiut

<img src="https://qtechtheory.org/wp-content/uploads/2018/02/github_circuit.png" alt="A quantum circuit" width=400px >

and after compiling (see section below), gives psuedo-random output

> ```
> Probability amplitude of |111>: 0.498751
> Probability of qubit 2 being in state 1: 0.749178
> Qubit 0 was measured in state 1
> Qubit 2 collapsed to 1 with probability 0.998752
> ```

> ```
> Probability amplitude of |111>: 0.498751
> Probability of qubit 2 being in state 1: 0.749178
> Qubit 0 was measured in state 0
> Qubit 2 collapsed to 1 with probability 0.499604
> ```

> standby for seeding instructions!

----------------------------

# Compiling

To compile, make sure there is a [makefile](makefile) in the same folder as your circuit code. An example makefile and the tutorial circuit file will be present in the root directory when you download the code.  

Edit the makefile, letting `MY_C_SOURCES` be a space-separated list of your source files, `EXE` be the name of the output executable, and `QUEST_DIR` point to the folder which contains `qubits.h`. 

```bash
# COMPILER options: {GNU, INTEL}
COMPILER = GNU

# EXECUTABLE TO GENERATE
EXE = myExecutable

# USER SOURCE FILES
MY_C_SOURCES = myCode1 myCode2 myCode3

# ENABLE DISTRIBUTED PROCESSING (ALLOWS MULTIPLE NODES)
USE_MPI=0

# ENABLE MULTIPLE THREADS PER PROCESS (RECOMMENDED)
USE_OPENMP=1

# PATH TO QUEST LIBRARY SOURCES FROM ROOT DIRECTORY
QUEST_DIR = QuEST
```
You can then compile your code by calling
```bash
make
```
at the terminal, in the directory of your code. For the above example, this performs
```bash
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c QuEST/qubits.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c QuEST/mt19937ar.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c myCode1.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c myCode2.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c QuEST/qubits_env_local.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -o myExecutable qubits.o mt19937ar.o myCode1.o myCode2.o qubits_env_local.o -lm
```
You can then call your code
```bash
./myExecutable
```
To control how many threads your code uses, modify `OMP_NUM_THREADS`. You should set it to the number of available cores on your machine
```bash
export OMP_NUM_THREADS=8
./myExecutable
```
QuEST will automatically allocate work between the given number of threads to speedup your simulation.

To compile for distributed simulation, set `USE_MPI=1`, and QuEST will spread the quantum state vector between available nodes. Launch with the appropriate mpi launcher for your system on a power of 2 number of processes, for example
```bash
mpirun -np 8 ./myExecutable
```

> For the moment, compiling for GPU use requires C++ source code (that your files are C++ compatible and have the `.cpp` extension). An example makefile for the GPU is provided [here](https://github.com/aniabrown/QuEST_GPU/blob/master/examples/makefile).

---------------------

# Running QuEST through a job submission system

There are no special requirements for running QuEST through job submission systems. Just call `./myExecutable` as you would any other binary.

Be sure to set `OMP_NUM_THREADS` appropriately, and that you target the hardware your job will ultimately run on when compiling (otherwise simply compile at runtime using the makefile, just as above).

For example, the [above code](tutorialExample.c) can be split over 4 MPI nodes (each with 8 cores) by setting `USE_MPI=1` (and `USE_OPENMP=1`) in the makefile, and with a SLURM submission script like
```bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1

module purge
module load mvapich2

make clean
make

export OMP_NUM_THREADS=8
mpirun ./myExecutable
```
or a PBS submission script like
```bash
#PBS -l select=4:ncpus=8

make clean
make

export OMP_NUM_THREADS=8
aprun -n 4 -d 8 -cc numa_node ./myExecutable
```

Running QuEST on a GPU partition is similarly easy (though currently requiring an [alternate makefile](https://github.com/aniabrown/QuEST_GPU/blob/master/examples/makefile)) in SLURM
```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1 

#SBATCH --partition=gpu    ## name may vary

module purge
module load cuda  ## name may vary

make clean
make

./myExecutable
```

On each platform, there is no change to our source code or our QuEST interface. We simply recompile, and QuEST will utilise the available hardware (a GPU, shared-memory or distributed CPUs) to speedup our code.

Note that parallelising with MPI will mean all code in your source file will be repeated on every node. To execute some code (e.g. printing) only on one node, do
```C
if (env.rank == 0)
    printf("Only one node executes this print!");
```
Such conditions are valid and always satisfied in code run on a single node.
