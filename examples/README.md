QuEST Tutorial
======

**Table of Contents**
- [Coding](#coding)
- [Compiling](#compiling)
- [Running](#running)


# Coding

QuEST can be used in your C or C++ code, simply by including
```C
#include <QuEST.h>
```

Independent of which platform you'll run your simulation on (multicore CPUS, a GPU, or over a network), your QuEST code will look the same, compile with the same [makefile](https://github.com/TysonRayJones/QuEST/blob/master/makefile), and use the same [API](https://tysonrayjones.github.io/QuEST/QuEST_8h.html).

Here's a simulation of a very simple circuit which measures  ![equation](https://latex.codecogs.com/gif.latex?C_0%28X_1%29%20H_0%20%7C00%5Crangle).
```C
#include <QuEST.h>

int main(int narg, char *varg[]) {

  // load QuEST
  QuESTEnv env;
  initQuESTEnv(&env);
  
  // create 2 qubits in the zero state
  MultiQubit qubits; 
  createMultiQubit(&qubits, 2, env);
  initStateZero(qubits);
	
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
Of course, this code doesn't output anything!


----------------------

Let's walk through a more sophisticated circuit.

We first construct a quest environment, which abstracts away any preparation of multithreading, distribution or GPU-acceleration strategies.
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
initStateZero(qubits);
```
We can create multiple `MultiQubit` instances, and QuEST will sort out allocating memory for the state-vectors, even over networks!

We're now ready to apply some gates to our qubits, which in this case have indices 0, 1 and 2.
When applying a gate, we pass along which quantum register to operate upon.
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
// sqrt(X) with a pi/4 global phase
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
even with multiple control qubits!
```C
multiControlledUnitary(qubits, (int []){0, 1}, 2, 2, u);
```

What has this done to the probability of the basis state |111> = |7>?
```C
REAL prob = getProbEl(qubits, 7);
printf("Probability amplitude of |111>: %lf\n", prob);
```
Here, `REAL` is a floating point number (e.g. `double`). The state-vector is stored as `REAL`s so that we can change its precision without any recoding, by configuring [QuEST_precision.h](../QuEST/QuEST_precision.h)

How probable is measuring our final qubit (2) in outcome `1`?
```C
prob = findProbabilityOfOutcome(qubits, 2, 1);
printf("Probability of qubit 2 being in state 1: %f\n", prob);
```

Let's measure the first qubit, randomly collapsing it to 0 or 1
```C
int outcome = measure(qubits, 0);
printf("Qubit 0 was measured in state %d\n", outcome);
```
and now measure our final qubit, while also learning of the probability of its outcome.
```
outcome = measureWithStats(qubits, 2, &prob);
printf("Qubit 2 collapsed to %d with probability %f\n", outcome, prob);
```

At the conclusion of our circuit, we should free up the memory used by our state-vector.
```C
destroyMultiQubit(qubits, env); 
closeQuESTEnv(env);
return 0;
```

The effect of the [code above](tutorial_example.c) is to simulate the below circuit

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

QuEST uses the [Mersenne Twister](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) algorithm to generate random numbers used for randomly collapsing the state-vector. The user can seed this RNG using `seedQuEST(arrayOfSeeds, arrayLength)`, otherwise QuEST will by default (through `seedQuESTDefault()`) create a seed from the current time, the process id, and the hostname.

----------------------------

# Compiling

To compile, copy the [makefile](../makefile) into the same folder as your circuit code. Adjust the *User Settings* section to configure compilation. You'll need to set
```bash
# name of the executable to create
EXE = myExecutable

# space-separated names (no file type) of all user source files (.c or .cpp) in the root directory
SOURCES = myCode1 myCode2

# path to QuEST library from root directory
QUEST_DIR = path/to/QuEST
```

Next, indicate which compiler you wish to use. For example, to use the default compiler on OSX:
```bash
# compiler to use, which should support both C and C++, to be wrapped by GPU/MPI compilers
COMPILER = clang

# type of above compiler, one of {GNU, INTEL, CLANG}, used for setting compiler flags
COMPILER_TYPE = CLANG
```

To compile your code to run on multicore/multi-CPU systems, and/or for distributed systems, or on GPUs, simply set the appropriate variables
```bash
# hardwares to target: 1 means use, 0 means don't use
MULTITHREADED = 0
DISTRIBUTED = 0
GPUACCELERATED = 0
```
Note that using multithreading requires an OpenMP compatible compiler (e.g. [GCC 4.9](https://gcc.gnu.org/gcc-4.9/changes.html)), using distribution requires an MPI compiler (`mpicc`)is installed on your system, and GPU acceleration requires a CUDA compiler (`nvcc`). We've made a comprehensive list of compatible compilers which you can view [here](../tests/compilers/compatibility.md). This does not change your `COMPILER` setting - the makefile will choose the appropriate MPI and CUDA wrappers automatically.

> Note also that GPU users must additionally specify the the *Compute Capability* of their GPU, which can be looked up at the [NVIDIA website](https://developer.nvidia.com/cuda-gpus)
> ```bash
> GPU_COMPUTE_CAPABILITY = 30
> ```
> An incorrect *Compute Capability* will lead to drastically incorrect computations. You can check if you've set the right *Compute Capability* by running the unit tests via `tests/runTests.c`.

You're now ready to compile your code by entering
```bash
make
```
at the terminal, in the directory of your code. For the above example, this performs
```bash
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c path/to/QuEST/CPU/QuEST.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c path/to/QuEST/mt19937ar.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c myCode1.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c myCode2.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -c path/to/QuEST/CPU/QuEST_env_local.c
gcc -O2 -std=c99 -mavx -Wall -fopenmp -o myExecutable QuEST.o mt19937ar.o myCode1.o myCode2.o QuEST_env_local.o -lm
```

----------------------------

# Running

## locally

You can then call your code
```bash
./myExecutable
```
If you enabled multithreading when compiling, you can control how many threads your code uses by setting `OMP_NUM_THREADS`, ideally to the number of available cores on your machine
```bash
export OMP_NUM_THREADS=8
./myExecutable
```
QuEST will automatically allocate work between the given number of threads to speedup your simulation.

If you compiled in distributed mode, your code can be run over a network (here, over 8 machines) using
```bash
mpirun -np 8 ./myExecutable
```
This will, if enabled, also utilise multithreading on each node with as many threads set in `OMP_NUM_THREADS`.

If you compiled for a GPU connected to your system, simply run
```bash
./myExecutable
```
as normal!

## through a job submission system

There are no special requirements for running QuEST through job submission systems. Just call `./myExecutable` as you would any other binary.

For example, the [above code](tutorial_example.c) can be split over 4 MPI nodes (each with 8 cores) by setting `DISTRIBUTED = 1` (and `MULTITHREADED = 1`) in the makefile, and writing a SLURM submission script
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

Running QuEST on a GPU partition is similarly easy in SLURM
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

Note that parallelising with MPI (`DISTRIBUTED = 1`) will mean all code in your source file will be repeated on every node. To execute some code (e.g. printing) only on one node, do
```C
if (env.rank == 0)
    printf("Only one node executes this print!");
```
Such conditions are valid and always satisfied in code run on a single node.
