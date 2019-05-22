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

Independent of which platform you'll run your simulation on (multicore CPUS, a GPU, or over a network), your QuEST code will look the same, compile with the same [makefile](https://github.com/quest-kit/QuEST/blob/master/makefile), and use the same [API](https://quest-kit.github.io/QuEST/QuEST_8h.html).

Here's a simulation of a very simple circuit which measures  ![equation](https://latex.codecogs.com/gif.latex?C_0%28X_1%29%20H_0%20%7C00%5Crangle).
```C
#include <QuEST.h>

int main(int narg, char *varg[]) {

  // load QuEST
  QuESTEnv env = createQuESTEnv();
  
  // create 2 qubits in the hadamard state
  Qureg qubits = createQureg(2, env);
  initPlusState(qubits);
	
  // apply circuit
  hadamard(qubits, 0);
  controlledNot(qubits, 0, 1);
  measure(qubits, 1);
	
  // unload QuEST
  destroyQureg(qubits, env); 
  destroyQuESTEnv(env);
  return 0;
}
```
Of course, this code doesn't output anything!


----------------------

Let's walk through a more sophisticated circuit.

We first construct a quest environment, which abstracts away any preparation of multithreading, distribution or GPU-acceleration strategies.
```C
QuESTEnv env = createQuESTEnv();
```

We then create a quantum register, in this case of 3 qubits.
```C
Qureg qubits = createQureg(3, env);
```
and set it to be in the zero state.
```C
initZeroState(qubits);
```
We can create multiple `Qureg` instances, and QuEST will sort out allocating memory for the state-vectors, even over networks! Note we can replace `createQureg` with `createDensityQureg`, a more powerful density matrix representation which can store mixed states!

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
qreal prob = getProbAmp(qubits, 7);
printf("Probability amplitude of |111>: %lf\n", prob);
```
Here, `qreal` is a floating point number (e.g. `double`). The state-vector is stored as `qreal`s so that we can change its precision without any recoding, by changing `PRECISION` in the [makefile](../makefile)

How probable is measuring our final qubit (2) in outcome `1`?
```C
prob = calcProbOfOutcome(qubits, 2, 1);
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
destroyQureg(qubits, env);
destroyQuESTEnv(env);
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

QuEST uses the [Mersenne Twister](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) algorithm to generate random numbers used for randomly collapsing the state-vector. The user can seed this RNG using `seedQuEST(arrayOfSeeds, arrayLength)`, otherwise QuEST will by default (through `seedQuESTDefault()`) create a seed from the current time and the process id.

----------------------------

# Compiling

QuEST uses CMake (3.1 or higher) as its build system.

To compile, make sure your circuit code is accessible from the root QuEST directory.
In the root directory, initially build using
```bash
mkdir build
cd build
cmake -DUSER_SOURCE="myCode1.c;myCode2.c" ..
make
```
Paths to target sources are set as a semi-colon separated list of paths to said sources relative to the root QuEST directory.

If you wish your executable to be named something other than `demo`, you can set this too by using:
```bash
cmake -DOUTPUT_EXE="myExecutable" ..
```

When using the cmake command as above, the -D[VAR=VALUE] option can be passed other options to further configure your build.

To compile your code to run on multi-CPU systems use
```bash
cmake -DDISTRIBUTED=1 ..
```

To compile for GPU, use
```bash
cmake -DGPUACCELERATED=1 -DGPU_COMPUTE_CAPABILITY=[COMPUTE_CAPABILITY] ..
```

Where COMPUTE_CAPABILITY is the compute cabability of your GPU, written without a decimal point. This can can be looked up at the [NVIDIA website](https://developer.nvidia.com/cuda-gpus). The default value is 30.

By default, QuEST will compile with OpenMP parallelism enabled if an OpenMP compatible compiler and OpenMP library can be found on your system (e.g. [GCC 4.9](https://gcc.gnu.org/gcc-4.9/changes.html)). Using distribution requires an MPI implementation is installed on your system, and GPU acceleration requires a CUDA compiler (`nvcc`). We've made a comprehensive list of compatible compilers which you can view [here](../tests/compilers/compatibility.md). This does not change your `COMPILER` setting - the makefile will choose the appropriate MPI and CUDA wrappers automatically.

You can additionally customise the precision with which the state-vector is stored.
```bash
cmake -DPRECISION=2 ..
```
Using greater precision means more precise computation but at the expense of additional memory requirements and runtime.
Checking results are unchanged when altaring the precision can be a great test that your calculations are sufficiently precise.


Please note that cmake caches these changes (per directory) so for any subsequent builds you should just type `make` from the build directory and the previously defined settings will be applied. If any parameters require changing, these can be redefined by:
```
cmake -D[VAR=VALUE] ..
```
as one would do in the initial configuration.

For a full list of available configuration parameters, use
```bash
cmake -LH ..
```

For manual configuration (not recommended) you can change the `CMakeLists.txt` in the root QuEST directory.

----------------------------

# Running unit tests

To confirm that QuEST has been compiled and is running correctly on your platform, unit tests can be run from the build folder with the command

```bash
make test
```

This will report whether the QuEST library has been built correctly and whether all unit tests have passed successfully. In case of failures, see utilities/QuESTLog.log for a detailed report. 

Tests will automatically run in distributed mode on four processes if -DDISTRIBUTED=1 is set at compile time, and on GPU if -DGPUACCELERATED=1 is set at compile time. In order to set the number of process on which the tests should be run, set:
```bash
cmake -DMPIEXEC_MAX_NUMPROCS=4 ..
```

Note, the most common reason for unit tests failing on a new platform is running on a GPU with the incorrect GPU_COMPUTE_CAPABILITY. Remember to specify this at compile time for [your device](https://developer.nvidia.com/cuda-gpus). Eg, for a P100, use


```bash
mkdir build
cd build
cmake -DGPUACCELERATED=1 -DGPU_COMPUTE_CAPABILIty=60 ..
make test
```

----------------------------

# Running

## locally

You can then call your code. From the build directory:
```bash
./myExecutable
```
If multithreading functionality was found when compiling, you can control how many threads your code uses by setting `OMP_NUM_THREADS`, ideally to the number of available cores on your machine
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

For example, the [above code](tutorial_example.c) can be split over 4 MPI nodes (each with 8 cores) on a SLURM system using the following SLURM submission script
```bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1

module purge
module load mvapich2

mkdir build
cd build
cmake -DDISTRIBUTED=1 ..
make

export OMP_NUM_THREADS=8
mpirun ./myExecutable
```
or a PBS submission script like
```bash
#PBS -l select=4:ncpus=8

module purge
module load mvapich2

mkdir build
cd build
cmake -DDISTRIBUTED=1 ..
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

mkdir build
cd build
cmake -DGPUACCELERATED=1 -DGPU_COMPUTE_CAPABILITY=[Compute capability] ..
make

./myExecutable
```

On each platform, there is no change to our source code or our QuEST interface. We simply recompile, and QuEST will utilise the available hardware (a GPU, shared-memory or distributed CPUs) to speedup our code.

Note that parallelising with MPI (`-DDISTRIBUTED=1`) will mean all code in your source file will be repeated on every node. To execute some code (e.g. printing) only on one node, do
```C
if (env.rank == 0)
    printf("Only one node executes this print!");
```
Such conditions are valid and always satisfied in code run on a single node.
