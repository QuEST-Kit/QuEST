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

Independent of which platform you'll run your simulation on (multicore CPUs, a GPU, or over a network), your QuEST code will look the same, compile with the same [makefile](https://github.com/quest-kit/QuEST/blob/master/makefile), and use the same [API](https://quest-kit.github.io/QuEST/QuEST_8h.html).

Here's a simulation of a very simple circuit which measures  ![equation](https://latex.codecogs.com/gif.latex?C_0%28X_1%29%20H_0%20%7C00%5Crangle).
```C
#include <QuEST.h>

int main() {

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

We first construct a QuEST environment, which abstracts away any preparation of multithreading, distribution or GPU-acceleration strategies.
```C
QuESTEnv env = createQuESTEnv();
```

We then create a quantum register, in this case of 3 qubits.
```C
Qureg qubits = createQureg(3, env);
```
and prepare it in the zero state.
```C
initZeroState(qubits);
```
We can create multiple `Qureg` instances, and QuEST will sort out allocating memory for the state-vectors, even over networks! If we wanted to simulate noise in our circuit, we can replace `createQureg` with `createDensityQureg` to create a more powerful density matrix capable of representing mixed states.

We're now ready to apply some gates to our qubits, which in this case have indices `0`, `1` and `2`.
When applying a gate, we pass along which quantum register to operate upon.
```C
hadamard(qubits, 0);
controlledNot(qubits, 0, 1);
rotateY(qubits, 2, .1);
```

Some gates allow us to specify a general number of control qubits
```C
multiControlledPhaseGate(qubits, (int[]) {0, 1, 2}, 3);
```

We can specify general single-qubit unitary operations as 2x2 matrices
```C
// sqrt(X) with a pi/4 global phase
ComplexMatrix2 u = {
    .real = {{.5, .5}, { .5,.5}},
    .imag = {{.5,-.5}, {-.5,.5}}};
unitary(qubits, 0, u);
```
or more compactly, foregoing the global phase factor,
```C
Complex a = {.real = .5, .imag = .5};
Complex b = {.real = .5, .imag =-.5};
compactUnitary(qubits, 1, a, b);
```
or even more compactly, as a rotation around an arbitrary axis on the Bloch-sphere
```C
Vector v = {.x=1, .y=0, .z=0};
rotateAroundAxis(qubits, 2, 3.14/2, v);
```

We can controlled-apply general unitaries
```C
controlledCompactUnitary(qubits, 0, 1, a, b);
```
even with multiple control qubits!
```C
multiControlledUnitary(qubits, (int[]) {0, 1}, 2, 2, u);
```

What has this done to the probability of the basis state `|111>` ?
```C
qreal prob = getProbAmp(qubits, 7);
printf("Probability amplitude of |111>: %lf\n", prob);
```
Here, `qreal` is an alias for a real floating point number, like `double`. This is to keep our code precision agnostic, so that we may change the numerical precision at compile time (by setting macro `PRECISION`) without any changes to our code. Changing the precision can be useful in verifying numerical convergences or studying rounding errors.

How probable is measuring our final qubit (with index `2`) in outcome `1`?
```C
prob = calcProbOfOutcome(qubits, 2, 1);
printf("Probability of qubit 2 being in state 1: %f\n", prob);
```

Let's measure the first qubit, randomly collapsing into outcome `0` or `1`
```C
int outcome = measure(qubits, 0);
printf("Qubit 0 was measured in state %d\n", outcome);
```
and now measure our final qubit, while also learning of the probability of its outcome.
```C
outcome = measureWithStats(qubits, 2, &prob);
printf("Qubit 2 collapsed to %d with probability %f\n", outcome, prob);
```

At the conclusion of our circuit, we should free up the memory used by our quantum registers.
```C
destroyQureg(qubits, env);
destroyQuESTEnv(env);
```

The effect of the [code above](tutorial_example.c) is to simulate the below circuit

<img src="https://github.com/QuEST-Kit/QuEST/raw/readme_update/examples/tutorial_circuit.png" width="50%"> <br>

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

QuEST uses the [Mersenne Twister](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) algorithm to generate random numbers used for randomly collapsing quantum states. The user can seed this RNG using `seedQuEST(arrayOfSeeds, arrayLength)`, otherwise QuEST will by default (through `seedQuESTDefault()`) create a seed from the current time and the process id.

----------------------------

# Compiling

QuEST uses CMake (`3.7` or higher) as its build system. Configure the build by supplying the below `-D[VAR=VALUE]` options after the `cmake ..` command. If you wish to avoid CMake, you can otherwise use GNUMake directly with the provided [makefile](makefile).

- To compile, make sure your circuit code is accessible from the root QuEST directory.
In the root directory, initially build using
```bash
mkdir build
cd build
cmake .. -DUSER_SOURCE="myCode1.c;myCode2.c"
make
```
Paths to target sources are set as a semi-colon separated list of paths to said sources relative to the root QuEST directory.

- To set the compilers used by cmake (to e.g. `gcc-6`), use
  ```bash 
   -DCMAKE_C_COMPILER=gcc-6
  ```
  and similarly to set the C++ compiler (as used in GPU mode), use
  ```bash 
   -DCMAKE_CXX_COMPILER=g++-6
  ```

- If you wish your executable to be named something other than `demo`, you can set this too by adding argument:
```bash
 -DOUTPUT_EXE="myExecutable" 
```



- To compile your code to use multithreading, for parallelism on multi-core or multi-CPU systems, use
```bash
 -DMULTITHREADED=1
```
Before launching your executable, set the number of participating threads using `OMP_NUM_THREADS`. For example,
```bash
export OMP_NUM_THREADS=16
./myExecutable
```

- To compile your code to run on distributed or networked systems use
```bash
 -DDISTRIBUTED=1
```
Depending on your MPI implementation, your executable can be launched via
```bash 
mpirun -np [NUM_NODES] [EXEC]
```
where `[NUM_NODES]` is the number of distributed compute nodes to use, and `[EXEC]` is the name of your executable. Note that QuEST *hybridises* multithreading and distribution. Hence you should set `[NUM_NODES]` to equal exactly the number of distinct compute nodes (which don't share memory), and set `OMP_NUM_THREADS` as above to assign the number of threads used on *each* compute node.

- To compile for GPU, use
```bash
 -DGPUACCELERATED=1 -DGPU_COMPUTE_CAPABILITY=[COMPUTE_CAPABILITY] ..
```
were `[COMPUTE_CAPABILITY]` is the compute cabability of your GPU, written without a decimal point. This can can be looked up at the [NVIDIA website](https://developer.nvidia.com/cuda-gpus).
> Note that CUDA is not compatible with all compilers. To force `cmake` to use a 
> compatible compiler, override `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`

- You can additionally customise the floating point precision used by QuEST's `qreal` type, via
```bash
 -DPRECISION=1
 -DPRECISION=2
 -DPRECISION=4
```
which uses single (`qreal = float`), double (`qreal = double`) and quad (`qreal = long double`) respectively.
Using greater precision means more precise computation but at the expense of additional memory requirements and runtime.
Checking results are unchanged when switching the precision can be a great test that your calculations are sufficiently precise.


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
