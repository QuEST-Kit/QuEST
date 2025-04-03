<!--
  Instructions for running QuEST with different parallelisations

  @author Tyson Jones
-->

# Run

Running your [compiled](compile.md) QuEST application can be as straightforward as running any other executable, though some additional steps are needed to make use of hardware acceleration. This page how to launch your own QuEST applications on different platforms, how to run the examples and unit tests, how to make use of multithreading, GPU-acceleration, distribution and supercomputer job schedulers, and monitor the hardware utilisation.

> [!NOTE]
> This page assumes you are working in a `build` directory into which all executables have been compiled.


---------------------


## Examples

> See [`compile.md`](compile.md#examples) for instructions on compiling the examples.

The example source codes are located in [`examples/`](/examples/) and are divided into subdirectories, e.g.
```
examples/
    krausmaps/
        initialisation.c
        initialisation.cpp
    reporters/
        env.c
        env.cpp
        matrices.c
        matrices.cpp
    ...
```
where `file.c` and `file.cpp` respectively demonstrate QuEST's `C11` and `C++14` interfaces.
These files are [compiled](compile.md#examples) into executables of the same name, respectively prefixed with `c_` or `cpp_`, and saved in subdirectories of `build` which mimic the structure of `examples/`. E.g.
```
build/
    examples/
        krausmaps/
            c_initialisation
            cpp_initialisation
        reporters/
            c_env
            cpp_env
            c_matrices
            cpp_matrices
    ...
```
Most of these executables can be run directly from within `build`, e.g.
```bash
./examples/reporters/cpp_paulis
```
while others require command-line arguments:
```bash
./examples/reporters/c_env

# output
Must pass single cmd-line argument:
  1 = serial
  2 = multithreaded
  3 = GPU-accelerated
  4 = distributed
  5 = all
  6 = auto
```


---------------------

## Tests

> See [`compile.md`](compile.md#tests) for instructions on compiling the `v4` and `v3` unit tests.

### v4

QuEST's unit and integration tests are compiled into executable `tests` within the `tests/` subdirectory, and can be directly run from within the `build` folder via
```bash
./tests/tests
```
which should, after some time, output something like
```
QuEST execution environment:
  precision:       2
  multithreaded:   1
  distributed:     1
  GPU-accelerated: 1
  cuQuantum:       1
  num nodes:       16
  num qubits:      6
  num qubit perms: 10

Tested Qureg deployments:
  GPU + MPI

Randomness seeded to: 144665856
===============================================================================
All tests passed (74214 assertions in 240 test cases)
```

This `tests` binary accepts all of the [Catch2 CLI arguments](https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md), for example to run specific tests
```bash
./tests/tests applyHadamard
```
or all tests within certain groups
```bash
./tests/tests "[qureg],[matrices]"
```
or specific test sections and subsections:
```bash
./tests/tests -c "validation" -c "matrix uninitialised"
```

If the tests were compiled with [distribution enabled](compile.md#distribution), they can distributed via
```bash
mpirun -np 8 ./tests/tests
```

Alternatively, the tests can be run through [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html) within the `build` directory via either
```bash
ctest
```
```bash
make test
```
which will log each passing test live, outputting something like
```
Test project /build
        Start   1: calcExpecPauliStr
  1/240 Test   #1: calcExpecPauliStr .............................   Passed   14.03 sec
        Start   2: calcExpecPauliStrSum
  2/240 Test   #2: calcExpecPauliStrSum ..........................   Passed   10.06 sec
        Start   3: calcExpecNonHermitianPauliStrSum
  3/240 Test   #3: calcExpecNonHermitianPauliStrSum ..............   Passed   10.34 sec
        Start   4: calcProbOfBasisState
  4/240 Test   #4: calcProbOfBasisState ..........................   Passed    0.33 sec
        Start   5: calcProbOfQubitOutcome
  5/240 Test   #5: calcProbOfQubitOutcome ........................   Passed    0.12 sec
        Start   6: calcProbOfMultiQubitOutcome
  6/240 Test   #6: calcProbOfMultiQubitOutcome ...................   Passed   15.07 sec
        Start   7: calcProbsOfAllMultiQubitOutcomes
...
```
Alas tests launched in this way cannot be deployed with distribution.


### v3

The deprecated tests, when [compiled](compile.md#v3), can be run from the `build` directory via
```bash
./tests/deprecated/dep_tests
```
which accepts the same [Catch2 CLI arguments](https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md) as the `v4` tests above, and can be distributed the same way.

To launch the tests with [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html), run
```bash
cd tests/deprecated
ctest
```

> [!CAUTION]
> The deprecated unit tests are non-comprehensive and the deprecated API should not be relied upon, for it may introduce
> undetected corner-case bugs. Please only use the deprecated API and tests for assistance porting your application
> from QuEST v3 to v4.





---------------------

## Multithreading

> [!NOTE]
> Parallelising QuEST over multiple cores and CPUs requires first compiling with 
> multithreading enabled, as detailed in [`compile.md`](compile.md#multithreading). 

### Choosing threads

The number of [threads](https://www.openmp.org/spec-html/5.0/openmpsu1.html) to use is decided before launching the compiled executable, using the [`OMP_NUM_THREADS`](https://www.openmp.org/spec-html/5.0/openmpse50.html) environment variable.

```bash
OMP_NUM_THREADS=32 ./myexec
```
```bash
export OMP_NUM_THREADS=32
./myexec
```

It is prudent to choose as many threads as your CPU(s) have total hardware threads or cores, which need not be a power of `2`. One can view this, and verify the number of available threads at runtime, by calling [`reportQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496) which outputs a subsection such as
```
  [cpu]
    numCpuCores.......10 per machine
    numOmpProcs.......10 per machine
    numOmpThrds.......32 per node
```
<!-- the doxygen-doc hyperlink above includes a hash of the function name which should be unchanging! -->

> [!NOTE]
> When running [distributed](#distribution), variable `OMP_NUM_THREADS` specifies the number of threads _per node_ and so should ordinarily be the number of hardware threads (or cores) **_per machine_**.


### Monitoring utilisation


The availability of multithreaded deployment can also be checked at runtime using [`reportQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496), which outputs something like:
```bash
  [compilation]
    isOmpCompiled...........1
  [deployment]
    isOmpEnabled............1
```
where `Omp` signifies OpenMP and the two `1` respectively indicate it has been compiled and runtime enabled.

Like all programs, the CPU utilisation of a running QuEST program can be viewed using

| OS    | Program | Method |
| -------- | ------- | ---- |
| Linux  | [HTOP](https://htop.dev/) | Run `htop` in terminal  |
| MacOS |  [Activity Monitor](https://support.apple.com/en-gb/guide/activity-monitor/welcome/mac) | Place on dock > right click icon > `Monitors` > `Show CPU usage` (see [here](https://stackoverflow.com/questions/50260592/getting-each-of-the-cpu-cores-usage-via-terminal-in-macos)) |
| Windows  |  [Task Manager](https://learn.microsoft.com/en-us/shows/inside/task-manager) |`Performance` > `CPU` > right click graph > `Change graph to` > `Logical processors` (see [here](https://superuser.com/questions/1398696/how-to-see-usage-of-each-core-in-windows-10)) |


Note however that QuEST will not leverage multithreading at runtime when either:
- Your `Qureg` created with [`createQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#gab3a231fba4fd34ed95a330c91fcb03b3) was too small to invoke automatic multithreading.
- You called [`initCustomQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga485268e52f838743357e7a4c8c241e57) and disabled multithreading for all subsequently-created `Qureg`.

Usage of multithreading can be (inadvisably) forced using [`createForcedQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga619bbba1cbc2f7f9bbf3d3b86b3f02be) or [`createCustomQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga849971f43e246d103da1731d0901f2e6).


### Improving performance

Performance may be improved by setting other [OpenMP variables](https://www.openmp.org/spec-html/5.0/openmpch6.html). Keep in mind that for large `Qureg`, QuEST's runtime is dominated by the costs of modifying large memory structures during long, uninterrupted loops: namely the updating of statevector amplitudes. Some sensible settings include

- [`OMP_DYNAMIC`](https://www.openmp.org/spec-html/5.0/openmpse51.html) `=false` to disable the costly runtime migration of threads between cores.
- [`OMP_PROC_BIND`](https://www.openmp.org/spec-html/5.0/openmpse52.html) `=spread` to (attemptedly) give threads their own caches (see [here](https://developer.arm.com/documentation/102580/0100/Control-the-placement-of-OpenMP-threads)).
  > Replace this with [`KMP_AFFINITY`](https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/developer-guide-reference/2023-0/thread-affinity-interface.html) on Intel compilers.
- [`OMP_PLACES`](https://www.openmp.org/spec-html/5.0/openmpse53.html) `=threads` to allocate each spawned thread to a CPU hardware thread.
  > Alternatively set `=cores` to assign one thread per core, helpful when the hardware threads interfere (e.g. due to caching conflicts).


OpenMP experts may further benefit from knowing that QuEST's multithreaded source code, confined to [`cpu_subroutines.cpp`](/quest/src/cpu/cpu_subroutines.cpp), is almost exclusively code similar to
```C++
#pragma omp parallel for if(qureg.isMultithreaded)
for (qindex n=0; n<numIts; n++)
```
```C++
#pragma omp parallel for reduction(+:val)
for (qindex n=0; n<numIts; n++)
    val += 
```
and never specifies [`schedule`](https://rookiehpc.org/openmp/docs/schedule/index.html) nor invokes setters in the [runtime library routines](https://www.openmp.org/spec-html/5.0/openmpch3.html). As such, all behaviour can be strongly controlled using environment variables, for example by:

- [`OMP_SCHEDULE`](https://www.openmp.org/spec-html/5.0/openmpse49.html)
- [`OMP_DEFAULT_DEVICE`](https://www.openmp.org/spec-html/5.0/openmpse63.html#x302-20790006.15)
- [`OMP_THREAD_LIMIT`](https://www.openmp.org/spec-html/5.0/openmpse58.html#x297-20700006.10)



> [!TIP]
> Sometimes the memory bandwidth between different sockets of a machine is poor, and it is substantially better to exchange memory in bulk between their NUMA nodes, rather than through repeated random access. In such settings, it can be worthwhile to hybridise multithreading and distribution, even upon a single machine, partitioning same-socket threads into their own MPI node. This forces inter-socket communication to happen in-batch, via message-passing, at the expense of using _double_ total memory (to store buffers). See the [distributed](#distribution) section.



---------------------


## GPU-acceleration

> [!NOTE]
> Using GPU-acceleration requires first compiling QuEST with `CUDA` or `HIP` enabled (to utilise NVIDIA and AMD GPUs respectively) as detailed in [`compile.md`](compile.md#gpu-acceleration).

### Launching

The compiled executable is launched like any other, via
```bash
./myexec
```

Using _multiple_ available GPUs, regardless of whether they are local or distributed, is done through additionally enabling [distribution](#multi-gpu).


### Monitoring utilisation


To runtime check whether GPU-acceleration was compiled and is being actively utilised, call [`reportQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496).
This will display a subsection like
```
  [compilation]
    isGpuCompiled...........1
  [deployment]
    isGpuEnabled............1
```
where the `1` indicate GPU-acceleration was respectively compiled and is runtime available (i.e. QuEST has found suitable GPUs). When this is the case, another section will be displayed detailing the discovered hardware properties, e.g.
```
  [gpu]
    numGpus...........2
    gpuDirect.........0
    gpuMemPools.......1
    gpuMemory.........15.9 GiB per gpu
    gpuMemoryFree.....15.6 GiB per gpu
    gpuCache..........0 bytes per gpu
```

Utilisation can also be externally monitored using third-party tools:


| GPU  | Type | Name |
| -------- | ------- | ---- |
| NVIDIA | CLI | [`nvidia-smi`](https://docs.nvidia.com/deploy/nvidia-smi/index.html) |
| NVIDIA | GUI | [Nsight](https://developer.nvidia.com/nsight-systems) |
| AMD | CLI, GUI | [`amdgpu_top`](https://github.com/Umio-Yasuno/amdgpu_top) |



Note however that GPU-acceleration might not be leveraged at runtime when either:
- Your `Qureg` created with [`createQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#gab3a231fba4fd34ed95a330c91fcb03b3) was too small to invoke automatic GPU-acceleration.
- You called [`initCustomQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga485268e52f838743357e7a4c8c241e57) and disabled GPU-acceleration for all subsequently-created `Qureg`.

Usage of GPU-acceleration can be (inadvisably) forced using [`createForcedQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga619bbba1cbc2f7f9bbf3d3b86b3f02be) or [`createCustomQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga849971f43e246d103da1731d0901f2e6).


### Configuring

There are a plethora of [environment variables](https://askubuntu.com/questions/58814/how-do-i-add-environment-variables) which be used to control the execution on [NVIDIA](https://docs.nvidia.com/cuda/cuda-c-programming-guide/#env-vars) and [AMD](https://rocm.docs.amd.com/projects/HIP/en/docs-develop/reference/env_variables.html) GPUs. We highlight only some below.

- Choose _which_ GPUs among multiple available to permit QuEST to utilise via [`CUDA_VISIBLE_DEVICES`](https://developer.nvidia.com/blog/cuda-pro-tip-control-gpu-visibility-cuda_visible_devices/) and [`ROCR_VISIBLE_DEVICES`](https://rocm.docs.amd.com/en/latest/conceptual/gpu-isolation.html).
- Alternatively, set the order of selected GPUs (`CUDA_DEVICE_ORDER`) to `FASTEST_FIRST` or `PCI_BUS_ID`.
  - In single-GPU mode, this informs which GPU QuEST will use (i.e. the first).
  - In multi-GPU mode, this informs which local GPUs are used.



### Benchmarking

Beware that the CPU dispatches tasks to the GPU _asynchronously_. Control flow returns immediately to the CPU, which will proceed to other duties (like dispatching the next several quantum operation's worth of instructions to the GPU) while the GPU undergoes independent computation (goes _brrrrr_).
This has no consequence to the user who uses only the QuEST API, which will automatically synchronise the CPU and GPU when necessary (like inside functions [`calcTotalProb()`](https://quest-kit.github.io/QuEST/group__calc__properties.html#gab082910d33473ec29e1d5852943de468)).

However, it _does_ mean codes which seeks to benchmark QuEST must be careful to _wait for the GPU to be ready_ before beginning the stopwatch, and _wait for the GPU to finish_ before stopping the stopwatch. This can be done with the [`syncQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#gaaa19c3112f1ecd80e3296df5c0ed058d), which incidentally also ensures nodes are synchronised when distributed.


---------------------


## Distribution


> [!NOTE]
> Distributing QuEST over multiple machines requires first compiling with 
> distribution enabled, as detailed in [`compile.md`](compile.md#distribution). 

> [!IMPORTANT]
> Simultaneously using distribution _and_ GPU-acceleration introduces additional considerations detailed in the [proceeding section](#multi-gpu).


### Launching

A distributed QuEST executable called `myexec` can be launched and distributed over (e.g.) `32` nodes using [`mpirun`](https://www.open-mpi.org/doc/v4.1/man1/mpirun.1.php) with the 
```bash
mpirun -np 32 ./myexec
```
or on some platforms (such as with Intel and Microsoft MPI):
```bash
mpiexec -n 32 myexec.exe
```

Some supercomputing facilities however may require custom or additional commands, like [SLURM](https://slurm.schedmd.com/documentation.html)'s [`srun`](https://slurm.schedmd.com/srun.html) command. See an excellent guide [here](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/distribution-binding/#distribution).
```bash
srun --nodes=8 --ntasks-per-node=4 --distribution=block:block
```


> [!IMPORTANT]
> QuEST can only be distributed with a _power of `2`_ number of nodes, i.e. `1`, `2`, `4`, `8`, `16`, ...


It is sometimes convenient (mostly for testing) to deploy QuEST across more nodes than there are available machines and sockets, inducing a gratuitous slowdown. Some MPI compilers like [OpenMPI](https://www.open-mpi.org/) forbid this by default, requiring additional commands to permit [oversubscription](https://docs.open-mpi.org/en/main/launching-apps/scheduling.html).
```bash
mpirun -np 1024 --oversubscribe ./mytests
```


### Configuring


> TODO:
> - detail environment variables


### Benchmarking

QuEST strives to reduce inter-node communication when performing distributed simulation, which can otherwise dominate runtime. Between these rare communications, nodes work in complete independence and are likely to desynchronise, especially when performing operations with non-uniform loads. In fact, many-controlled quantum gates are skipped by non-participating nodes which would otherwise wait idly!

Nodes will only synchronise when forced by the user (with [`syncQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#gaaa19c3112f1ecd80e3296df5c0ed058d)), or when awaiting necessary communication (due to functions like [`calcTotalProb()`](https://quest-kit.github.io/QuEST/group__calc__properties.html#gab082910d33473ec29e1d5852943de468)). Furthermore, `Qureg` created with [`createQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#gab3a231fba4fd34ed95a330c91fcb03b3) will automatically disable distribution (and be harmlessly cloned upon every node) when they are too small to outweigh the performance overheads.

This can make monitoring difficult; CPU loads on different nodes can correspond to different stages of execution, and memory loads may fail to distinguish whether a large `Qureg` is distributed or a small `Qureg` is duplicated! Further, a node reaching the end of the program and terminating does not indicate the simulation has finished - other desynchronised nodes may still be working.

It is ergo always prudent to explicitly call [`syncQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#gaaa19c3112f1ecd80e3296df5c0ed058d) immediately before starting and ending a performance timer. This way, the recorded runtime should reflect that of the slowest node (and ergo, the full calculation) rather than that of the node which happened to have its timer output logged.


---------------------


## Multi-GPU


> TODO:
> - explain usecases (multi local GPU, multi remote GPU, hybrid)
> - explain GPUDirect
> - explain CUDA-aware MPI
> - explain UCX
> - detail environment variables
> - detail controlling local vs distributed gpus with device visibility


---------------------

## Supercomputer job schedulers

> TODO:
> - slurm examples

