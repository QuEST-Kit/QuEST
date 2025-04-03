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
build /
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
This binary accepts all of the [Catch2 CLI arguments](https://github.com/catchorg/Catch2/blob/devel/docs/command-line.md), for example to run specific tests
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
which will log each passing test alive, but alas cannot be deployed with distribution.


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

The number of [threads](https://www.openmp.org/spec-html/5.0/openmpsu1.html) to use is decided immediately before launching the compiled executable, using the [`OMP_NUM_THREADS`](https://www.openmp.org/spec-html/5.0/openmpse50.html) environment variable.

```bash
OMP_NUM_THREADS=32 ./myexec
```
```bash
export OMP_NUM_THREADS=32
./myexec
```

It is prudent to choose as many threads as your CPU(s) have total hardware threads or cores. One can view this, and verify the number of available threads at runtime by calling [`reportQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496) which outputs a subsection such as
```
  [cpu]
    numCpuCores.......10 per machine
    numOmpProcs.......10 per machine
    numOmpThrds.......32 per node
```
<!-- the doxygen-doc hyperlink above includes a hash of the function name which should be unchanging! -->

> [!NOTE]
> When running [distributed](#distribution), variable `OMP_NUM_THREADS` specifies the number of threads _per node_ and so should usually be the total number of hardware threads (or cores) _per machine_.


### Monitoring utilisation


The availability of multithreaded deployment can also be checked at runtime using [`reportQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496), which outputs something like:
```bash
  [compilation]
    isOmpCompiled...........1
  [deployment]
    isOmpEnabled............1
```
where `Omp` signifies OpenMP and `1` indicates it is respectively compiled and runtime enabled.

Like all programs, the CPU utilisation of a running QuEST program can be viewed using

| OS    | Program | Method |
| -------- | ------- | ---- |
| Linux  | [`HTOP`](https://htop.dev/) | Run `htop` in terminal  |
| MacOS |  [Activity Monitor](https://support.apple.com/en-gb/guide/activity-monitor/welcome/mac) | Place on dock > right click icon > `Monitors` > `Show CPU usage` [see [here](https://stackoverflow.com/questions/50260592/getting-each-of-the-cpu-cores-usage-via-terminal-in-macos)] |
| Windows  |  [Task Manager](https://learn.microsoft.com/en-us/shows/inside/task-manager) |`Performance` > `CPU` > right click graph > `Change graph to` > `Logical processors` [see [here](https://superuser.com/questions/1398696/how-to-see-usage-of-each-core-in-windows-10)] |


Note however that QuEST will not leverage multithreading at runtime when either:
- Your `Qureg` created with [`createQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#gab3a231fba4fd34ed95a330c91fcb03b3) was too small to invoke automatic multithreading.
- You called [`initCustomQuESTEnv()`](https://quest-kit.github.io/QuEST/group__environment.html#ga485268e52f838743357e7a4c8c241e57) and disabled multithreading for all subsequently-created `Qureg`.

Usage of multithreading can be (inadvisably) forced using [`createForcedQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga619bbba1cbc2f7f9bbf3d3b86b3f02be) or [`createCustomQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga849971f43e246d103da1731d0901f2e6).


### Improving performance

Performance may be improved by setting other [OpenMP variables](https://www.openmp.org/spec-html/5.0/openmpch6.html). Keep in mind that for large `Qureg`s, QuEST's runtime is dominated by the costs of modifying large memory structures during long, uninterrupted loops: namely the updating of statevector amplitudes. For example

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




---------------------


## GPU-acceleration

> [!NOTE]
> Using GPU-acceleration requires first compiling QuEST with `CUDA` or `HIP` enabled (to utilise NVIDIA and AMD GPUs respectively), as detailed in [`compile.md`](compile.md#gpu-acceleration).

The compiled executable is launched like any other, via
```bash
./myexec
```

Using _multiple_ available GPUs, regardless of whether they are local or distributed, is done through additionally enabling [distribution](#distributed-gpu-acceleration).


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




---------------------


## Distribution



> TODO:
> - launching
> - power-of-2 nodes
> - memory move implications
> - distributing over sockets trick
> - oversubscribing
> - env vars



---------------------


## Distributed GPU-acceleration


> TODO:
> - env vars etc
> - ucx
> - controlling local vs distributed gpus with device visibility etc


---------------------

## Supercomputer job schedulers

> TODO:
> - slurm examples




