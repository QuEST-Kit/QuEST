<!--
  Examples and tutorials
  
  @author Tyson Jones
-->

# 🔖  Examples

The above folders contain example `C` and `C++` files which use QuEST's [API](https://quest-kit.github.io/QuEST/group__api.html), helping illustrate how to use specific functions. Instructions for compiling and running them are given in [`compile.md`](/docs/compile.md#tests) and [`run.md`](/docs/run.md#tests) respectively.


# 🎓  Tutorial

QuEST is included into a `C` or `C++` project via
```C++
#include "quest.h"
```

> [!TIP]
> Some of QuEST's deprecated `v3` API can be accessed by specifying `ENABLE_DEPRECATED_API` when [compiling](/docs/compile.md), or defining it before import, i.e. 
> ```C++
> #define ENABLE_DEPRECATED_API 1
> #include "quest.h"
> ```
> We recommend migrating to the latest `v4` API however, demonstrated below.

Simulation typically proceeds as:
1. [Initialise](https://quest-kit.github.io/QuEST/group__environment.html#gab89cfc1bf94265f4503d504b02cf54d4) the QuEST [environment](https://quest-kit.github.io/QuEST/group__environment.html), preparing available GPUs and networks.
2. [Configure](https://quest-kit.github.io/QuEST/group__debug.html) the environment, such as through [seeding](https://quest-kit.github.io/QuEST/group__debug__seed.html).
3. [Create](https://quest-kit.github.io/QuEST/group__qureg__create.html) a [`Qureg`](https://quest-kit.github.io/QuEST/structQureg.html), allocating memory for its amplitudes.
4. Prepare its [initial state](https://quest-kit.github.io/QuEST/group__initialisations.html), overwriting its amplitudes.
5. Apply [operators](https://quest-kit.github.io/QuEST/group__operations.html) and [decoherence](https://quest-kit.github.io/QuEST/group__decoherence.html), expressed as [matrices](https://quest-kit.github.io/QuEST/group__matrices.html) and [channels](https://quest-kit.github.io/QuEST/group__channels.html).
6. Perform [calculations](https://quest-kit.github.io/QuEST/group__calculations.html), potentially using [Pauli](https://quest-kit.github.io/QuEST/group__paulis.html) observables.
7. [Report](https://quest-kit.github.io/QuEST/group__types.html) or log the results to file.
8. Destroy any heap-allocated [`Qureg`](https://quest-kit.github.io/QuEST/group__qureg__destroy.html) or [matrices](https://quest-kit.github.io/QuEST/group__matrices__destroy.html).
8. [Finalise](https://quest-kit.github.io/QuEST/group__environment.html#ga428faad4d68abab20f662273fff27e39) the QuEST environment.

Of course, the procedure is limited only by the programmers imagination `¯\_(ツ)_/¯` Let's see an example of these steps below.


## 1. Initialise the environment

Before calling any other QuEST functions, we must [_initialise_](https://quest-kit.github.io/QuEST/group__environment.html#gab89cfc1bf94265f4503d504b02cf54d4) the QuEST [_environment_](https://quest-kit.github.io/QuEST/group__environment.html).
```C++
initQuESTEnv();
```
This does several things, such as
- assessing which hardware accelerations (multithreading, GPU-acceleration, distribution, cuQuantum) were compiled and are currently available to use.
- initialising any external libraries as needed, like MPI, CUDA and cuQuantum.
- seeding the random number generators (informing measurements and random states), using a [CSPRNG](https://en.wikipedia.org/wiki/Cryptographically_secure_pseudorandom_number_generator) if available.

We could instead forcefully [disable](https://quest-kit.github.io/QuEST/group__environment.html#ga485268e52f838743357e7a4c8c241e57) certain hardware accelerations
```C++
int useMPI = 0;
int useGPU = 0;
int useOMP = 0;
initCustomQuESTEnv(useMPI, useGPU, useOMP);
```

> [!TIP]
> We recommend enabling _all_ deployments, as automated by `initQuESTEnv()`, which
> permits QuEST to choose how to best accelerate subsequently created `Qureg`.

We can [view](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496) the environment configuration at runtime, via
```C++
reportQuESTEnv();
```
which might output something like
```
QuEST execution environment:
  [precision]
    qreal.................double (8 bytes)
    qcomp.................std::__1::complex<double> (16 bytes)
    qindex................long long int (8 bytes)
    validationEpsilon.....1e-12
  [compilation]
    isMpiCompiled...........1
    isGpuCompiled...........1
    isOmpCompiled...........1
    isCuQuantumCompiled.....0
  [deployment]
    isMpiEnabled.....0
    isGpuEnabled.....1
    isOmpEnabled.....1
  [cpu]
    numCpuCores.......10 per machine
    numOmpProcs.......10 per machine
    numOmpThrds.......8 per node
    cpuMemory.........32 GiB per node
    cpuMemoryFree.....7.1 GiB per node
  [gpu]
    numGpus...........1
    gpuDirect.........1
    gpuMemPools.......1
    gpuMemory.........15.9 GiB per node
    gpuMemoryFree.....15.2 GiB per node
    gpuCache..........1 GiB
  [distribution]
    isMpiGpuAware.....0
    numMpiNodes.......8
  [statevector limits]
    minQubitsForMpi.............3
    maxQubitsForCpu.............30
    maxQubitsForGpu.............29
    maxQubitsForMpiCpu..........35
    maxQubitsForMpiGpu..........34
    maxQubitsForMemOverflow.....59
    maxQubitsForIndOverflow.....63
  [density matrix limits]
    minQubitsForMpi.............2
    maxQubitsForCpu.............15
    maxQubitsForGpu.............14
    maxQubitsForMpiCpu..........17
    maxQubitsForMpiGpu..........16
    maxQubitsForMemOverflow.....29
    maxQubitsForIndOverflow.....31
  [statevector autodeployment]
    8 qubits.....[omp]
    12 qubits....[gpu]
    29 qubits....[gpu] [mpi]
  [density matrix autodeployment]
    4 qubits.....[omp]
    6 qubits.....[gpu]
    15 qubits....[gpu] [mpi]
```

We can also [obtain](https://quest-kit.github.io/QuEST/group__environment.html#ga6b9e84b462a999a1fbb9a372f990c491) some of the environment information [programmatically](https://quest-kit.github.io/QuEST/structQuESTEnv.html)
```C++
QuESTEnv env = getQuESTEnv();

if (env.isGpuAccelerated)
   printf("vroom vroom");
```


## 2. Configure the environment

Configuring the environment is ordinarily not necessary, but convenient in certain applications.

For example, we may wish our simulations to deterministically obtain the same measurement outcomes and random states as a previous or future run, and ergo choose to [override](https://quest-kit.github.io/QuEST/group__debug__seed.html#ga9e3a6de413901afbf50690573add1587) the default seeds.
```C++
unsigned seeds[] = {123u, 1u << 10};
setSeeds(seeds, 2);
```

We may wish further to [adjust](https://quest-kit.github.io/QuEST/group__debug__reporting.html) how subsequent functions will display information to the screen
```C++
int maxRows = 8;
int maxCols = 4;
setMaxNumReportedItems(maxRows, maxCols);
setMaxNumReportedSigFigs(3);
```
or [add](https://quest-kit.github.io/QuEST/group__debug__reporting.html#ga29413703d609254244d6b13c663e6e06) extra spacing between QuEST's printed outputs
```C++
setNumReportedNewlines(3);
```

Perhaps we also wish to relax the [precision](https://quest-kit.github.io/QuEST/group__debug__validation.html#gae395568df6def76045ec1881fcb4e6d1) with which our future inputs will be asserted unitary or Hermitian
```C++
setValidationEpsilon(0.001);
```
but when unitarity _is_ violated, or we otherwise pass an invalid input, we wish to execute a [custom function](https://quest-kit.github.io/QuEST/group__debug__validation.html#ga14b6e7ce08465e36750da3acbc41062f) before exiting.
```C++
#include <stdlib.h>

void myErrorHandler(const char *func, const char *msg) {
    printf("QuEST function '%s' encountered error '%s'\n", func, msg);
    printf("Exiting...\n");
    exit(1);
}

setInputErrorHandler(myErrorHandler);
```
`C++` users may prefer to throw an exception which can be caught, safely permitting execution to continue. In such cases, the erroneous function will _never_ corrupt any passed inputs like `Qureg` nor matrices, nor cause leaks.
```C++
#include <stdexcept>
#include <string>

void myErrorHandlerA(const char* errFunc, const char* errMsg) {
    std::string func(errFunc);
    std::string msg(errMsg);
    throw std::runtime_error(func + ": " + msg);
}

setInputErrorHandler(myErrorHandler);
```

## 3. Create a `Qureg`

To [create](https://quest-kit.github.io/QuEST/group__qureg__create.html) a statevector of `10` qubits, we call
```C++
Qureg qureg = createQureg(10);
```
which we can [verify](https://quest-kit.github.io/QuEST/group__qureg__report.html#ga2a9df2538e537332b1aef8596ce337b2) has begun in the very boring zero state.
```C++
reportQureg(qureg);
```
```
Qureg (10 qubit statevector, 1024 qcomps, 16.1 KiB):
    1  |0⟩
    0  |1⟩
    0  |2⟩
    0  |3⟩
    ⋮
    0  |1020⟩
    0  |1021⟩
    0  |1022⟩
    0  |1023⟩
```
> This printed only `8` amplitudes as per our setting of [`setMaxNumReportedItems()`](https://quest-kit.github.io/QuEST/group__debug__reporting.html#ga093c985b1970a0fd8616c01b9825979a) above.

Behind the scenes, the function `createQureg` did something clever; it consulted the compiled deployments and available hardware to decide whether to distribute `qureg`, or dedicate it persistent GPU memory, and marked whether or not to multithread its subsequent modification. It attempts to choose _optimally_, avoiding gratuitous parallelisation if the overheads outweigh the benefits, or if the hardware devices have insufficient memory.

We call this **_auto-deployment_**, and the chosen configuration can be [previewed](https://quest-kit.github.io/QuEST/group__qureg__report.html#ga97d96af7c7ea7b31e32cbe3b25377e09) via
```C++
reportQuregParams(qureg);
```
```
Qureg:
  [deployment]
    isMpiEnabled.....0
    isGpuEnabled.....0
    isOmpEnabled.....1
  [dimension]
    isDensMatr.....0
    numQubits......10
    numCols........N/A
    numAmps........2^10 = 1024
  [distribution]
    numNodes.....N/A
    numCols......N/A
    numAmps......N/A
  [memory]
    cpuAmps...........16 KiB
    gpuAmps...........N/A
    cpuCommBuffer.....N/A
    gpuCommBuffer.....N/A
    globalTotal.......16 KiB
```
The above output informs us that the `qureg` has not been distributed nor GPU-accelerated, but _will_ be multithreaded.

If we so wished, we could [_force_](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga619bbba1cbc2f7f9bbf3d3b86b3f02be) the use of all deployments available to the environment
```C++
Qureg qureg = createForcedQureg(10);
reportQuregParams(qureg);
```
```
Qureg:
  [deployment]
    isMpiEnabled.....1
    isGpuEnabled.....1
    isOmpEnabled.....1
  [dimension]
    isDensMatr.....0
    numQubits......10
    numCols........N/A
    numAmps........2^10 = 1024
  [distribution]
    numNodes.....2^3 = 8
    numCols......N/A
    numAmps......2^7 = 128 per node
  [memory]
    cpuAmps...........2 KiB per node
    gpuAmps...........2 KiB per node
    cpuCommBuffer.....2 KiB per node
    gpuCommBuffer.....2 KiB per node
    globalTotal.......64 KiB
```
or [select](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga849971f43e246d103da1731d0901f2e6) specific deployments
```C++
int useMPI = 1;
int useGPU = 0;
int useOMP = 0;
Qureg qureg = createCustomQureg(10, 0, useMPI, useGPU, useOMP);
```

In lieu of a statevector, we could create a [density matrix](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga1470424b0836ae18b5baab210aedf5d9)
```C++
Qureg qureg = createDensityQureg(10);
```
which is also auto-deployed. Note this contains _square_ as many amplitudes as the equal-dimensin statevector, and ergo requires _square_ as much memory.
```C++
reportQureg(qureg);
reportQuregParams(qureg);
```
```
Qureg (10 qubit density matrix, 1024x1024 qcomps, 16 MiB):
    1  0  …  0  0
    0  0  …  0  0
    0  0  …  0  0
    0  0  …  0  0
    ⋮
    0  0  …  0  0
    0  0  …  0  0
    0  0  …  0  0
    0  0  …  0  0


Qureg:
  ...
  [dimension]
    isDensMatr.....1
    numQubits......10
    numCols........2^10 = 1024
    numAmps........2^20 = 1048576
  ...
  [memory]
    cpuAmps...........16 MiB
    ...
    globalTotal.......16 MiB
```

> The spacing between the outputs of those two consecutive QuEST functions was determined by our earlier call to [`setMaxNumReportedSigFigs()`](https://quest-kit.github.io/QuEST/group__debug__reporting.html#ga29413703d609254244d6b13c663e6e06).


A density matrix `Qureg` can model classical uncertainty as results from [decoherence](https://quest-kit.github.io/QuEST/group__decoherence.html), and proves useful when simulating quantum operations on a noisy quantum computer.


## 4. Prepare an initial state

In lieu of manually [modifying](https://quest-kit.github.io/QuEST/group__init__amps.html) the state amplitudes, QuEST includes functions to prepare a `Qureg` in some common [initial states](https://quest-kit.github.io/QuEST/group__init__states.html)

```C++
initZeroState(qureg);         // |0> or |0><0|
initPlusState(qureg);         // |+> or |+><+|
initClassicalState(qureg, i); // |i> or |i><i|
initPureState(rho, psi);      // rho = |psi><psi|
```
or random states
```C++
initRandomPureState(psi);

int numPureStates = 15;
initRandomMixedState(rho, numPureStates);

reportQureg(psi);
reportQureg(rho);
```
```
Qureg (5 qubit statevector, 32 qcomps, 616 bytes):
    0.0884-0.164i     |0⟩
    0.149+0.207i      |1⟩
    0.232+0.0656i     |2⟩
    -0.0435+0.0332i   |3⟩
            ⋮
    -0.108-0.0431i    |28⟩
    -0.0161-0.121i    |29⟩
    -0.0463+0.00341i  |30⟩
    -0.0491-0.186i    |31⟩


Qureg (5 qubit density matrix, 32x32 qcomps, 16.1 KiB):
    0.0256+(1.08e-19)i  -0.000876+0.00412i  …  0.000912+0.00869i   -0.00597+0.00615i
    -0.000876-0.00412i  0.033-(6.78e-20)i   …  0.000223+0.00369i   -0.00207+0.00451i
    -0.00443-0.00871i   0.0155-0.000843i    …  0.00375+0.00669i    (8.5e-5)-0.000851i
    0.00287-0.00397i    0.00637-0.000315i   …  0.00486+0.00218i    0.00268+0.0053i
             ⋮
    -0.00385-0.000732i  0.00965+0.00542i    …  0.00162-0.0112i     0.00404+0.00685i
    0.00491+0.00245i    -0.000319+0.0021i   …  -0.00902-0.00312i   -0.00465+0.00275i
    0.000912-0.00869i   0.000223-0.00369i   …  0.0183+(1.32e-19)i  0.000509+0.00401i
    -0.00597-0.00615i   -0.00207-0.00451i   …  0.000509-0.00401i   0.0173+(3.12e-19)i
```

> The number of printed significant figures above results from our earlier calling of [`setMaxNumReportedSigFigs()`](https://quest-kit.github.io/QuEST/group__debug__reporting.html#ga15d46e5d813f70b587762814964e1994).


## 5. Apply operators

> TODO