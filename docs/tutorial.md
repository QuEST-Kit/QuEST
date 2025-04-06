# 🎓  Tutorial

<!--
  Tutorial
  (this comment must be under the title for valid doxygen rendering)
  
  @author Tyson Jones
-->

QuEST is included into a `C` or `C++` project via
```cpp
#include "quest.h"
```

<!-- @todo the below link fails in Doxygen; it's too stupid to recognise the section ref -->
> [!TIP]
> Some of QuEST's deprecated `v3` API can be accessed by specifying `ENABLE_DEPRECATED_API` when [compiling](/docs/compile.md#v3), or defining it before import, i.e. 
> ```cpp
> #define ENABLE_DEPRECATED_API 1
> #include "quest.h"
> ```
> We recommend migrating to the latest `v4` API however as will be showcased below.

Simulation typically proceeds as:
1. [Initialise](https://quest-kit.github.io/QuEST/group__environment.html#gab89cfc1bf94265f4503d504b02cf54d4) the QuEST [environment](https://quest-kit.github.io/QuEST/group__environment.html), preparing available GPUs and networks.
2. [Configure](https://quest-kit.github.io/QuEST/group__debug.html) the environment, such as through [seeding](https://quest-kit.github.io/QuEST/group__debug__seed.html).
3. [Create](https://quest-kit.github.io/QuEST/group__qureg__create.html) a [`Qureg`](https://quest-kit.github.io/QuEST/structQureg.html), allocating memory for its amplitudes.
4. Prepare its [initial state](https://quest-kit.github.io/QuEST/group__initialisations.html), overwriting its amplitudes.
5. Apply [operators](https://quest-kit.github.io/QuEST/group__operations.html) and [decoherence](https://quest-kit.github.io/QuEST/group__decoherence.html), expressed as [matrices](https://quest-kit.github.io/QuEST/group__matrices.html) and [channels](https://quest-kit.github.io/QuEST/group__channels.html).
6. Perform [calculations](https://quest-kit.github.io/QuEST/group__calculations.html), potentially using [Pauli](https://quest-kit.github.io/QuEST/group__paulis.html) observables.
7. [Report](https://quest-kit.github.io/QuEST/group__types.html) or log the results to file.
8. Destroy any heap-allocated [`Qureg`](https://quest-kit.github.io/QuEST/group__qureg__destroy.html) or [matrices](https://quest-kit.github.io/QuEST/group__matrices__destroy.html).
9. [Finalise](https://quest-kit.github.io/QuEST/group__environment.html#ga428faad4d68abab20f662273fff27e39) the QuEST environment.

Of course, the procedure is limited only by the programmers imagination `¯\_(ツ)_/¯` Let's see an example of these steps below.


<!-- 
    we are using explicit <a>, rather than markdown links,
    for Doxygen compatibility. It cannot handle [](#sec)
    links, and its <a> anchors are not scoped to files, so
    we here prefix each name with the filename. Grr!
-->

> **TOC**:
> - <a href="#tutorial_initialise-the-environment">Initialise the environment</a>
> - <a href="#tutorial_configure-the-environment">Configure the environment</a>
> - <a href="#tutorial_create-a-qureg">Create a `Qureg`</a>
> - <a href="#tutorial_prepare-an-initial-state">Prepare an initial state</a>
> - <a href="#tutorial_apply-operators">Apply operators</a>
>   * <a href="#tutorial_controls">controls</a>
>   * <a href="#tutorial_paulis">paulis</a>
>   * <a href="#tutorial_matrices">matrices</a>
>   * <a href="#tutorial_circuits">circuits</a>
>   * <a href="#tutorial_measurements">measurements</a>
>   * <a href="#tutorial_decoherence">decoherence</a>
> - <a href="#tutorial_perform-calculations">Perform calculations</a>
> - <a href="#tutorial_report-the-results">Report the results</a>
> - <a href="#tutorial_cleanup">Cleanup</a>
> - <a href="#tutorial_finalise-quest">Finalise QuEST</a>



--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_initialise-the-environment"></a>

## Initialise the environment


Before calling any other QuEST functions, we must [_initialise_](https://quest-kit.github.io/QuEST/group__environment.html#gab89cfc1bf94265f4503d504b02cf54d4) the QuEST [_environment_](https://quest-kit.github.io/QuEST/group__environment.html).
```cpp
initQuESTEnv();
```
This does several things, such as
- assessing which hardware accelerations (multithreading, GPU-acceleration, distribution, cuQuantum) were compiled and are currently available to use.
- initialising any external libraries as needed, like MPI, CUDA and cuQuantum.
- seeding the random number generators (informing measurements and random states), using a [CSPRNG](https://en.wikipedia.org/wiki/Cryptographically_secure_pseudorandom_number_generator) if available.

We could instead forcefully [disable](https://quest-kit.github.io/QuEST/group__environment.html#ga485268e52f838743357e7a4c8c241e57) certain hardware accelerations
```cpp
int useMPI = 0;
int useGPU = 0;
int useOMP = 0;
initCustomQuESTEnv(useMPI, useGPU, useOMP);
```

> [!TIP]
> We recommend enabling _all_ deployments, as automated by `initQuESTEnv()`, which
> permits QuEST to choose how to best accelerate subsequently created `Qureg`.

We can [view](https://quest-kit.github.io/QuEST/group__environment.html#ga08bf98478c4bf21b0759fa7cd4a97496) the environment configuration at runtime, via
```cpp
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
```cpp
QuESTEnv env = getQuESTEnv();

if (env.isGpuAccelerated)
    printf("vroom vroom");
```



--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_configure-the-environment"></a>

## Configure the environment


Configuring the environment is ordinarily not necessary, but convenient in certain applications.

For example, we may wish our simulations to deterministically obtain the same measurement outcomes and random states as a previous or future run, and ergo choose to [override](https://quest-kit.github.io/QuEST/group__debug__seed.html#ga9e3a6de413901afbf50690573add1587) the default seeds.
```cpp
unsigned seeds[] = {123u, 1u << 10};
setSeeds(seeds, 2);
```

We may wish further to [adjust](https://quest-kit.github.io/QuEST/group__debug__reporting.html) how subsequent functions will display information to the screen
```cpp
int maxRows = 8;
int maxCols = 4;
setMaxNumReportedItems(maxRows, maxCols);
setMaxNumReportedSigFigs(3);
```
or [add](https://quest-kit.github.io/QuEST/group__debug__reporting.html#ga29413703d609254244d6b13c663e6e06) extra spacing between QuEST's printed outputs
```cpp
setNumReportedNewlines(3);
```

Perhaps we also wish to relax the [precision](https://quest-kit.github.io/QuEST/group__debug__validation.html#gae395568df6def76045ec1881fcb4e6d1) with which our future inputs will be asserted unitary or Hermitian
```cpp
setValidationEpsilon(0.001);
```
but when unitarity _is_ violated, or we otherwise pass an invalid input, we wish to execute a [custom function](https://quest-kit.github.io/QuEST/group__debug__validation.html#ga14b6e7ce08465e36750da3acbc41062f) before exiting.
```cpp
#include <stdlib.h>

void myErrorHandler(const char *func, const char *msg) {
    printf("QuEST function '%s' encountered error '%s'\n", func, msg);
    printf("Exiting...\n");
    exit(1);
}

setInputErrorHandler(myErrorHandler);
```

> [!TIP]
> `C++` users may prefer to throw an exception which can be caught, safely permitting execution to continue. In such cases, the erroneous function will _never_ corrupt any passed inputs like `Qureg` nor matrices, nor cause leaks.
> ```cpp
> #include <stdexcept>
> #include <string>
> void myErrorHandlerA(const char* errFunc, const char* errMsg) {
>     std::string func(errFunc);
>     std::string msg(errMsg);
>     throw std::runtime_error(func + ": " + msg);
> }
> setInputErrorHandler(myErrorHandler);
> ```
<!-- newlines removed above because doxygen renders them as <br> text, how stupid! -->




--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_create-a-qureg"></a>

## Create a Qureg


To [create](https://quest-kit.github.io/QuEST/group__qureg__create.html) a statevector of `10` qubits, we call
```cpp
Qureg qureg = createQureg(10);
```
which we can [verify](https://quest-kit.github.io/QuEST/group__qureg__report.html#ga2a9df2538e537332b1aef8596ce337b2) has begun in the very boring zero state.
```cpp
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

We call this _auto-deployment_, and the chosen configuration can be [previewed](https://quest-kit.github.io/QuEST/group__qureg__report.html#ga97d96af7c7ea7b31e32cbe3b25377e09) via
```cpp
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
```cpp
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
```cpp
int useMPI = 1;
int useGPU = 0;
int useOMP = 0;
Qureg qureg = createCustomQureg(10, 0, useMPI, useGPU, useOMP);
```

In lieu of a statevector, we could create a [density matrix](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga1470424b0836ae18b5baab210aedf5d9)
```cpp
Qureg qureg = createDensityQureg(10);
```
which is also auto-deployed. Note this contains _square_ as many amplitudes as the equal-dimension statevector and ergo requires _square_ as much memory.
```cpp
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




--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_prepare-an-initial-state"></a>

## Prepare an initial state


In lieu of manually [modifying](https://quest-kit.github.io/QuEST/group__init__amps.html) the state amplitudes, QuEST includes functions to prepare a `Qureg` in some common [initial states](https://quest-kit.github.io/QuEST/group__init__states.html)

```cpp
initZeroState(qureg);         // |0> or |0><0|
initPlusState(qureg);         // |+> or |+><+|
initClassicalState(qureg, i); // |i> or |i><i|
initPureState(rho, psi);      // rho = |psi><psi|
```
or random states
```cpp
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



--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_apply-operators"></a>

## Apply operators


QuEST supports an extensive set of [operators](https://quest-kit.github.io/QuEST/group__operations.html) to effect upon a `Qureg`. 
```cpp
int target = 2;
applyHadamard(qureg, target);

qreal angle = 3.14 / 5;
int targets[]  = {4,5,6};
applyPhaseGadget(qureg, targets, 3, angle);
```

> [!IMPORTANT]  
> Notice the type of `angle` is [`qreal`](https://quest-kit.github.io/QuEST/group__types.html#ga2d479c159621c76ca6f96abe66f2e69e) rather than the expected `double`. This is a precision agnostic alias for a floating-point, real scalar which allows you to recompile QuEST with a varying [precision](/docs/compile.md#precision) with no modifications to your code. 
<!-- @todo the above link fails in Doxygen; it's too stupid to recognise the section ref -->


<!-- permit doxygen to reference section -->
<a id="tutorial_controls"></a>

### controls


All unitary operations accept any number of control qubits
```cpp
int controls[] = {0,1,2,3,7,8,9};
applyMultiControlledSqrtSwap(qureg, controls, 7, targets[0], targets[1]);
```
and even _control states_ which specify the bits (`0` or `1`) that the respective controls must be in to effect the non-identity operation.
```cpp
int states[] = {0,0,0,1,1,1,0};
applyMultiStateControlledRotateX(qureg, controls, states, 7, target, angle);
```

> [!TIP]
> `C` users can pass inline list arguments using [compound literals](https://en.cppreference.com/w/c/language/compound_literal)
> ```C
> applyMultiControlledMultiQubitNot(qureg, (int[]) {0,1,2}, 3, (int[]) {4,5}, 2);
> ```
> while `C++` users can pass [vector](https://en.cppreference.com/w/cpp/container/vector) literals or [initializer lists](https://en.cppreference.com/w/cpp/utility/initializer_list), alleviating the need to specify the list lengths.
> ```cpp
> applyMultiControlledMultiQubitNot(qureg, {0,1,2}, {4,5});
> ```


<!-- permit doxygen to reference section -->
<a id="tutorial_paulis"></a>

### paulis


Some operators accept [`PauliStr`](https://quest-kit.github.io/QuEST/structPauliStr.html) which can be [constructed](https://quest-kit.github.io/QuEST/group__paulis__create.html) all sorts of ways - even inline!
```cpp
applyPauliGadget(qureg, getPauliStr("XYZ"), angle);
```

> [!TIP]
> Using _one_ QuEST function is _always_ faster than using an equivalent sequence. So
> ```cpp
> applyPauliStr(qureg, getPauliStr("YYYYYYY"));
> ```
> is _much_ faster than
> ```cpp
> for (int i=0; i<7; i++)
>     applyPauliY(qureg, i);
> ```


<!-- permit doxygen to reference section -->
<a id="tutorial_matrices"></a>

### matrices


<!-- giving no hyperlink -->

#### `CompMatr1`

Don't see your operation in the API? You can specify it as a general [matrix](https://quest-kit.github.io/QuEST/group__matrices.html).
```cpp
qcomp x = 1i/sqrt(2);
CompMatr1 matrix = getInlineCompMatr({{-x,x},{-x,-x}});
applyCompMatr1(qureg, target, matrix);
```

> [!IMPORTANT]  
> The type [`qcomp`](https://quest-kit.github.io/QuEST/group__types.html#ga4971f489e74bb185b9b2672c14301983) above is a precision agnostic complex scalar, and has beautiful arithmetic overloads!
> ```cpp
> qcomp x = 1.5 + 3.14i;
> qcomp *= 1E3i - 1E-5i;
> ```
> Beware that in `C++`, `1i` is a _double precision_ literal, so `C++` users should instead
> use the custom precision-agnostic literal `1_i`.
> ```cpp
> qcomp x = 1.5 + 3.14_i;
> ```


<!-- giving no hyperlink -->

#### `CompMatr`

Want a bigger matrix? No problem - they can be [any size](https://quest-kit.github.io/QuEST/group__matrices__create.html#ga634309472d1edf400174680af0685b89), with many ways to [initialise](https://quest-kit.github.io/QuEST/group__matrices__setters.html) them.
```cpp
CompMatr bigmatrix = createCompMatr(8);
setCompMatr(bigmatrix, {{1,2,3,...}});
applyCompMatr(qureg, ..., bigmatrix);
```
Matrix elements can be manually modified, though this requires we [synchronise](https://quest-kit.github.io/QuEST/group__matrices__sync.html) them with GPU memory once finished.
```cpp
qindex dim = bigmatrix.numRows;

// initialise random diagonal unitary
for (qindex r=0; r<dim; r++)
    for (qindex c=0; c<dim; c++)
        bigmatrix.cpuElems[r][c] = exp(rand() * 1i) * (r==c);

// update the GPU copy 
syncCompMatr(bigmatrix);
```

> [!IMPORTANT]  
> The created `CompMatr` is a [heap object](https://craftofcoding.wordpress.com/2015/12/07/memory-in-c-the-stack-the-heap-and-static/) and must be [destroyed](https://quest-kit.github.io/QuEST/group__matrices__destroy.html) when we are finished with it, to free up its memory and avoid leaks.
> ```cpp
> destroyCompMatr(bigmatrix);
> ```
> This is true of any QuEST structure returned by a `create*()` function. It is _not_ true of functions prefixed with `get*()` with are always [stack variables](https://craftofcoding.wordpress.com/2015/12/07/memory-in-c-the-stack-the-heap-and-static/), hence why functions like `getCompMatr1()` can be called inline!


<!-- giving no hyperlink -->

#### `FullStateDiagMatr`

Above, we initialised [`CompMatr`](https://quest-kit.github.io/QuEST/structCompMatr.html) to a diagonal unitary. This is incredibly wasteful; only `256` of its `65536` elements are non-zero! We should instead use [`DiagMatr`](https://quest-kit.github.io/QuEST/structDiagMatr.html) or [`FullStateDiagMatr`](https://quest-kit.github.io/QuEST/structFullStateDiagMatr.html). The latter is even distributed (if chosen by the autodeployer), permitting it to be as large as a `Qureg` itself!
```cpp
FullStateDiagMatr fullmatrix = createFullStateDiagMatr(qureg.numQubits);
```
and can be [initialised](https://quest-kit.github.io/QuEST/group__matrices__setters.html) in many ways, including from all-`Z` pauli sums!
```cpp
PauliStrSum sum = createInlinePauliStrSum(R"(
    1   II
    1i  ZI
    1i  IZ
    -1  ZZ
)");

setFullStateDiagMatrFromPauliStrSum(fullmatrix, sum);
```
> [!IMPORTANT]  
> The argument to `createInlinePauliStrSum` is a multiline string for which the syntax differs between `C` and `C++`; we used the latter above. See examples [`initialisation.c`](/examples/paulis/initialisation.c) and [`initialisation.cpp`](/paulis/matrices/initialisation.cpp) for clarity.

> [!CAUTION]
> Beware that in distributed settings, because `fullmatrix` _may_ be distributed, we should must exercise extreme caution when modifying its `fullmatrix.cpuElems` directly. 


A `FullStateDiagMatr` acts upon all qubits of a qureg
```cpp
applyFullStateDiagMatr(qureg, fullmatrix);
```
and can be raised to an arbitrary power, helpful for example in simulating [quantum spectral methods](https://www.science.org/doi/10.1126/sciadv.abo7484).
```cpp
qcomp exponent = 3.5;
applyFullStateDiagMatrPower(qureg, fullmatrix, exponent);
```

Notice the `exponent` is a `qcomp` and ergo permitted to be a complex number. Unitarity requires `exponent` is strictly real, but we can always relax the unitarity validation...


<!-- giving no hyperlink -->

#### validation


Our example above initialised `CompMatr` to a diagonal because it is tricky to generate random non-diagonal _unitary_ matrices - and QuEST checks for unitarity!
```cpp
// m * dagger(m) != identity
CompMatr1 m = getCompMatr1({{.1,.2},{.3,.4}});
applyCompMatr1(qureg, 0, m);
```
```
QuEST encountered a validation error during function 'applyCompMatr1':
The given matrix was not (approximately) unitary.
Exiting...
```
If we're satisfied our matrix _is_ sufficiently approximately unitary, we can [adjust](https://quest-kit.github.io/QuEST/group__debug__validation.html#gae395568df6def76045ec1881fcb4e6d1) or [disable](https://quest-kit.github.io/QuEST/group__debug__validation.html#ga5999824df0785ea88fb2d5b5582f2b46) the validation.
```cpp
// max(norm(m * dagger(m) - identity)) = 0.9025
setValidationEpsilon(0.903);
applyCompMatr1(qureg, 0, m);
```


<!-- permit doxygen to reference section -->
<a id="tutorial_circuits"></a>

### circuits


QuEST includes a few convenience functions for effecting [QFT](https://quest-kit.github.io/QuEST/group__op__qft.html) and [Trotter](https://quest-kit.github.io/QuEST/group__op__paulistrsum.html) circuits.

```cpp
applyQuantumFourierTransform(qureg, targets, 3);

qreal time = .3;
int order = 4;
int reps = 10;
applyTrotterizedPauliStrSumGadget(qureg, sum, time, order, reps);
```


<!-- permit doxygen to reference section -->
<a id="tutorial_measurements"></a>

### measurements


We can also effect a wide range of non-unitary operations, such as destructive [measurements](https://quest-kit.github.io/QuEST/group__op__measurement.html)
```cpp
int outcome1 = applyQubitMeasurement(qureg, 0);

qreal prob;
qindex outcome2 = applyMultiQubitMeasurementAndGetProb(qureg, targets, 3, &prob);
```
and conveniently [report](https://quest-kit.github.io/QuEST/group__types.html#ga2be8a4433585a8d737c02128b4754a03) their outcome.
```cpp
reportScalar("one qubit outcome", outcome1);
reportScalar("three qubit outcome", outcome2);
```

> [!IMPORTANT]  
> Notice the type of `outcome2` is a [`qindex`](https://quest-kit.github.io/QuEST/group__types.html#ga6017090d3ed4063ee7233e20c213424b) rather than an `int`. This is a larger type which can store much larger numbers without overflow - up to `2^63` - and is always used by the API for many-qubit indices.

Should we wish to leave the state unnormalised, we can instead use [projectors](https://quest-kit.github.io/QuEST/group__op__projectors.html).



<!-- permit doxygen to reference section -->
<a id="tutorial_decoherence"></a>

### decoherence


Density matrices created with [`createDensityQureg()`](https://quest-kit.github.io/QuEST/group__qureg__create.html#ga1470424b0836ae18b5baab210aedf5d9) can undergo [decoherence](https://quest-kit.github.io/QuEST/group__decoherence.html) channels.

```cpp
qreal prob = 0.1;
mixDamping(rho, target, prob);
mixDephasing(rho, target, prob);
mixTwoQubitDepolarising(rho, targets[0], targets[1], prob);
```
which we can specify as inhomogeneous Pauli channels
```cpp
// passing probabilities of X, Y, Z errors respectively
mixPaulis(Qureg qureg, target, .05, .10, .15);
```
or completely generally as [Kraus maps](https://quest-kit.github.io/QuEST/group__channels.html) and [superoperators](https://quest-kit.github.io/QuEST/group__channels.html)!
```cpp
int numTargets = 1;
int numOperators = 4;

qreal p = 0.1;
qreal l = 0.3;

// generalised amplitude damping
KrausMap map = createInlineKrausMap(numTargets, numOperators, {
    {
        {sqrt(p), 0},
        {0, sqrt(p*(1-l))}
    }, {
        {0, sqrt(p*l)}, 
        {0, 0}
    }, {
        {sqrt((1-p)*(1-l)), 0},
        {0, sqrt(1-p)}
    }, {
        {0, 0},
        {sqrt((1-p)*l), 0}
    }
});

int victims[] = {2};
mixKrausMap(rho, victims, 1, map);
```
We can even directy mix density matrices together
```cpp
mixQureg(rho1, rho2, prob);
```

Sometimes we wish to left-multiply general operators upon density matrices without also right-multiplying their adjoint - i.e. our operators should _not_ be effected as unitaries. We can do this with the `multiply*()` functions.
```cpp
multiplyDiagMatrPower(rho, fullmatrix, 0.5);
```




--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_perform-calculations"></a>

## Perform calculations


After so much modification to our state, we will find that its amplitudes have differed substantially. But it's impractical to observe the exponentially-many amplitudes with [`reportQureg()`](https://quest-kit.github.io/QuEST/group__qureg__report.html#ga2a9df2538e537332b1aef8596ce337b2). We can instead give QuEST the [questions](https://quest-kit.github.io/QuEST/group__calculations.html) we wish to answer about the resulting state.

For example, we can find the [probability](https://quest-kit.github.io/QuEST/group__calc__prob.html) of measurement outcomes _without_ modifying the state.
```cpp
int outcome = 1;
qreal prob1 = calcProbOfQubitOutcome(qureg, target, outcome);

int qubits[]   = {2,3,4};
int outcomes[] = {0,1,1};
qreal prob2 = calcProbOfMultiQubitOutcome(qureg, qubits, outcomes, 3);
```
We can obtain _all_ outcome probabilities in one swoop:
```cpp
qreal probs[8];
calcProbsOfAllMultiQubitOutcomes(probs, qureg, qubits, 3);
```

> [!TIP]
> `C++` users can also obtain the result as a natural `std::vector<qreal>`.
> ```cpp
> auto probs = calcProbsOfAllMultiQubitOutcomes(qureg, {2,3,4});
> ```

It is similarly trivial to find [expectation values](https://quest-kit.github.io/QuEST/group__calc__expec.html)
```cpp
qreal expec1 = calcExpecPauliStr(qureg, getPauliStr("XYZIII"));
qreal expec2 = calcExpecPauliStrSum(qureg, sum);
qreal expec3 = calcExpecFullStateDiagMatr(qureg, fullmatrix);
```
or [distance measures](https://quest-kit.github.io/QuEST/group__calc__comparisons.html) between states, including between statevectors and density matrices.
```cpp
qreal pur = calcPurity(rho);
qreal fid = calcFidelity(rho, psi);
qreal dist = calcDistance(rho, psi);
```

We can even find reduced density matrices resulting from [partially tracing](https://quest-kit.github.io/QuEST/group__calc__partialtrace.html) out qubits.
```cpp
Qureg reduced = calcPartialTrace(qureg, targets, 3);

reportScalar("entanglement", calcPurity(reduced));
```



--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_report-the-results"></a>

## Report the results


We've seen above that [scalars](https://quest-kit.github.io/QuEST/group__types.html) can be reported, handling the pretty formatting of real and complex numbers, controlled by settings like [`setMaxNumReportedSigFigs()`](https://quest-kit.github.io/QuEST/group__debug__reporting.html#ga15d46e5d813f70b587762814964e1994). But we can also report every data structure in the QuEST API, such as Pauli strings
```cpp
reportPauliStr(
    getInlinePauliStr("XXYYZZ", {5,50, 10,60, 30,40})
);
```
```
YIIIIIIIIIXIIIIIIIIIZIIIIIIIIIZIIIIIIIIIIIIIIIIIIIYIIIIXIIIII
```
and their weighted sums
```cpp
reportPauliStrSum(sum);
```
```
PauliStrSum (4 terms, 160 bytes):
    1   II
    i   ZI
    i   IZ
    -1  ZZ
```
All outputs are affected by the [reporter settings](https://quest-kit.github.io/QuEST/group__debug__reporting.html).
```cpp
setMaxNumReportedItems(4,4);
setMaxNumReportedSigFigs(1);
reportCompMatr(bigmatrix);
```
```
CompMatr (8 qubits, 256x256 qcomps, 1 MiB):
    0.9-0.5i  0         …  0          0
    0         0.8-0.6i  …  0          0
        ⋮
    0         0         …  -0.5-0.9i  0
    0         0         …  0          0.4+0.9i
```



> [!NOTE]  
> Facilities for automatically logging to file are coming soon!



--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_cleanup"></a>

## Cleanup


While not strictly necessary before the program ends, it is a good habit to destroy data structures as soon as you are finished with them, freeing their memory.

```cpp
destroyQureg(qureg);
destroyCompMatr(bigmatrix);
destroyFullStateDiagMatr(fullmatrix);
destroyPauliStrSum(sum);
destroyKrausMap(map);
```



--------------------------------------------

<!-- permit doxygen to reference section -->
<a id="tutorial_finalise-quest"></a>

## Finalise QuEST


The _final_ [step](https://quest-kit.github.io/QuEST/group__environment.html#ga428faad4d68abab20f662273fff27e39) of our program should be to call
```cpp
finalizeQuESTEnv();
```
which ensures everything is synchronised, frees accelerator resources, and finalises MPI.
This is important because it ensures:
-  _everything is done_, and that distributed nodes that are still working (e.g. haven't yet logged to their own file) are not interrupted by early termination of another node.
- the MPI process ends gracefully, and doesn't spew out messy errors!
- our GPU processes are killed quickly, freeing resources for other processes.

> [!CAUTION]
> After calling `finalizeQuESTEnv()`, MPI will close and each if being accessed directly by the user, will enter an undefined state. Subsequent calls to MPI routines may return gibberish, and distributed machines will have lost their ability to communicate. It is recommended to call `finalizeQuESTEnv()` immediately before exiting.

You are now a QuEST expert 🎉 though there are _many_ more functions in the [API](https://quest-kit.github.io/QuEST/group__api.html) not covered here. Go forth and simulate!