<!-- 
  QuEST README page and doxygen-doc mainpage.

  Because this file doubles as the doxygen mainpage, it must use absolute paths.

  @author Tyson Jones
-->


<!-- banner and top badges (centered) -->
<div align="center">

  <!-- banner -->
  <a href="https://quest.qtechtheory.org">
    <img src="https://raw.githubusercontent.com/QuEST-Kit/QuEST/refs/heads/main/utils/docs/logos/banner.png" alt="The QuEST logo" width=400>
  </a>

  <!-- TODO: restore CI 'compilation/test pass' badge! -->
  [![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41598--019--47174--9-yellow.svg)](https://doi.org/10.1038/s41598-019-47174-9)
  <br>
  [![GitHub release](https://img.shields.io/github/release/QuEST-Kit/QuEST)](https://GitHub.com/QuEST-Kit/QuEST/releases/) 
  [![Doc](https://img.shields.io/badge/doc-Github.io-orange.svg)](https://quest-kit.github.io/QuEST/modules.html)
  [![MIT license](https://img.shields.io/badge/license-MIT-lightgrey.svg)](LICENCE.txt)


</div>



<!-- intro -->

> [!NOTE]
> QuEST `v4` has been released which re-designed QuEST from the ground up. Read about the exciting new features [here](docs/v4.md).

The **Quantum Exact Simulation Toolkit** (QuEST) is a high-performance simulator of quantum statevectors and density matrices.
It hybridises **multithreading**, **GPU acceleration** and **distribution** to run lightning fast on laptops, desktops and 
supercomputers, parallelising over multiple cores, CPUs and GPUs. Behind the scenes, QuEST leverages [OpenMP](https://www.openmp.org/),
[MPI](https://www.mpi-forum.org/), [CUDA](https://developer.nvidia.com/cuda-zone), [HIP](https://rocm.docs.amd.com/projects/HIP/en/docs-develop/what_is_hip.html),
[Thrust](https://developer.nvidia.com/thrust), [cuQuantum](https://developer.nvidia.com/cuquantum-sdk) and [GPUDirect](https://developer.nvidia.com/gpudirect)
for cutting-edge performance on modern multi-GPU clusters, and compatibility with older NVIDIA and AMD GPUs. These deployments can 
be combined in any combination, or automatically decided at runtime, yet are abstracted behind a single, seamless interface, accessible 
by both `C` and `C++` and all the major compilers (detailed [here](docs/compilers.md)).


<!-- detail badgets (centered) -->

<div align="center">

  [![Languages](https://img.shields.io/badge/C-11-ff69b4.svg)](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3631.pdf)
  [![Languages](https://img.shields.io/badge/C++-14-ff69b4.svg)](https://isocpp.org/wiki/faq/cpp14)
  ![OS](https://img.shields.io/badge/os-MacOS-9cbd3c.svg)
  ![OS](https://img.shields.io/badge/os-Linux-9cbd3c.svg)
  ![OS](https://img.shields.io/badge/os-Windows-9cbd3c.svg) <br>
  [![Platforms](https://img.shields.io/badge/multithreaded-OpenMP-6699ff.svg)](https://www.openmp.org/)
  [![Platforms](https://img.shields.io/badge/distributed-MPI-6699ff.svg)](https://www.mpi-forum.org/) 
  [![Platforms](https://img.shields.io/badge/GPU-CUDA-6699ff.svg)](https://developer.nvidia.com/cuda-zone)
  [![Platforms](https://img.shields.io/badge/GPU-AMD-6699ff.svg)](https://docs.amd.com/bundle/HIP-Programming-Guide-v5.3)

</div>


QuEST development is led by the [QTechTheory](http://qtechtheory.org/) group at the University of Oxford, with active contributions from the [EPCC](https://www.epcc.ed.ac.uk/) team at the University of Edinburgh, and support from the below organisations.
In particular, QuEST `v4` was made possible through the support of the UK National Quantum Computing centre (_NQCC200921_) and the [UKRI SEEQA](https://gtr.ukri.org/projects?ref=EP%2FY004655%2F1#/tabOverview) project.

<div align="center">

  <img src="https://raw.githubusercontent.com/QuEST-Kit/QuEST/refs/heads/main/utils/docs/logos/nqcc.png" alt="NQCC" height=30> &nbsp;
  <img src="https://raw.githubusercontent.com/QuEST-Kit/QuEST/refs/heads/main/utils/docs/logos/amd.png" alt="AMD" height=25> &nbsp;
  <img src="https://raw.githubusercontent.com/QuEST-Kit/QuEST/refs/heads/main/utils/docs/logos/nvidia.png" alt="NVIDIA" height=25> &nbsp;
  <img src="https://raw.githubusercontent.com/QuEST-Kit/QuEST/refs/heads/main/utils/docs/logos/qmt.png" alt="Quantum Motion" height=25> &nbsp;
  <img src="https://raw.githubusercontent.com/QuEST-Kit/QuEST/refs/heads/main/utils/docs/logos/edinburgh.png" alt="University of Edinburgh" height=25> &nbsp;
  <img src="https://raw.githubusercontent.com/QuEST-Kit/QuEST/refs/heads/main/utils/docs/logos/oxford.png" alt="University of Oxford" height=28> &nbsp;

</div>


<!-- <a> used below for doxygen compatibility -->

To learn more:
- view the <a href="#main_documentation">documentation</a>
- visit the [website](https://quest.qtechtheory.org/)
- read the [whitepaper](https://www.nature.com/articles/s41598-019-47174-9), which featured in Scientific Report's [Top 100 in Physics](https://www.nature.com/collections/ecehgdfcba/) :trophy:



---------------------------------


<!-- BEWARE that we use two non-breaking spaces after each emoji in
     a section title, to add spacing between emoji and text -->

## 🎉  Introduction

QuEST has a simple interface which is agnostic to whether it's running on CPUs, GPUs or a networked supercomputer.
```cpp
Qureg qureg = createQureg(30);
initRandomPureState(qureg);

applyHadamard(qureg, 0);
applyControlledRotateX(qureg, 0, 1, angle);
applyFullQuantumFourierTransform(qureg);

reportQureg(qureg);

qreal prob  = calcProbOfQubitOutcome(qureg, 0, outcome);
qreal expec = calcExpecPauliStr(qureg, getPauliStr("XYZ"));
```
Yet, it is flexible
```cpp
mixDepolarising(qureg, targ, prob);
mixKrausMap(qureg, targs, ntargs, krausmap);

applyQubitMeasurement(qureg, targ);
applyMultiQubitProjector(qureg, targs, outcomes, ntargs);

applyControlledPauliGadget(qureg, ctrl, paulistr, angle);
applyMultiStateControlledSwap(qureg, ctrls, states, nctrls, targ1, targ2);

multiplyCompMatr1(qureg, targ, getInlineCompMatr1( {{1,2i},{3i,4}} ));
multiplyDiagMatrPower(qureg, targs, ntargs, diagmatr, exponent);
```
and extremely powerful
```cpp
setFullStateDiagMatrFromMultiVarFunc(fullmatr, myfunc, ntargsPerVar, nvars);
applyFullStateDiagMatrPower(qureg, fullmatr, exponent);

CompMatr matr = createCompMatr(10);
setCompMatr(matr, ...);
applyCompMatr(qureg, targs, 10, matr);

PauliStrSum observ = createInlinePauliStrSum(R"(
    0.123 XXIXX
    1.23i XYZXZ
    -1-6i IIZII
)");
applyTrotterizedPauliStrSumGadget(qureg, observ, time, order, nreps);

Qureg reduce = calcPartialTrace(qureg, targs, ntargs);
qreal expec1 = calcExpecPauliStrSum(reduce, observ);
qreal expec2 = calcExpecFullStateDiagMatr(qureg, fullmatr);
```

---------------------------------

## ✅  Features 

<!-- BEWARE that a bug in Doxygen v1.13.2 (github.com/doxygen/doxygen/issues/11515)
     means we cannot immediately follow a non-breaking space (inserted below after
     each emoji to effect spacing) with markdown syntax like **. Instead, we insert
     one final regular/non-breaking space before ** which isn't rendered, and which
     works around the bug -->

QuEST supports:  
- ☑️   **density matrices** for precise simulation of noisy quantum computers  
- ☑️   **general unitaries** with any number of control, control-states, and target qubits  
- ☑️   **general decoherence channels** of any dimension  
- ☑️   **general observables** in the Pauli or diagonal-Z bases  
- ☑️   **many *many* operators**, including Pauli gadgets, trotterised time evolutions, and projectors
- ☑️   **many tools to analyse** quantum states, such as calculations of probability, fidelity, expectation value, distances and partial traces
- ☑️   **variable precision** through `qreal` and `qcomp` numerical types which can use single, double or quad precision  
- ☑️   **direct access to amplitudes** for rapid custom modification of the quantum state 
- ☑️   **native compilation** on MacOS, Linux and Windows, through Clang, GNU, Intel, and MSVC compilers
- ☑️   **hybridisation** of multithreading, GPU-acceleration, distribution and GPU-distribution
- ☑️   **optimisation** using NVLink'd GPUs, cuQuantum, and CUDA-aware MPI
- ☑️   **automatic deployment** of a `Qureg` to the optimal hardware at runtime
- ☑️   **hardware probing** to determine how many qubits can be simulated at runtime
- ☑️   **bespoke algorithms** to optimally simulate a wide variety of esoteric operations

---------------------------------

<!-- permit doxygen to reference section -->
<a id="main_documentation"></a>

## 📖  Documentation

> [!IMPORTANT]
> QuEST v4's documentation is still under construction!

Visit the [docs](docs/README.md) for guides and tutorials, or the [API](https://quest-kit.github.io/QuEST/group__api.html) which is divided into:
  - [calculations](https://quest-kit.github.io/QuEST/group__calculations.html)
  - [channels](https://quest-kit.github.io/QuEST/group__channels.html)
  - [debug](https://quest-kit.github.io/QuEST/group__debug.html)
  - [decoherence](https://quest-kit.github.io/QuEST/group__decoherence.html)
  - [environment](https://quest-kit.github.io/QuEST/group__environment.html)
  - [initialisations](https://quest-kit.github.io/QuEST/group__initialisations.html)
  - [matrices](https://quest-kit.github.io/QuEST/group__matrices.html)
  - [modes](https://quest-kit.github.io/QuEST/group__modes.html)
  - [operations](https://quest-kit.github.io/QuEST/group__operations.html)
  - [paulis](https://quest-kit.github.io/QuEST/group__paulis.html)
  - [precision](https://quest-kit.github.io/QuEST/group__precision.html)
  - [qureg](https://quest-kit.github.io/QuEST/group__qureg.html)
  - [types](https://quest-kit.github.io/QuEST/group__types.html)
  - [tests](https://quest-kit.github.io/QuEST/group__tests.html)

<!-- hiding test doc since too large, but preserved here for useful links -->
<!--
You can also browse QuEST's extensive [tests](https://quest-kit.github.io/QuEST/group__tests.html):
  - [test utilities](https://quest-kit.github.io/QuEST/group__testutils.html)
    - [cache](https://quest-kit.github.io/QuEST/group__testutilscache.html)
    - [compare](https://quest-kit.github.io/QuEST/group__testutilscompare.html)
    - [convert](https://quest-kit.github.io/QuEST/group__testutilsconvert.html)
    - [evolve](https://quest-kit.github.io/QuEST/group__testutilsevolve.html)
    - [linalg](https://quest-kit.github.io/QuEST/group__testutilslinalg.html)
    - [lists](https://quest-kit.github.io/QuEST/group__testutilslists.html)
    - [macros](https://quest-kit.github.io/QuEST/group__testutilsmacros.html)
    - [measure](https://quest-kit.github.io/QuEST/group__testutilsmeasure.html)
    - [qmatrix](https://quest-kit.github.io/QuEST/group__testutilsqmatrix.html)
    - [qvector](https://quest-kit.github.io/QuEST/group__testutilsqvector.html)
    - [random](https://quest-kit.github.io/QuEST/group__testutilsrandom.html)
  - [unit tests](https://quest-kit.github.io/QuEST/group__unittests.html)
    - [calculations](https://quest-kit.github.io/QuEST/group__unitcalcs.html)
    - [channels](https://quest-kit.github.io/QuEST/group__unitchannels.html)
    - [debug](https://quest-kit.github.io/QuEST/group__unitdebug.html)
    - [decoherence](https://quest-kit.github.io/QuEST/group__unitdeco.html)
    - [environment](https://quest-kit.github.io/QuEST/group__unitenv.html)
    - [initialisations](https://quest-kit.github.io/QuEST/group__unitinit.html)
    - [matrices](https://quest-kit.github.io/QuEST/group__unitmatr.html)
    - [operations](https://quest-kit.github.io/QuEST/group__unitops.html)
    - [paulis](https://quest-kit.github.io/QuEST/group__unitpaulis.html)
    - [qureg](https://quest-kit.github.io/QuEST/group__unitqureg.html)
    - [types](https://quest-kit.github.io/QuEST/group__unittypes.html)
  - [integration tests](https://quest-kit.github.io/QuEST/group__integrationtests.html)
  - [deprecated tests](https://quest-kit.github.io/QuEST/group__deprecatedtests.html)
  - [deprecated test utilities](https://quest-kit.github.io/QuEST/group__deprecatedutils.html)
-->

---------------------------------

## 🚀  Getting started 

To rocket right in, download QuEST with [git](https://git-scm.com/) at the terminal
```bash
git clone https://github.com/quest-kit/QuEST.git
cd QuEST
```
We recommend working in a `build` directory:
```bash
mkdir build
cd build
```

Compile the [minimum example](/examples/tutorials/min_example.cpp) using [cmake](https://cmake.org/):
```bash
cmake .. 
make
```
then run it with
```bash
./min_example
```

See the [docs](docs/README.md) for enabling acceleration and running the unit tests.

---------------------------------

## ❤  Acknowledgements

In addition to QuEST's [authors](AUTHORS.txt), we sincerely thank the following external contributors to QuEST.

- [Jakub Adamski](https://github.com/jjacobx) for optimising distributed communication of max-size messages.
- [Bruno Villasenor Alvarez](https://github.com/bvillasen) of [AMD](https://www.amd.com/en.html) for porting the v3 GPU backend to [HIP](https://github.com/ROCm-Developer-Tools/HIP), for compatibility with AMD GPUs.
- [HQS Quantum simulations](https://quantumsimulations.de/) for prototyping v3 `mixDamping`.
- [Kshitij Chhabra](https://github.com/kshitijc) for patching some v3 validation bugs.
- [Drew Silcock](https://github.com/drewsilcock) for patching the v3 multithreaded build on MacOS.
- [Zach van Rijn](https://github.com/zv-io) for patching the v3 multithreading code for GCC-9 OpenMP-5 compatibility.
- [SchineCompton](https://github.com/SchineCompton) for patching the v3 GPU CMake release build.
- [Christopher J. Anders](https://github.com/chr5tphr) for patching the v3 multithreading (when default off) and GPU builds (revising min cmake).
- [Gleb Struchalin](https://github.com/glebx-f) for patching the v3 cmake standalone build.
- [Milos Prokop](https://github.com/Milos9304) for serial prototyping of v3's `initDiagonalOpFromPauliHamil`.



---------------------------------

## 📰  Related projects

- [QuESTlink](https://questlink.qtechtheory.org)   <br>
  a Mathematica package enabling symbolic circuit manipulation, analytic simulation, visualisation and high performance simulation with remote accelerated hardware.
  
- [pyQuEST](https://github.com/rrmeister/pyQuEST)   <br>
  a python interface to QuEST, based on Cython, developed by [Richard Meister](https://github.com/rrmeister) within the [QTechTheory](https://qtechtheory.org) group.

- [QuEST.jl](https://github.com/ediparquantum/QuEST.jl) <br>
  a Julia interface to QuEST, developed by [Jonathan Miller](https://github.com/fieldofnodes).

- [qoqo-quest](https://github.com/HQSquantumsimulations/qoqo-quest) <br>
  a Rust interface to QuEST, developed by [HQS Quantum Simulations](https://quantumsimulations.de/).
   
- [PyQuEST-cffi](https://github.com/HQSquantumsimulations/PyQuEST-cffi)   <br>
  a python interface to QuEST developed by [HQS Quantum Simulations](https://quantumsimulations.de/). 
