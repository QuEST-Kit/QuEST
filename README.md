
<!-- We avoid a relative path of 'doxyconfig/logo.png' here, because this README
     will also serve as Doxygen's mainpage in the doc, where the root directory 
     will differ -->
<img align="left" src="https://github.com/QuEST-Kit/QuEST/blob/master/doxyconfig/logo.png?raw=true" alt="The QuEST logo">

# [QuEST](https://quest.qtechtheory.org)


<!-- temporarily hiding Github Actions badges (pending aesthetic customisation)
    [![Ubuntu unit](https://github.com/QuEST-Kit/QuEST/workflows/Ubuntu%20unit/badge.svg?branch=develop)](https://github.com/QuEST-Kit/QuEST/actions)
    [![macOS unit](https://github.com/QuEST-Kit/QuEST/workflows/macOS%20unit/badge.svg)](https://github.com/QuEST-Kit/QuEST/actions)
    [![LLVM](https://github.com/QuEST-Kit/QuEST/workflows/LLVM%20asan/badge.svg)](https://github.com/QuEST-Kit/QuEST/actions)
 -->
<!-- temporarily hiding incorrect coverage statistics 
(currently only considers serial CPU; needs also GPU and distributed test contributions)
    [![codecov](https://codecov.io/gh/QuEST-Kit/QuEST/branch/develop/graph/badge.svg)](https://codecov.io/gh/QuEST-Kit/QuEST)
-->

[![GitHub release](https://img.shields.io/github/release/QuEST-Kit/QuEST)](https://GitHub.com/QuEST-Kit/QuEST/releases/) 
[![Doc](https://img.shields.io/badge/doc-Github.io-orange.svg)](https://quest-kit.github.io/QuEST/modules.html)
[![unit tests](https://action-badges.now.sh/QuEST-Kit/QuEST)](https://github.com/QuEST-Kit/QuEST/actions)
[![MIT license](https://img.shields.io/badge/license-MIT-lightgrey.svg)](LICENCE.txt)


The **Quantum Exact Simulation Toolkit** is a high performance simulator of quantum circuits, state-vectors and density matrices. QuEST uses **multithreading**, **GPU acceleration** and **distribution** to run lightning first on laptops, desktops and networked supercomputers. QuEST *just works*; it is stand-alone, and trivial to compile and getting running.

[![Languages](https://img.shields.io/badge/C-99-ff69b4.svg)](http://www.open-std.org/jtc1/sc22/wg14/www/standards.html#9899)
[![Languages](https://img.shields.io/badge/C++-11-ff69b4.svg)](http://www.open-std.org/jtc1/sc22/wg14/www/standards.html#9899)
[![Platforms](https://img.shields.io/badge/multithreaded-OpenMP-6699ff.svg)](https://www.openmp.org/)
[![Platforms](https://img.shields.io/badge/GPU-CUDA-6699ff.svg)](https://developer.nvidia.com/cuda-zone)
[![Platforms](https://img.shields.io/badge/distributed-MPI-6699ff.svg)](https://www.mpi-forum.org/) 

QuEST is developed by the [QTechTheory](http://qtechtheory.org/) group at the University of Oxford, and [these authors](https://github.com/QuEST-Kit/QuEST/blob/master/AUTHORS.txt). To learn more:
- see the [tutorial](https://github.com/QuEST-Kit/QuEST/blob/master/examples/README.md)
- view the [documentation](https://quest-kit.github.io/QuEST/modules.html)
- visit our [website](https://quest.qtechtheory.org/)
- read our [whitepaper](https://www.nature.com/articles/s41598-019-47174-9), which featured in Scientific Report's [Top 100 in Physics](https://www.nature.com/collections/ecehgdfcba/)

[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41598--019--47174--9-yellow.svg)](https://doi.org/10.1038/s41598-019-47174-9)
[![Email](https://img.shields.io/badge/email-quest@materials.ox.ac.uk-red.svg)](mailto:quest@materials.ox.ac.uk)

## Introduction

QuEST has a simple interface, which is agnostic to its runtime environment, between CPUs, GPUs and over networks.
```C
hadamard(qubits, 0);

controlledRotateX(qubits, 0, 1, angle);

double prob = calcProbOfOutcome(qubits, 0, outcome);
```
Yet, it is flexible
```C
Vector v;
v.x = 1; v.y = .5; v.z = 0;
rotateAroundAxis(qubits, 0, angle, v);

ComplexMatrix2 u = {
    .real = {{.5, .5}, { .5,.5}},
    .imag = {{.5,-.5}, {-.5,.5}}};
unitary(qubits, 0, u);

mixDepolarising(qubits, 0, prob);
```
and extremely powerful
```C
ComplexMatrixN u = createComplexMatrixN(5);
int ctrls[] = {0, 1, 2};
int targs[] = {5, 20, 15, 10, 25};
multiControlledMultiQubitUnitary(qubits, ctrls, 3, targs, 5, u);

ComplexMatrixN k1, k2, k3 = ...
mixMultiQubitKrausMap(qubits, targs, 5, {k1, k2, k3}, 3);

double val = calcExpecPauliHamil(qubits, hamiltonian, workspace);

applyTrotterCircuit(qubits, hamiltonian, time, order, repetitions);
```

## Features 
QuEST supports:  
- :white_check_mark: &nbsp; **density matrices** for precise simulation of noisy quantum computers  
- :white_check_mark: &nbsp; **general unitaries** with any number of control and target qubits  
- :white_check_mark: &nbsp; **general decoherence channels** of any dimension  
- :white_check_mark: &nbsp; **general Hermitian operators** in the Pauli basis  
- :white_check_mark: &nbsp; **many *many* operators**, including even [Pauli gadgets](https://quest-kit.github.io/QuEST-develop-doc/group__unitary.html#ga34aa4865c92f9aa5d898c91286c9eca5), [analytic phase functions](https://quest-kit.github.io/QuEST-develop-doc/group__operator.html#ga467f517abd18dbc3d6fced84c6589161) and [Trotter circuits](https://quest-kit.github.io/QuEST-develop-doc/group__operator.html#ga35b6321c578a8c69470132b5ee95f930)  
- :white_check_mark: &nbsp; **many tools to analyse** quantum states, such as calculations of [probability](https://quest-kit.github.io/QuEST-develop-doc/group__calc.html#gad0cc08d52cad5062553d6f78126780cc), [fidelity](https://quest-kit.github.io/QuEST-develop-doc/group__calc.html#gaa266ed6c8ae5d0d0f49e1ac50819cffc), and [expected value](https://quest-kit.github.io/QuEST-develop-doc/group__calc.html#ga82f17e96a4cb7612fb9c6ef856df3810)  
- :white_check_mark: &nbsp; **variable precision** through a `qreal` numerical type which can use single, double or quad precision  
- :white_check_mark: &nbsp; **QASM output** to verify simulated circuits  
- :white_check_mark: &nbsp; **direct access to amplitudes** for rapid custom modification of the quantum state  

---------------------------------

## Documentation

- The [tutorial](/examples/README.md) includes instructions for
  - [compiling](/examples/README.md#compiling) QuEST
  - [running](/examples/README.md#running) QuEST locally and on supercomputers
  - [testing](/examples/README.md#testing) QuEST using the comprehensive [unit tests](https://quest-kit.github.io/QuEST/group__unittest.html)
  <br><br>
- The [documentation](https://quest-kit.github.io/QuEST/modules.html) is divided into the following modules (collated [here](https://quest-kit.github.io/QuEST/QuEST_8h.html)):
  - [data structures](https://quest-kit.github.io/QuEST/group__type.html)
  - [initialisations](https://quest-kit.github.io/QuEST/group__init.html)
  - [unitaries](https://quest-kit.github.io/QuEST/group__unitary.html)
  - [gates](https://quest-kit.github.io/QuEST/group__normgate.html)
  - [decoherence](https://quest-kit.github.io/QuEST/group__decoherence.html)
  - [operators](https://quest-kit.github.io/QuEST/group__operator.html)
  - [calculations](https://quest-kit.github.io/QuEST/group__calc.html)
  - [QASM](https://quest-kit.github.io/QuEST/group__qasm.html)
  <br><br>
- Additionally, 
  - [debugging](https://quest-kit.github.io/QuEST/group__debug.html)
  - [unit tests](https://quest-kit.github.io/QuEST/group__unittest.html)
  - [testing utilities](https://quest-kit.github.io/QuEST/group__testutilities.html)


> **For developers**: To regenerate the doc after making changes to the code, run `doxygen doxyconfig/config` in the root directory. This will generate documentation in `Doxygen_doc/html`, the contents of which should be copied into [`docs/`](/docs/)) 

---------------------------------

## Getting started 

To dive right in, download QuEST with [git](https://git-scm.com/) at the terminal
```bash
git clone https://github.com/quest-kit/QuEST.git
cd QuEST
```
Compile the [tutorial](/examples/README.md) example ([source](/examples/tutorial_example.c)) using [cmake](https://cmake.org/) and [make](https://www.gnu.org/software/make/)
```bash
mkdir build
cd build
cmake ..
make
```
then run it with
```bash
./demo
```


OLD GETTING STARTED BELOW


#### Windows 

Windows users should install [Build Tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019) for Visual Studio, [CMake](https://cmake.org/download/) and [MinGW-w64](https://sourceforge.net/projects/mingw-w64/). 
Then, in a *Developer Command Prompt for VS*, run
```bash
git clone "https://github.com/quest-kit/QuEST.git"
cd QuEST
```
```bash
mkdir build
cd build
cmake .. -G "MinGW Makefiles"
make
```
```bash
demo.exe
```



---------------------------------

## Acknowledgements

We sincerely thank the following external contributors to QuEST.

- [HQS Quantum simulations](https://quantumsimulations.de/) for contributing `mixDamping` on CPU.
- [Kshitij Chhabra](https://github.com/kshitijc) for patching some validation bugs.
- [Drew Silcock](https://github.com/drewsilcock) for patching the multithreaded build on MacOS.
- [Zach van Rijn](https://github.com/zv-io) for patching the multithreading code for GCC-9 OpenMP-5 compatibility.
- [SchineCompton](https://github.com/SchineCompton) for patching the GPU CMake release build.
- [Christopher J. Anders](https://github.com/chr5tphr) for patching the multithreading (when default off) and GPU builds (revising min cmake).
- [Gleb Struchalin](https://github.com/glebx-f) for patching the cmake standalone build

QuEST uses the [mt19937ar](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) Mersenne Twister algorithm for random number generation, under the BSD licence. QuEST optionally (by additionally importing `QuEST_complex.h`) integrates the [language agnostic complex type](http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/2003/0303/cuj0303meyers/index.htm) by Randy Meyers and Dr. Thomas Plum

---------------------------------

## Related projects

* [QuESTlink](https://questlink.qtechtheory.org): a Mathematica package allowing symbolic circuit manipulation and high performance simulation with remote accelerated hardware.
* [PyQuEST-cffi](https://github.com/HQSquantumsimulations/PyQuEST-cffi): a python interface to QuEST based on cffi developed by HQS Quantum Simulations. Please note, PyQuEST-cffi is currently in the alpha stage and not an official QuEST project.
