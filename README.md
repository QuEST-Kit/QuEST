
<!-- This README avoids relative paths since it is also used as Doxygen's mainpage -->

<!-- banner -->
[<img src="https://github.com/QuEST-Kit/QuEST/blob/master/doxyconfig/banner.png?raw=true" alt="The QuEST logo" width=500>](https://quest.qtechtheory.org)


<!-- temporarily hiding Github Actions badges (pending aesthetic customisation)
    [![Ubuntu unit](https://github.com/QuEST-Kit/QuEST/workflows/Ubuntu%20unit/badge.svg?branch=develop)](https://github.com/QuEST-Kit/QuEST/actions)
    [![macOS unit](https://github.com/QuEST-Kit/QuEST/workflows/macOS%20unit/badge.svg)](https://github.com/QuEST-Kit/QuEST/actions)
    [![LLVM](https://github.com/QuEST-Kit/QuEST/workflows/LLVM%20asan/badge.svg)](https://github.com/QuEST-Kit/QuEST/actions)
 -->
<!-- temporarily hiding incorrect coverage statistics 
(currently only considers serial CPU; needs also GPU and distributed test contributions)
    [![codecov](https://codecov.io/gh/QuEST-Kit/QuEST/branch/develop/graph/badge.svg)](https://codecov.io/gh/QuEST-Kit/QuEST)
-->



<!-- action-badges has broken and is incorrectly showing build failure. 
     temporarily forcing badge to show the correct pass while replacement is sought 
     (that isn't Github's hideous default CI badge)
     [![unit tests](https://action-badges.now.sh/QuEST-Kit/QuEST)](https://github.com/QuEST-Kit/QuEST/actions) 
--> 
[![GitHub release](https://img.shields.io/github/release/QuEST-Kit/QuEST)](https://GitHub.com/QuEST-Kit/QuEST/releases/) 
[![Doc](https://img.shields.io/badge/doc-Github.io-orange.svg)](https://quest-kit.github.io/QuEST/modules.html)
[![unit tests](https://img.shields.io/badge/build-passing-green.svg)](https://github.com/QuEST-Kit/QuEST/actions)  <!-- forgive my sins (see above) -->
[![MIT license](https://img.shields.io/badge/license-MIT-lightgrey.svg)](LICENCE.txt)


The **Quantum Exact Simulation Toolkit** is a high performance simulator of quantum circuits, state-vectors and density matrices. QuEST uses **multithreading**, **GPU acceleration** and **distribution** to run lightning first on laptops, desktops and networked supercomputers. 
QuEST *just works*; it is stand-alone, requires no installation, and is trivial to compile and run. 
QuEST hybridises [OpenMP](https://www.openmp.org/) and [MPI](https://www.mpi-forum.org/) with _huge_ compiler support to run on all sorts of multicore, multi-CPU and distributed hardware, uses [HIP](https://docs.amd.com/bundle/HIP-Programming-Guide-v5.3) to run on AMD GPUs, integrates [cuQuantum](https://developer.nvidia.com/cuquantum-sdk) and [Thrust](https://developer.nvidia.com/thrust) for cutting-edge performance on modern NVIDIA GPUs, and has a custom kernel backend to run on older CUDA-compatible GPUs. 
And it hides these deployment modes behind a single, seamless interface.

[![Languages](https://img.shields.io/badge/C-99-ff69b4.svg)](http://www.open-std.org/jtc1/sc22/wg14/www/standards.html#9899)
[![Languages](https://img.shields.io/badge/C++-11-ff69b4.svg)](https://isocpp.org/wiki/faq/cpp11)
![OS](https://img.shields.io/badge/os-MacOS-9cbd3c.svg)
![OS](https://img.shields.io/badge/os-Linux-9cbd3c.svg)
![OS](https://img.shields.io/badge/os-Windows-9cbd3c.svg)
[![Platforms](https://img.shields.io/badge/multithreaded-OpenMP-6699ff.svg)](https://www.openmp.org/)
[![Platforms](https://img.shields.io/badge/distributed-MPI-6699ff.svg)](https://www.mpi-forum.org/) 
[![Platforms](https://img.shields.io/badge/GPU-CUDA-6699ff.svg)](https://developer.nvidia.com/cuda-zone)
[![Platforms](https://img.shields.io/badge/GPU-AMD-6699ff.svg)](https://docs.amd.com/bundle/HIP-Programming-Guide-v5.3)
[![Platforms](https://img.shields.io/badge/GPU-cuQuantum-6699ff.svg)](https://developer.nvidia.com/cuquantum-sdk)

QuEST is developed by the [QTechTheory](http://qtechtheory.org/) group at the University of Oxford, and [these authors](https://github.com/QuEST-Kit/QuEST/blob/master/AUTHORS.txt). To learn more:
- see the [tutorial](https://github.com/QuEST-Kit/QuEST/blob/master/examples/README.md)
- view the [documentation](https://quest-kit.github.io/QuEST/modules.html)
- visit the [website](https://quest.qtechtheory.org/)
- see some [examples](https://github.com/QuEST-Kit/QuEST/blob/master/examples/)
- read the [whitepaper](https://www.nature.com/articles/s41598-019-47174-9), which featured in Scientific Report's [Top 100 in Physics](https://www.nature.com/collections/ecehgdfcba/) :trophy:

[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41598--019--47174--9-yellow.svg)](https://doi.org/10.1038/s41598-019-47174-9)
[![Email](https://img.shields.io/badge/email-quest@materials.ox.ac.uk-red.svg)](mailto:quest@materials.ox.ac.uk)

---------------------------------


## :tada:&nbsp; Introduction

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

---------------------------------

## :white_check_mark:&nbsp; Features 
QuEST supports:  
- :ballot_box_with_check: &nbsp; **density matrices** for precise simulation of noisy quantum computers  
- :ballot_box_with_check: &nbsp; **general unitaries** with any number of control and target qubits  
- :ballot_box_with_check: &nbsp; **general decoherence channels** of any dimension  
- :ballot_box_with_check: &nbsp; **general Hermitian operators** in the Pauli basis  
- :ballot_box_with_check: &nbsp; **many *many* operators**, including even [Pauli gadgets](https://quest-kit.github.io/QuEST-develop-doc/group__unitary.html#ga34aa4865c92f9aa5d898c91286c9eca5), [analytic phase functions](https://quest-kit.github.io/QuEST-develop-doc/group__operator.html#ga467f517abd18dbc3d6fced84c6589161) and [Trotter circuits](https://quest-kit.github.io/QuEST-develop-doc/group__operator.html#ga35b6321c578a8c69470132b5ee95f930)  
- :ballot_box_with_check: &nbsp; **many tools to analyse** quantum states, such as calculations of [probability](https://quest-kit.github.io/QuEST-develop-doc/group__calc.html#gad0cc08d52cad5062553d6f78126780cc), [fidelity](https://quest-kit.github.io/QuEST-develop-doc/group__calc.html#gaa266ed6c8ae5d0d0f49e1ac50819cffc), and [expected value](https://quest-kit.github.io/QuEST-develop-doc/group__calc.html#ga82f17e96a4cb7612fb9c6ef856df3810)  
- :ballot_box_with_check: &nbsp; **variable precision** through a `qreal` numerical type which can use single, double or quad precision  
- :ballot_box_with_check: &nbsp; **QASM output** to verify simulated circuits  
- :ballot_box_with_check: &nbsp; **direct access to amplitudes** for rapid custom modification of the quantum state 
- :ballot_box_with_check: &nbsp; **native compilation** on MacOS, Linux and Windows, through Clang, GNU, Intel, and MSVC compilers

---------------------------------

## :book:&nbsp; Documentation

- The [tutorial](/examples/README.md) includes instructions for
  - [compiling](/examples/README.md#compiling) QuEST
  - [running](/examples/README.md#running) QuEST locally and on supercomputers
  - [testing](/examples/README.md#testing) QuEST using the comprehensive [unit tests](https://quest-kit.github.io/QuEST/group__unittest.html)
  <br><br>
- The [documentation](https://quest-kit.github.io/QuEST/modules.html) is divided into the following modules (collated [here](https://quest-kit.github.io/QuEST/QuEST_8h.html))
  - [data structures](https://quest-kit.github.io/QuEST/group__type.html)
  - [initialisations](https://quest-kit.github.io/QuEST/group__init.html)
  - [unitaries](https://quest-kit.github.io/QuEST/group__unitary.html)
  - [gates](https://quest-kit.github.io/QuEST/group__normgate.html)
  - [decoherence](https://quest-kit.github.io/QuEST/group__decoherence.html)
  - [operators](https://quest-kit.github.io/QuEST/group__operator.html)
  - [calculations](https://quest-kit.github.io/QuEST/group__calc.html)
  - [QASM](https://quest-kit.github.io/QuEST/group__qasm.html)
  <br><br>
- Additional utilities for debugging and testing are documented below
  - [debugging](https://quest-kit.github.io/QuEST/group__debug.html)
  - [unit tests](https://quest-kit.github.io/QuEST/group__unittest.html)
  - [testing utilities](https://quest-kit.github.io/QuEST/group__testutilities.html)


> **For developers**: QuEST's doc is automatically regenerated when the `master` branch is updated via [Github Actions](https://github.com/features/actions). To locally regenerate the doc, run `doxygen doxyconfig/config` in the root directory, which generates html documentation in `Doxygen_doc/html`.

---------------------------------

## :rocket:&nbsp; Getting started 

To rocket right in, download QuEST with [git](https://git-scm.com/) at the terminal
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
<br>

> **Windows** users should install [Build Tools](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019) for Visual Studio, and [CMake](https://cmake.org/download/), and run the above commmands in the *Developer Command Prompt for VS*, though using build commands
> ```bash 
> cmake .. -G "NMake Makefiles"
> nmake
> ```
> If using MSVC and NMake in this way fails, users can forego GPU acceleration, download
> [MinGW-w64](https://sourceforge.net/projects/mingw-w64/), and compile via 
> ```bash 
> cmake .. -G "MinGW Makefiles"
> make
> ```




---------------------------------

## :heart:&nbsp; Acknowledgements

We sincerely thank the following external contributors to QuEST.

- [Jakub Adamski](https://github.com/jjacobx) for optimising distributed communication of max-size messages.
- [Bruno Villasenor Alvarez](https://github.com/bvillasen) of [AMD](https://www.amd.com/en.html) for porting the GPU backend to [HIP](https://github.com/ROCm-Developer-Tools/HIP), for compatibility with AMD GPUs.
- [HQS Quantum simulations](https://quantumsimulations.de/) for contributing `mixDamping` on CPU.
- [Kshitij Chhabra](https://github.com/kshitijc) for patching some validation bugs.
- [Drew Silcock](https://github.com/drewsilcock) for patching the multithreaded build on MacOS.
- [Zach van Rijn](https://github.com/zv-io) for patching the multithreading code for GCC-9 OpenMP-5 compatibility.
- [SchineCompton](https://github.com/SchineCompton) for patching the GPU CMake release build.
- [Christopher J. Anders](https://github.com/chr5tphr) for patching the multithreading (when default off) and GPU builds (revising min cmake).
- [Gleb Struchalin](https://github.com/glebx-f) for patching the cmake standalone build.
- [Milos Prokop](https://github.com/Milos9304) for serial prototyping of `initDiagonalOpFromPauliHamil`.

QuEST uses the [mt19937ar](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) Mersenne Twister algorithm for random number generation, under the BSD licence. QuEST optionally (by additionally importing `QuEST_complex.h`) integrates the [language agnostic complex type](http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/2003/0303/cuj0303meyers/index.htm) by Randy Meyers and Dr. Thomas Plum

---------------------------------

## :newspaper:&nbsp; Related projects

- [QuESTlink](https://questlink.qtechtheory.org)   <br>
  a Mathematica package enabling symbolic circuit manipulation, analytic simulation, visualisation and high performance simulation with remote accelerated hardware.
  
- [pyQuEST](https://github.com/rrmeister/pyQuEST/tree/master)   <br>
  a python interface to QuEST, based on Cython, developed within the [QTechTheory](https://qtechtheory.org) group. Please note, pyQuEST is currently in the alpha stage.
   
- [PyQuEST-cffi](https://github.com/HQSquantumsimulations/PyQuEST-cffi)   <br>
  a python interface to QuEST based on cffi developed by HQS Quantum Simulations. Please note, PyQuEST-cffi is currently in the alpha stage and not an official QuEST project.
