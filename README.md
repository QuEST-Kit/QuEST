<img align="left" src="doxyconfig/logo.png" alt="The QuEST logo">

# [QuEST](https://quest.qtechtheory.org)

[![Ubuntu unit](https://github.com/QuEST-Kit/QuEST/workflows/Ubuntu%20unit/badge.svg)](https://github.com/QuEST-Kit/QuEST/actions)
[![macOS unit](https://github.com/QuEST-Kit/QuEST/workflows/macOS%20unit/badge.svg)](https://github.com/QuEST-Kit/QuEST/actions)
[![LLVM](https://github.com/QuEST-Kit/QuEST/workflows/LLVM%20asan/badge.svg)](https://github.com/QuEST-Kit/QuEST/actions)

<!--- 
temporarily hiding incorrect coverage statistics 
(currently only considers serial CPU; needs also GPU and distributed test contributions)
[![codecov](https://codecov.io/gh/QuEST-Kit/QuEST/branch/develop/graph/badge.svg)](https://codecov.io/gh/QuEST-Kit/QuEST)
--->

## Introduction

The **Quantum Exact Simulation Toolkit** is a high performance simulator of universal quantum circuits, state-vectors and density matrices. QuEST is written in C, hybridises OpenMP and MPI, and can run on a GPU. Needing only compilation, QuEST is easy to run both on laptops and supercomputers (in both C and C++), where it can take advantage of multicore, GPU-accelerated and networked machines to quickly simulate circuits on many qubits.

QuEST has a simple interface, independent of its run environment (on CPUs, GPUs or over networks),
```C
hadamard(qubits, 0);

controlledNot(qubits, 0, 1);

rotateY(qubits, 0, .1);
```
though is flexible
```C
Vector v;
v.x = 1; v.y = .5; v.z = 0;
rotateAroundAxis(qubits, 0, 3.14/2, v);
```
and powerful
```C
// sqrt(X) with pi/4 global phase
ComplexMatrix2 u = {
    .real = {{.5, .5}, { .5,.5}},
    .imag = {{.5,-.5}, {-.5,.5}}};
unitary(qubits, 0, u);

int controls[] = {1, 2, 3, 4, 5};
multiControlledUnitary(qureg, controls, 5, 0, u);
```

QuEST can simulate decoherence on mixed states, output [QASM](https://arxiv.org/abs/1707.03429), perform measurements, apply general unitaries with any number of control and target qubits, and boasts cheap/fast access to the underlying numerical representation of the state. QuEST offers precision-agnostic real and imaginary (additionally include `QuEST_complex.h`) number types, the precision of which can be modified at compile-time, as can the target hardware.

Learn more about QuEST at [quest.qtechtheory.org](https://quest.qtechtheory.org), or read the [whitepaper](https://www.nature.com/articles/s41598-019-47174-9). If you find QuEST useful, feel free to cite
```
Jones, T., Brown, A., Bush, I. et al. 
QuEST and High Performance Simulation of Quantum Computers. 
Sci Rep 9, 10736 (2019) doi:10.1038/s41598-019-47174-9
```
```
@article{Jones2019,
  title={{QuEST} and high performance simulation of quantum computers},
  author={Jones, Tyson and Brown, Anna and Bush, Ian and Benjamin, Simon C},
  journal={Scientific reports},
  volume={9},
  number={1},
  pages={1--11},
  year={2019},
  publisher={Nature Publishing Group}
}
```

---------------------------------

## Documentation

Full documentation is available at [quest.qtechtheory.org/docs](https://quest.qtechtheory.org/docs/), and the API is available [here](https://quest-kit.github.io/QuEST/modules.html) (all functions listed [here](https://quest-kit.github.io/QuEST/QuEST_8h.html)). See also the [tutorial](/examples/README.md).

> **For developers**: To regenerate the API doc after making changes to the code, run `doxygen doxyconfig/config` in the root directory. This will generate documentation in `Doxygen_doc/html`, the contents of which should be copied into [`docs/`](/docs/)). Make sure that `PROJECT_NUMBER` in `doxyconfig/config` is up to date!

---------------------------------

## Getting started

QuEST is contained entirely in the files in the `QuEST/` folder. To use QuEST, copy this folder to your computer and include `QuEST.h` in your `C` or `C++` code, and compile using cmake with the provided [CMakeLists.txt file](/CMakeLists.txt). See the [tutorial](/examples/README.md) for an introduction. We also include example [submission scripts](examples/submissionScripts/) for using QuEST with SLURM and PBS. 

### Quick Start

#### MacOS & Linux

MacOS and Linux users can clone this repository to your machine through the terminal:
```bash
git clone https://github.com/quest-kit/QuEST.git
cd QuEST
```
Compile the [example](examples/tutorial_example.c) using
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


#### Tests
Additionally, you can run QuEST's rigorous unit tests in your own environment, 
which should take no longer than ten minutes.
```bash
mkdir build
cd build
cmake .. -DTESTING=ON
make 
make test
```

---------------------------------

## Contact

To file a bug report or feature request, [raise a github issue](https://github.com/QuEST-Kit/QuEST/issues). For additional support, email quest@materials.ox.ac.uk. You can view the list of contributors to QuEST in [`AUTHORS.txt`](AUTHORS.txt).

---------------------------------

## Acknowledgements

We sincerely thank the following external contributors to QuEST.

- [HQS Quantum simulations](https://quantumsimulations.de/) for contributing `mixDamping` on CPU.
- [Kshitij Chhabra](https://github.com/kshitijc) for patching some validation bugs.
- [Drew Silcock](https://github.com/drewsilcock) for patching the multithreaded build on MacOS.
- [Zach van Rijn](https://github.com/zv-io) for patching the multithreading code for GCC-9 OpenMP-5 compatibility
- @SachinCompton for patching the GPU CMake Release build.
- [Christopher J. Anders](https://github.com/chr5tphr) for patching the multithreading (when default off) and GPU builds (revising min cmake)

QuEST uses the [mt19937ar](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) Mersenne Twister algorithm for random number generation, under the BSD licence. QuEST optionally (by additionally importing `QuEST_complex.h`) integrates the [language agnostic complex type](http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/2003/0303/cuj0303meyers/index.htm) by Randy Meyers and Dr. Thomas Plum

---------------------------------

## Licence

QuEST is released under a [MIT Licence](LICENCE.txt)

---------------------------------

## Related projects -- QuEST utilities and extensions

* [QuESTlink](https://questlink.qtechtheory.org): a Mathematica package allowing symbolic circuit manipulation and high performance simulation with remote accelerated hardware.
* [PyQuEST-cffi](https://github.com/HQSquantumsimulations/PyQuEST-cffi): a python interface to QuEST based on cffi developed by HQS Quantum Simulations. Please note, PyQuEST-cffi is currently in the alpha stage and not an official QuEST project.
