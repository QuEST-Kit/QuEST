# [QuEST](https://quest.qtechtheory.org)

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
ComplexMatrix2 u;
u.r0c0 = (Complex) {.real=.5, .imag= .5};
u.r0c1 = (Complex) {.real=.5, .imag=-.5}; 
u.r1c0 = (Complex) {.real=.5, .imag=-.5};
u.r1c1 = (Complex) {.real=.5, .imag= .5};
unitary(qubits, 0, u);

int[] controls = {1, 2, 3, 4, 5};
multiControlledUnitary(qureg, controls, 5, 0, u);
```

QuEST can simulate decoherence on mixed states, output [QASM](https://arxiv.org/abs/1707.03429), perform measurements, apply gates with any number of control qubits, and provides cheap/fast access to the underlying statevector. QuEST offers precision-agnostic real and imaginary (additionally include `QuEST_complex.h`) number types, the precision of which can be modified at compile-time, as can the target hardware.

Learn more about QuEST at [quest.qtechtheory.org](https://quest.qtechtheory.org).

## Getting started

QuEST is contained entirely in the files in the `QuEST/` folder. To use QuEST, copy this folder to your computer and include `QuEST.h` in your `C` or `C++` code, and compile using cmake with the provided [CMakeLists.txt](CMakeLists.txt file). See the [tutorial](/examples/README.md) for an introduction, and view the full API [here](https://quest-kit.github.io/QuEST/QuEST_8h.html).

We also include example [submission scripts](examples/submissionScripts/) for using QuEST with SLURM and PBS. 

### Quick Start

Copy or clone this repository to your machine. E.g. in the desired directory, enter
```bash
git clone https://github.com/quest-kit/QuEST.git
cd QuEST
```
at terminal. You can then compile the [example](examples/tutorial_example.c) using
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
and afterward, clean up with
```bash
make clean
````

or, to remove the build directory entirely, from the root directory
```bash
rm -r build
```

The program will print information about your execution environment and some simple operations on a three qubit system. See the [tutorial](examples/README.md) for a better introduction. 

Additionally, you can run unit tests to see if QuEST runs correctly in your environment, using
```bash
make test
```

This requires Python 3.4+. 

## Documentation

View the API [here](https://quest-kit.github.io/QuEST/QuEST_8h.html), and check compatible compiler versions [here](tests/compilers/compatibility.md).

> For developers: To recreate the full documentation after making changes to the code, run doxygen doxyconf in the root directory. This will generate documentation in Doxygen_doc/html, and can be accessed through index.html in that folder. 

## Acknowledgements

QuEST uses the [mt19937ar](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) Mersenne Twister algorithm for random number generation, under the BSD licence. QuEST optionally (by additionally importing `QuEST_complex.h`) integrates the [language agnostic complex type](http://collaboration.cmc.ec.gc.ca/science/rpn/biblio/ddj/Website/articles/CUJ/2003/0303/cuj0303meyers/index.htm) by Randy Meyers and Dr. Thomas Plum

Thanks to [HQS Quantum simulations](https://quantumsimulations.de/) for contributing the applyOneQubitDampingError function.

## Licence

QuEST is released under a [MIT Licence](https://github.com/quest-kit/QuEST/blob/master/LICENCE.txt)


## Related projects -- QuEST utilities and extensions

* [PyQuEST-cffi](https://github.com/HQSquantumsimulations/PyQuEST-cffi): a python interface to QuEST based on cffi developed by HQS Quantum Simulations. Please note, PyQuEST-cffi is currently in the alpha stage and not an official QuEST project.


