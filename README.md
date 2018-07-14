# [QuEST](https://quest.qtechtheory.org)

## Versions


Latest version: [1.1.0](https://github.com/aniabrown/QuEST/releases/tag/v1.1.0) 

Please report errors or feedback to anna.brown@oerc.ox.ac.uk 

## Quick Start

Copy or clone this repository to your machine. E.g. in the desired directory, enter
```bash
git clone https://github.com/aniabrown/QuEST.git
```
at terminal.

In the root directory, compile the [example](examples/tutorial_example.c) using
```bash
cp examples/tutorial_example.c tutorial_example.c
make
```
then run it with
```bash
./demo
```

The program will print information about your execution environment and some simple operations on a three qubit system. See the [tutorial](examples/README.md) for a better introduction. Additionally, run `tests/runTests.sh` to test QuEST runs correctly in your environment.


## Introduction

The **Quantum Exact Simulation Toolkit** is a high performance simulator of universal quantum circuits. QuEST is written in C, hybridises OpenMP and MPI, and can run on a GPU. Needing only compilation, QuEST is easy to run both on laptops and supercomputers, where it can take advantage of multicore, GPU-accelerated and networked machines to quickly simulate circuits on many qubits.

QuEST has a simple interface, independent of its run environment (on CPUs, GPUs or over networks),
```C
hadamard(qubits, 0);

controlledNot(qubits, 0, 1);

rotateY(qubits, 0, .1);
```
though is flexible
```C
Vector v;
v.x = 1; v.y = 0; v.z = 0;
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

## Getting started

QuEST is contained entirely in the files in the `QuEST/` folder. To use QuEST, copy this folder to your computer and include `QuEST.h` in your `C` or `C++` code, and compile using the [makefile](makefile). We include [submission scripts](examples/submissionScripts/) for using QuEST with SLURM and PBS. See the [tutorial](/examples/README.md) for an introduction.

Explicit instructions to download and run QuEST from the command line can be found at [quest.qtechtheory.org/download](https://quest.qtechtheory.org/download/).

## API Documentation

View the API [here](https://aniabrown.github.io/QuEST/QuEST_8h.html), and check compatible compiler versions [here](tests/compilers/compatibility.md).

> For developers: To recreate the full documentation after making changes to the code, run doxygen doxyconf in the root directory. This will generate documentation in Doxygen_doc/html, and can be accessed through index.html in that folder. 

## Acknowledgements

QuEST uses the [mt19937ar](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) Mersenne Twister algorithm for random number generation, under the BSD licence. 

## Licence

QuEST is released under a [MIT Licence](https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt)



