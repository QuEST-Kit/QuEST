# QuEST

## Versions

QuEST is currently in prerelease and may be unstable.  

Latest version: [0.10.0](https://github.com/aniabrown/QuEST/releases/tag/v0.10.0) 

Please report errors or feedback to anna.brown@oerc.ox.ac.uk 

## Introduction

The **Quantum Exact Simulation Toolkit** is a high performance simulator of universal quantum circuits. QuEST is written in C, hybridises OpenMP and MPI, and can run on a GPU. Needing only compilation, QuEST is easy to run both on laptops and supercomputers, where it can take advantage of multicore and networked machines to quickly simulate circuits on many qubits.

> The GPU version of QuEST is currently located in another repository [here](https://github.com/aniabrown/QuEST_GPU), though the introduction in [examples/tutorial.md](/examples/tutorial.md) is still relevant.

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

QuEST is contained entirely in the `.c` and `.h` files in the `QuEST/` folder. To use QuEST, copy these files to your computer and include `qubits.h` in your C code. We include make files for compiling QuEST, and submission scripts for using QuEST with SLURM and PBS. See [examples/tutorial.md](/examples/tutorial.md) for an introduction. Clone or download this entire repository to include all examples as well as tests and documentation. 

## API Documentation

View the API [here](https://aniabrown.github.io/QuEST/qubits_8h.html), and the full documentation at https://aniabrown.github.io/QuEST/

## Licence

QuEST is released under a [MIT Licence](https://github.com/aniabrown/QuEST/blob/master/LICENCE.txt)



