# QuEST

The **Quantum Exact Simulation Toolkit** is a high performance simulator of universal quantum circuits. QuEST is written in C, hybridises OpenMP and MPI, and can run on a GPU. Needing only compilation, QuEST is easy to run both on laptops and supercomputers, where it can take advantage of multicore and networked machines to quickly simulate circuits on many qubits.

> The GPU version of QuEST is currently located in another repository [here](https://github.com/aniabrown/QuEST_GPU).

## Getting started

QuEST is contained entirely in the `.c` and `.h` files in the `QuEST/` folder. To use QuEST, copy these files to your computer and include `qubits.h` in your C code. We include make files for compiling QuEST, and submission scripts for using QuEST with SLURM and PBS. See [examples/tutorial.md](/examples/tutorial.md) for an introduction.

## API Documentation

View the API [here](https://aniabrown.github.io/QuEST/qubits_8h.html), and the full documentation at https://aniabrown.github.io/QuEST/
