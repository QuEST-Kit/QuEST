# 🔖🤖  Automated examples

<!--
  CI-triggered examples
  (this comment must be under the title for valid doxygen rendering)
  
  @author Tyson Jones
-->

This folder contains examples which are compiled and executed by QuEST's `compile` Github Action, on [free runners](https://docs.github.com/en/actions/using-github-hosted-runners/using-github-hosted-runners/about-github-hosted-runners).
Contributors can include rudimentary testing or demo code in this folder (as standalone `.c` or `.cpp` files)
in their pull requests and preview its output on Github, where it is run with every combination of operating system, compiler and precision. 
This is intended for debugging and/or logging purposes, rather than for presenting example codes for users.

> [!IMPORTANT] 
> Beware that the configurations include use of multithreading, distribution and GPU-acceleration (as `OMP`, `MPI` and `CUDA` respectively) though this reflects only their _compilation_. All files are executed serially and locally on the CPU.

The output can be viewed at the `compile.yml`
[workflow tab](https://github.com/QuEST-Kit/QuEST/actions/workflows/compile.yml), clicking on the PR name, selecting a configuration (such as `Windows [1] OMP`) and expanding the `run automated examples` section. All files within `automated/` will be run in-turn and their outputs sequentially presented.

> [!IMPORTANT] 
> Since these examples are executed at every invocation of the CI, they should _not_ perform intensive, long computations. Such examples are better suited for the [`extended/`](../extended/) folder.

> [!TIP] 
> Files with extensions `.c` or `.cpp` will be respectively compiled in `C11` and `C++17`. An example which uses a facility sensitive to the language can be duplicated as both `file.c` and `file.cpp` so that both versions are tested.
