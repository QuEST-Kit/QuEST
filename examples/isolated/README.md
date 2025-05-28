# 🔖🔍  Isolated examples

<!--
  Isolated examples
  (this comment must be under the title for valid doxygen rendering)
  
  @author Tyson Jones
-->

This folder contains _isolated_ examples which demonstrate the possible uses of a specific QuEST function or facility. They mostly serve to illustrate the various ways that a single task (such as initialising a `KrausMap`) can be performed with the QuEST API. We provide both `C` and `C++` versions (compiled against standards `C11` and `C++17`) of each example to highlight the interfaces available to the respective languages.

Every `.c` and `.cpp` in this folder is compiled _and_ executed by the [compile](https://github.com/QuEST-Kit/QuEST/actions/workflows/compile.yml) Github Action, to automatically check for compilation and link-time issues.
