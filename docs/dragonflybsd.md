# Compiling and Running QuEST on DragonFlyBSD

## Step-by-step tutorial for QuEST 4.x, GCC 13, and OpenMP

This tutorial documents the installation path successfully tested with
**QuEST 4.3.0 on DragonFlyBSD**. It includes the DragonFlyBSD-specific
fixes needed for the C example linker and the GCC 13 C++ runtime.

## 1. Install the required build tools

Update the package database and install Git, CMake, GNU Make, GCC 13,
and `pkgconf`:

``` csh
sudo pkg update
sudo pkg install git cmake gmake gcc13 pkgconf
rehash
```

Verify the tools:

``` csh
git --version
cmake --version
gmake --version
gcc13 --version
g++13 --version
```

## 2. Download QuEST

Clone the official QuEST repository:

``` csh
cd ~
git clone https://github.com/QuEST-Kit/QuEST.git
cd QuEST
```

For reproducibility, record the exact revision:

``` csh
git rev-parse HEAD
git describe --tags --always
```

## 3. Create a clean build directory

``` csh
cd ~/QuEST
rm -rf build
mkdir build
cd build
```

## 4. Configure QuEST for GCC 13 and OpenMP

Use QuEST's current CMake option names:

``` csh
cmake .. \
 -DCMAKE_C_COMPILER=/usr/local/bin/gcc13 \
 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++13 \
 -DQUEST_ENABLE_OMP=ON \
 -DQUEST_ENABLE_MPI=OFF \
 -DQUEST_ENABLE_CUDA=OFF \
 -DQUEST_BUILD_EXAMPLES=ON \
 -DCMAKE_BUILD_TYPE=Release
```

A successful configuration should report multithreading as ON and
MPI/distribution and CUDA as OFF.

You may see:

``` text
libnuma not found, QuEST will not be aware of numa locality
```

This was only a warning in the tested DragonFlyBSD configuration and did
not prevent QuEST from compiling or running.

## 5. Determine the CPU count

DragonFlyBSD commonly uses `csh`/`tcsh`. Check the number of visible
CPUs:

``` csh
sysctl -n hw.ncpu
```

If it reports `2`:

``` csh
gmake -j2
```

If it reports `8`:

``` csh
gmake -j8
```

A csh-compatible automatic form is:

``` csh
gmake -j`sysctl -n hw.ncpu`
```

Avoid Bash-style `$(...)` command substitution in `csh`/`tcsh`.

## 6. Compile QuEST

For the tested 2-vCPU system:

``` csh
gmake -j2
```

The QuEST library successfully reached:

``` text
[94%] Linking CXX shared library libQuEST.so
[94%] Built target QuEST
```

## 7. Handle the DragonFlyBSD C-example linker issue

On the tested system, CMake attempted to link `min_example.c` with
`gcc13`, but `libQuEST.so` contains C++ code and requires the C++
runtime.

The build failed with a GLIBCXX-related unresolved symbol even though
the QuEST library itself had already compiled.

From `~/QuEST/build`, manually link the generated C object with `g++13`:

``` csh
/usr/local/bin/g++13 \
 CMakeFiles/min_example.dir/examples/tutorials/min_example.c.o \
 -o min_example \
 ./libQuEST.so.4.3.0
```

Confirm that the executable exists:

``` csh
ls -l min_example
```

## 8. Verify the GCC 13 C++ runtime

Find GCC 13's `libstdc++`:

``` csh
/usr/local/bin/g++13 -print-file-name=libstdc++.so.6
```

Verify that the tested library contains the required symbol version:

``` csh
strings /usr/local/lib/gcc13/libstdc++.so.6 | grep GLIBCXX_3.4.32
```

The successful test returned:

``` text
GLIBCXX_3.4.32
```

## 9. Set the runtime library path

The tested DragonFlyBSD environment initially loaded GCC 11's
`libstdc++`, even though QuEST had been compiled with GCC 13.

Set the runtime path to GCC 13 and the QuEST build directory:

``` csh
setenv LD_LIBRARY_PATH /usr/local/lib/gcc13:/home/user/QuEST/build
```

Replace `/home/user` if your home directory is different.

## 10. Verify dynamic library resolution

Run:

``` csh
ldd ./min_example
```

The important entries should resolve similarly to:

``` text
libQuEST.so.4  => /home/user/QuEST/build/libQuEST.so.4
libstdc++.so.6 => /usr/local/lib/gcc13/libstdc++.so.6
libgomp.so.1   => /usr/local/lib/gcc13/libgomp.so.1
```

`libgomp` is the GNU OpenMP runtime used for CPU multithreading.

## 11. Run QuEST

``` csh
./min_example
```

In the successful DragonFlyBSD test, QuEST created a 20-qubit
state-vector simulation with 1,048,576 amplitudes using 16 MiB.

The output ended with:

``` text
Total probability: 1
```

This confirms successful execution and a normalized quantum state.

## 12. Test OpenMP multithreading

One thread:

``` csh
setenv OMP_NUM_THREADS 1
./min_example
```

Two threads:

``` csh
setenv OMP_NUM_THREADS 2
./min_example
```

Check the CPU count before choosing larger thread counts:

``` csh
sysctl -n hw.ncpu
```

## Troubleshooting summary

  -----------------------------------------------------------------------
  Symptom                             Action
  ----------------------------------- -----------------------------------
  `Illegal variable name` with        Use an explicit `-jN` or csh
  `gmake -j$(...)`                    backticks.

  `libnuma not found`                 Warning only in the tested CPU
                                      build.

  Undefined C++ / GLIBCXX symbol      Link the generated object with
  while linking `min_example`         `g++13`.

  `GLIBCXX_3.4.32 not found` at       Ensure `/usr/local/lib/gcc13` is
  runtime                             first in `LD_LIBRARY_PATH`.

  `libQuEST.so` cannot be found       Add `~/QuEST/build` to
                                      `LD_LIBRARY_PATH` and verify with
                                      `ldd`.
  -----------------------------------------------------------------------

## Quick command reference

``` csh
sudo pkg update
sudo pkg install git cmake gmake gcc13 pkgconf
rehash

cd ~
git clone https://github.com/QuEST-Kit/QuEST.git
cd QuEST

rm -rf build
mkdir build
cd build

cmake .. \
 -DCMAKE_C_COMPILER=/usr/local/bin/gcc13 \
 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++13 \
 -DQUEST_ENABLE_OMP=ON \
 -DQUEST_ENABLE_MPI=OFF \
 -DQUEST_ENABLE_CUDA=OFF \
 -DQUEST_BUILD_EXAMPLES=ON \
 -DCMAKE_BUILD_TYPE=Release

sysctl -n hw.ncpu
gmake -j2

/usr/local/bin/g++13 \
 CMakeFiles/min_example.dir/examples/tutorials/min_example.c.o \
 -o min_example \
 ./libQuEST.so.4.3.0

setenv LD_LIBRARY_PATH /usr/local/lib/gcc13:/home/user/QuEST/build
ldd ./min_example
./min_example
```

## Result

**Successfully tested configuration:** DragonFlyBSD + GCC/G++ 13.3 +
QuEST 4.3.0 + OpenMP/libgomp + CPU state-vector simulation.
