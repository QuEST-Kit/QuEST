/** @file
 * CUDA GPU-accelerated definitions of the subroutines called by
 * accelerator.cpp. This file contains the host definitions and
 * associated memory and thread management, and invocations of
 * Thrust and cuQuantum subroutines (if the latter is enabled).
 * All custom kernels are defined in kernels.hpp, which is never
 * parsed by non-CUDA compilers and ignored when not compiling GPU.
 * When compiling for AMD GPUs, the CUDA symbols invoked herein are
 * mapped to HIP symbols by cuda_to_hip.h
 */