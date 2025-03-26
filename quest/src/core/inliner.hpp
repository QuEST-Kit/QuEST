/** @file
 * This file defines INLINE to aggressively force inlining
 * of crucial functions (typically bitwise functions which
 * get invoked exponentially-many times in hot loops) in a
 * compiler agnostic way. Inlined functions may be called
 * invoked by OpenMP threads, or CUDA kernels. They never
 * cause symbol duplication when multiply-imported.
 * 
 * @author Tyson Jones
 * @author Oliver Thomson Brown (patched HIP guards)
 */

#ifndef INLINER_HPP
#define INLINER_HPP



/*
 * We must choose the right INLINE keyword for the user's compiler.
 * Note the below logic means all INLINE functions through the entire
 * QuEST source will declared __device__ when the user opts to compile
 * everything with NVCC. It is ergo important that all INLINE functions
 * are CUDA-kernel-compatible (e.g. don't make use of std::vector).
 */


#if defined(__NVCC__)

    // CUDA compilers define '__forceinline__'. We declare
    // all inlined functions as __device__ so they can be
    // invoked by CUDA kernels. Note this means that if the
    // user opts to use nvcc to compile the entire QuEST
    // source, then every INLINE function will be twice-
    // compiled; as a host and device function. That's fine,
    // though it constrains us to never use non-kernel
    // compatible code inside an INLINE function.
    #define INLINE __forceinline__ __device__ __host__

#elif defined(__HIP__)

    // HIP has the same syntax as GNU below, but we must also
    // declare that all INLINE functions can be invoked by kernels
    #define INLINE inline __attribute__((always_inline)) __device__ __host__

#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)

    // MSVC and Intel compilers define '__forceinline'
    #define INLINE __forceinline

#elif defined(__GNUC__)
    #define INLINE inline __attribute__((always_inline))

#else

    // #warning command safe in non-MSVC compiler
    #warning "Could not ascertain compiler type in order to choose the correct inline attribute. Assuming GNU and proceeding..."
    #define INLINE inline __attribute__((always_inline))

#endif



#endif // INLINER_HPP