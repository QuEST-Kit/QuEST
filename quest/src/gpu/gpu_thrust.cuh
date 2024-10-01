/** @file
 * Subroutines which invoke Thrust. This file is only ever included
 * when COMPILE_CUDA=1 so it can safely invoke CUDA
 * signatures without guards.
 */

#ifndef GPU_THRUST_HPP
#define GPU_THRUST_HPP

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_thrust.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#include "quest/include/modes.h"
#include "quest/include/types.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/gpu/gpu_types.cuh"



void thrust_setElemsToConjugate(cu_qcomp* matrElemsPtr, qindex matrElemsLen) {

    // TODO
    error_functionNotImplemented();
}



#endif // GPU_THRUST_HPP