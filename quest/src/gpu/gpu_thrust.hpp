/** @file
 * Subroutines which invoke Thrust. This file is only ever included
 * when ENABLE_GPU_ACCELERATION=1 so it can safely invoke CUDA
 * signatures without guards.
 */

#ifndef GPU_THRUST_HPP
#define GPU_THRUST_HPP

#if ! ENABLE_GPU_ACCELERATION
    #error "A file being compiled somehow included gpu_thrust.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif




#endif // GPU_THRUST_HPP