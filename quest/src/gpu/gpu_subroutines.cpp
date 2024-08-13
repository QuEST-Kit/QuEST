/** @file
 * CUDA GPU-accelerated definitions of the subroutines called by
 * accelerator.cpp. This file is always compiled, even when GPU
 * acceleration is disabled and when parsed by a non-CUDA compiler,
 * and so uses precompiler guards to disable CUDA-only code.
 * 
 * This file contains host definitions and associated memory and
 * thread management, and invokes custom kernels defined in
 * kernels.hpp which is never parsed by non-CUDA compilers. This
 * file also invokes Thrust and cuQuantum routines, defined in
 * gpu_thrust.hpp and gpu_cuquantum.hpp respectively, which are
 * also never parsed by non-CUDA compilers.
 * 
 * Note that some custom kernels are templated in order to apply
 * compile-time optimisations like automatic loop unrolling. So
 * too are their calling host definitions in this file, which are
 * called by accelerator.cpp which chooses the template parameter.
 * This unnecessarily duplicates other parts of the host functions
 * responsible for dispatching to Thrust or cuQuantum, bloating the
 * compiled binary size; we accept this design wart over having
 * non-templated host functions because this requires duplicating
 * the template-dispatch logic (which would then also be defined in
 * cpu_subroutines.cpp) and moving it out of the aptly-named
 * accelerator.cpp file.
 *
 * When compiling for AMD GPUs, the CUDA symbols invoked herein are
 * mapped to HIP symbols by cuda_to_hip.h 
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/accelerator.hpp"

#if COMPILE_CUDA
    #include "quest/src/gpu/gpu_types.hpp"
    #include "quest/src/gpu/gpu_kernels.hpp"
    #include "quest/src/gpu/gpu_thrust.hpp"
#endif

#if COMPILE_CUQUANTUM
    #include "quest/src/gpu/gpu_cuquantum.hpp"
#endif

#include <vector>

using std::vector;



/*
 * BUFFER PACKING
 *
 * which are templated and require explicit instantiation below
 */


template <int NumCtrls>
void gpu_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {
    
    // there is no cuQuantum facility for packing

#if COMPILE_CUDA


    // TODO:
    //      I think we have to device-malloc space for ctrls.
    //      or like, copy it into a single compile-time array



    auto [qubits, stateMask] = getSortedQubitsAndMask(ctrls, ctrlStates, {}, {});
    qindex numThreads = qureg.numAmpsPerNode / powerOf2(NumCtrls);

    // cast qcomp to cu_qcomp
    qindex numBlocks = getNumBlocks(numThreads);
    cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);
    cu_qcomp* buff = toCuQcomps(qureg.gpuCommBuffer);

    kernel_statevec_packAmpsIntoBuffer 
        <NumCtrls> 
        <<<numBlocks, NUM_THREADS_PER_BLOCK>>> 
        (amps, buff, qubits.data(), stateMask, numThreads);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_packAmpsIntoBuffer, (Qureg, vector<int>, vector<int>) )






/*
 * TODO SECTION
 */


template <int NumCtrls>
void gpu_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
#if COMPILE_CUQUANTUM

    // ignore ctrlFlag
    cuquantum_statevec_anyCtrlAnyTargDenseMatrix_subA(qureg, ctrls, ctrlStates, {targ}, unpackMatrixToCuQcomps(matr).data());

#elif COMPILE_CUDA

    // prepare parameters needed for optimal indexing logic with the given ctrls

    qindex numBlocks = getNumBlocks(params.numInds);
    cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);

    // there are sufficiently few matrix elements (<256 bytes) to pass directly to kernel
    cu_qcomp m00, m01, m10, m11;
    unpackMatrixToCuQcomps(matr, m00,m01,m10,m11);

    // use ctrlFlag to dispatch to optimised kernel
    kernel_statevec_anyCtrlOneTargDenseMatr_subA <ctrlFlag> <<<numBlocks,NUM_THREADS_PER_BLOCK>>> (amps, params, targ, m00,m01,m10,m11);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlOneTargDenseMatr_subA, (Qureg, vector<int>, vector<int>, int, CompMatr1) )



template <int NumCtrls>
void gpu_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {
    
    // there is no cuQuantum facility for subB

#if COMPILE_CUDA

    // prepare parameters needed for optimal indexing logic with the given ctrls
    CtrlIndParams params = getParamsInformingIndsWhereCtrlsAreActive(qureg.numAmpsPerNode, ctrls, ctrlStates);

    // cast qcomp to cu_qcomp
    qindex numBlocks = getNumBlocks(params.numInds);
    cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);
    cu_qcomp* buff = toCuQcomps(qureg.gpuCommBuffer);
    cu_qcomp f0 = toCuQcomp(fac0);
    cu_qcomp f1 = toCuQcomp(fac1);

    // use ctrlFlag to dispatch to optimised kernel
    kernel_statevec_anyCtrlOneTargDenseMatr_subB <ctrlFlag> <<<numBlocks,NUM_THREADS_PER_BLOCK>>> (amps, buff, params, f0, f1);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}


INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS( void, gpu_statevec_anyCtrlOneTargDenseMatr_subB, (Qureg, vector<int>, vector<int>, qcomp, qcomp) )







template <int NumCtrls, int NumTargs>
void gpu_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {
#if COMPILE_CUQUANTUM

    cuquantum_statevec_anyCtrlAnyTargDiagMatr_subA(qureg, ctrls, ctrlStates, targs, toCuQcomps(matr.elems));

#elif COMPILE_CUDA

    qindex numThreads = qureg.numAmpsPerNode / powerOf2(NumCtrls);
    qindex numBlocks = getNumBlocks(numThreads);
    cu_qcomp* amps = toCuQcomps(qureg.gpuAmps);
    cu_qcomp d0 = toCuQcomp(matr.elems[0]);
    cu_qcomp d1 = toCuQcomp(matr.elems[1]);

    // TODO: fix kernel args etc

    // use ctrlFlag to dispatch to optimised kernel
    kernel_statevec_anyCtrlOneTargDiagMatr_sub <ctrlFlag> <<<numBlocks, NUM_THREADS_PER_BLOCK>>> (amps, params, targ, d0, d1);

#else
    error_gpuSimButGpuNotCompiled();
#endif
}

INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( void, gpu_statevec_anyCtrlAnyTargDiagMatr_sub, (Qureg, vector<int>, vector<int>, vector<int>, DiagMatr) )





