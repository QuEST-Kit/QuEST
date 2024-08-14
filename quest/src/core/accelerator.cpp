/** @file
 * Internal functions for choosing which accelerator backend
 * (CPU or GPU) to dispatch to, and which preconditions the 
 * qubit indices satisfy in order (informing which compile-time
 * optimisations to use), to effect local simulation subroutines 
 * upon Quregs.
 * 
 * These routines are called by localiser.cpp and are embarrassingly 
 * parallel, so are always called before/after any necessary
 * communication has happened. The data they need must already be
 * localised into the appropriate memory (RAM or VRAM) and location
 * (qureg's amplitudes or buffer space).
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/accelerator.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/cpu/cpu_subroutines.hpp"
#include "quest/src/gpu/gpu_subroutines.hpp"

#include <vector>
#include <algorithm>

using std::vector;



/*
 * MACROS
 *
 * which automate the choosing of the appropriate backend template function,
 * optimised for the given configuration of qubit indices, for example through
 * automatic unrolling of loops with bounds known at compile-time. When the 
 * number of controls or targets exceeds that which have optimised compilations, 
 * we fall back to using a generic implementation, indicated by <-1>. In essence,
 * these macros simply call func<ctrls.size()> albeit without illegally passing
 * a runtime variable as a template parameter.
 */


#if (MAX_OPTIMISED_NUM_CTRLS != 5) || (MAX_OPTIMISED_NUM_TARGS != 5)
    #error "The number of optimised, templated functions was inconsistent between accelerator's source and header."
#endif


#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS(func, numctrls) \
    (vector{func<0>, func<1>, func<2>, func<3>, func<4>, func<5>, func<-1>}) \
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS - 1)]


#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(func, numctrls, numtargs) \
    (vector{ \
        vector{func<0,0>,  func<0,1>,  func<0,2>,  func<0,3>,  func<0,4>,  func<0,5>,  func<0,-1>}, \
        vector{func<1,0>,  func<1,1>,  func<1,2>,  func<1,3>,  func<1,4>,  func<1,5>,  func<1,-1>}, \
        vector{func<2,0>,  func<2,1>,  func<2,2>,  func<2,3>,  func<2,4>,  func<2,5>,  func<2,-1>}, \
        vector{func<3,0>,  func<3,1>,  func<3,2>,  func<3,3>,  func<3,4>,  func<3,5>,  func<3,-1>}, \
        vector{func<4,0>,  func<4,1>,  func<4,2>,  func<4,3>,  func<4,4>,  func<4,5>,  func<4,-1>}, \
        vector{func<5,0>,  func<5,1>,  func<5,2>,  func<5,3>,  func<5,4>,  func<5,5>,  func<5,-1>}, \
        vector{func<-1,0>, func<-1,1>, func<-1,2>, func<-1,3>, func<-1,4>, func<-1,5>, func<-1,-1>}}) \
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS - 1)] \
    [std::min((int) numtargs, MAX_OPTIMISED_NUM_TARGS - 1)]




/*
 * COMMUNICATION BUFFER PACKING
 */


void accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    // we can never pack and swap buffers when there are no ctrl qubits, because we'd fill the entire buffer
    // andhave no room to receive the other node's buffer; caller would instead send amps straight to buffer
    if (ctrls.empty())
        error_noCtrlsGivenToBufferPacker();

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_statevec_packAmpsIntoBuffer, ctrls.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_statevec_packAmpsIntoBuffer, ctrls.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates);
}



/*
 * DENSE MATRIX
 */


void accel_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_statevec_anyCtrlOneTargDenseMatr_subA, ctrls.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_statevec_anyCtrlOneTargDenseMatr_subA, ctrls.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates, targ, matr);
}

void accel_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_statevec_anyCtrlOneTargDenseMatr_subB, ctrls.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_statevec_anyCtrlOneTargDenseMatr_subB, ctrls.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates, fac0, fac1);
}

void accel_statevec_anyCtrlAnyTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_statevec_anyCtrlAnyTargDenseMatr_subA, ctrls.size(), targs.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_statevec_anyCtrlAnyTargDenseMatr_subA, ctrls.size(), targs.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates, targs, matr);
}



/*
 * DIAGONAL MATRIX
 */


void accel_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_statevec_anyCtrlAnyTargDiagMatr_sub, ctrls.size(), targs.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_statevec_anyCtrlAnyTargDiagMatr_sub, ctrls.size(), targs.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates, targs, matr);
}