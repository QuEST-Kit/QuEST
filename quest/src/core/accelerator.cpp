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
#include "quest/src/core/bitwise.hpp"
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
 * a runtime variable as a template parameter. Note an awkward use of decltype()
 * is to workaround a GCC <12 bug with implicitly-typed vector initialisations.
 */


#if (MAX_OPTIMISED_NUM_CTRLS != 5) || (MAX_OPTIMISED_NUM_TARGS != 5)
    #error "The number of optimised, templated functions was inconsistent between accelerator's source and header."
#endif


#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS(f, numctrls) \
    (vector <decltype(&f<0>)> {&f<0>, &f<1>, &f<2>, &f<3>, &f<4>, &f<5>, &f<-1>}) \
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS - 1)]


#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(f, numctrls, numtargs) \
    (vector <ARR(f)> { \
        ARR(f) {&f<0,0>,  &f<0,1>,  &f<0,2>,  &f<0,3>,  &f<0,4>,  &f<0,5>,  &f<0,-1>}, \
        ARR(f) {&f<1,0>,  &f<1,1>,  &f<1,2>,  &f<1,3>,  &f<1,4>,  &f<1,5>,  &f<1,-1>}, \
        ARR(f) {&f<2,0>,  &f<2,1>,  &f<2,2>,  &f<2,3>,  &f<2,4>,  &f<2,5>,  &f<2,-1>}, \
        ARR(f) {&f<3,0>,  &f<3,1>,  &f<3,2>,  &f<3,3>,  &f<3,4>,  &f<3,5>,  &f<3,-1>}, \
        ARR(f) {&f<4,0>,  &f<4,1>,  &f<4,2>,  &f<4,3>,  &f<4,4>,  &f<4,5>,  &f<4,-1>}, \
        ARR(f) {&f<5,0>,  &f<5,1>,  &f<5,2>,  &f<5,3>,  &f<5,4>,  &f<5,5>,  &f<5,-1>}, \
        ARR(f) {&f<-1,0>, &f<-1,1>, &f<-1,2>, &f<-1,3>, &f<-1,4>, &f<-1,5>, &f<-1,-1>}}) \
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS - 1)] \
    [std::min((int) numtargs, MAX_OPTIMISED_NUM_TARGS - 1)]


#define ARR(f) vector<decltype(&f<0,0>)>



/*
 * COMMUNICATION BUFFER PACKING
 */


qindex accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubits, vector<int> qubitStates) {

    // we can never pack and swap buffers when there are no constrained qubit states, because we'd 
    // then fill the entire buffer andhave no room to receive the other node's buffer; caller would 
    // instead send amps straight to buffer
    if (qubitStates.empty())
        error_noCtrlsGivenToBufferPacker();

    // packing treats qubits as if they were ctrl qubits
    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_statevec_packAmpsIntoBuffer, qubits.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_statevec_packAmpsIntoBuffer, qubits.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, qubits, qubitStates);

    // return the number of packed amps, for caller convenience
    return qureg.numAmpsPerNode / powerOf2(qubits.size());
}



/*
 * SWAPS
 */


void accel_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_statevec_anyCtrlSwap_subA, ctrls.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_statevec_anyCtrlSwap_subA, ctrls.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates, targ1, targ2);
}

void accel_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_statevec_anyCtrlSwap_subB, ctrls.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_statevec_anyCtrlSwap_subB, ctrls.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates);
}

void accel_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_statevec_anyCtrlSwap_subC, ctrls.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_statevec_anyCtrlSwap_subC, ctrls.size() );
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates, targ, targState);
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

void accel_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_statevec_anyCtrlAnyTargDenseMatr_sub, ctrls.size(), targs.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_statevec_anyCtrlAnyTargDenseMatr_sub, ctrls.size(), targs.size() );
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

