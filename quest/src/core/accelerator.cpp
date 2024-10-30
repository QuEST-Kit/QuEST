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
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"
#include "quest/src/cpu/cpu_subroutines.hpp"
#include "quest/src/gpu/gpu_subroutines.hpp"

#include <vector>
#include <algorithm>

using std::vector;
using std::min;



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
 * 
 * BEWARE that these macros are single-line expressions, so they can be used in
 * braceless if/else or ternary operators - but stay vigilant!
 */


#define GET_FUNC_OPTIMISED_FOR_BOOL(funcname, value) \
    ((value)? funcname<true> : funcname<false>)


#define GET_FUNC_OPTIMISED_FOR_TWO_BOOLS(funcname, b1, b2) \
    ((b1)? \
        ((b2)? funcname<true, true> : funcname<true, false>) : \
        ((b2)? funcname<false,true> : funcname<false,false>))


#if (MAX_OPTIMISED_NUM_CTRLS != 5) || (MAX_OPTIMISED_NUM_TARGS != 5)
    #error "The number of optimised, templated QuEST functions was inconsistent between accelerator's source and header."
#endif


#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS(f, numctrls) \
    (vector <decltype(&f<0>)> {&f<0>, &f<1>, &f<2>, &f<3>, &f<4>, &f<5>, &f<-1>}) \
    [min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)]

#define GET_FUNC_OPTIMISED_FOR_NUM_TARGS(f, numtargs) \
    (vector <decltype(&f<0>)> {&f<0>, &f<1>, &f<2>, &f<3>, &f<4>, &f<5>, &f<-1>}) \
    [min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(f, numctrls, numtargs) \
    (vector <ARR(f)> { \
        ARR(f) {&f<0,0>,  &f<0,1>,  &f<0,2>,  &f<0,3>,  &f<0,4>,  &f<0,5>,  &f<0,-1>}, \
        ARR(f) {&f<1,0>,  &f<1,1>,  &f<1,2>,  &f<1,3>,  &f<1,4>,  &f<1,5>,  &f<1,-1>}, \
        ARR(f) {&f<2,0>,  &f<2,1>,  &f<2,2>,  &f<2,3>,  &f<2,4>,  &f<2,5>,  &f<2,-1>}, \
        ARR(f) {&f<3,0>,  &f<3,1>,  &f<3,2>,  &f<3,3>,  &f<3,4>,  &f<3,5>,  &f<3,-1>}, \
        ARR(f) {&f<4,0>,  &f<4,1>,  &f<4,2>,  &f<4,3>,  &f<4,4>,  &f<4,5>,  &f<4,-1>}, \
        ARR(f) {&f<5,0>,  &f<5,1>,  &f<5,2>,  &f<5,3>,  &f<5,4>,  &f<5,5>,  &f<5,-1>}, \
        ARR(f) {&f<-1,0>, &f<-1,1>, &f<-1,2>, &f<-1,3>, &f<-1,4>, &f<-1,5>, &f<-1,-1>}}) \
    [min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)] \
    [min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

#define ARR(f) vector<decltype(&f<0,0>)>


#define GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS(funcsuffix, qureg, numctrls) \
    ((qureg.isGpuAccelerated)? \
        GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( gpu_##funcsuffix, numctrls ) : \
        GET_FUNC_OPTIMISED_FOR_NUM_CTRLS( cpu_##funcsuffix, numctrls ))

#define GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS(funcsuffix, qureg, numtargs) \
    ((qureg.isGpuAccelerated)? \
        GET_FUNC_OPTIMISED_FOR_NUM_TARGS( gpu_##funcsuffix, numtargs ) : \
        GET_FUNC_OPTIMISED_FOR_NUM_TARGS( cpu_##funcsuffix, numtargs ))

#define GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(funcsuffix, qureg, numctrls, numtargs) \
    ((qureg.isGpuAccelerated)? \
        GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs ) : \
        GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs ))


#define GET_BOOLEAN_FUNC_OPTIMISED_FOR_NUM_TARGS(f, numtargs, boolval) \
    ((boolval)? \
        GET_BOOLEAN_FUNC_LIST( f, numtargs, true ) : \
        GET_BOOLEAN_FUNC_LIST( f, numtargs, false ))

#define GET_BOOLEAN_FUNC_LIST(f, numtargs, b) \
    (vector <decltype(&f<0,b>)> {&f<0,b>, &f<1,b>, &f<2,b>, &f<3,b>, &f<4,b>, &f<5,b>, &f<-1,b>}) \
    [min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

#define GET_CPU_OR_GPU_BOOLEAN_FUNC_OPTIMISED_FOR_NUM_TARGS(funcsuffix, qureg, numtargs, boolval) \
    ((qureg.isGpuAccelerated)? \
        GET_BOOLEAN_FUNC_OPTIMISED_FOR_NUM_TARGS( gpu_##funcsuffix, numtargs, boolval ) : \
        GET_BOOLEAN_FUNC_OPTIMISED_FOR_NUM_TARGS( cpu_##funcsuffix, numtargs, boolval ))


#define GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(f, numctrls, numtargs, c) \
    (vector <CONJ_ARR(f)> { \
        CONJ_ARR(f) {&f<0,0,c>,  &f<0,1,c>,  &f<0,2,c>,  &f<0,3,c>,  &f<0,4,c>,  &f<0,5,c>,  &f<0,-1,c>}, \
        CONJ_ARR(f) {&f<1,0,c>,  &f<1,1,c>,  &f<1,2,c>,  &f<1,3,c>,  &f<1,4,c>,  &f<1,5,c>,  &f<1,-1,c>}, \
        CONJ_ARR(f) {&f<2,0,c>,  &f<2,1,c>,  &f<2,2,c>,  &f<2,3,c>,  &f<2,4,c>,  &f<2,5,c>,  &f<2,-1,c>}, \
        CONJ_ARR(f) {&f<3,0,c>,  &f<3,1,c>,  &f<3,2,c>,  &f<3,3,c>,  &f<3,4,c>,  &f<3,5,c>,  &f<3,-1,c>}, \
        CONJ_ARR(f) {&f<4,0,c>,  &f<4,1,c>,  &f<4,2,c>,  &f<4,3,c>,  &f<4,4,c>,  &f<4,5,c>,  &f<4,-1,c>}, \
        CONJ_ARR(f) {&f<5,0,c>,  &f<5,1,c>,  &f<5,2,c>,  &f<5,3,c>,  &f<5,4,c>,  &f<5,5,c>,  &f<5,-1,c>}, \
        CONJ_ARR(f) {&f<-1,0,c>, &f<-1,1,c>, &f<-1,2,c>, &f<-1,3,c>, &f<-1,4,c>, &f<-1,5,c>, &f<-1,-1,c>}}) \
    [min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)] \
    [min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

#define CONJ_ARR(f) vector<decltype(&f<0,0,false>)>

#define GET_CPU_OR_GPU_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(funcsuffix, qureg, numctrls, numtargs, conj) \
    ((qureg.isGpuAccelerated)? \
        ((conj)? \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, true ) : \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, false ) ) : \
        ((conj)? \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, true ) : \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, false ) ) )


// TODO:
// This has gotten a bit ridiculous. Is there a way to use (likely)
// more abominable pre-processor mischief which negates the need
// to repeat the entire macro(s) when the number of templated
// parameters grows?


#define GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(f, numctrls, numtargs, c, h) \
    (vector <POWER_CONJ_ARR(f)> { \
        POWER_CONJ_ARR(f) {&f<0,0,c,h>,  &f<0,1,c,h>,  &f<0,2,c,h>,  &f<0,3,c,h>,  &f<0,4,c,h>,  &f<0,5,c,h>,  &f<0,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<1,0,c,h>,  &f<1,1,c,h>,  &f<1,2,c,h>,  &f<1,3,c,h>,  &f<1,4,c,h>,  &f<1,5,c,h>,  &f<1,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<2,0,c,h>,  &f<2,1,c,h>,  &f<2,2,c,h>,  &f<2,3,c,h>,  &f<2,4,c,h>,  &f<2,5,c,h>,  &f<2,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<3,0,c,h>,  &f<3,1,c,h>,  &f<3,2,c,h>,  &f<3,3,c,h>,  &f<3,4,c,h>,  &f<3,5,c,h>,  &f<3,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<4,0,c,h>,  &f<4,1,c,h>,  &f<4,2,c,h>,  &f<4,3,c,h>,  &f<4,4,c,h>,  &f<4,5,c,h>,  &f<4,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<5,0,c,h>,  &f<5,1,c,h>,  &f<5,2,c,h>,  &f<5,3,c,h>,  &f<5,4,c,h>,  &f<5,5,c,h>,  &f<5,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<-1,0,c,h>, &f<-1,1,c,h>, &f<-1,2,c,h>, &f<-1,3,c,h>, &f<-1,4,c,h>, &f<-1,5,c,h>, &f<-1,-1,c,h>}}) \
    [min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)] \
    [min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

#define POWER_CONJ_ARR(f) vector<decltype(&f<0,0,false,false>)>

#define GET_CPU_OR_GPU_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(funcsuffix, qureg, numctrls, numtargs, conj, haspower) \
    ((qureg.isGpuAccelerated)? \
        ((conj)? \
            ((haspower)? \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, true, true ) : \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, true, false ) ) : \
            ((haspower)? \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, false, true ) : \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, false, false ) ) ) : \
        ((conj)? \
            ((haspower)? \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, true, true ) : \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, true, false ) ) : \
            ((haspower)? \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, false, true ) : \
                GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, false, false ) ) ) )



/*
 * GETTERS 
 */


qcomp accel_statevec_getAmp_sub(Qureg qureg, qindex ind) {

    return (qureg.isGpuAccelerated)?
        gpu_statevec_getAmp_sub(qureg, ind):
        cpu_statevec_getAmp_sub(qureg, ind);
}



/*
 * COMMUNICATION BUFFER PACKING
 */


qindex accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubits, vector<int> qubitStates) {

    // we can never pack and swap buffers when there are no constrained qubit states, because we'd 
    // then fill the entire buffer andhave no room to receive the other node's buffer; caller would 
    // instead send amps straight to buffer
    if (qubitStates.empty())
        error_noCtrlsGivenToBufferPacker();

    // note qubits may incidentally be ctrls or targs; it doesn't matter
    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( statevec_packAmpsIntoBuffer, qureg, qubits.size() );
    
    // return the number of packed amps, for caller convenience
    return func(qureg, qubits, qubitStates);
}


qindex accel_statevec_packPairSummedAmpsIntoBuffer(Qureg qureg, int qubit1, int qubit2, int qubit3, int bit2) {

    return (qureg.isGpuAccelerated)?
        gpu_statevec_packPairSummedAmpsIntoBuffer(qureg, qubit1, qubit2, qubit3, bit2):
        cpu_statevec_packPairSummedAmpsIntoBuffer(qureg, qubit1, qubit2, qubit3, bit2);
}



/*
 * SWAPS
 */


void accel_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlSwap_subA, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, targ1, targ2);
}
void accel_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlSwap_subB, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates);
}
void accel_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlSwap_subC, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, targ, targState);
}



/*
 * DENSE MATRIX
 */


void accel_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlOneTargDenseMatr_subA, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, targ, matr);
}
void accel_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlOneTargDenseMatr_subB, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, fac0, fac1);
}


void accel_statevec_anyCtrlTwoTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, CompMatr2 matr) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlTwoTargDenseMatr_sub, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, targ1, targ2, matr);
}


void accel_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr, bool conj) {

    auto func = GET_CPU_OR_GPU_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( statevec_anyCtrlAnyTargDenseMatr_sub, qureg, ctrls.size(), targs.size(), conj );
    func(qureg, ctrls, ctrlStates, targs, matr);
}



/*
 * ANY-TARG DIAGONAL MATRIX
 */


void accel_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlOneTargDiagMatr_sub, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, targ, matr);
}


void accel_statevec_anyCtrlTwoTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, DiagMatr2 matr) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevec_anyCtrlTwoTargDiagMatr_sub, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, targ1, targ2, matr);
}


void accel_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, qcomp exponent, bool conj) {

    bool hasPower = exponent != qcomp(1, 0);

    auto func = GET_CPU_OR_GPU_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( statevec_anyCtrlAnyTargDiagMatr_sub, qureg, ctrls.size(), targs.size(), conj, hasPower );
    func(qureg, ctrls, ctrlStates, targs, matr, exponent);
}



/*
 * ALL-TARGS DIAGONAL MATRIX
 */


void accel_statevec_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent) {

    // qureg and matr are equal size and identically distributed...
    assert_quregAndFullStateDiagMatrAreBothOrNeitherDistrib(qureg, matr);
    
    // but they may have differing GPU deployments
    bool quregGPU = qureg.isGpuAccelerated;
    bool matrGPU = matr.isGpuAccelerated;

    // our chosen function must be correctly templated
    bool hasPower = exponent != qcomp(1, 0);
    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_BOOL( cpu_statevec_allTargDiagMatr_sub, hasPower );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_BOOL( gpu_statevec_allTargDiagMatr_sub, hasPower );

    // when deployments match, we trivially call the common backend
    if ( quregGPU &&  matrGPU) gpuFunc(qureg, matr, exponent);
    if (!quregGPU && !matrGPU) cpuFunc(qureg, matr, exponent);

    // deployments differing is a strange and expectedly rare scenario;
    // why use GPU-acceleration for a Qureg but not the equally-sized
    // matrix? We provide the below fallbacks for defensive design.

    // When GPU-accel differs, we fall-back to copying memory to RAM and
    // using the CPU backend. In theory, we could leverage exsting GPU 
    // memory of the Qureg's communication buffer (if it existed), but
    // this is an even rarer situation and is hacky. We could also 
    // create new, temporary GPU memory and graft it to the non-
    // accelerated object, but the new allocation would be the same
    // size as the objects and ergo be dangerously large.

    if (!quregGPU && matrGPU) {

        // copying matr GPU memory to CPU is unnecessary,
        // because it should never have diverged
        cpuFunc(qureg, matr, exponent);
    }

    if (quregGPU && !matrGPU) {
        gpu_copyGpuToCpu(qureg);
        cpuFunc(qureg, matr, exponent);
        gpu_copyCpuToGpu(qureg);
    }
}


void accel_densmatr_allTargDiagMatr_subA(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool multiplyOnly) {

    // matr is always local, qureg can be local or distributed...
    assert_fullStateDiagMatrIsLocal(matr);

    // and their GPU deployments can differ
    bool quregGPU = qureg.isGpuAccelerated;
    bool matrGPU = matr.isGpuAccelerated;

    // our chosen CPU or GPU dispatched function must be correctly templated
    bool hasPower = exponent != qcomp(1, 0);
    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( cpu_densmatr_allTargDiagMatr_sub, hasPower, multiplyOnly );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( cpu_densmatr_allTargDiagMatr_sub, hasPower, multiplyOnly );

    // when deployments match, we trivially call the common backend
    if ( quregGPU &&  matrGPU) gpuFunc(qureg, matr, exponent);
    if (!quregGPU && !matrGPU) cpuFunc(qureg, matr, exponent);

    // when only the matr is GPU-accelerated (which is strange, but
    // supported for defensive design), we must use CPU simulation. 
    // No need to copy memory; matr's CPU copy should be unchanged
    if (!quregGPU && matrGPU)
        cpuFunc(qureg, matr, exponent);

    // the most common scenario is that qureg (which is quadratically 
    // larger than matr) is GPU-accelerated, while matr is not. In that
    // case, we graft GPU memory onto matr and call gpuFunc(). If
    // qureg is distributed, we can re-use its existing GPU communication
    // buffer memory, otherwise we will have to allocate temporary memory; 
    // not a big deal given it is quadratically smaller than Qureg's memory
    if (quregGPU && !matrGPU) {

        // binding qureg's GPU communication buffer to matrix is safe,
        // even when subB() below (which itself grafts qureg's buffer to matr)
        // calls this function; that scenario never triggers condition (GPU
        // deployments will match) and instead calls the both-gpu function above.
        assert_quregGpuBufferIsNotGraftedToMatrix(qureg, matr);

        // spoof a GPU-accelerated matrix, grafting buffer or new memory
        // (we use a paranoid copy of matr, even though matr is already a 
        // mere copy of the user's matrix, in case this code changes to
        // accept a reference. Still, beware addressing temp's ptr fields!)
        FullStateDiagMatr temp = matr;
        temp.isGpuAccelerated = 1;
        temp.gpuElems = (qureg.isDistributed)?
            qureg.gpuCommBuffer : 
            gpu_allocArray(temp.numElems);

        // error if that (relatively) small allocation failed (always succeeds if buffer)
        assert_applyFullStateDiagMatrTempGpuAllocSucceeded(temp.gpuElems);

        // harmlessly overwrite new memory or qureg's buffer, and call GPU routine
        gpu_copyCpuToGpu(temp);
        gpuFunc(qureg, temp, exponent);

        // free new GPU memory, but do NOT free qureg's communication buffer
        if (!qureg.isDistributed)
            gpu_deallocArray(temp.gpuElems);
    }
}


void accel_densmatr_allTargDiagMatr_subB(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool multiplyOnly) {

    assert_fullStateDiagMatrIsDistributed(matr);
    assert_acceleratorQuregIsDistributed(qureg);

    // qureg's communication buffer (matching its own CPU or GPU deployment) 
    // already contains all elements of matr; so we simply spoof matr having
    // its own full-size local memory (matching qureg's GPU/CPU), by grafting 
    // qureg's buffer to it, and call _subA() above. It's ergo crucial _subA()
    // does not try to access qureg's communication buffer, which it safely
    // does not in this "qureg deployment = matr deployment" scenario.

    // we use a paranoid copy of matr, even though it is already a mere copy
    // of the user's matr, in case this one day changes to a reference
    FullStateDiagMatr temp = matr;

    // which is non-distributed
    temp.isDistributed = 0;
    temp.numElemsPerNode = temp.numElems;

    // and matches qureg's CPU vs GPU deployment
    temp.isGpuAccelerated = qureg.isGpuAccelerated;
    temp.cpuElems = qureg.cpuCommBuffer;
    temp.gpuElems = qureg.gpuCommBuffer;

    accel_densmatr_allTargDiagMatr_subA(qureg, temp, exponent, multiplyOnly);
}



/*
 * PAULI TENSOR AND GADGET
 */


void accel_statevector_anyCtrlPauliTensorOrGadget_subA(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> suffixTargsXY, 
    qindex suffixMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
) {
    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( statevector_anyCtrlPauliTensorOrGadget_subA, qureg, ctrls.size(), suffixTargsXY.size() );
    func(qureg, ctrls, ctrlStates, suffixTargsXY, suffixMaskXY, allMaskYZ, powI, fac0, fac1);
}
void accel_statevector_anyCtrlPauliTensorOrGadget_subB(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates,
    qindex suffixMaskXY, qindex bufferMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
) {
    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevector_anyCtrlPauliTensorOrGadget_subB, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, suffixMaskXY, bufferMaskXY, allMaskYZ, powI, fac0, fac1);
}


void accel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, 
    qcomp fac0, qcomp fac1
) {
    // no template nor compile-time optimisation necessary for the number of targs
    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevector_anyCtrlAnyTargZOrPhaseGadget_sub, qureg, ctrls.size() );
    func(qureg, ctrls, ctrlStates, targs, fac0, fac1);
}



/*
 * QUREG COMBINATION
 */


void accel_densmatr_mixQureg_subA(qreal outProb, Qureg out, qreal inProb, Qureg in) {

    // quregs are equally-sized density matrices and are equally-distributed... 
    assert_mixedQuregIsDensityMatrix(out);
    assert_mixedQuregIsDensityMatrix(in);
    assert_mixedQuregsAreBothOrNeitherDistributed(out, in);

    // but may differ in GPU accel
    bool outGPU = out.isGpuAccelerated;
    bool inGPU = in.isGpuAccelerated;

    // when deployments match, we trivially call the common backend
    if (outGPU && inGPU)
        gpu_densmatr_mixQureg_subA(outProb, out, inProb, in);
    if (!outGPU && !inGPU)
        cpu_densmatr_mixQureg_subA(outProb, out, inProb, in);

    // deployments differing is a strange and expectedly rare scenario;
    // why use GPU-acceleration for one Qureg but not the equally-sized
    // other qureg? We provide the below fallbacks for defensive design.

    // When GPU-accel differs, we fall-back to copying memory to RAM and
    // using the CPU backend. In theory, we could instead copy
    // the non-GPU qureg into the VRAM buffer of the GPU qureg and
    // always use the GPU backend, but this is only possible when
    // the buffer exists (GPU qureg is distributed), and complicates
    // the backend; unworthwhile for such a rare scenario. We could
    // also create new, temporary GPU memory and attach it to the
    // non-GPU Qureg, but the new allocation is the same size of
    // the Quregs so is dangerously large, and may fail.

    if (!outGPU && inGPU) {
        gpu_copyGpuToCpu(in);
        cpu_densmatr_mixQureg_subA(outProb, out, inProb, in);
    }

    if (outGPU && !inGPU) {
        gpu_copyGpuToCpu(out);
        cpu_densmatr_mixQureg_subA(outProb, out, inProb, in);
        gpu_copyCpuToGpu(out);
    }
}


void accel_densmatr_mixQureg_subB(qreal outProb, Qureg out, qreal inProb, Qureg in) {

    // quregs are densmatr and statevec, and are both non-distributed...
    assert_mixedQuregIsDensityMatrix(out);
    assert_mixedQuregIsStatevector(in);
    assert_mixedQuregIsLocal(out);
    assert_mixedQuregIsLocal(in);

    // but may differ in GPU accel
    bool outGPU = out.isGpuAccelerated;
    bool inGPU = in.isGpuAccelerated;

    // when deployments match, we trivially call the common backend
    if (outGPU && inGPU)
        gpu_densmatr_mixQureg_subB(outProb, out, inProb, in);
    if (!outGPU && !inGPU)
        cpu_densmatr_mixQureg_subB(outProb, out, inProb, in);

    // GPU-accelarated smaller register defaults to CPU
    if (!outGPU && inGPU) {
        gpu_copyGpuToCpu(in);
        cpu_densmatr_mixQureg_subB(outProb, out, inProb, in);
    }
    
    // GPU-accelerated larger register is a very common scenario,
    // but is irksome because without communication buffers, there
    // is no existing GPU memory to copy CPU-only small register to.
    // So we regrettably create temporary GPU memory, which will 
    // thankfully be very small; quadratically smaller than 'out').
    // Because quregs are local, there are no buffers to re-use
    if (outGPU && !inGPU) {

        // make a cheap copy of 'in' but with GPU memory
        // (we use a paranoid copy of 'in', even though 'in'' is already a 
        // mere copy of the user's qureg, in case this code changes to
        // accept a reference. Still, beware addressing in's ptr fields!)
        Qureg temp = in;
        temp.isGpuAccelerated = 1;
        temp.gpuAmps = gpu_allocArray(temp.numAmpsPerNode);
        assert_mixQuregTempGpuAllocSucceeded(temp.gpuAmps);

        // clone in's CPU memory to copy's new GPU memory, simulate, then free
        gpu_copyCpuToGpu(temp);
        gpu_densmatr_mixQureg_subB(outProb, out, inProb, temp);
        gpu_deallocArray(temp.gpuAmps);
    }
}


void accel_densmatr_mixQureg_subC(qreal outProb, Qureg out, qreal inProb) {

    // statevector has been copied to out's GPU or CPU buffer
    assert_mixedQuregIsDensityMatrix(out);
    assert_mixedQuregIsDistributed(out);

    (out.isGpuAccelerated)?
        gpu_densmatr_mixQureg_subC(outProb, out, inProb):
        cpu_densmatr_mixQureg_subC(outProb, out, inProb);
}


void accel_densmatr_mixQureg_subD(qreal outProb, Qureg out, qreal inProb, Qureg in) {

    // 'in' is local statevec and 'out' is a distributed density matrix...
    assert_mixedQuregIsDensityMatrix(out);
    assert_mixedQuregIsStatevector(in);
    assert_mixedQuregIsDistributed(out);
    assert_mixedQuregIsLocal(in);

    // but they may differ in GPU deployment
    bool outGPU = out.isGpuAccelerated;
    bool inGPU = in.isGpuAccelerated;

    // we copy 'in' into 'out's communication buffer and invoke subC;
    // the choice of buffer (CPU or GPU) depends on 'out's deployment
    qindex len = in.numAmps;

    if (outGPU && !inGPU)
        gpu_copyCpuToGpu(in.cpuAmps, out.gpuCommBuffer, len);
    if (!outGPU && inGPU)
        gpu_copyGpuToCpu(in.gpuAmps, out.cpuCommBuffer, len);

    // when 'in' and 'out' are identically deployed, we can
    // avoid copies by temporarily re-assigning pointers
    qcomp* cpuPtr = out.cpuCommBuffer;
    qcomp* gpuPtr = out.gpuCommBuffer; // may be nullptr

    if ( outGPU &&  inGPU) out.gpuCommBuffer = in.gpuAmps;
    if (!outGPU && !inGPU) out.cpuCommBuffer = in.cpuAmps;

    accel_densmatr_mixQureg_subC(outProb, out, inProb);

    // restore pointers in case they were modified
    out.cpuCommBuffer = cpuPtr;
    out.gpuCommBuffer = gpuPtr;
}



/*
 * DEPHASING
 */


void accel_densmatr_oneQubitDephasing_subA(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDephasing_subA(qureg, qubit, prob):
        cpu_densmatr_oneQubitDephasing_subA(qureg, qubit, prob);
}
void accel_densmatr_oneQubitDephasing_subB(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDephasing_subB(qureg, qubit, prob):
        cpu_densmatr_oneQubitDephasing_subB(qureg, qubit, prob);
}


void accel_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDephasing_subA(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDephasing_subA(qureg, qubit1, qubit2, prob);
}
void accel_densmatr_twoQubitDephasing_subB(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDephasing_subB(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDephasing_subB(qureg, qubit1, qubit2, prob);
}



/*
 * DEPOLARISING
 */


void accel_densmatr_oneQubitDepolarising_subA(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDepolarising_subA(qureg, qubit, prob):
        cpu_densmatr_oneQubitDepolarising_subA(qureg, qubit, prob);
}
void accel_densmatr_oneQubitDepolarising_subB(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDepolarising_subB(qureg, qubit, prob):
        cpu_densmatr_oneQubitDepolarising_subB(qureg, qubit, prob);
}


void accel_densmatr_twoQubitDepolarising_subA(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDepolarising_subA(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDepolarising_subA(qureg, qubit1, qubit2, prob);
}
void accel_densmatr_twoQubitDepolarising_subB(Qureg qureg, int qubit1, int qubit2, qreal prob) {
    
    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDepolarising_subB(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDepolarising_subB(qureg, qubit1, qubit2, prob);
}
void accel_densmatr_twoQubitDepolarising_subC(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDepolarising_subC(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDepolarising_subC(qureg, qubit1, qubit2, prob);
}
void accel_densmatr_twoQubitDepolarising_subD(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDepolarising_subD(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDepolarising_subD(qureg, qubit1, qubit2, prob);
}
void accel_densmatr_twoQubitDepolarising_subE(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDepolarising_subE(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDepolarising_subE(qureg, qubit1, qubit2, prob);
}
void accel_densmatr_twoQubitDepolarising_subF(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDepolarising_subF(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDepolarising_subF(qureg, qubit1, qubit2, prob);
}
void accel_densmatr_twoQubitDepolarising_subG(Qureg qureg, int qubit1, int qubit2, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_twoQubitDepolarising_subG(qureg, qubit1, qubit2, prob):
        cpu_densmatr_twoQubitDepolarising_subG(qureg, qubit1, qubit2, prob);
}



/*
 * PAULI CHANNEL
 */


void accel_densmatr_oneQubitPauliChannel_subA(Qureg qureg, int qubit, qreal pI, qreal pX, qreal pY, qreal pZ) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitPauliChannel_subA(qureg, qubit, pI, pX, pY, pZ):
        cpu_densmatr_oneQubitPauliChannel_subA(qureg, qubit, pI, pX, pY, pZ);
}
void accel_densmatr_oneQubitPauliChannel_subB(Qureg qureg, int ketQubit, qreal pI, qreal pX, qreal pY, qreal pZ) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitPauliChannel_subB(qureg, ketQubit, pI, pX, pY, pZ):
        cpu_densmatr_oneQubitPauliChannel_subB(qureg, ketQubit, pI, pX, pY, pZ);
}



/*
 * AMPLITUDE DAMPING CHANNEL
 */


void accel_densmatr_oneQubitDamping_subA(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDamping_subA(qureg, qubit, prob):
        cpu_densmatr_oneQubitDamping_subA(qureg, qubit, prob);
}
void accel_densmatr_oneQubitDamping_subB(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDamping_subB(qureg, qubit, prob):
        cpu_densmatr_oneQubitDamping_subB(qureg, qubit, prob);
}
void accel_densmatr_oneQubitDamping_subC(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDamping_subC(qureg, qubit, prob):
        cpu_densmatr_oneQubitDamping_subC(qureg, qubit, prob);
}
void accel_densmatr_oneQubitDamping_subD(Qureg qureg, int qubit, qreal prob) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_oneQubitDamping_subD(qureg, qubit, prob):
        cpu_densmatr_oneQubitDamping_subD(qureg, qubit, prob);
}



/*
 * PARTIAL TRACE
 */


void accel_densmatr_partialTrace_sub(Qureg inQureg, Qureg outQureg, vector<int> targs, vector<int> pairTargs) {

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_TARGS( cpu_densmatr_partialTrace_sub, targs.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_TARGS( gpu_densmatr_partialTrace_sub, targs.size() );
    
    // GPU-acceleration only possible if both Quregs are GPU-enabled
    auto useFunc = (inQureg.isGpuAccelerated && outQureg.isGpuAccelerated)? gpuFunc : cpuFunc;
    useFunc(inQureg, outQureg, targs, pairTargs);
}



/*
 * PROBABILITIES
 */


qreal accel_statevec_calcTotalProb_sub(Qureg qureg) {

    return (qureg.isGpuAccelerated)?
        gpu_statevec_calcTotalProb_sub(qureg):
        cpu_statevec_calcTotalProb_sub(qureg);
}
qreal accel_densmatr_calcTotalProb_sub(Qureg qureg) {

    return (qureg.isGpuAccelerated)?
        gpu_densmatr_calcTotalProb_sub(qureg):
        cpu_densmatr_calcTotalProb_sub(qureg);
}


qreal accel_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, bool realOnly) {

    auto func = GET_CPU_OR_GPU_BOOLEAN_FUNC_OPTIMISED_FOR_NUM_TARGS( statevec_calcProbOfMultiQubitOutcome_sub, qureg, qubits.size(), realOnly );
    return func(qureg, qubits, outcomes);
}


void accel_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( statevec_calcProbsOfAllMultiQubitOutcomes_sub, qureg, qubits.size() );
    return func(outProbs, qureg, qubits);
}
void accel_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( densmatr_calcProbsOfAllMultiQubitOutcomes_sub, qureg, qubits.size() );
    return func(outProbs, qureg, qubits);
}
