/** @file
 * Internal functions for choosing which accelerator backend
 * (CPU or GPU) to dispatch to, and which preconditions the 
 * qubit indices satisfy (informing which compile-time
 * optimisations to use) in order to effect local simulation 
 * subroutines upon Quregs.
 * 
 * These routines are called by localiser.cpp and are embarrassingly 
 * parallel, so are always called before/after any necessary
 * communication has happened. The data they need must already be
 * localised into the appropriate memory (RAM or VRAM) and location
 * (qureg's amplitudes or buffer space).
 * 
 * @author Tyson Jones
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/core/accelerator.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/memory.hpp"
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
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)]

#define GET_FUNC_OPTIMISED_FOR_NUM_TARGS(f, numtargs) \
    (vector <decltype(&f<0>)> {&f<0>, &f<1>, &f<2>, &f<3>, &f<4>, &f<5>, &f<-1>}) \
    [std::min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(f, numctrls, numtargs) \
    (vector <ARR(f)> { \
        ARR(f) {&f<0,0>,  &f<0,1>,  &f<0,2>,  &f<0,3>,  &f<0,4>,  &f<0,5>,  &f<0,-1>}, \
        ARR(f) {&f<1,0>,  &f<1,1>,  &f<1,2>,  &f<1,3>,  &f<1,4>,  &f<1,5>,  &f<1,-1>}, \
        ARR(f) {&f<2,0>,  &f<2,1>,  &f<2,2>,  &f<2,3>,  &f<2,4>,  &f<2,5>,  &f<2,-1>}, \
        ARR(f) {&f<3,0>,  &f<3,1>,  &f<3,2>,  &f<3,3>,  &f<3,4>,  &f<3,5>,  &f<3,-1>}, \
        ARR(f) {&f<4,0>,  &f<4,1>,  &f<4,2>,  &f<4,3>,  &f<4,4>,  &f<4,5>,  &f<4,-1>}, \
        ARR(f) {&f<5,0>,  &f<5,1>,  &f<5,2>,  &f<5,3>,  &f<5,4>,  &f<5,5>,  &f<5,-1>}, \
        ARR(f) {&f<-1,0>, &f<-1,1>, &f<-1,2>, &f<-1,3>, &f<-1,4>, &f<-1,5>, &f<-1,-1>}}) \
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)] \
    [std::min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

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


/// @todo
/// GET_CPU_OR_GPU_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS,
/// as defined below, is only ever called by used by anyCtrlAnyTargDenseMatr,
/// which only ever receives numTargs>=3 (due to accelerator redirecting 
/// fewer targets to faster bespoke functions which e.g. avoid global GPU
/// cache emory access). This means its instantiation with numTargs=0,1,2
/// is useless, though contributes to 42% of the function's compilation
/// time which is large because of the 7*7*2=98 unique instantiations. We
/// can ergo non-negligibly speed up compilation by avoiding these redundant 
/// instances at the cost of increased code complexity/asymmetry. Consider!

#define GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(f, numctrls, numtargs, c) \
    (vector <CONJ_ARR(f)> { \
        CONJ_ARR(f) {&f<0,0,c>,  &f<0,1,c>,  &f<0,2,c>,  &f<0,3,c>,  &f<0,4,c>,  &f<0,5,c>,  &f<0,-1,c>}, \
        CONJ_ARR(f) {&f<1,0,c>,  &f<1,1,c>,  &f<1,2,c>,  &f<1,3,c>,  &f<1,4,c>,  &f<1,5,c>,  &f<1,-1,c>}, \
        CONJ_ARR(f) {&f<2,0,c>,  &f<2,1,c>,  &f<2,2,c>,  &f<2,3,c>,  &f<2,4,c>,  &f<2,5,c>,  &f<2,-1,c>}, \
        CONJ_ARR(f) {&f<3,0,c>,  &f<3,1,c>,  &f<3,2,c>,  &f<3,3,c>,  &f<3,4,c>,  &f<3,5,c>,  &f<3,-1,c>}, \
        CONJ_ARR(f) {&f<4,0,c>,  &f<4,1,c>,  &f<4,2,c>,  &f<4,3,c>,  &f<4,4,c>,  &f<4,5,c>,  &f<4,-1,c>}, \
        CONJ_ARR(f) {&f<5,0,c>,  &f<5,1,c>,  &f<5,2,c>,  &f<5,3,c>,  &f<5,4,c>,  &f<5,5,c>,  &f<5,-1,c>}, \
        CONJ_ARR(f) {&f<-1,0,c>, &f<-1,1,c>, &f<-1,2,c>, &f<-1,3,c>, &f<-1,4,c>, &f<-1,5,c>, &f<-1,-1,c>}}) \
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)] \
    [std::min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

#define CONJ_ARR(f) vector<decltype(&f<0,0,false>)>

#define GET_CPU_OR_GPU_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(funcsuffix, qureg, numctrls, numtargs, conj) \
    ((qureg.isGpuAccelerated)? \
        ((conj)? \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, true ) : \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( gpu_##funcsuffix, numctrls, numtargs, false ) ) : \
        ((conj)? \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, true ) : \
            GET_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( cpu_##funcsuffix, numctrls, numtargs, false ) ) )


/// @todo
/// This has gotten a bit ridiculous. Is there a way to use (likely)
/// more abominable pre-processor mischief which negates the need
/// to repeat the entire macro(s) when the number of templated
/// parameters grows?


#define GET_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(f, numctrls, numtargs, c, h) \
    (vector <POWER_CONJ_ARR(f)> { \
        POWER_CONJ_ARR(f) {&f<0,0,c,h>,  &f<0,1,c,h>,  &f<0,2,c,h>,  &f<0,3,c,h>,  &f<0,4,c,h>,  &f<0,5,c,h>,  &f<0,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<1,0,c,h>,  &f<1,1,c,h>,  &f<1,2,c,h>,  &f<1,3,c,h>,  &f<1,4,c,h>,  &f<1,5,c,h>,  &f<1,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<2,0,c,h>,  &f<2,1,c,h>,  &f<2,2,c,h>,  &f<2,3,c,h>,  &f<2,4,c,h>,  &f<2,5,c,h>,  &f<2,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<3,0,c,h>,  &f<3,1,c,h>,  &f<3,2,c,h>,  &f<3,3,c,h>,  &f<3,4,c,h>,  &f<3,5,c,h>,  &f<3,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<4,0,c,h>,  &f<4,1,c,h>,  &f<4,2,c,h>,  &f<4,3,c,h>,  &f<4,4,c,h>,  &f<4,5,c,h>,  &f<4,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<5,0,c,h>,  &f<5,1,c,h>,  &f<5,2,c,h>,  &f<5,3,c,h>,  &f<5,4,c,h>,  &f<5,5,c,h>,  &f<5,-1,c,h>}, \
        POWER_CONJ_ARR(f) {&f<-1,0,c,h>, &f<-1,1,c,h>, &f<-1,2,c,h>, &f<-1,3,c,h>, &f<-1,4,c,h>, &f<-1,5,c,h>, &f<-1,-1,c,h>}}) \
    [std::min((int) numctrls, MAX_OPTIMISED_NUM_CTRLS + 1)] \
    [std::min((int) numtargs, MAX_OPTIMISED_NUM_TARGS + 1)]

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


qcomp accel_statevec_getAmp_sub(Qureg qureg, qindex localInd) {

    // we use a bespoke function, rather than merely invoking
    // getAmps() below, so that the CPU implementation can
    // make use of the faster array access, rather than memcpy,
    // and we keep the bespoke GPU function for symmetry/consistency 

    return (qureg.isGpuAccelerated)?
        gpu_statevec_getAmp_sub(qureg, localInd):
        cpu_statevec_getAmp_sub(qureg, localInd);
}


void accel_statevec_getAmps_sub(qcomp* outAmps, Qureg qureg, qindex localStartInd, qindex numLocalAmps) {

    // copy directly from GPU/CPU to outAmps
    (qureg.isGpuAccelerated)?
        gpu_copyGpuToCpu(&qureg.gpuAmps[localStartInd], outAmps, numLocalAmps): // (src, dest) = (gpu, cpu)
        cpu_copyArray(   outAmps, &qureg.cpuAmps[localStartInd], numLocalAmps); // (dest, src)
}



/*
 * SETTERS 
 */


void accel_statevec_setAmps_sub(qcomp* inAmps, Qureg qureg, qindex localStartInd, qindex numLocalAmps) {

    // in CPU settings, we use memory-copying rather than OpenMP
    // loop updating, because the latter is only faster when carefully
    // optimising parallelisation granularity with the memory
    // architecture, which we cannot reliably do in a platform
    // agnostic way (except via hwloc or something)

    // copy directly from inAmps to GPU/CPU
    (qureg.isGpuAccelerated)?
        gpu_copyCpuToGpu(inAmps, &qureg.gpuAmps[localStartInd], numLocalAmps): // (src, dest) = (cpu, gpu)
        cpu_copyArray(   &qureg.cpuAmps[localStartInd], inAmps, numLocalAmps); // (dest, src)
}


void accel_densmatr_setAmpsToPauliStrSum_sub(Qureg qureg, PauliStrSum sum) {

    (qureg.isGpuAccelerated)?
        gpu_densmatr_setAmpsToPauliStrSum_sub(qureg, sum):
        cpu_densmatr_setAmpsToPauliStrSum_sub(qureg, sum);
}


void accel_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in) {

    // use GPU to populate FullStateDiagMatr if available
    (out.isGpuAccelerated)?
        gpu_fullstatediagmatr_setElemsToPauliStrSum(out, in):
        cpu_fullstatediagmatr_setElemsToPauliStrSum(out, in);

    // but thereafter copy to CPU, to keep GPU and CPU consistent
    if (out.isGpuAccelerated)
        gpu_copyGpuToCpu(out.gpuElems, out.cpuElems, out.numElemsPerNode);
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

    bool hasPower = exponent != qcomp(1, 0);
    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_BOOL( cpu_statevec_allTargDiagMatr_sub, hasPower );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_BOOL( gpu_statevec_allTargDiagMatr_sub, hasPower );

    // when deployments match, we trivially call the common backend
    if ( quregGPU &&  matrGPU) gpuFunc(qureg, matr, exponent);
    if (!quregGPU && !matrGPU) cpuFunc(qureg, matr, exponent);

    // deployments differing is a strange and expectedly rare scenario;
    // why use GPU-acceleration for a Qureg but not the equally-sized
    // matrix? We provide the below fallbacks for defensive design, and
    // fall-back to copying memory to RAM and using the CPU backend. 
    // In theory, we could leverage exsting GPU memory of the Qureg's 
    // communication buffer (if it existed), but this is an even rarer 
    // situation and is hacky. We could also create new, temporary GPU 
    // memory and graft it to the non-accelerated object, but the new 
    // allocation would be the same size as the objects and ergo be 
    // dangerously large.

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

    bool hasPower = exponent != qcomp(1, 0);
    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( cpu_densmatr_allTargDiagMatr_sub, hasPower, multiplyOnly );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( gpu_densmatr_allTargDiagMatr_sub, hasPower, multiplyOnly );

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
            gpu_allocArray(temp.numElemsPerNode);

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


void accel_statevector_anyCtrlPauliTensorOrGadget_subA(Qureg qureg, vector<int> ctrls, vector<int> states, vector<int> x, vector<int> y, vector<int> z, qcomp f0, qcomp f1) {

    // only X and Y constitute target qubits (Z merely induces a phase)
    int numTargs = x.size() + y.size();

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS( statevector_anyCtrlPauliTensorOrGadget_subA, qureg, ctrls.size(), numTargs );
    func(qureg, ctrls, states, x, y, z, f0, f1);
}
void accel_statevector_anyCtrlPauliTensorOrGadget_subB(Qureg qureg, vector<int> ctrls, vector<int> states, vector<int> x, vector<int> y, vector<int> z, qcomp f0, qcomp f1, qindex mask) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevector_anyCtrlPauliTensorOrGadget_subB, qureg, ctrls.size() );
    func(qureg, ctrls, states, x, y, z, f0, f1, mask);
}


void accel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(Qureg qureg, vector<int> ctrls, vector<int> states, vector<int> targs, qcomp f0, qcomp f1) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_CTRLS( statevector_anyCtrlAnyTargZOrPhaseGadget_sub, qureg, ctrls.size() );
    func(qureg, ctrls, states, targs, f0, f1);
}



/*
 * QUREG COMBINATION
 */


void accel_statevec_setQuregToSuperposition_sub(qcomp facOut, Qureg outQureg, qcomp fac1, Qureg inQureg1, qcomp fac2, Qureg inQureg2) {

    // consult outQureg's deployment (other quregs should match, though we dangerously do not assert this post-validation)
    (outQureg.isGpuAccelerated)?
        gpu_statevec_setQuregToSuperposition_sub(facOut, outQureg, fac1, inQureg1, fac2, inQureg2):
        cpu_statevec_setQuregToSuperposition_sub(facOut, outQureg, fac1, inQureg1, fac2, inQureg2); 
}


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
    // why use GPU-acceleration for a Qureg but not the equally-sized
    // matrix? We provide the below fallbacks for defensive design, and
    // fall-back to copying memory to RAM and using the CPU backend. 
    // In theory, we could leverage exsting GPU memory of the Qureg's 
    // communication buffer (if it existed), but this is an even rarer 
    // situation and is hacky. We could also create new, temporary GPU 
    // memory and graft it to the non-accelerated object, but the new 
    // allocation would be the same size as the objects and ergo be 
    // dangerously large.

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
    assert_partialTraceQuregsAreIdenticallyDeployed(inQureg, outQureg);

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_TARGS( cpu_densmatr_partialTrace_sub, targs.size() );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_TARGS( gpu_densmatr_partialTrace_sub, targs.size() );

    // inQureg == outQureg except for dimension, so use common backend
    auto useFunc = (inQureg.isGpuAccelerated)? gpuFunc : cpuFunc;
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


qreal accel_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( statevec_calcProbOfMultiQubitOutcome_sub, qureg, qubits.size() );
    return func(qureg, qubits, outcomes);
}
qreal accel_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( densmatr_calcProbOfMultiQubitOutcome_sub, qureg, qubits.size() );
    return func(qureg, qubits, outcomes);
}


void accel_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( statevec_calcProbsOfAllMultiQubitOutcomes_sub, qureg, qubits.size() );
    func(outProbs, qureg, qubits);
}
void accel_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( densmatr_calcProbsOfAllMultiQubitOutcomes_sub, qureg, qubits.size() );
    func(outProbs, qureg, qubits);
}


qreal accel_densmatr_calcHilbertSchmidtDistance_sub(Qureg quregA, Qureg quregB) {

    // quregs are gauranteed to be identically deployed
    return (quregA.isGpuAccelerated)?
        gpu_densmatr_calcHilbertSchmidtDistance_sub(quregA, quregB):
        cpu_densmatr_calcHilbertSchmidtDistance_sub(quregA, quregB);
}



/*
 * INNER PRODUCTS
 */


qcomp accel_statevec_calcInnerProduct_sub(Qureg quregA, Qureg quregB) {
    assert_innerProductedSameDimQuregsHaveSameGpuAccel(quregA, quregB);

    // in theory, we could permit them to differ in GPU-acceleration
    // if one (or both) is distributed; we could then hijack the
    // GPU communication buffer and copy over the CPU-only Qureg's
    // amps. But this is a nonsensical and inefficient scenario to support.

    return (quregA.isGpuAccelerated)?
        gpu_statevec_calcInnerProduct_sub(quregA, quregB):
        cpu_statevec_calcInnerProduct_sub(quregA, quregB);
}


qcomp accel_densmatr_calcFidelityWithPureState_sub(Qureg rho, Qureg psi, bool conj) {
    assert_calcFidStateVecIsLocal(psi);

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_BOOL( cpu_densmatr_calcFidelityWithPureState_sub, conj );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_BOOL( gpu_densmatr_calcFidelityWithPureState_sub, conj );

    // quregs may differ in their GPU vs CPU deployments
    bool rhoGpu = rho.isGpuAccelerated;
    bool psiGpu = psi.isGpuAccelerated;

    // if deployments agree, trivially call the common backend
    if (rhoGpu == psiGpu)
        return (rhoGpu)? gpuFunc(rho,psi) : cpuFunc(rho,psi);

    // if only the smaller psi is GPU-accel (which is sensible when the larger
    // rho is distributed and/or exceeds the GPU memory capacity), copy psi's
    // GPU memory to CPU and proceed with CPU calculation
    if (!rhoGpu && psiGpu) {
        gpu_copyGpuToCpu(psi);
        return cpuFunc(rho, psi);
    }

    // it is also possible/sensible that rho is GPU-accelerated while the quadratically-smaller
    // psi is not. In that case, we spoof a GPU-accelerated psi which re-uses rho's
    // GPU communication buffer if it exists, else creates temporary memory (not so big).
    Qureg temp = psi;
    temp.isGpuAccelerated = 1;
    temp.gpuAmps = (rho.isDistributed)?
        rho.gpuCommBuffer :
        gpu_allocArray(temp.numAmpsPerNode);

    // error if that (relatively) small allocation failed (always succeeds if buffer)
    assert_calcFidTempGpuAllocSucceeded(temp.gpuAmps);

    // harmlessly overwrite new memory or rho's buffer, and call GPU routine
    gpu_copyCpuToGpu(temp);
    qcomp prod = gpuFunc(rho, temp);
    
    // free new GPU memory, but do NOT free rho's communication buffer
    if (!rho.isDistributed)
        gpu_deallocArray(temp.gpuAmps);

    return prod;
}



/*
 * EXPECTATION VALUES
 */


qreal accel_statevec_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs) {

    return (qureg.isGpuAccelerated)?
        gpu_statevec_calcExpecAnyTargZ_sub(qureg, targs):
        cpu_statevec_calcExpecAnyTargZ_sub(qureg, targs);
}
qcomp accel_densmatr_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs) {

    return (qureg.isGpuAccelerated)?
        gpu_densmatr_calcExpecAnyTargZ_sub(qureg, targs):
        cpu_densmatr_calcExpecAnyTargZ_sub(qureg, targs);
}


qcomp accel_statevec_calcExpecPauliStr_subA(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

    return (qureg.isGpuAccelerated)?
        gpu_statevec_calcExpecPauliStr_subA(qureg, x, y, z):
        cpu_statevec_calcExpecPauliStr_subA(qureg, x, y, z);
}
qcomp accel_statevec_calcExpecPauliStr_subB(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

    return (qureg.isGpuAccelerated)?
        gpu_statevec_calcExpecPauliStr_subB(qureg, x, y, z):
        cpu_statevec_calcExpecPauliStr_subB(qureg, x, y, z);
}
qcomp accel_densmatr_calcExpecPauliStr_sub(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

    return (qureg.isGpuAccelerated)?
        gpu_densmatr_calcExpecPauliStr_sub(qureg, x, y, z):
        cpu_densmatr_calcExpecPauliStr_sub(qureg, x, y, z);
}


qcomp accel_statevec_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool useRealPow) {

     // qureg and matr are identically distributed (caller may have spoofed)...
    assert_quregAndFullStateDiagMatrAreBothOrNeitherDistrib(qureg, matr);
    
    // but they may have differing GPU deployments
    bool quregGPU = qureg.isGpuAccelerated;
    bool matrGPU = matr.isGpuAccelerated;

    // disable useRealPow when exponent==1 which never invokes pow()
    bool hasPower = exponent != qcomp(1, 0);
    useRealPow = useRealPow && hasPower; 
    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( cpu_statevec_calcExpecFullStateDiagMatr_sub, hasPower, useRealPow );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( gpu_statevec_calcExpecFullStateDiagMatr_sub, hasPower, useRealPow );

    // when both are GPU-accelerated, we trivially use the GPU backend
    if (quregGPU && matrGPU)
        return gpuFunc(qureg, matr, exponent);

    // but otherwise, we must fall back to a CPU function, copying
    // Qureg's memory (if GPU-accelerated) into RAM. This is because
    // we lack existing GPU memory to copy the non-GPU object into
    // (communication buffers might not exist; note even a distributed
    // qureg might merely be spoofed). This is anyway a silly scenario
    // only supported for defensive design

    // GPU qureg may have diverged from CPU (GPU matrix never diverges)
    if (quregGPU)
        gpu_copyGpuToCpu(qureg);

    return cpuFunc(qureg, matr, exponent);
}

qcomp accel_densmatr_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool useRealPow) {

     // qureg and matr are identically distributed (caller may have spoofed)...
    assert_quregAndFullStateDiagMatrAreBothOrNeitherDistrib(qureg, matr);

    // but they may have differing GPU deployments
    bool quregGPU = qureg.isGpuAccelerated;
    bool matrGPU = matr.isGpuAccelerated;

    // disable useRealPow when exponent==1 which never invokes pow()
    bool hasPower = exponent != qcomp(1, 0);
    useRealPow = useRealPow && hasPower;
    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( cpu_densmatr_calcExpecFullStateDiagMatr_sub, hasPower, useRealPow );
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_TWO_BOOLS( gpu_densmatr_calcExpecFullStateDiagMatr_sub, hasPower, useRealPow );

    // if deployments agree, trivially call the common backend
    if (quregGPU == matrGPU)
        return (quregGPU)? 
            gpuFunc(qureg, matr, exponent): 
            cpuFunc(qureg, matr, exponent);

    // if only the smaller matr is GPU-accel (which is sensible when the larger
    // qureg is distributed and/or exceeds the GPU memory capacity), use the CPU
    // function; matr's CPU memory should be unchanged from its GPU memory
    if (!quregGPU && matrGPU)
        return cpuFunc(qureg, matr, exponent);

    // it is also possible/sensible that the larger qureg is GPU-accelerated while 
    // the quadratically-smaller matr is not. In that case, we spoof a GPU-accelerated
    // matr which re-uses quregs's GPU communication buffer if it exists, else creates 
    // temporary memory (not so big). Note that we cannot qureg.isDistributed to
    // check whether the buffer exists, since qureg may be merely spoofed by localiser,
    // so we must explicitly use mem_isAllocated() 
    
    FullStateDiagMatr temp = matr;
    temp.isGpuAccelerated = 1;
    temp.gpuElems = mem_isAllocated(qureg.gpuCommBuffer)?
        qureg.gpuCommBuffer :
        gpu_allocArray(matr.numElemsPerNode);

    // error if that (relatively) small allocation failed (always succeeds if buffer)
    assert_calcExpecDiagTempGpuAllocSucceeded(temp.gpuElems);

    // harmlessly overwrite new memory or qureg's buffer, and call GPU routine
    gpu_copyCpuToGpu(temp);
    qcomp value = gpuFunc(qureg, temp, exponent);
    
    // free new GPU memory, but do NOT free qureg's communication buffer
    if (!mem_isAllocated(qureg.gpuCommBuffer))
        gpu_deallocArray(temp.gpuElems);

    return value;
}



/*
 * PROJECTORS 
 */


void accel_statevec_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( statevec_multiQubitProjector_sub, qureg, qubits.size() );
    func(qureg, qubits, outcomes, prob);
}
void accel_densmatr_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob) {

    auto func = GET_CPU_OR_GPU_FUNC_OPTIMISED_FOR_NUM_TARGS( densmatr_multiQubitProjector_sub, qureg, qubits.size() );
    func(qureg, qubits, outcomes, prob);
}



/*
 * STATE INITIALISATION
 */


void accel_statevec_initUniformState_sub(Qureg qureg, qcomp amp) {

    (qureg.isGpuAccelerated)?
        gpu_statevec_initUniformState_sub(qureg, amp):
        cpu_statevec_initUniformState_sub(qureg, amp);
}


void accel_statevec_initDebugState_sub(Qureg qureg) {

    (qureg.isGpuAccelerated)?
        gpu_statevec_initDebugState_sub(qureg):
        cpu_statevec_initDebugState_sub(qureg);
}


void accel_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(Qureg qureg) {

    (qureg.isGpuAccelerated)?
        gpu_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(qureg):
        cpu_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(qureg);
}
