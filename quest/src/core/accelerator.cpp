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
#include "quest/src/core/indexer.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/cpu/cpu_subroutines.hpp"
#include "quest/src/gpu/gpu_subroutines.hpp"

#include <vector>
#include <algorithm>

using std::vector;
using std::min;



/*
 * MACROS
 *
 * which automate the choosing of the appropriate backend template function
 * (optimised for the given configuration of qubit indices), and whether to
 * dispatch to the OpenMP-accelerated backend (cpu_subroutines.cpp) or the
 * CUDA-accelerated backend (gpu_subroutines.cpp). The arguments to wrapped
 * function calls are given as variadic arguments.
 */




// TODO: assert this switch matches accelerator size

// switch (num) {
//     case 0:  cpu_statevec_packAmpsIntoBuffer<0>(qureg, ctrls, ctrlStates); break;
//     case 1:  cpu_statevec_packAmpsIntoBuffer<1>(qureg, ctrls, ctrlStates);
//     case 2:  cpu_statevec_packAmpsIntoBuffer<2>(qureg, ctrls, ctrlStates);
//     case 3:  cpu_statevec_packAmpsIntoBuffer<3>(qureg, ctrls, ctrlStates);
//     case 4:  cpu_statevec_packAmpsIntoBuffer<4>(qureg, ctrls, ctrlStates);
//     case 5:  cpu_statevec_packAmpsIntoBuffer<5>(qureg, ctrls, ctrlStates);
//     default: cpu_statevec_packAmpsIntoBuffer<-1>(qureg, ctrls, ctrlStates);
// }


// #define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS(funcname, numctrls) \
//     [](int arg) { \
//         switch (arg) { \
//             case 0: return funcname<0>; \
//             case 1: return funcname<1>; \
//             case 2: return funcname<2>; \
//             case 3: return funcname<3>; \
//             case 4: return funcname<4>; \
//             case 5: return funcname<5>; \
//         } \
//         if (arg <= MAX_OPTIMISED_NUM_CTRLS) \
//             /* TODO: error */ \
//             ; \
//         return funcname<-1>; \
//     }(numctrls)



// #define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(funcname, numctrls, numtargs) \
//     [](int arg) { \
//         switch (arg) { \
//             case 0: return funcname<0>; \
//             case 1: return funcname<1>; \
//             case 2: return funcname<2>; \
//             case 3: return funcname<3>; \
//             case 4: return funcname<4>; \
//             case 5: return funcname<5>; \
//         } \
//         if (arg <= MAX_OPTIMISED_NUM_CTRLS) \
//             /* TODO: error */ \
//             ; \
//         return funcname<-1>; \
//     }(numctrls)




#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS(func, numctrls)         \
    [](size_t arg) {                                                    \
        auto funcs = vector{ func<0>, func<1>, func<2>, func<3>, func<4>, func<5>, func<-1> }; \
        auto num = min(arg, funcs.size() - 1);                   \
        return funcs[num];                                        \
    }(numctrls)







#define GET_ARRAY_OF_FUNCS_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(func) \
    vector{ INNER(func,0), INNER(func,1), INNER(func,2), INNER(func,3), INNER(func,4), INNER(func,5), INNER(func,-1) }

#define INNER(func, n) \
    vector{ func<n,0>, func<n,1>, func<n,2>, func<n,3>, func<n,4>, func<n,5>, func<n,-1> }

#define GET_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(funcname, numctrls, numtargs) \
    [](size_t nctrls, size_t ntargs) {                                                    \
        auto arr = GET_ARRAY_OF_FUNCS_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(funcname); \
        auto i = min(nctrls, arr   .size() - 1); \
        auto j = min(ntargs, arr[0].size() - 1); \
        return arr[i][j];                                  \
    }(numctrls, numtargs)










// cleanup
#define CALL_FUNC_OPTIMISED_FOR_NUM_CTRLS(...) ;



/*
 * COMMUNICATION BUFFER PACKING
 */


constexpr int MAX_NUM_CTRLS = 63;


void accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {

    // we can never pack and swap buffers when there are no ctrl qubits, because we'd fill the entire buffer
    // andhave no room to receive the other node's buffer; caller would instead send amps straight to buffer
    if (ctrls.empty())
        error_noCtrlsGivenToBufferPacker();

    auto cpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS(cpu_statevec_packAmpsIntoBuffer, ctrls.size());
    auto gpuFunc = GET_FUNC_OPTIMISED_FOR_NUM_CTRLS(gpu_statevec_packAmpsIntoBuffer, ctrls.size());
    auto useFunc = (qureg.isGpuAccelerated)? gpuFunc : cpuFunc;

    useFunc(qureg, ctrls, ctrlStates);
}



/*
 * MATRICES
 */


void accel_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {

    if (qureg.isGpuAccelerated) {
        CALL_FUNC_OPTIMISED_FOR_NUM_CTRLS( ctrls.size(), gpu_statevec_anyCtrlOneTargDenseMatr_subA, qureg, ctrls, ctrlStates, targ, matr );
    } else {
        CALL_FUNC_OPTIMISED_FOR_NUM_CTRLS( ctrls.size(), cpu_statevec_anyCtrlOneTargDenseMatr_subA, qureg, ctrls, ctrlStates, targ, matr );
    }
}

void accel_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {

    if (qureg.isGpuAccelerated) {
        CALL_FUNC_OPTIMISED_FOR_NUM_CTRLS( ctrls.size(), gpu_statevec_anyCtrlOneTargDenseMatr_subB, qureg, ctrls, ctrlStates, fac0, fac1 );
    } else {
        CALL_FUNC_OPTIMISED_FOR_NUM_CTRLS( ctrls.size(), cpu_statevec_anyCtrlOneTargDenseMatr_subB, qureg, ctrls, ctrlStates, fac0, fac1 );
    }
}


void accel_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr) {



    // if (qureg.isGpuAccelerated) {
    //     _CALL_FUNC_OPTIMISED_FOR_NUM_CTRLS_THEN_RETURN(ctrls.size(), gpu_statevec_anyCtrlAnyTargDiagMatr_sub, qureg, ctrls, ctrlStates, targ, matr);
    // } else {
    //    CALL_FUNC_OPTIMISED_FOR_NUM_CTRLS_THEN_RETURN(ctrls.size(), cpu_statevec_anyCtrlAnyTargDiagMatr_sub, qureg, ctrls, ctrlStates, targ, matr);
    // }
}