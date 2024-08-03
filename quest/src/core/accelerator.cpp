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

#include "quest/src/core/indexer.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/cpu/cpu_subroutines.hpp"
#include "quest/src/gpu/gpu_subroutines.hpp"

#include <vector>

using std::vector;
using namespace index_flags;



/*
 * MACROS
 *
 * which automate the choosing of the appropriate backend template function
 * (optimised for the given configuration of qubit indices), and whether to
 * dispatch to the OpenMP-accelerated backend (cpu_subroutines.cpp) or the
 * CUDA-accelerated backend (gpu_subroutines.cpp). The arguments to wrapped
 * function calls are given as variadic arguments.
 */


#define CALL_FUNC_WITH_CTRL_FLAG(flag, func, ...) \
    switch (flag) { \
        case NO_CTRLS          : func <NO_CTRLS>          (__VA_ARGS__); break; \
        case ONE_CTRL          : func <ONE_CTRL>          (__VA_ARGS__); break; \
        case ONE_STATE_CTRL    : func <ONE_STATE_CTRL>    (__VA_ARGS__); break; \
        case MULTI_CTRLS       : func <MULTI_CTRLS>       (__VA_ARGS__); break; \
        case MULTI_STATE_CTRLS : func <MULTI_STATE_CTRLS> (__VA_ARGS__); break; \
    }


#define CALL_CPU_OR_GPU_FUNC_WITH_CTRL_FLAG(isGpuAccel, flag, funcSuffix, ...) \
    if (isGpuAccel) { \
        CALL_FUNC_WITH_CTRL_FLAG(flag, gpu_##funcSuffix, __VA_ARGS__); \
    } else { \
        CALL_FUNC_WITH_CTRL_FLAG(flag, cpu_##funcSuffix, __VA_ARGS__); \
    }
        


/*
 * ANY-CTRL ONE-TARG MATRIX
 */


template <class MatrType> 
void inner_statevector_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, MatrType matr) {

    // determine ctrls preconditions to optimise indexing and memory access in backend function
    CtrlFlag ctrlFlag = indexer_getCtrlFlag(ctrls, ctrlStates);

    CALL_CPU_OR_GPU_FUNC_WITH_CTRL_FLAG(
        qureg.isGpuAccelerated, 
        ctrlFlag, 
        statevector_anyCtrlOneTargMatrix_subA, 
        qureg, ctrls, ctrlStates, targ, matr
    );
}

void statevector_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
    inner_statevector_anyCtrlOneTargMatrix_subA(qureg, ctrls, ctrlStates, targ, matr);
}

void statevector_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {
    inner_statevector_anyCtrlOneTargMatrix_subA(qureg, ctrls, ctrlStates, targ, matr);
}
