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


#define CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG(flag, func, ...) \
    switch (flag) { \
        case NO_CTRLS          : return func <NO_CTRLS>          (__VA_ARGS__); \
        case ONE_CTRL          : return func <ONE_CTRL>          (__VA_ARGS__); \
        case ONE_STATE_CTRL    : return func <ONE_STATE_CTRL>    (__VA_ARGS__); \
        case MULTI_CTRLS       : return func <MULTI_CTRLS>       (__VA_ARGS__); \
        case MULTI_STATE_CTRLS : return func <MULTI_STATE_CTRLS> (__VA_ARGS__); \
    }



/*
 * COMMUNICATION BUFFER PACKING
 */


qindex accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates) {
    indexer_assertValidCtrls(ctrls, ctrlStates);

    // we can never pack and swap buffers when there are no ctrl qubits, because we'd fill the entire buffer
    // and would ergo have no room to receive the other node's buffer; we'd send amps straight to buffer instead
    if (ctrls.empty())
        error_noCtrlsGivenToBufferPacker();

    // packing can be optimised depending on the precondition of control qubits
    CtrlFlag flag = indexer_getCtrlFlag(ctrls, ctrlStates);

    // return the number of packed amplitudes, as returned by backend function
    if (qureg.isGpuAccelerated) {
        CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG( flag, gpu_statevec_packAmpsIntoBuffer, qureg, ctrls, ctrlStates )
    } else {
        CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG( flag, cpu_statevec_packAmpsIntoBuffer, qureg, ctrls, ctrlStates )
    }

    // unreachable; above if/else contains switches where each case returns
    return 0;
}



/*
 * ANY-CTRL ONE-TARG MATRIX
 */


template <class MatrType> 
void inner_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, MatrType matr) {
    indexer_assertValidCtrls(ctrls, ctrlStates);

    // simulation is optimised depending on the precondition of control qubits
    CtrlFlag flag = indexer_getCtrlFlag(ctrls, ctrlStates);

    if (qureg.isGpuAccelerated) {
        CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG( flag, gpu_statevec_anyCtrlOneTargMatrix_subA, qureg, ctrls, ctrlStates, targ, matr )
    } else {
        CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG( flag, cpu_statevec_anyCtrlOneTargMatrix_subA, qureg, ctrls, ctrlStates, targ, matr )
    }
}

void accel_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
    inner_statevec_anyCtrlOneTargMatrix_subA(qureg, ctrls, ctrlStates, targ, matr);
}

void accel_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr) {
    inner_statevec_anyCtrlOneTargMatrix_subA(qureg, ctrls, ctrlStates, targ, matr);
}

void accel_statevec_anyCtrlOneTargDenseMatrix_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1) {
    indexer_assertValidCtrls(ctrls, ctrlStates);

    // simulation is optimised depending on the precondition of control qubits
    CtrlFlag flag = indexer_getCtrlFlag(ctrls, ctrlStates);

    if (qureg.isGpuAccelerated) {
        CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG( flag, gpu_statevec_anyCtrlOneTargDenseMatrix_subB, qureg, ctrls, ctrlStates, fac0, fac1 )
    } else {
        CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG( flag, cpu_statevec_anyCtrlOneTargDenseMatrix_subB, qureg, ctrls, ctrlStates, fac0, fac1 )
    }
}





// NEW STUFF

void NEW_accel_statevec_anyCtrlOneTargMatrix_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr) {
    indexer_assertValidCtrls(ctrls, ctrlStates);


    if (qureg.isGpuAccelerated) {
        
        // EH DOING NOTHIN
    } else {


        int numQubits = 1 + ctrls.size();

        // ETC
        if (numQubits == 1)
            NEW_cpu_statevec_anyCtrlOneTargMatrix_subA<1>(qureg, ctrls, ctrlStates, targ, matr);
        if (numQubits == 2)
            NEW_cpu_statevec_anyCtrlOneTargMatrix_subA<2>(qureg, ctrls, ctrlStates, targ, matr);
        if (numQubits == 3)
            NEW_cpu_statevec_anyCtrlOneTargMatrix_subA<3>(qureg, ctrls, ctrlStates, targ, matr);
        

       // CALL_AND_RETURN_FUNC_WITH_CTRL_FLAG( flag, cpu_statevec_anyCtrlOneTargMatrix_subA, qureg, ctrls, ctrlStates, targ, matr )
    }
}