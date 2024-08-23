/** @file
 * Internal functions for choosing which accelerator backend
 * (CPU or GPU) to call, and which qubit preconditions are
 * satisfied in order to accelerate simulation.
 */

#ifndef ACCELERATOR_HPP
#define ACCELERATOR_HPP

#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include <vector>

using std::vector;


/*
 * TEMPLATE INSTANTIATION MACROS
 *
 * used by cpu_subroutines.cpp and gpu_subroutines to force the compiler
 * to instantiate and compile their template definitions with the given
 * explicit parameters below. Notice the final parameter is always -1, 
 * to handle when the number of controls or targets is not known at 
 * compile-time (it is larger than a bespoke, optimised instantiations), 
 * causing the optimised function to fallback to a suboptimal but general 
 * implementation.
 */

// must match the macros below, and those in accelerator.cpp
#define MAX_OPTIMISED_NUM_CTRLS 5
#define MAX_OPTIMISED_NUM_TARGS 5

#define INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(returntype, funcname, args) \
    template returntype funcname <0> args; \
    template returntype funcname <1> args; \
    template returntype funcname <2> args; \
    template returntype funcname <3> args; \
    template returntype funcname <4> args; \
    template returntype funcname <5> args; \
    template returntype funcname<-1> args;

#define INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(returntype, funcname, args) \
    private_INSTANTIATE(returntype, funcname, 0, args); \
    private_INSTANTIATE(returntype, funcname, 1, args); \
    private_INSTANTIATE(returntype, funcname, 2, args); \
    private_INSTANTIATE(returntype, funcname, 3, args); \
    private_INSTANTIATE(returntype, funcname, 4, args); \
    private_INSTANTIATE(returntype, funcname, 5, args); \
    private_INSTANTIATE(returntype, funcname,-1, args);

#define private_INSTANTIATE(returntype, funcname, numtargs, args) \
    template returntype funcname <0, numtargs> args; \
    template returntype funcname <1, numtargs> args; \
    template returntype funcname <2, numtargs> args; \
    template returntype funcname <3, numtargs> args; \
    template returntype funcname <4, numtargs> args; \
    template returntype funcname <5, numtargs> args; \
    template returntype funcname <-1,numtargs> args;


/*
 * COMPILE-TIME VARIABLE MACROS
 *
 * used by cpu_subroutines.cpp and gpu_subroutines to attemptedly set
 * a variable to a value known at compile-time (like a templated function's
 * parameter), enabling compile-time optimisations of subsequent code which 
 * uses the variable such a loop unrolling. If the value is not known at
 * compile-time (compileval==-1 which indicates a templated function has
 * been called with more qubits than it has been explicitly instantiated and
 * optimised for), the runtime value is used, precluding optimisations.
 */

#define SET_VAR_AT_COMPILE_TIME(type, name, compileval, runtimeval) \
    type name; \
    if constexpr (compileval == -1) \
        name = runtimeval; \
    else \
        name = compileval;


/*
 * COMMUNICATION BUFFER PACKING
 */

qindex accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubits, vector<int> qubitStates);


/*
 * SWAPS
 */

void accel_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2);
void accel_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates);
void accel_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState);


/*
 * ONE-TARGET MATRIX
 */

void accel_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
void accel_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);


/*
 * ANY-TARGET MATRIX
 */

void accel_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr);

void accel_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr);


/*
 * PAULI TENSOR AND GADGET
 */

void accel_statevector_anyCtrlPauliTensorOrGadget_subA(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> suffixTargsXY, 
    qindex suffixMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1);

void accel_statevector_anyCtrlPauliTensorOrGadget_subB(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates,
    qindex suffixMaskXY, qindex bufferMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1);

void accel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs,
    qcomp fac0, qcomp fac1);


/*
 * DECOHERENCE
 */

void accel_densmatr_oneQubitDephasing_subA(Qureg qureg, int qubit, qreal prob);
void accel_densmatr_oneQubitDephasing_subB(Qureg qureg, int qubit, qreal prob);

void accel_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubitA, int qubitB, qreal prob);
void accel_densmatr_twoQubitDephasing_subB(Qureg qureg, int qubitA, int qubitB, qreal prob);


#endif // ACCELERATOR_HPP