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


#define INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS(returntype, funcname, args) \
    template returntype funcname <0> args; \
    template returntype funcname <1> args; \
    template returntype funcname <2> args; \
    template returntype funcname <3> args; \
    template returntype funcname <4> args; \
    template returntype funcname <5> args; \
    template returntype funcname<-1> args;

#define INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_CTRLS(returntype, funcname, args) \
    INSTANTIATE_FUNC_OPTIMISED_FOR_NUM_TARGS(returntype, funcname, args)


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


#define INSTANTIATE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(returntype, funcname, args) \
    private_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, true,  args) \
    private_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, false, args)

#define private_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, conj, args) \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 0, conj, args); \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 1, conj, args); \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 2, conj, args); \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 3, conj, args); \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 4, conj, args); \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 5, conj, args); \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname,-1, conj, args);

#define private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, numtargs, conj, args) \
    template returntype funcname <0, numtargs, conj>  args; \
    template returntype funcname <1, numtargs, conj>  args; \
    template returntype funcname <2, numtargs, conj>  args; \
    template returntype funcname <3, numtargs, conj>  args; \
    template returntype funcname <4, numtargs, conj>  args; \
    template returntype funcname <5, numtargs, conj>  args; \
    template returntype funcname <-1,numtargs, conj>  args;



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

#define SET_CONJ_AT_COMPILE_TIME(type, name, runtimeval, compileconj) \
    type name = runtimeval; \
    if constexpr (compileconj) \
        name = conj(name);



/*
 * COMMUNICATION BUFFER PACKING
 */

qindex accel_statevec_packAmpsIntoBuffer(Qureg qureg, vector<int> qubits, vector<int> qubitStates);

qindex accel_statevec_packPairSummedAmpsIntoBuffer(Qureg qureg, int qubit1, int qubit2, int qubit3, int bit2);


/*
 * SWAPS
 */

void accel_statevec_anyCtrlSwap_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2);
void accel_statevec_anyCtrlSwap_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates);
void accel_statevec_anyCtrlSwap_subC(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, int targState);


/*
 * MATRICES
 */

void accel_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
void accel_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);

void accel_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr, bool conj);

void accel_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, bool conj);


/*
 * PAULI TENSOR AND GADGET
 */

void accel_statevector_anyCtrlPauliTensorOrGadget_subA(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> suffixTargsXY, 
    qindex suffixMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
);
void accel_statevector_anyCtrlPauliTensorOrGadget_subB(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates,
    qindex suffixMaskXY, qindex bufferMaskXY, qindex allMaskYZ, qcomp powI, qcomp fac0, qcomp fac1
);

void accel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(
    Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs,
    qcomp fac0, qcomp fac1);


/*
 * DECOHERENCE
 */

void accel_densmatr_oneQubitDephasing_subA(Qureg qureg, int qubit, qreal prob);
void accel_densmatr_oneQubitDephasing_subB(Qureg qureg, int qubit, qreal prob);

void accel_densmatr_twoQubitDephasing_subA(Qureg qureg, int qubit1, int qubit2, qreal prob);
void accel_densmatr_twoQubitDephasing_subB(Qureg qureg, int qubit1, int qubit2, qreal prob);

void accel_densmatr_oneQubitDepolarising_subA(Qureg qureg, int qubit, qreal prob);
void accel_densmatr_oneQubitDepolarising_subB(Qureg qureg, int qubit, qreal prob);

void accel_densmatr_twoQubitDepolarising_subA(Qureg qureg, int qubit1, int qubit2, qreal prob);
void accel_densmatr_twoQubitDepolarising_subB(Qureg qureg, int qubit1, int qubit2, qreal prob);
void accel_densmatr_twoQubitDepolarising_subC(Qureg qureg, int qubit1, int qubit2, qreal prob);
void accel_densmatr_twoQubitDepolarising_subD(Qureg qureg, int qubit1, int qubit2, qreal prob);
void accel_densmatr_twoQubitDepolarising_subE(Qureg qureg, int qubit1, int qubit2, qreal prob);
void accel_densmatr_twoQubitDepolarising_subF(Qureg qureg, int qubit1, int qubit2, qreal prob);
void accel_densmatr_twoQubitDepolarising_subG(Qureg qureg, int qubit1, int qubit2, qreal prob);

void accel_densmatr_oneQubitPauliChannel_subA(Qureg qureg, int qubit, qreal pI, qreal pX, qreal pY, qreal pZ);
void accel_densmatr_oneQubitPauliChannel_subB(Qureg qureg, int qubit, qreal pI, qreal pX, qreal pY, qreal pZ);

void accel_densmatr_oneQubitDamping_subA(Qureg qureg, int qubit, qreal prob);
void accel_densmatr_oneQubitDamping_subB(Qureg qureg, int qubit, qreal prob);
void accel_densmatr_oneQubitDamping_subC(Qureg qureg, int qubit, qreal prob);
void accel_densmatr_oneQubitDamping_subD(Qureg qureg, int qubit, qreal prob);



/*
 * PARTIAL TRACE
 */

void accel_densmatr_partialTrace_sub(Qureg inQureg, Qureg outQureg, vector<int> targs, vector<int> pairTargs);



#endif // ACCELERATOR_HPP