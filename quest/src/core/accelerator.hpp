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
 * This header additionally defines macros for cpu_subroutines.cpp 
 * and gpu_subroutines.cpp to use to concisely instantiate definitions
 * for their templated functions.
 * 
 * @author Tyson Jones
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
    private_INSTANTIATE(returntype, funcname, 0, args) \
    private_INSTANTIATE(returntype, funcname, 1, args) \
    private_INSTANTIATE(returntype, funcname, 2, args) \
    private_INSTANTIATE(returntype, funcname, 3, args) \
    private_INSTANTIATE(returntype, funcname, 4, args) \
    private_INSTANTIATE(returntype, funcname, 5, args) \
    private_INSTANTIATE(returntype, funcname,-1, args)

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
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 0, conj, args) \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 1, conj, args) \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 2, conj, args) \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 3, conj, args) \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 4, conj, args) \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 5, conj, args) \
    private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname,-1, conj, args)

#define private_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, numtargs, conj, args) \
    template returntype funcname <0, numtargs, conj>  args; \
    template returntype funcname <1, numtargs, conj>  args; \
    template returntype funcname <2, numtargs, conj>  args; \
    template returntype funcname <3, numtargs, conj>  args; \
    template returntype funcname <4, numtargs, conj>  args; \
    template returntype funcname <5, numtargs, conj>  args; \
    template returntype funcname <-1,numtargs, conj>  args;


#define INSTANTIATE_EXPONENTIABLE_CONJUGABLE_FUNC_OPTIMISED_FOR_NUM_CTRLS_AND_TARGS(returntype, funcname, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, true,  true,  args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, true,  false, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, false, true,  args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, false, false, args)

#define private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_outer(returntype, funcname, conj, haspower, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 0, conj, haspower, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 1, conj, haspower, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 2, conj, haspower, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 3, conj, haspower, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 4, conj, haspower, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, 5, conj, haspower, args) \
    private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname,-1, conj, haspower, args)

#define private_EXPONENTIABLE_CONJUGABLE_INSTANTIATE_inner(returntype, funcname, numtargs, conj, haspower, args) \
    template returntype funcname <0, numtargs, conj, haspower>  args; \
    template returntype funcname <1, numtargs, conj, haspower>  args; \
    template returntype funcname <2, numtargs, conj, haspower>  args; \
    template returntype funcname <3, numtargs, conj, haspower>  args; \
    template returntype funcname <4, numtargs, conj, haspower>  args; \
    template returntype funcname <5, numtargs, conj, haspower>  args; \
    template returntype funcname <-1,numtargs, conj, haspower>  args;



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
 * GETTERS 
 */

qcomp accel_statevec_getAmp_sub(Qureg qureg, qindex localInd);

void accel_statevec_getAmps_sub(qcomp* outAmps, Qureg qureg, qindex localStartInd, qindex numLocalAmps);


/*
 * SETTERS 
 */

void accel_statevec_setAmps_sub(qcomp* inAmps, Qureg qureg, qindex localStartInd, qindex numLocalAmps);

void accel_densmatr_setAmpsToPauliStrSum_sub(Qureg qureg, PauliStrSum sum);

void accel_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in);


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
 * DENSE MATRICES
 */

void accel_statevec_anyCtrlOneTargDenseMatr_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, CompMatr1 matr);
void accel_statevec_anyCtrlOneTargDenseMatr_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, qcomp fac0, qcomp fac1);

void accel_statevec_anyCtrlTwoTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, CompMatr2 matr);

void accel_statevec_anyCtrlAnyTargDenseMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, CompMatr matr, bool conj);


/*
 * DIAGONAL MATRICES
 */

void accel_statevec_anyCtrlOneTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ, DiagMatr1 matr);

void accel_statevec_anyCtrlTwoTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, int targ1, int targ2, DiagMatr2 matr);

void accel_statevec_anyCtrlAnyTargDiagMatr_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> targs, DiagMatr matr, qcomp exponent, bool conj);

void accel_statevec_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent);

void accel_densmatr_allTargDiagMatr_subA(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool multiplyOnly);
void accel_densmatr_allTargDiagMatr_subB(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool multiplyOnly);



/*
 * PAULI TENSOR AND GADGET
 */

void accel_statevector_anyCtrlAnyTargZOrPhaseGadget_sub(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> z, qcomp ampFac, qcomp pairAmpFac);

void accel_statevector_anyCtrlPauliTensorOrGadget_subA(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> x, vector<int> y, vector<int> z, qcomp ampFac, qcomp pairAmpFac);
void accel_statevector_anyCtrlPauliTensorOrGadget_subB(Qureg qureg, vector<int> ctrls, vector<int> ctrlStates, vector<int> x, vector<int> y, vector<int> z, qcomp ampFac, qcomp pairAmpFac, qindex bufferMaskXY);


/*
 * QUREG COMBINATION
 */

void accel_statevec_setQuregToSuperposition_sub(qcomp facOut, Qureg outQureg, qcomp fac1, Qureg inQureg1, qcomp fac2, Qureg inQureg2);

void accel_densmatr_mixQureg_subA(qreal outProb, Qureg out, qreal inProb, Qureg in);
void accel_densmatr_mixQureg_subB(qreal outProb, Qureg out, qreal inProb, Qureg in);
void accel_densmatr_mixQureg_subC(qreal outProb, Qureg out, qreal inProb);
void accel_densmatr_mixQureg_subD(qreal outProb, Qureg out, qreal inProb, Qureg in);


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


/*
 * PROBABILITIES
 */

qreal accel_statevec_calcTotalProb_sub(Qureg qureg);
qreal accel_densmatr_calcTotalProb_sub(Qureg qureg);

qreal accel_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes);
qreal accel_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes);

void accel_statevec_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits);
void accel_densmatr_calcProbsOfAllMultiQubitOutcomes_sub(qreal* outProbs, Qureg qureg, vector<int> qubits);


/*
 * INNER PRODUCTS
 */

qcomp accel_statevec_calcInnerProduct_sub(Qureg quregA, Qureg quregB);

qcomp accel_densmatr_calcFidelityWithPureState_sub(Qureg rho, Qureg psi, bool conj);

qreal accel_densmatr_calcHilbertSchmidtDistance_sub(Qureg quregA, Qureg quregB);


/*
 * EXPECTATION VALUES
 */

qreal accel_statevec_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> sufTargs);
qcomp accel_densmatr_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> allTargs);;

qcomp accel_statevec_calcExpecPauliStr_subA(Qureg qureg, vector<int> x, vector<int> y, vector<int> z);
qcomp accel_statevec_calcExpecPauliStr_subB(Qureg qureg, vector<int> x, vector<int> y, vector<int> z);
qcomp accel_densmatr_calcExpecPauliStr_sub (Qureg qureg, vector<int> x, vector<int> y, vector<int> z);

qcomp accel_statevec_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool useRealPow);
qcomp accel_densmatr_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, qcomp exponent, bool useRealPow);


/*
 * PROJECTORS 
 */

void accel_statevec_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob);
void accel_densmatr_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal prob);


/*
 * STATE INITIALISATION
 */

void accel_statevec_initUniformState_sub(Qureg qureg, qcomp amp);

void accel_statevec_initDebugState_sub(Qureg qureg);

void accel_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(Qureg qureg);


#endif // ACCELERATOR_HPP