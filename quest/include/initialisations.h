/** @file
 * API signatures for initialising Quregs into 
 * particular states. Note when a Qureg is GPU-
 * accelerated, these functions only update the
 * state in GPU memory; the CPU amps are unchanged.
 * 
 * @author Tyson Jones
 * 
 * @defgroup initialisations Initialisations
 * @ingroup api
 * @brief Functions for preparing Quregs in particular states.
 * @{
 */

#ifndef INITIALISATIONS_H
#define INITIALISATIONS_H

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/** 
 * @defgroup init_states States
 * @brief Functions for initialising Qureg into physical states.
 * @{
 */


/** Initialises @p qureg to the unnormalised all-zero-amplitude state.
 *
 * Every statevector amplitude, or every density-matrix element, is set to zero.
 * This is not a physical quantum state, but is useful as a blank
 * workspace before manually setting amplitudes.
 * 
 * > [!CAUTION]
 * > When @p qureg is GPU-accelerated, this function modifies only its GPU
 * > amplitudes (Qureg::gpuAmps), leaving its CPU amps (Qureg::cpuAmps)
 * > unchanged (like almost all QuEST operations). It is therefore necessary
 * > to follow this function with syncQuregFromGpu() in order to make 
 * > further, manual changes from the host side.
 * 
 * @equivalences
 * 
 * - This function is equivalent to (but much faster than) overwriting every
 *   amplitude to zero, _except_ that it does not modify Qureg::cpuAmps.
 *   ```cpp
     for (int i=0; i<qureg.numAmpsPerNode; i++)
         qureg.cpuAmps[i] = 0;
     syncQuregToGpu(qureg);

     // restore qureg.cpuAmps when !qureg.isGpuAccelerated
 *   ```
 *
 * @myexample
 * 
 * This function is useful for preparing sparse states, noting we must
 * explicitly copy the newly-zeroed amplitudes from GPU memory, when
 * @p qureg is GPU-accelerated (though such functions are always safe
 * to call).
 * 
 * ```cpp
   initBlankState(qureg);
   syncQuregFromGpu(qureg);
   
   // manually modify qureg.cpuAmps in some manner that wouldn't
   // be more sensible to perform with setQuregAmps, such as to a
   // uniform superposition of random basis states
   for (qindex i=0; i<500; i++) {
       qindex j = rand % qureg.numAmpsPerNode;
       qureg.cpuAmps[j] = 1;
   }

   syncQuregToGpu(qureg);
   setQuregToRenormalized();
 * ```
 *
 * @param[in,out] qureg  the Qureg to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * @see
 * - setQuregAmps()
 * - initZeroState()
 * - initClassicalState()
 * - initPlusState()
 * - initRandomPureState()
 * @author Tyson Jones
 */
void initBlankState(Qureg qureg);


/** Initialises @p qureg to the zero computational basis state.
 * 
 * > [!NOTE]
 * > Like most of QuEST's API, this function leaves Qureg::cpuAmps unchanged when
 * > @p qureg is GPU-accelerated, overwriting only the relevant GPU buffer Qureg::gpuAmps.
 * > See initBlankState() for more information.
 * 
 * @formulae
 *
 * Let @f$N@f$ be the number of qubits in @p qureg. 
 * 
 * - If @p qureg is a statevector, it is initialised to @f$\ket{0}^{\otimes N}@f$.
 * - If @p qureg is a density matrix, it is initialised to @f$\ket{0}\bra{0}^{\otimes N}@f$.'
 *
 * @equivalences
 * 
 * - The zero state is the first enumerated classical state.
 *   ```cpp
     initClassicalState(qureg, 0);
 *   ```
 * - The zero state has a zero amplitude everywhere except at the first index, which has one.
 *   The code below is equivalent to this function, _except_ Qureg::cpuAmps are also modified
 *   below, whereas initZeroState() leaves them unchanged when @p qureg is GPU-accelerated.
 *   ```cpp
     initBlankState(qureg);
     if (qureg.rank == 0)
         qureg.cpuAmps[0] = 1;
     syncQureg(qureg);
     
     // restore qureg.cpuAmps when !qureg.isGpuAccelerated
 *   ```
 *
 * @param[in,out] qureg  the Qureg to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * @author Tyson Jones
 */
void initZeroState(Qureg qureg);


/** Initialises @p qureg to the uniform plus state.
 * 
 * > [!NOTE]
 * > Like most of QuEST's API, this function leaves Qureg::cpuAmps unchanged when
 * > @p qureg is GPU-accelerated, overwriting only the relevant GPU buffer Qureg::gpuAmps.
 * > See initBlankState() for more information.
 * 
 * @formulae
 *
 * Let @f$N@f$ be the number of qubits in @p qureg. 
 * 
 * - If @p qureg is a statevector, it is initialised to
 *   @f[
        \begin{aligned}
        \ket{+}^{\otimes N} &= \left( \frac{1}{\sqrt{2}} \ket{0} + \frac{1}{\sqrt{2}} \ket{1} \right)^{\otimes N} \\
                            &= \frac{1}{\sqrt{2^N}} \{ 1, 1, \dots, 1 \}
        \end{aligned}
 *   @f]
 * - If @p qureg is a density matrix, it is initialised to 
 *   @f[
        \ket{+}\bra{+}^{\otimes N} = \frac{1}{2^N}
            \begin{pmatrix} 
            1 & 1 & \dots \\ 1 & \ddots \\ \vdots 
            \end{pmatrix}
 *   @f]
 *
 * @equivalences
 * 
 * - The plus state can also be produced by applying a Hadamard gate upon every zero-state qubit.
 *   ```cpp
     initZeroState(qureg);
     for (int i=0; i<qureg.numQubits; i++)
        applyHadamard(qureg, i);
 *   ``` 
 *
 * @param[in,out] qureg  the Qureg to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * @author Tyson Jones
 */
void initPlusState(Qureg qureg);


/** Initialises @p qureg to the state in statevector @p pure.
 * 
 * > [!NOTE]
 * > Like most of QuEST's API, this function leaves Qureg::cpuAmps unchanged when
 * > @p qureg is GPU-accelerated, overwriting only the relevant GPU buffer Qureg::gpuAmps.
 * > See initBlankState() for more information.
 * 
 * @formulae
 * 
 * Let @f$N@f$ be the number of qubits in @p qureg or @p pure, and let @f$\ket{\psi} = @f$ @p pure,
 * with amplitudes @f$\ket{\psi} = \sum_i \alpha_i \ket{i}@f$.
 * 
 * - If @p qureg is a statevector, it is overwritten by the state in @p pure.
 * - If @p qureg is a density matrix, it is initialised to 
 *   @f[
        \ket{\psi}\bra{\psi} = 
            \sum\limits_i\sum\limits_j \alpha_i \,\alpha_j^* \, \ket{i}\bra{j}
 *   @f]
 *
 * @equivalences
 * 
 * - When @p qureg is a statevector, this function is entirely equivalent to
 *   ```cpp
     setQuregToClone(qureg, pure);
 *   ``` 
 * - When @p qureg is a density matrix, this function is equivalent to
 *   ```cpp
     double prob = 1;
     setQuregToMixture(qureg, &prob, &pure, 1);
 *   ```
 *
 * @param[in,out] qureg  the Qureg to overwrite.
 * @param[in]     pure   the statevector pure state to copy.
 * @throws @validationerror
 * - if @p qureg or @p pure are uninitialised.
 * - if @p pure is not a statevector.
 * - if @p qureg and @p pure have incompatible dimensions or deployments.
 * @author Tyson Jones
 */
void initPureState(Qureg qureg, Qureg pure);


/** Initialises @p qureg to a computational basis state.
 * 
 * > [!NOTE]
 * > Like most of QuEST's API, this function leaves Qureg::cpuAmps unchanged when
 * > @p qureg is GPU-accelerated, overwriting only the relevant GPU buffer Qureg::gpuAmps.
 * > See initBlankState() for more information.
 * 
 * @formulae
 * 
 * Let @f$N@f$ be the number of qubits in @p qureg, and let @f$i=@f$ @p stateInd.
 * 
 * States are enumerated from @f$0@f$ to @f$2^N-1@f$, such that the bits of the indices
 * match the qubits of the corresponding basis states.
 * 
 * - If @p qureg is a statevector, it is initialised to @f$\ket{i}@f$.
 * - If @p qureg is a density matrix, it is initialised to @f$\ket{i}\bra{i}@f$.
 *
 * The bits of @f$i@f$ will match the qubit values of the resulting state in @p qureg, 
 *  where the zero-th qubit is the rightmost bit.
 *
 * @equivalences
 * 
 * - The resulting state contains zero for all amplitudes except that at global index @f$i@f$
 *   (when @p qureg is a statevector) or the @f$i@f$-th diagonal (when @p qureg is a
 *   density matrix).
 *   ```cpp
     initBlankState(qureg);

     // determine where the single, global amp to modify is located
     qindex numNewAmps = 1;
     qindex densityDim = 1 + (1 << qureg.numQubits);
     qindex globalAmpInd = stateInd * (qureg.isDensityMatrix? densityDim : 1);
     qindex localAmpInd = globalAmpInd % qureg.numAmpsPerNode;
     int rankContainingAmp = i / qureg.numAmpsPerNode;
     bool isAmpInThisNode = (rankContainingAmp == qureg.rank);

     // one node modifies 1 CPU amp
     if (isAmpInThisNode)
         qureg.cpuAmps[localAmpInd] = 1;
     
     // all nodes sync to GPU but only one node specifies a more than zero amps 
     syncSubQuregToGpu(qureg, i, numNewAmps * isAmpInThisNode);
 *   ```
 * - The resulting state can be (pointlessly slowly) produced by qubit flips from the
 *   zero state, according to the bits in @p stateInd.
 *   ```cpp
     initZeroState(qureg);
     for (int i=0; i<qureg.numQubits; i++)
         if ((stateInd >> i) & 1)
            applyPauliX(qureg, i);
 *   ```
 *
 * @param[in,out] qureg     the Qureg to overwrite.
 * @param[in]     stateInd  the computational basis-state index.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p stateInd is outside the computational basis of @p qureg.
 * @author Tyson Jones
 */
void initClassicalState(Qureg qureg, qindex stateInd);


/** Initialises @p qureg to the debug state.
 *
 * This is a non-physical, deterministic pattern useful for debugging.
 * The @f$j@f$-th local amplitude becomes 
 * @f[
    2j/10 + \iu(2j+1)/10,
 * @f]
 * even if @p qureg is a density matrix, in which case it is enumerated
 * column-major.
 * 
 * > [!CAUTION]
 * > When @p qureg is GPU-accelerated, this function modifies only its GPU
 * > amplitudes (Qureg::gpuAmps), leaving its CPU amps (Qureg::cpuAmps)
 * > unchanged (like almost all QuEST operations). It is therefore necessary
 * > to follow this function with syncQuregFromGpu() in order to make 
 * > further, manual changes from the host side. See initBlankState() for more
 * > information.
 *
 * @myexample
 * 
 * ```cpp
   Qureg qureg = createQureg(3);
   initDebugState(qureg);
   reportQureg(qureg);
 * ```
 * ```text
    Qureg (3 qubit statevector, 8 qcomps, 232 bytes):
        0.1i      |0⟩
        0.2+0.3i  |1⟩
        0.4+0.5i  |2⟩
        0.6+0.7i  |3⟩
        0.8+0.9i  |4⟩
        1+1.1i    |5⟩
        1.2+1.3i  |6⟩
        1.4+1.5i  |7⟩
 * ```
 *
 * @param[in,out] qureg  the Qureg to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * @author Tyson Jones
 */
void initDebugState(Qureg qureg);


/** Initialises @p qureg from the statevector amplitudes in @p amps.
 * 
 * @formulae
 * 
 * Let @f$N@f$ be the number of qubits in @p qureg, and let @f$\alpha_i@f$ be
 * the amplitude `amps[i]`. Array @p amps must be length @f$2^N@f$.
 * 
 * - If @p qureg is a statevector, its amplitudes are overwritten by @p amps, to become
 *   @f[
        \sum\limits_i \alpha_i \ket{i}.
 *   @f]
 * 
 * - If @p qureg is a density matrix, it is initialised to the pure state @f$\ket{\psi}@f$ 
 *   encoded by @p amps, i.e.
 *   @f[
        \ket{\psi}\bra{\psi} = 
            \sum\limits_i\sum\limits_j \alpha_i \,\alpha_j^* \, \ket{i}\bra{j}
 *   @f]
 *
 * There is no need for @p amps to be normalised, although @p qureg will otherwise be left
 * in an unnormalised, non-physical state.
 *
 * @param[in,out] qureg  the Qureg to overwrite.
 * @param[in]     amps   an array of @f$2^N@f$ pure-state amplitudes.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * @throws seg-fault
 * - if @p amps has fewer than @f$2^N@f$ elements.
 * @author Tyson Jones
 */
void initArbitraryPureState(Qureg qureg, qcomp* amps);


/** Initialises @p qureg (a statevector or density matrix) to a pure state with 
 * uniformly random amplitudes.
 * 
 * The resulting state is normalised, with basis state probabilities sampled
 * from a chi-squared variate, as described
 * [here](https://sumeetkhatri.com/wp-content/uploads/2020/05/random_pure_states.pdf). 
 * 
 * @param[in,out] qureg  the Qureg to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * @see 
 * - initRandomMixedState()
 * @author Tyson Jones
 */
void initRandomPureState(Qureg qureg);


/** Initialises a density matrix to a mixture of uniformly random pure states.
 *
 * The resulting density matrix is the equally weighted mixture of @p numPureStates
 * independently sampled random pure states, each sampled as per initRandomPureState().
 *
 * @formulae
 * 
 * Let @f$n=@f$ @p numPureStates, and let @f$\ket{\psi_i}@f$ be a random pure
 * state with number of qubits as @p qureg.
 * 
 * This function overwrites @p qureg to
 * @f[
 *      \sum\limits_i^n \frac{1}{n} \ket{\psi_i}\bra{\psi_i}.
 * @f]
 *
 * @param[in,out] qureg          the density matrix to overwrite.
 * @param[in]     numPureStates  the number of random pure states in the mixture.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p qureg is not a density matrix.
 * - if @p numPureStates is invalid.
 * @see 
 * - initRandomPureState()
 * @author Tyson Jones
 */
void initRandomMixedState(Qureg qureg, qindex numPureStates);


/** @} */



/** 
 * @defgroup init_amps Amplitudes
 * @brief Functions for overwriting Qureg amplitudes.
 * @{
 */


/** Overwrites a contiguous range of statevector amplitudes.
 * 
 * - Amplitudes outside the given range are unchanged.
 * - There is no validation nor requirement that the new amplitudes,
 *   together with the remaining original amplitudes, produce a validly
 *   normalised state. Normalization can be re-established with a subsequent
 *   call to setQuregToRenormalized().
 * - When @p qureg is distributed, @p startInd and @p numAmps are treated
 *   _globally_ and use of this function is ergo agnostic to distribution.
 *   Therefore, every process should contain identical @p amps, although
 *   only elements which fall within a process' statevector partition will
 *   be consulted by a particular process.
 * - When @p qureg is GPU-accelerated, only its GPU amplitudes are updated.
 *  
 * The equivalent function for a density matrix is setDensityQuregAmps().
 * 
 * @formulae
 * 
 * Let @f$\svpsi=@f$ @p qureg with @f$N@f$ qubits, and with @f$i@f$-th global amplitude @f$\alpha_i@f$.
 * Let @f$s=@f$ @p startInd, @f$n=@f$ @p numAmps, and let @f$\beta_j@f$ be the @f$j@f$-th element of @p amps. 
 * 
 * This function overwrites @p qureg from @f$\svpsi=\sum_{i=0}^{2^N-1} \alpha_i \ket{i}@f$ to
 * @f[
        \svpsi \rightarrow 
            \sum\limits_{i=0}^{s-1} \alpha_i \ket{i} +
            \sum\limits_{j=0}^{n-1} \beta_j \ket{j + s} +
            \sum\limits_{i=s+n}^{2^N-1} \alpha_i \ket{i}
 * @f]
 * where amplitudes at global indices in @f$[s,s+n)@f$ have been modified. Expressed as a row-vector,
 * @f[
        \svpsi = \begin{pmatrix} \alpha_0 & \alpha_1 &  \dots & \alpha_{2^N-1} \end{pmatrix}
 * @f]
 * is modified to become
 * @f[
        \svpsi \rightarrow \begin{pmatrix} 
            \alpha_0 & \alpha_1 & \dots & \alpha_{s-1} & 
            \beta_0 & \beta_1 & \dots & \beta_{n-1} &
            \alpha_{s + n} & \dots & \alpha_{2^N-1}
        \end{pmatrix}.
 * @f]
 *
 * @constraints
 * 
 * - Argument @p qureg must be a statevector, and ergo compatible with a 1D range.
 *   Density matrices can be overwritten at a 2D range with setDensityQuregAmps(),
 *   or with a 1D contiguous range when flattening the density matrix column-major
 *   with setDensityQuregFlatAmps(). Alternatively, setQuregAmps() can be called
 *   with validation disabled via setQuESTValidationOff(), accepting density matrices,
 *   and behaving identically to setDensityQuregFlatAmps().
 * 
 * @equivalences
 * 
 * - When @p qureg is **_not_** distributed, this function is equivalent to (but
 *   much faster than) manual modification of the CPU elements, followed by a copy
 *   to GPU (_except_ that this function does not modify Qureg::cpuAmps when @p qureg
 *   is not GPU-accelerated).
 *   ```cpp
     for (qindex i=0; i<numAmps; i++)
         qureg.cpuAmps[i + startInd] = amps[i];
     syncSubQuregToGpu(qureg, startInd, numAmps); 
     // beware, syncQuregToGpu() would copy over stale, unmodified CPU amps

     // restore qureg.cpuAmps when !qureg.isGpuAccelerated
 *   ```
 * - When @p qureg _is_ distributed, the logic is complicated by the specified global
 *   range of amplitudes overlapping some, none or all of a node's partition.
 * 
 *   ```cpp
     qcomp* amps[numAmps] = // global

     qindex dim = qureg.numAmpsPerNode;
     qindex endInd = startInd + numAmps;

     qindex nodeStartInd = (qureg.rank    ) * dim;
     qindex nodeEndInd   = (qureg.rank + 1) * dim;
     bool nodeContainsAmps = (startInd < nodeEndInd) && (endInd > nodeStartInd);

     qindex localStartInd = (startInd < nodeStartInd)? 0 : startInd % dim;
     qindex localEndInd = (endInd > nodeEndInd)? dim : endInd % dim;
     qindex numLocalAmps = nodeContainsAmps * (localEndInd - localStartInd);

     qindex nodeOffset = nodeStartInd + localStartInd - startInd;
     for (qindex i=0; i<numLocalAmps; i++)
         qureg.cpuAmps[localStartInd + i] = amps[nodeOffset + i]

     syncSubQuregToGpu(qureg, localStartInd, numLocalAmps);

     // restore qureg.cpuAmps when !qureg.isGpuAccelerated
 *   ```
 * 
 * @myexample
 * 
 * - When @p numAmps is sufficiently small such that the array @p amps can
 *   fit onto every distributed node, this function can be used in a manner
 *   totally agnostic to distribution and/or @p qureg deployments.
 * 
 *   ```cpp
     Qureg qureg = createQureg(35);
     initBlankState(qureg);

     qcomp amps[1000] = { ... };
     setQuregAmps(qureg, 300000000, amps, 1000);
 *   ```
 * - When @p numAmps is large, one can avoid the superfluous storing of all @p amps
 *   simultaneously, by repeatedly calling setQuregAmps(), each time passing a tractable
 *   sub-range, regardless of how @p qureg is distributed.
 *   ```cpp
     Qureg qureg = createQureg(35);
     initBlankState(qureg);

     // global range
     const qindex startInd = 1234567;
     const qindex totalNumAmps = 1000000000; // 16  GB worth of double-prec qcomp

     // local memory budget
     const qindex batchSize = 10000000;   // 160 MB worth
     qcomp amps[batchSize];

     int numBatches = totalNumAmps / batchSize; // divides evenly here for simplicity
     
     for (int batchInd=0; batchInd<numBatches; batchInd++) {

         // update amps, such that amps[i] is the desired amplitude
         // for global index (startInd + batchInd * batchSize)
         ...

         setQuregAmps(qureg, startInd + batchInd * batchSize, amps, batchSize);
     }
 *   ```
 *
 * @param[in,out] qureg     the statevector to modify.
 * @param[in]     startInd  the first global computational-basis index to overwrite.
 * @param[in]     amps      an array of @p numAmps amplitudes.
 * @param[in]     numAmps   the total number of amplitudes to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p qureg is not a statevector.
 * - if @p startInd and @p numAmps describes a range outside @p qureg.
 * @see
 * - setDensityQuregAmps()
 * - setDensityQuregFlatAmps()
 * - setQuregToWeightedSum()
 * - setQuregToRenormalized()
 * @author Tyson Jones
 */
void setQuregAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);


/** Overwrites a rectangular block of density-matrix amplitudes.
 *
 * - Amplitudes outside the given block are unchanged.
 * - There is no validation nor requirement that the new amplitudes,
 *   together with the remaining original amplitudes, produce a validly
 *   normalised density matrix. Normalization can be re-established with a
 *   subsequent call to setQuregToRenormalized().
 * - When @p qureg is distributed, @p startRow, @p startCol, @p numRows and
 *   @p numCols are treated _globally_ and use of this function is ergo
 *   agnostic to distribution. Therefore, every process should contain
 *   identical @p amps.
 * - When @p qureg is GPU-accelerated, only its GPU amplitudes are updated.
 *
 * The equivalent function for a statevector is setQuregAmps().
 *
 * @formulae
 *
 * Let @f$\dmrho=@f$ @p qureg with @f$N@f$ qubits, and with @f$(r,c)@f$-th global amplitude
 * @f$\alpha_{r,c}@f$. Let @f$s_r=@f$ @p startRow, @f$s_c=@f$ @p startCol, @f$n_r=@f$
 * @p numRows and @f$n_c=@f$ @p numCols, and let @f$\beta_{j,k}@f$ be the @f$(j,k)@f$-th
 * element of @p amps, i.e. `amps[j][k]`.
 *
 * This function overwrites @p qureg from
 * @f[
        \dmrho = \sum\limits_{r=0}^{2^N-1} \sum\limits_{c=0}^{2^N-1}
            \alpha_{r,c} \ket{r}\bra{c}
 * @f]
 * by modifying only the global rows @f$[s_r,s_r+n_r)@f$ and columns @f$[s_c,s_c+n_c)@f$,
 * such that
 * @f[
        \alpha_{s_r+j,\,s_c+k} \rightarrow \beta_{j,k}
        \quad\quad
        \forall \; j \in [0,n_r), \; k \in [0,n_c).
 * @f]
 * Expressed as a matrix,
 * @f[
        \dmrho = 
        \begin{pmatrix}
            \alpha_{0,0} & \alpha_{0,1} & \dots & \alpha_{0,2^N-1} \\
            \alpha_{1,0} & \alpha_{1,1} &  \\
            \vdots & & \ddots \\
            \alpha_{2^N-1,0} & & & \alpha_{2^N-1,2^N-1}
        \end{pmatrix},
 * @f]
 * the state is modified to contain the sub-matrix below. Grey dots indicate @f$\alpha_{ij}@f$ above.
 * @f[
        \def\x{{\color{gray}\circ}}
        \dmrho \rightarrow
        \begin{array}{c@{\;}c}
            & \hspace{0.5em}
              \overset{\scriptstyle [s_c,\,s_c+n_c)}{\overline{\hspace{10.5em}}}
              \hspace{-0.5em} \\[1ex]
            \lower1.0em\hbox{$
                \scriptstyle [s_r,\,s_r+n_r) \quad
                \left\{\vphantom{\begin{matrix}
                    \beta_{0,0} \\
                    \beta_{1,0} \\
                    \vdots \\
                    \beta_{n_r-1,0}
                \end{matrix}}\right.
            $}
            &
            \begin{pmatrix}
                \x & \x & \x & \x & \x & \x & \x \\
                \x & \x & \x & \x & \x & \x & \x \\
                \x & \x & \beta_{0,0} & \beta_{0,1} & \cdots & \beta_{0,n_c-1} & \x \\
                \x & \x & \beta_{1,0} & \beta_{1,1} & \cdots & \beta_{1,n_c-1} & \x \\
                \x & \x & \vdots & \vdots & \ddots & \vdots & \x \\
                \x & \x & \beta_{n_r-1,0} & \beta_{n_r-1,1} & \cdots & \beta_{n_r-1,n_c-1} & \x \\
                \x & \x & \x & \x & \x & \x & \x
            \end{pmatrix}
        \end{array}
 * @f]
 *
 * @constraints
 *
 * - Argument @p qureg must be a density matrix, and ergo compatible with a 2D range.
 *   Statevectors can be overwritten at a 1D range with setQuregAmps(). Density
 *   matrices can also be overwritten with a 1D contiguous range when flattening
 *   the density matrix column-major with setDensityQuregFlatAmps().
 *
 * @equivalences
 *
 * - When @p qureg is **_not_** distributed, this function is equivalent to manual
 *   modification of the CPU elements, followed by copies to GPU (_except_ that this
 *   function does not modify Qureg::cpuAmps when @p qureg is not GPU-accelerated).
 *   ```cpp
     for (qindex c=0; c<numCols; c++) {
         qindex flatInd = (startCol + c) * (1LL << qureg.numQubits) + startRow;
         for (qindex r=0; r<numRows; r++)
             qureg.cpuAmps[flatInd + r] = amps[r][c];
         syncSubQuregToGpu(qureg, flatInd, numRows);
     }
     // beware, syncQuregToGpu() would copy over stale, unmodified CPU amps

     // restore qureg.cpuAmps when !qureg.isGpuAccelerated
 *   ```
 * - When @p qureg _is_ distributed, the logic is complicated by each specified
 *   global column range overlapping some, none or all of a node's partition.
 *   It follows the same pattern as demonstrated in setQuregAmps(), through
 *   column-wise linearisation of the density matrix.
 * - When @p amps are within a single column (`numCols==1`), or span multiple
 *   _full_ columns (`numRows==(1<<qureg.numQubits)`), this function becomes
 *   equivalent to calling setDensityQuregFlatAmps(), passing @p amps as a
 *   column-flattened 1D array.
 *
 * @myexample
 * 
 * - When @p numRows and @p numCols are sufficiently small such that the matrix
 *   @p amps can fit onto every distributed node, this function can be used in a
 *   manner totally agnostic to distribution and/or @p qureg deployments.
 *   In C++, an overload accepts @p amps as nested `std::vector<qcomp>`:
 *   ```cpp
     Qureg qureg = createDensityQureg(35);

     std::vector<std::vector<qcomp>> amps = {
         {1, 2, 3},
         {4, 5, 6}};

     setDensityQuregAmps(qureg, startRow, startCol, amps, 2, 3);
 *   ```
 *
 *   In C, @p amps must be a double pointer (and alas not an array):
 *   ```cpp
     Qureg qureg = createDensityQureg(35);
     initBlankState(qureg);

     qcomp ampsArr[2][3] = {
         {1, 2, 3},
         {4, 5, 6}};
     qcomp* ampsPtr[] = {ampsArr[0], ampsArr[1]};

     setDensityQuregAmps(qureg, startRow, startCol, ampsPtr, 2, 3);
 *   ```
 * - When @p numRows or @p numCols is large, one can avoid the superfluous storing
 *   of all @p amps simultaneously, by repeatedly calling setDensityQuregAmps(),
 *   each time passing a tractable sub-block, regardless of how @p qureg is
 *   distributed.
 *   ```cpp
     Qureg qureg = createDensityQureg(35);
     initBlankState(qureg);

     // global block
     const qindex startRow = 1234567;
     const qindex startCol = 7654321;
     const qindex totalNumRows = 1000000; // divides evenly into batches below
     const qindex totalNumCols = 1000000;

     // local memory budget
     const qindex batchNumRows = 100;
     const qindex batchNumCols = 200;

     // memory for a single batch
     qcomp** ampsBatch = ... // 2D malloc of size (batchNumRows, batchNumCols)

     for (qindex row=0; row<totalNumRows; row+=batchNumRows) {
         for (qindex col=0; col<totalNumCols; col+=batchNumCols) {

             // populate this two-dimensional batch
             for (qindex r=0; r<batchNumRows; r++)
                 for (qindex c=0; c<batchNumCols; c++)
                     amps[r][c] = ... 

             setDensityQuregAmps(
                 qureg, startRow + row, startCol + col,
                 ampsBatch, batchNumRows, batchNumCols);
         }
     }
 *   ```
 *
 * @param[in,out] qureg     the density matrix to modify.
 * @param[in]     startRow  the first global density-matrix row index to overwrite.
 * @param[in]     startCol  the first global density-matrix column index to overwrite.
 * @param[in]     amps      a @p numRows by @p numCols matrix of amplitudes, as row-major nested pointers.
 * @param[in]     numRows   the total number of rows to overwrite.
 * @param[in]     numCols   the total number of columns to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if @p qureg is not a density matrix.
 * - if @p startRow, @p startCol, @p numRows and @p numCols describe a block outside @p qureg.
 * @see
 * - setQuregAmps()
 * - setDensityQuregFlatAmps()
 * - setQuregToWeightedSum()
 * - setQuregToRenormalized()
 * @notyetvalidated
 * @author Tyson Jones
 */
void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols);


/** Overwrites a contiguous range of density-matrix amplitudes, indexing @p qureg
 * as a column-linearised form.
 *
 * @formulae
 * 
 * Let @f$ \dmrho = @f$ @p qureg contain @f$N@f$ qubits, with amplitudes @f$ \alpha_{ij} @f$.
 * @f[
        \dmrho = \sum\limits_i^{2^N} \sum\limits_j^{2^N} \alpha_{ij} \ket{i}\bra{j}.
 * @f] 
 * Internally, this matrix of dimension @f$ 2^N \times 2^N @f$ is stored in a vectorised
 * form @f$ \ket{\rho} @f$ of dimension @f$ 2^{2N} \times 1 @f$, which concatenates the columns
 * of @f$ \dmrho @f$.
 * @f[  
        \begin{aligned}
        \ket{\rho} &= \sum\limits_i^{2^N} \sum\limits_j^{2^N} \alpha_{ij} \ket{j} \ket{i} \\
                   &= \sum\limits_k^{2^{2N}} \beta_k \ket{k}
        \end{aligned}
 * @f]
 * Let @f$s = @f$ @p startInd, @f$n = @f$ @p numAmps, and @f$\gamma_i = @f$ `amps[i]`.
 * This function overwrites amplitudes @f$\beta_k : s \le k < s + n @f$ with @f$\gamma_k@f$.
 * 
 * @equivalences
 * 
 * - This is equivalent to calling setQuregAmps() upon @p qureg, with identical parameters,
 *   when treating the column-wise linearised @p qureg as a statevector.
 * 
 * @param[in,out] qureg     the density matrix to modify.
 * @param[in]     startInd  the first flattened, global index to overwrite.
 * @param[in]     amps      an array of @p numAmps amplitudes.
 * @param[in]     numAmps   the number of flattened amplitudes to overwrite.
 * @throws @validationerror
 * - if @p qureg is uninitialised or is not a density matrix.
 * - if @p startInd or @p numAmps describes a range outside the flattened matrix.
 * @notyetvalidated
 * @author Tyson Jones
 */
void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, qcomp* amps, qindex numAmps);


/// @notyetdoced
/// @notyettested
void setQuregToClone(Qureg outQureg, Qureg inQureg);


/// @notyetdoced
/// @notyettested
void setQuregToWeightedSum(Qureg out, qcomp* coeffs, Qureg* in, int numIn);


/// @notyetdoced
/// @notyettested
void setQuregToMixture(Qureg out, qreal* probs, Qureg* in, int numIn);


/// @notyetdoced
/// @notyetvalidated
qreal setQuregToRenormalized(Qureg qureg);


/// @notyetdoced
/// @notyetvalidated
void setQuregToPauliStrSum(Qureg qureg, PauliStrSum sum);


/// @notyetdoced
/// @notyettested
void setQuregToPartialTrace(Qureg out, Qureg in, int* traceOutQubits, int numTraceQubits);


/// @notyetdoced
/// @notyettested
void setQuregToReducedDensityMatrix(Qureg out, Qureg in, int* retainQubits, int numRetainQubits);


/** @} */



// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We 
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

#include <vector>


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
/// @see setQuregAmps()
void setQuregAmps(Qureg qureg, qindex startInd, std::vector<qcomp> amps);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
/// @see setDensityQuregAmps()
void setDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, std::vector<std::vector<qcomp>> amps);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
/// @see setDensityQuregFlatAmps()
void setDensityQuregFlatAmps(Qureg qureg, qindex startInd, std::vector<qcomp> amps);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
/// @see setQuregToPartialTrace()
void setQuregToPartialTrace(Qureg out, Qureg in, std::vector<int> traceOutQubits);


/// @ingroup init_amps
/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cpponly
/// @see setQuregToReducedDensityMatrix()
void setQuregToReducedDensityMatrix(Qureg out, Qureg in, std::vector<int> retainQubits);


/// @ingroup init_amps
/// @notyetdoced
/// @cpponly
/// @see setQuregToWeightedSum()
void setQuregToWeightedSum(Qureg out, std::vector<qcomp> coeffs, std::vector<Qureg> in);


/// @ingroup init_amps
/// @notyetdoced
/// @cpponly
/// @see setQuregToMixture()
void setQuregToMixture(Qureg out, std::vector<qreal> probs, std::vector<Qureg> in);


#endif // __cplusplus



#endif // INITIALISATIONS_H

/** @} */ // (end file-wide doxygen defgroup)
