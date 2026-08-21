/** @file
 * Definitions of all dense and diagonal matrices, their getters and setters,
 * as well as their reporting utilities. Note that Kraus maps are treated in
 * a bespoke file (channels.h).
 * 
 * This file uses extensive preprocessor trickery to achieve overloaded,
 * platform agnostic, C and C++ compatible, precision agnostic, getters 
 * and setters of complex matrices. All macros herein expand to single-line 
 * definitions, for safety. Some intendedly private functions are necessarily
 * exposed here to the user, and are prefixed with an underscore.
 * 
 * Note too that designated initialisers are generally avoided, since while
 * C99/C11 compatible, they are not explicitly supported by C++14 (only from
 * C++20) except through compiler extensions
 * 
 * @author Tyson Jones
 * @author Richard Meister (aided in design)
 * @author Erich Essmann (aided in design, patched on MSVC)
 * 
 * @defgroup matrices Matrices
 * @ingroup api
 * @brief Data structures for representing operator matrices.
 * @{
 */

#ifndef MATRICES_H
#define MATRICES_H

#include "quest/include/types.h"
#include "quest/include/paulis.h"

// C++ gets vector initialiser overloads, whereas C gets a macro
#ifdef __cplusplus
    #include <vector>
#endif



/*
 * unlike some other headers, we here intermix the C and C++-only
 * signatures, grouping them semantically & by their doc groups
 */



/** 
 * @defgroup matrices_structs Structs
 * @brief Data structures for representing operator matrices.
 * @{
 */


/*
 * DENSE MATRIX STRUCTS
 *
 * which are visible to both C and C++, where qcomp resolves
 * to the native complex type. These are not de-mangled because
 * C++ structs are already C compatible.
 * 
 * The compile-time sized structs have field 'elems', while
 * dynamic-sized structs have separate 'cpuElems' and 
 * 'gpuElemsFlat', for persistent GPU allocation, and ergo need 
 * syncing. Note 'gpuElemsFlat' is always 1D (hence the name).
 */


/** A one-qubit, @f$2\times 2@f$, dense complex matrix.
 * 
 * > [!TIP]
 * > This matrix is stored in stack memory, and so never needs explicit destruction,
 * > and like all QuEST types, is safe to pass and return to/from functions by-value.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.cpp) 
 * examples of initialising CompMatr1.
 * 
 * @see
 * - getCompMatr1()
 * - getInlineCompMatr1()
 * - reportCompMatr1()
 * - applyCompMatr1()
 */
typedef struct {

    /** The number of qubits upon which the matrix can act, which is always @c =1.
     */
    int numQubits;

    /** The matrix dimension, which is always @c =2.
     */
    qindex numRows;

    /** The 2D matrix elements. These are safe to update at any time.
     */
    qcomp elems[2][2];

} CompMatr1;


/** A two-qubit, @f$4\times 4@f$, dense complex matrix.
 * 
 * > [!TIP]
 * > This matrix is stored in stack memory, and so never needs explicit destruction,
 * > and like all QuEST types, is safe to pass and return to/from functions by-value.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.cpp) 
 * examples of initialising CompMatr2.
 * 
 * @see
 * - getCompMatr2()
 * - getInlineCompMatr2()
 * - reportCompMatr2()
 * - applyCompMatr2()
 */
typedef struct {

    /** The number of qubits upon which the matrix can act, which is always @c =2.
     */
    int numQubits;

    /** The matrix dimension, which is always @c =4.
     */
    qindex numRows;

    /** The 2D matrix elements. These are safe to update at any time.
     */
    qcomp elems[4][4];

} CompMatr2;


/** A variable-size, @f$n@f$-qubit, @f$2^n \times 2^n@f$, dense complex matrix.
 * 
 * Unlike CompMatr1 and CompMatr2, the elements of this matrix are stored in
 * heap memory. But like all QuEST structs, it is safe to copy, and ergo to pass
 * to functions by-value, or be returned from them. However, its destructor
 * must only ever be called upon one such copy.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.cpp) 
 * examples of initialising CompMatr.
 * 
 * @see
 * - createCompMatr()
 * - [createInlineCompMatr()](https://quest-kit.github.io/QuEST/group__matrices__create.html#gad0ac51cb08a5d0b4e0878a5fde7ca893)
 * - [setInlineCompMatr()](https://quest-kit.github.io/QuEST/group__matrices__setters.html#ga00993b6b571985fd04697e487a561517)
 * - setCompMatr()
 * - syncCompMatr()
 * - reportCompMatr()
 * - applyCompMatr()
 * - destroyCompMatr()
 */
typedef struct {

    // beware that CompMatr instances are sometimes 'spoofed' inside localiser.cpp,
    // which will set the fields from other object instances (like a SuperOp). As
    // such, adding new fields to this struct may require updating these spoofers.

    /** The number of qubits upon which the matrix can act.
     * 
     * The total memory costs of the matrix scale exponentially with the
     * number of qubits; an @f$n@f$-qubit matrix contains @f$4^n@f$
     * complex elements.
     */
    int numQubits;

    /** The dimension of the matrix, which is a square matrix, and so
     * equals both the number of rows and columns.
     * 
     * Letting @f$n=@f$ #numQubits, then #numRows @f$=2^n@f$.
     */
    qindex numRows;

    /** Whether the matrix is (within validation epsilon tolerance) unitary (@c =1), or known to be
     * non-unitary (@c =0), or whether it is unknown (@c =-1).
     * 
     * This is a heap pointer to a persistent flag which is initially @c =-1 at matrix creation, and
     * is only ever consulted and/or updated by input validation within functions like applyCompMatr(),
     * when such validation is enabled. Calling setCompMatr() or syncCompMatr() restores #isApproxUnitary 
     * to @c -1.
     * 
     * Unitarity is measured approximately, with reference to the validation
     * epsilon as modified with setQuESTValidationEpsilon(). This flag will never be updated from
     * @c =-1 when validation is disabled. Like all epsilon-dependent fields, it is restored
     * to @c =-1 automatically whenever setQuESTValidationEpsilon() or setQuESTValidationEpsilonToDefault()
     * are called.
     * 
     * @see
     * - setQuESTValidationEpsilon()
     */
    int* isApproxUnitary;

    /** Whether the matrix is (within validation epsilon tolerance) Hermitian (@c =1), or known to
     * be non-Hermitian (@c =0), or whether it is unknown (@c =-1).
     * 
     * This is a heap pointer to a persistent flag which is initially @c =-1 at matrix creation, and
     * is only ever consulted and/or updated by input validation within functions which expect
     * Hermiticity. Calling setCompMatr() or syncCompMatr() restores 
     * #isApproxHermitian to @c -1.
     * 
     * > [!NOTE]
     * > There are presently _no_ QuEST functions which expect Hermiticity of a CompMatr, 
     * > and so this field is presently never modified from its initial value of @c =-1.
     * 
     * Hermiticity is measured approximately, with reference to the validation
     * epsilon as modified with setQuESTValidationEpsilon(). This flag will never be updated from
     * @c =-1 when validation is disabled. Like all epsilon-dependent fields, it is restored
     * to @c =-1 automatically whenever setQuESTValidationEpsilon() or setQuESTValidationEpsilonToDefault()
     * are called.
     * 
     * @see
     * - setQuESTValidationEpsilon()
     */
    int* isApproxHermitian; // @todo currently unused (relevant to not-yet-implemented calc-expec-val)

    /** Whether the matrix elements were ever synchronised.
     * 
     * This is a heap pointer to a persistent flag which is initially @c 0 at matrix creation,
     * but which is permanently overwritten to @c 1 when synchronisation is performed, such as
     * via syncCompMatr() or setCompMatr(). The flag indicates whether the matrix
     * elements have been initialised (and when QuEST is GPU-accelerated, whether they have been
     * copied to GPU device memory), and ergo whether it is valid to pass the matrix to a
     * simulation function like applyCompMatr().
     * 
     * Note this flag can only indicate whether the matrix has _ever_ been synced; it cannot be
     * used to detect whether manual modification of #cpuElems made after an initial sync have been
     * re-synced, as required for correct behaviour in GPU mode.
     * 
     * @see
     * - syncCompMatr()
     */
    int* wasGpuSynced;

    /** The 2D matrix elements, stored in CPU host memory.
     * 
     * It is safest to modify this matrix through setCompMatr(), but direct modification
     * is possible; the matrix element of the `r`-th row and `c`-th column is stored
     * at `cpuElems[r][c]`.
     * 
     * > [!IMPORTANT]
     * > It is _critical_ to call syncCompMatr() after direct modification of
     * > #cpuElems in order to update persistent matrix properties,
     * > such as its data in GPU device memory (even when not running in
     * > GPU-accelerated mode).
     * 
     * The field #cpuElems merely aliases the 1D #cpuElemsFlat field, such that
     * modifications of #cpuElems also updates #cpuElemsFlat.
     * 
     * @see
     * - setCompMatr()
     * - [setInlineCompMatr()](https://quest-kit.github.io/QuEST/group__matrices__setters.html#ga00993b6b571985fd04697e487a561517)
     * - syncCompMatr()
     */
    qcomp** cpuElems;

    /** A 1D row-major form of #cpuElems.
     * 
     * Modification of the matrix elements should be done through #cpuElems
     * for mathematical clarity, though this 1D contiguous form may be convenient
     * when performing copying.
     * 
     * > [!IMPORTANT]
     * > It is _critical_ to call syncCompMatr() after direct modification of
     * > #cpuElemsFlat in order to update persistent matrix properties,
     * > such as its data in GPU device memory (even when not running in
     * > GPU-accelerated mode).
     * 
     * @see
     * - setCompMatr()
     * - [setInlineCompMatr()](https://quest-kit.github.io/QuEST/group__matrices__setters.html#ga00993b6b571985fd04697e487a561517)
     * - syncCompMatr()
     */
    qcomp* cpuElemsFlat;

    /** The matrix elements, stored in GPU device memory, in a 1D row-major form.
     * 
     * This is a copy of #cpuElemsFlat consulted by QuEST's GPU backend and should
     * _never_ be modified directly. Instead, it is updated by calling syncCompMatr()
     * after modifying #cpuElems or #cpuElemsFlat, which is performed automatically by
     * API functions like setCompMatr(). In this way, #gpuElemsFlat and #cpuElemsFlat
     * should never be out of sync with one another.
     * 
     * Within a GPU-enabled QuEST environment, every CompMatr allocates #gpuElemsFlat
     * in GPU memory, even if never ultimately consulted from the GPU backend.
     */
    qcomp* gpuElemsFlat;

} CompMatr;


/*
 * DIAGONAL MATRIX STRUCTS
 *
 * with all the same nuances as the CompMatr structs described above. 
 */


/** A one-qubit, @f$2\times 2@f$, diagonal complex matrix, stored as a length-2 array.
 * 
 * > [!TIP]
 * > This matrix is stored in stack memory, and so never needs explicit destruction,
 * > and like all QuEST types, is safe to pass and return to/from functions by-value.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.cpp) 
 * examples of initialising DiagMatr1.
 * 
 * @see
 * - getDiagMatr1()
 * - getInlineDiagMatr1()
 * - reportDiagMatr1()
 * - applyDiagMatr1()
 */
typedef struct {

    /** The number of qubits upon which the matrix can act, which is always @c =1.
     */
    int numQubits;

    /** The matrix dimension, equal to the number of contained diagonal elements, which is always @c =2.
     */
    qindex numElems;

    /** The diagonal elements of the matrix, stored as a flat array, in CPU stack memory.
     * These are safe to directly modify at any time, to update the matrix.
     */
    qcomp elems[2];

} DiagMatr1;


/** A two-qubit, @f$4\times 4@f$, diagonal complex matrix, stored as a length-4 array.
 * 
 * > [!TIP]
 * > This matrix is stored in stack memory, and so never needs explicit destruction,
 * > and like all QuEST types, is safe to pass and return to/from functions by-value.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.cpp) 
 * examples of initialising DiagMatr2.
 * 
 * @see
 * - getDiagMatr2()
 * - getInlineDiagMatr2()
 * - reportDiagMatr2()
 * - applyDiagMatr2()
 */
typedef struct {

    /** The number of qubits upon which the matrix can act, which is always @c =2.
     */
    int numQubits;

    /** The matrix dimension, equal to the number of contained diagonal elements, which is always @c =4.
     */
    qindex numElems;

    /** The diagonal elements of the matrix, stored as a flat array, in CPU stack memory.
     * These are safe to directly modify at any time, to update the matrix.
     */
    qcomp elems[4];

} DiagMatr2;


/** A variable-size, @f$n@f$-qubit, @f$2^n \times 2^n@f$, diagonal complex matrix,
 * stored as a length-@f$2^n@f$ array.
 * 
 * This is intended for _local_ diagonal matrices, smaller than the Qureg upon
 * which they operate. For a full-state diagonal matrix, see FullStateDiagMatr.
 * 
 * Unlike DiagMatr1 and DiagMatr2, the diagonal elements of this matrix are stored in
 * heap memory. But like all QuEST structs, it is safe to copy, and ergo to pass
 * to functions by-value, or be returned from them. However, its destructor
 * must only ever be called upon one such copy.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.cpp) 
 * examples of initialising DiagMatr.
 * 
 * @see
 * - createDiagMatr()
 * - [createInlineDiagMatr()](https://quest-kit.github.io/QuEST/group__matrices__create.html#gad0ac51cb08a5d0b4e0878a5fde7ca893)
 * - [setInlineDiagMatr()](https://quest-kit.github.io/QuEST/group__matrices__setters.html#ga7b2c2e134546df6101c6001727bc6d69)
 * - setDiagMatr()
 * - syncDiagMatr()
 * - reportDiagMatr()
 * - applyDiagMatr()
 * - destroyDiagMatr()
 */
typedef struct {

    /** The number of qubits upon which the matrix can act.
     * 
     * The total memory costs of the matrix scale exponentially with the
     * number of qubits; an @f$n@f$-qubit diagonal matrix contains @f$2^n@f$
     * complex elements.
     */
    int numQubits;

   /** The dimension of the matrix, which is a square matrix, and so
     * equals both the number of rows and columns, and the number of diagonal terms.
     * 
     * Letting @f$n=@f$ #numQubits, then #numElems @f$=2^n@f$.
     */
    qindex numElems;

    /** Whether the matrix is (within validation epsilon tolerance) unitary (@c =1), or known to be
     * non-unitary (@c =0), or whether it is unknown (@c =-1).
     * 
     * This is a heap pointer to a persistent flag which is initially @c =-1 at matrix creation, and
     * is only ever consulted and/or updated by input validation within functions like applyDiagMatr(),
     * when such validation is enabled. Calling setDiagMatr() or syncDiagMatr() restores #isApproxUnitary 
     * to @c -1.
     * 
     * Unitarity is measured approximately, with reference to the validation
     * epsilon as modified with setQuESTValidationEpsilon(). This flag will never be updated from
     * @c =-1 when validation is disabled. Like all epsilon-dependent fields, it is restored
     * to @c =-1 automatically whenever setQuESTValidationEpsilon() or setQuESTValidationEpsilonToDefault()
     * are called.
     * 
     * @see
     * - setQuESTValidationEpsilon()
     */
    int* isApproxUnitary;

    /** Whether the matrix is (within validation epsilon tolerance) Hermitian (@c =1), or known to be
     * non-Hermitian (@c =0), or whether it is unknown (@c =-1).
     * 
     * > [!NOTE]
     * > There are presently _no_ QuEST functions which expect Hermiticity of a DiagMatr, 
     * > and so this field is presently never modified from its initial value of @c =-1.
     * 
     * This field has all the same nuances and mechanisms as #isApproxUnitary.
     */
    int* isApproxHermitian;  // @todo currently unused (relevant to not-yet-implemented calc-expec-val)

    /** Whether all diagonal elements of the matrix are "approximately" non-zero (have a magnitude
     * exceeding the validation epsilon) in which case this field @c =1, or whether the matrix is
     * known to contain an approximate zero (@c =0), or whether the property is unknown (@c =-1).
     * 
     * This property is relevant when attempting to raise the DiagMatr to a negative power,
     * such as through applyDiagMatrPower().
     * 
     * This field has all the same nuances and mechanisms as #isApproxUnitary.
     */
    int* isApproxNonZero;

    /** Whether all diagonal elements of the matrix are strictly zero or positive (@c =1), or
     * whether one or more elements are known to be negative (@c =0), or whether the property is
     * unknown (@c =-1).
     * 
     * This property is relevant when attempting to raise the DiagMatr to non-integer power,
     * such as through applyDiagMatrPower().
     * 
     * This field has similar nuances and mechanics as #isApproxUnitary, but is _not_ dependent
     * upon the validation epsilon (as controlled with setQuESTValidationEpsilon()).
     */
    int* isStrictlyNonNegative;
    
    /** Whether the matrix elements were ever synchronised.
     * 
     * This is a heap pointer to a persistent flag which is initially @c 0 at matrix creation,
     * but which is permanently overwritten to @c 1 when synchronisation is performed, such as
     * via syncDiagMatr() or setDiagMatr(). The flag indicates whether the matrix
     * elements have been initialised (and when QuEST is GPU-accelerated, whether they have been
     * copied to GPU device memory), and ergo whether it is valid to pass the matrix to a
     * simulation function like applyDiagMatr().
     * 
     * Note this flag can only indicate whether the matrix has _ever_ been synced; it cannot be
     * used to detect whether manual modification of #cpuElems made after an initial sync have been
     * re-synced, as required for correct behaviour in GPU mode.
     * 
     * @see
     * - syncDiagMatr()
     */
    int* wasGpuSynced;

    /** The matrix diagonal elements, stored in CPU host memory.
     * 
     * It is safest to modify this matrix through setDiagMatr(), but direct modification
     * is possible.
     * 
     * > [!IMPORTANT]
     * > It is _critical_ to call syncDiagMatr() after direct modification of
     * > #cpuElems in order to update persistent matrix properties,
     * > such as its data in GPU device memory (even when not running in
     * > GPU-accelerated mode).
     * 
     * @see
     * - syncDiagMatr()
     */
    qcomp* cpuElems;

    /** The matrix diagonal elements, stored in GPU device memory.
     * 
     * This is a copy of #cpuElems consulted by QuEST's GPU backend and should
     * _never_ be modified directly. Instead, it is updated by calling syncDiagMatr()
     * after modifying #cpuElems, which is performed automatically by
     * API functions like setDiagMatr(). In this way, #gpuElems and #cpuElems
     * should never be out of sync with one another.
     * 
     * Within a GPU-enabled QuEST environment, every DiagMatr allocates #gpuElems
     * in GPU memory, even if never ultimately consulted from the GPU backend.
     */
    qcomp* gpuElems;

} DiagMatr;


/*
 * DISTRIBUTED MATRIX STRUCTS
 */


/** A variable-size but full-state diagonal complex matrix, which can only act upon a
 * Qureg of the same dimension.
 * 
 * Unlike DiagMatr, a FullStateDiagMatr is intendedly non-local and can be distributed,
 * just like a Qureg of the same dimension.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_matrices.cpp) 
 * examples of initialising FullStateDiagMatr.
 * 
 * @see
 * - createFullStateDiagMatr()
 * - [setInlineFullStateDiagMatr()](fhttps://quest-kit.github.io/QuEST/group__matrices__setters.html#gaa27aaa1a001a41f36a60847a19d824d4)
 * - setFullStateDiagMatr()
 * - setFullStateDiagMatrFromMultiVarFunc()
 * - setFullStateDiagMatrFromMultiDimLists()
 * - setFullStateDiagMatrFromPauliStrSum()
 * - syncFullStateDiagMatr()
 * - reportFullStateDiagMatr()
 * - applyFullStateDiagMatr()
 * - calcExpecFullStateDiagMatr()
 * - destroyFullStateDiagMatr()
 */
typedef struct {

    /** The number of qubits in the Qureg upon which the matrix can act.
     * 
     * The total memory costs of the matrix scale exponentially with the
     * number of qubits; an @f$n@f$-qubit full-state diagonal matrix contains @f$2^n@f$
     * complex elements, which may be distributed across MPI nodes.
     */
    int numQubits;

   /** The total dimension of the matrix, which is a square matrix, and so
     * equals both the number of rows and columns, and the number of diagonal terms.
     * 
     * Letting @f$n=@f$ #numQubits, then #numElems @f$=2^n@f$. This is the total number
     * of elements across all nodes, and so is independent of whether the matrix is
     * distributed.
     */
    qindex numElems;

    /** Whether the matrix has allocated GPU memory, and so can be used by QuEST's
     * GPU-accelerated backend. 
     * 
     * Unlike CompMatr and DiagMatr, a FullStateDiagMatr can opt _not_ to allocate GPU memory when
     * QuEST is running in GPU-accelerated mode.
     */
    int isGpuAccelerated;

    /** Whether multithreading will be employed when initialising or simulating the matrix.
     * 
     * When a QuEST function receives both a Qureg and a FullStateDiagMatr, multithreading
     * is performed when _either_ or _both_ structs have multithreading enabled.
     */
    int isMultithreaded;

    /** Whether the matrix elements are distributed, informing the size of #cpuElems on
     * every node.
     * 
     * When distributed, each node contains only a subset of the full matrix's diagonal
     * elements, in an identical fashion to how statevector Qureg are distributed.
     * Unlike a Qureg however, a FullStateDiagMatr has no attached communication buffer,
     * so have _half_ the total memory costs of a same-size distributed statevector Qureg.
     * 
     * When not distributed, every node contains a duplicate copy of all elements in the matrix.
     */
    int isDistributed;

    /** The number of diagonal matrix elements stored on each distributed node, and the
     * length of #cpuElems and #gpuElems.
     * 
     * This differs from #numElems only when the matrix is distributed, in which case it
     * equals #numElems divided by the number of nodes.
     */
    qindex numElemsPerNode;

    /** Whether the matrix is (within validation epsilon tolerance) unitary (@c =1), or known to be
     * non-unitary (@c =0), or whether it is unknown (@c =-1).
     *
     * This is a heap pointer to a persistent flag which is initially @c =-1 at matrix creation, and
     * is only ever consulted and/or updated by input validation within functions like
     * applyFullStateDiagMatr(), when such validation is enabled. Calling functions like
     * setFullStateDiagMatr() or syncFullStateDiagMatr() restores #isApproxUnitary to @c -1.
     * 
     * Unitarity is measured approximately, with reference to the validation
     * epsilon as modified with setQuESTValidationEpsilon(). This flag will never be updated from
     * @c =-1 when validation is disabled. Like all epsilon-dependent fields, it is restored
     * to @c =-1 automatically whenever setQuESTValidationEpsilon() or setQuESTValidationEpsilonToDefault()
     * are called.
     * 
     * When distributed, this property reflects the _global_ state of the entire matrix across
     * all nodes, and so should always be identical across nodes.
     * 
     * @see
     * - setQuESTValidationEpsilon()
     */
    int* isApproxUnitary;

    /** Whether the matrix is (within validation epsilon tolerance) Hermitian (@c =1), or known to be
     * non-Hermitian (@c =0), or whether it is unknown (@c =-1).
     *
     * This field is relevant to functions like calcExpecFullStateDiagMatr(), and
     * otherwise has all the same nuances and mechanisms as #isApproxUnitary.
     */
    int* isApproxHermitian;

    /** Whether all diagonal elements of the matrix are "approximately" non-zero (have a magnitude
     * exceeding the validation epsilon) in which case this field @c =1, or whether the matrix is
     * known to contain an approximate zero (@c =0), or whether the property is unknown (@c =-1).
     * 
     * This property is relevant when attempting to raise the matrix to a negative power,
     * such as through applyFullStateDiagMatrPower().
     * 
     * This field has all the same nuances and mechanisms as #isApproxUnitary.
     */
    int* isApproxNonZero;

    /** Whether all diagonal elements of the matrix are strictly zero or positive (@c =1), or
     * whether one or more elements are known to be negative (@c =0), or whether the property is
     * unknown (@c =-1).
     * 
     * This property is relevant when attempting to raise the matrix to non-integer power,
     * such as through applyFullStateDiagMatrPower().
     * 
     * This field has similar nuances and mechanics as #isApproxUnitary, but is _not_ dependent
     * upon the validation epsilon (as controlled with setQuESTValidationEpsilon()).
     */
    int* isStrictlyNonNegative;

    /** Whether the matrix elements were ever synchronised.
     * 
     * This is a heap pointer to a persistent flag which is initially @c 0 at matrix creation,
     * but which is permanently overwritten to @c 1 when synchronisation is performed, such as
     * via syncFullStateDiagMatr() or setFullStateDiagMatr(). The flag indicates whether the matrix
     * elements have been initialised (and when QuEST is GPU-accelerated, whether they have been
     * copied to GPU device memory), and ergo whether it is valid to pass the matrix to a
     * simulation function like applyFullStateDiagMatr().
     * 
     * Note this flag can only indicate whether the matrix has _ever_ been synced; it cannot be
     * used to detect whether manual modification of #cpuElems made after an initial sync have been
     * re-synced, as required for correct behaviour in GPU mode.
     * 
     * @see
     * - syncFullStateDiagMatr()
     */
    int* wasGpuSynced;

    /** This node's matrix diagonal elements, stored in CPU host memory.
     * 
     * This is an array of length #numElemsPerNode, containing the diagonal elements of
     * the matrix. It is safest to modify these elements through functions like 
     * setFullStateDiagMatr(). Direct modification is possible, but in distributed
     * simulation, must be done with care. The @f$i@f$-th local element on the node
     * with rank @f$r@f$ corresponds to the full matrix element @f$(m r + i, m r + i)@f$,
     * where @f$m=@f$ #numElemsPerNode.
     * 
     * > [!IMPORTANT]
     * > It is _critical_ to call syncFullStateDiagMatr() after direct modification of
     * > #cpuElems on any node, in order to update persistent matrix properties,
     * > such as its data in GPU device memory (even when not running in
     * > GPU-accelerated mode).
     * 
     * @see
     * - setFullStateDiagMatr()
     * - setFullStateDiagMatrFromMultiVarFunc()
     * - setFullStateDiagMatrFromMultiDimLists()
     * - setFullStateDiagMatrFromPauliStrSum()
     * - syncFullStateDiagMatr()
     */
    qcomp* cpuElems;

    /** The matrix diagonal elements, stored in GPU device memory.
     * 
     * This is a copy of #cpuElems consulted by QuEST's GPU backend and should
     * _never_ be modified directly. Instead, it is updated by calling
     * syncFullStateDiagMatr()
     * after modifying #cpuElems, which is performed automatically by
     * API functions like setFullStateDiagMatr(). In this way, #gpuElems and #cpuElems
     * should never be out of sync with one another.
     * 
     * The #gpuElems memory is only allocated when #isGpuAccelerated @c =1.
     */
    qcomp* gpuElems;

} FullStateDiagMatr;


/** @} */



// we define the remaining doc groups in advance, since their signatures are
// more naturally grouped in an implementation-specific way below. Note the
// above structs were not doc'd this way (which would be more consistent)
// because it inexplicably causes Doxygen to duplicate their section at the
// top-level under Matrices (rather than under Structs). Bizarre! The order
// of declaration below will match the order shown in the html doc.
/** 
 * @defgroup matrices_getters Getters
 * @brief Functions for obtaining fixed-size matrices.
 * 
 * @defgroup matrices_create Constructors
 * @brief Functions for creating variable-size matrices.
 * 
 * @defgroup matrices_destroy Destructors
 * @brief Functions for destroying existing matrices.
 * 
 * @defgroup matrices_reporters Reporters
 * @brief Functions for printing matrices.
 * 
 * @defgroup matrices_setters Setters
 * @brief Functions for overwriting the elements of matrices.
 * 
 * @defgroup matrices_sync Synchronisation
 * @brief Functions for overwriting a matrix's GPU (VRAM) memory with its CPU (RAM) contents.
 * @details These functions are only necessary when the user wishes to manually modify the
 *          elements of a matrix (in lieu of using the @ref matrices_setters "Setters"), to
 *          thereafter synchronise the changes to the GPU copy of the channel. These functions 
 *          have no effect when running without GPU-acceleration, but remain legal and harmless 
 *          to call (to achieve platform agnosticism).
 */



/*
 * FIZED-SIZE MATRIX GETTERS VIA POINTERS
 *
 * which are defined here in the header because the 'qcomp' type is interpreted
 * distinctly by C++ (the backend) and C (user code). The C and C++ ABIs do not
 * agree on a complex type, so a qcomp (to C; a _Complex, and to C++; a std::complex)
 * cannot be directly passed between C and C++ compiled binaries; nor can a CompMatr1
 * struct which unwraps the qcomp[][] array. However, the C and C++ complex types have 
 * identical memory layouts, so pointers to qcomp types can safely be passed between
 * C and C++ binaries. Ordinarily we leverage this by defining all qcomp-handling API
 * functions in C++, and defining additional C-only wrappers in wrappers.h, which pass 
 * only pointers to qcomp. Alas, we cannot use this trick here, because the CompMatr1/2 
 * fields are declared 'const'; we cannot modify them through pointers, nor should we 
 * try to address them. Ergo we directly define these functions below (static inline to 
 * avoid symbol duplication), initializing the struct in one line. These functions will be 
 * separately interpreted by the C and C++ compilers, resolving qcomp to their individual 
 * native complex types.
 * 
 * These functions permit users to pass heap and stack pointers:
 *   - qcomp** ptr = malloc(...); getCompMatr1(ptr);
 *   - qcomp* ptrs[2]; getCompMatr1(ptrs);
 * in both C and C++. Because of 1D pointer decay, they also permit:
 *   - qcomp* ptr = malloc(...); getDiagMatr1(ptr);
 *   - qcomp arr[2]; getDiagMatr1(arr);
 */


// private validators
#ifdef __cplusplus
extern "C" {
#endif

/// @private
extern void _validateNewNestedElemsPtrNotNull(qcomp** ptrs, int numQubits, const char* caller);

/// @private
extern void _validateNewElemsPtrNotNull(qcomp* ptr, const char* caller);

#ifdef __cplusplus
}
#endif


/** @ingroup matrices_getters
 * @notyetdoced
 * 
 * @see
 * - reportCompMatr1()
 * - getInlineCompMatr2()
 * - getDiagMatr1()
 * - getCompMatr2()
 * - createCompMatr()
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
 */
static inline CompMatr1 getCompMatr1(qcomp** in) {
    _validateNewNestedElemsPtrNotNull(in, 1, __func__);

    CompMatr1 out;

    out.numQubits = 1;
    out.numRows = 2;
    for (int r=0; r<2; r++)
        for (int c=0; c<2; c++)
            out.elems[r][c] = in[r][c];

    return out;
}


/** @ingroup matrices_getters
 * @notyetdoced
 * 
 * @see
 * - reportCompMatr2()
 * - getInlineCompMatr2()
 * - getDiagMatr2()
 * - createCompMatr()
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
 */
static inline CompMatr2 getCompMatr2(qcomp** in) {
    _validateNewNestedElemsPtrNotNull(in, 2, __func__);

    CompMatr2 out;

    out.numQubits = 2;
    out.numRows = 4;
    for (int r=0; r<4; r++)
        for (int c=0; c<4; c++)
            out.elems[r][c] = in[r][c];
    
    return out;
}


/** @ingroup matrices_getters
 * @notyetdoced
 * 
 * @see
 * - reportDiagMatr1()
 * - getInlineDiagMatr1()
 * - getCompMatr1()
 * - getDiagMatr2()
 * - createDiagMatr()
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
 */
static inline DiagMatr1 getDiagMatr1(qcomp* in) {
    _validateNewElemsPtrNotNull(in, __func__);

    DiagMatr1 out;

    out.numQubits = 1;
    out.numElems = 2;
    for (int i=0; i<2; i++)
        out.elems[i] = in[i];

    return out;
}


/** @ingroup matrices_getters
 * @notyetdoced
 * 
 * @see
 * - reportDiagMatr2()
 * - getInlineDiagMatr2()
 * - getCompMatr2()
 * - createDiagMatr()
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
 */
static inline DiagMatr2 getDiagMatr2(qcomp* in) {
    _validateNewElemsPtrNotNull(in, __func__);

    DiagMatr2 out;
    
    out.numQubits = 2;
    out.numElems = 4;
    for (int i=0; i<4; i++)
        out.elems[i] = in[i];

    return out;
}



/*
 * FIZED-SIZE MATRIX GETTERS VIA ARRAYS & VECTORS
 *
 * which define additional overloads for arrays, VLAs, C99 temporary arrays,
 * vectors and vector initialisation lists. This empowers C users to call:
 *   - qcomp arr[2][2]; getCompMatr1(arr);
 *   - int n=2; qcomp arr[n][n]; getCompMatr1(arr);
 *   - getCompMatr1( (qcomp[2][2]) {...} );
 * and C++ users call:
 *   - qcomp arr[2][2]; getCompMatr1(arr);
 *   - std::vector vec(2); getCompMatr1(vec);
 *   - getCompMatr1( {...} );
 */


// define the array overloads with a distinct name from the base
// C function - we will alias it with getCompMatr() using Generics

/// @private
static inline CompMatr1 _getCompMatr1FromArr(qcomp in[2][2]) {

    qcomp* rowPtrs[] = {in[0], in[1]};
    return getCompMatr1(rowPtrs);
}

/// @private
static inline CompMatr2 _getCompMatr2FromArr(qcomp in[4][4]) {

    qcomp* rowPtrs[] = {in[0], in[1], in[2], in[3]};
    return getCompMatr2(rowPtrs);
}


// no array overloads are necessary for getDiagMatr(), because
// a 1D array automatically decays to a pointer


#ifdef __cplusplus

    // C++ defines overloads which merely wrap _getCompMatr1FromArr(), for fixed-size arrays.
    // C++ also defines additional std::vector overloads (for convenience, and for inline initialisation).
    // these are defined in matrices.cpp because they invoke validation (checking vector sizes)


    /// @ingroup matrices_getters
    /// @notyetdoced
    static inline CompMatr1 getCompMatr1(qcomp in[2][2]) { return _getCompMatr1FromArr(in); }


    /// @ingroup matrices_getters
    /// @notyetdoced
    static inline CompMatr2 getCompMatr2(qcomp in[4][4]) { return _getCompMatr2FromArr(in); }


    /** @ingroup matrices_getters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportCompMatr1()
     * - getInlineCompMatr1()
     * - getDiagMatr1()
     * - getCompMatr2()
     * - createCompMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    CompMatr1 getCompMatr1(std::vector<std::vector<qcomp>> in);


    /** @ingroup matrices_getters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportCompMatr2()
     * - getInlineCompMatr2()
     * - getDiagMatr2()
     * - createCompMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    CompMatr2 getCompMatr2(std::vector<std::vector<qcomp>> in);


    /** @ingroup matrices_getters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportDiagMatr1()
     * - getInlineDiagMatr1()
     * - getCompMatr1()
     * - getDiagMatr2()
     * - createDiagMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    DiagMatr1 getDiagMatr1(std::vector<qcomp> in);


    /** @ingroup matrices_getters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportDiagMatr2()
     * - getInlineDiagMatr2()
     * - getCompMatr2()
     * - createDiagMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    DiagMatr2 getDiagMatr2(std::vector<qcomp> in);


#else

    // C uses C11 Generics to effectively overload getCompMatr1/2 to accept both
    // pointers (as prior defined) and arrays (wrapping _getCompMatr1FromArr()). Note:
    // - our macros below accept C99 variadic arguments so that users pass C99
    //   compound literals (e.g. (qcomp[]) {1,2}) in addition to existing ptrs.
    //   they cannot however exclude the (qcomp[]) syntax like C++ users enjoy, 
    //   which is why we will subsequently define a getInlineCompMatr1()
    // - Generics evaluate at compile-time (AFTER preprocessing) so their RHS
    //   expressions are limited; because of this, it is impossible to avoid
    //   defining the _getCompMatr1FromArr() inner functions to avoid exposing them.
    // - our Generics explicitly check for pointer types (qcomp**), but we use default 
    //   to catch all array types (qcomp[][n], or qcomp(*)[] due to automatic Generic 
    //   pointer decay in GCC). This makes the code more consistent with our variable-size 
    //   CompMatr macros later in this  file, which cannot use VLA in Generics at all. It
    //   also avoids the user having to see a Generic compilation error message when they 
    //   pass an invalid type. 
    // - Generic expansion does not recurse, hence our macro safely has the same name
    //   (e.g. getCompMatr1) as the inner function, defining a true overload 
    // - we could not have _Generic's 'default' to catch unrecognised types at compile
    //   time to issue a custom message, because we must expand _Generic to a function 
    //   rather than a macro; preprocessing is finished by the time _Generic evaluates,
    //   so a macro would always be substituted before compilation and if it contained
    //   a compile-time error, it will always be triggered. A function error however
    //   would compile fine, but the error message would only be triggered at runtime
    //   when the user actually calls getCompMatr1() which is much worse than a slightly
    //   less clear compile-time error! A non-portable solution to this is to use
    //   _Pragma() in the RHS which is evaluated at compile-time (NOT pre-procesing),
    //   e.g. default: _Pragma("GCC error \"arg not allowed\"").
    

    /// @neverdoced
    #define getCompMatr1(...) \
        _Generic((__VA_ARGS__), \
            qcomp** : getCompMatr1, \
            default : _getCompMatr1FromArr \
        )((__VA_ARGS__))


    /// @neverdoced
    #define getCompMatr2(...) \
        _Generic((__VA_ARGS__), \
            qcomp** : getCompMatr2, \
            default : _getCompMatr2FromArr \
        )((__VA_ARGS__))


    // note the above macros do not need explicit, separate doxygen
    // doc because the C++ overloads above it have identical signatures

#endif



/*
 * FIXED-SIZE MATRIX GETTERS VIA LITERALS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * compound literal syntax. We expose these macros to C++ too for API consistency.
 * although C++'s getCompMatr1 vector overload achieves the same thing, and cannot
 * use C-style temporary arrays.
 * 
 * These empower C and C++ users to call
 *   - getInlineCompMatr1( {{1,2},{3,4}} )
 */


#ifdef __cplusplus

    // C++ merely invokes the std::vector initialiser overload

    /// @neverdoced
    #define getInlineCompMatr1(...) \
        getCompMatr1(__VA_ARGS__)

    /// @neverdoced
    #define getInlineCompMatr2(...) \
        getCompMatr2(__VA_ARGS__)

    /// @neverdoced
    #define getInlineDiagMatr1(...) \
        getDiagMatr1(__VA_ARGS__)

    /// @neverdoced
    #define getInlineDiagMatr2(...) \
        getDiagMatr2(__VA_ARGS__)

#else

    // C adds compound literal syntax to make a temporary array. Helpfully, 
    // explicitly specifying the DiagMatr dimension enables defaulting-to-zero

    /// @neverdoced
    #define getInlineCompMatr1(...) \
        _getCompMatr1FromArr((qcomp[2][2]) __VA_ARGS__)

    /// @neverdoced
    #define getInlineCompMatr2(...) \
        _getCompMatr2FromArr((qcomp[4][4]) __VA_ARGS__)

    /// @neverdoced
    #define getInlineDiagMatr1(...) \
        getDiagMatr1((qcomp[2]) __VA_ARGS__)

    /// @neverdoced
    #define getInlineDiagMatr2(...) \
        getDiagMatr2((qcomp[4]) __VA_ARGS__)

#endif

// spoofing above macros as functions to doc
#if 0

    /** @ingroup matrices_getters
     * @notyetdoced
     * @macrodoc
     * 
     * @see 
     * - reportCompMatr1()
     * - getInlineDiagMatr1()
     * - getInlineCompMatr2()
     * - createCompMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    CompMatr1 getInlineCompMatr1({{ matrix }});

    /** @ingroup matrices_getters
     * @notyetdoced
     * @macrodoc
     * 
     * @see 
     * - reportCompMatr2()
     * - getInlineDiagMatr2()
     * - createCompMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    CompMatr2 getInlineCompMatr2({{ matrix }});

    /** @ingroup matrices_getters
     * @notyetdoced
     * @macrodoc
     * 
     * @see 
     * - reportDiagMatr1()
     * - getInlineCompMatr1()
     * - getInlineDiagMatr2()
     * - createDiagMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    DiagMatr1 getInlineDiagMatr1({ list });

    /** @ingroup matrices_getters
     * @notyetdoced
     * @macrodoc
     * 
     * @see 
     * - reportDiagMatr2()
     * - getInlineCompMatr2()
     * - createDiagMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    DiagMatr2 getInlineDiagMatr2({ list });

#endif



/*
 * VARIABLE-SIZE MATRIX CONSTRUCTORS, DESTRUCTORS, SYNC
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif


    /** @ingroup matrices_create
     * @notyetdoced
     * 
     * @see
     * - setCompMatr()
     * - syncCompMatr()
     * - createInlineCompMatr()
     * - createDiagMatr()
     * - getCompMatr1()
     * - getCompMatr2()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    CompMatr createCompMatr(int numQubits);


    /** @ingroup matrices_create
     * @notyetdoced
     * 
     * @see
     * - setDiagMatr()
     * - syncDiagMatr()
     * - createInlineDiagMatr()
     * - createCompMatr()
     * - getDiagMatr1()
     * - getDiagMatr2()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    DiagMatr createDiagMatr(int numQubits);


    /** @ingroup matrices_create
     * @notyetdoced
     * 
     * @see
     * - createCustomFullStateDiagMatr()
     * - setFullStateDiagMatr()
     * - syncFullStateDiagMatr()
     * - createDiagMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    FullStateDiagMatr createFullStateDiagMatr(int numQubits);


    /** @ingroup matrices_create
     * @notyetdoced
     * 
     * @see
     * - createFullStateDiagMatr()
     * - setFullStateDiagMatr()
     * - syncFullStateDiagMatr()
     * - createDiagMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    FullStateDiagMatr createCustomFullStateDiagMatr(int numQubits, int useDistrib, int useGpuAccel, int useMultithread);


    /// @ingroup matrices_destroy
    /// @notyetdoced
    void destroyCompMatr(CompMatr matrix);


    /// @ingroup matrices_destroy
    /// @notyetdoced
    void destroyDiagMatr(DiagMatr matrix);


    /// @ingroup matrices_destroy
    /// @notyetdoced
    void destroyFullStateDiagMatr(FullStateDiagMatr matrix);


    /// @ingroup matrices_sync
    /// @notyetdoced
    void syncCompMatr(CompMatr matr);


    /// @ingroup matrices_sync
    /// @notyetdoced
    void syncDiagMatr(DiagMatr matr);


    /// @ingroup matrices_sync
    /// @notyetdoced
    void syncFullStateDiagMatr(FullStateDiagMatr matr);


#ifdef __cplusplus
}
#endif



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA POINTERS
 *
 * These functions permit users to pass heap and stack pointers:
 *   - qcomp** ptr = malloc(...); setCompMatr(m, ptr);
 *   - qcomp* ptrs[8]; setCompMatr(m, ptrs);
 * in both C and C++. By decay, they also permit arrays to diagonals:
 *   - qcomp* ptr = malloc(...); setDiagMatr(m, ptr);
 *   - qcomp arr[8]; setDiagMatr(m, arr); 
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif


    /** @ingroup matrices_setters
     * @notyetdoced
     * 
     * @see
     * - setInlineCompMatr()
     * - reportCompMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setCompMatr(CompMatr matr, qcomp** vals);


    /** @ingroup matrices_setters
     * @notyetdoced
     * 
     * @see
     * - setInlineDiagMatr()
     * - reportDiagMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setDiagMatr(DiagMatr out, qcomp* in);


    /** @ingroup matrices_setters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [setInlineFullStateDiagMatr()](fhttps://quest-kit.github.io/QuEST/group__matrices__setters.html#gaa27aaa1a001a41f36a60847a19d824d4)
     * - reportFullStateDiagMatr()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, qcomp* in, qindex numElems);


#ifdef __cplusplus
}
#endif



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA ARRAYS & VECTORS
 *
 * which define additional overloads for arrays, VLAs, vectors and vector initialisation lists.
 * C users can call:
 *   - qcomp arr[8][8]; setCompMatr(m, arr);
 *   - int n=8; qcomp arr[n][n]; setCompMatr(m, arr);
 *   - setCompMatr(m, (qcomp[8][8]) {{...}});
 *   - inline temporary VLA remains impossible even in C99, however
 * and C++ users can call:
 *   - int n=8; std::vector vec(n); setCompMatr(m, vec);
 *   - setCompMatr(m, {{...}});
 */


#if defined(__cplusplus)

    // C++ defines vector overloads, permitting inline initialisation


    /** @ingroup matrices_setters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - setInlineCompMatr()
     * - reportCompMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setCompMatr(CompMatr out, std::vector<std::vector<qcomp>> in);


    /** @ingroup matrices_setters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - setInlineDiagMatr()
     * - reportDiagMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setDiagMatr(DiagMatr out, std::vector<qcomp> in);


    /** @ingroup matrices_setters
     * @notyetdoced
     * @notyettested
     * @cpponly
     * 
     * @see
     * - [setInlineFullStateDiagMatr()](fhttps://quest-kit.github.io/QuEST/group__matrices__setters.html#gaa27aaa1a001a41f36a60847a19d824d4)
     * - reportFullStateDiagMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setFullStateDiagMatr(FullStateDiagMatr out, qindex startInd, std::vector<qcomp> in);


    // C++ cannot accept 2D arrays at all, because it does not support C99 VLA. 
    // It can however accept 1D arrays (which decay to pointers) already to setDiagMatr()

#elif !defined(_MSC_VER)

    // C first defines a bespoke functions receiving C99 VLAs, which we have to define here in
    // the header becauses the C++ source cannot use VLA, nor should we pass a 2D qcomp array
    // directly between C and C++ binaries (due to limited interoperability)


    // C must validate struct fields before accessing passed 2D arrays to avoid seg-faults
    /// @private
    extern void _validateParamsToSetCompMatrFromArr(CompMatr matr);


    // static inline to avoid header-symbol duplication
    /// @private
    static inline void _setCompMatrFromArr(CompMatr matr, qcomp arr[matr.numRows][matr.numRows]) {
        _validateParamsToSetCompMatrFromArr(matr);

        // new ptrs array safely fits in stack, since it's sqrt-smaller than user's passed stack array
        qcomp* ptrs[matr.numRows];

        // collect pointers to each row of arr
        for (qindex r=0; r<matr.numRows; r++)
            ptrs[r] = arr[r];

        // array decays to qcomp**, and *FromPtr function re-performs validation (eh)
        setCompMatr(matr, ptrs); // validation gauranteed to pass
    }


    // C then overloads setCompMatr() to call the above VLA when given arrays, using C11 Generics.
    // See the doc of getCompMatr1() above for an explanation of Generic, and its nuances


    /// @neverdoced
    #define setCompMatr(matr, ...) \
        _Generic((__VA_ARGS__), \
            qcomp** : setCompMatr, \
            default : _setCompMatrFromArr \
        )((matr), (__VA_ARGS__))

    // spoofing above macro as functions to doc
    #if 0

        /** @ingroup matrices_setters
         * @notyetdoced
         * @macrodoc
         * @conly
         * 
         * @see
         * - setInlineCompMatr()
         * - reportCompMatr()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) examples
         */
        void setCompMatr(CompMatr matr, qcomp arr[matr.numRows][matr.numRows]);

    #endif


    // no need to define bespoke overload for diagonal matrices, because 1D arrays decay to pointers

#else

    // MSVC's C11 does not support C99 VLAs (which the standard left optional, grr!), so
    // we cannot support 2D-array initialisation of CompMatr at all. This means only the 
    // existing setCompMatr(qcomp**) declared previously is usable by MSVC C users

#endif



/*
 * VARIABLE-SIZE MATRIX SETTERS VIA LITERALS
 *
 * which enable C users to give inline 2D array literals without having to use the
 * VLA compound literal syntax. We expose these macros to C++ too for API consistency,
 * although C++'s vector overloads achieve the same thing.
 * 
 * These empower C and C++ users to call e.g.
 *   - setInlineCompMatr(m, 1, {{1,2},{3,4}})
 */


#if defined(__cplusplus)

    // C++ redirects to vector overloads, passing initialiser lists.  The args like 'numQb'
    // are superfluous, but needed for consistency with the C API, so we additionally
    // validate that they match the struct dimensions (which requires validating the structs).


    /** @ingroup matrices_setters
     * @notyetdoced
     * @cpponly
     *
     * @see
     * - reportCompMatr()
     * - createInlineCompMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setInlineCompMatr(CompMatr matr, int numQb, std::vector<std::vector<qcomp>> in);


    /** @ingroup matrices_setters
     * @notyetdoced
     * @cpponly
     *
     * @see
     * - reportDiagMatr()
     * - createInlineDiagMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setInlineDiagMatr(DiagMatr matr, int numQb, std::vector<qcomp> in);


    /** @ingroup matrices_setters
     * @notyetdoced
     * @notyettested
     * @cpponly
     *
     * @see
     * - reportFullStateDiagMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    void setInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems, std::vector<qcomp> in);


#elif !defined(_MSC_VER)

    // C defines macros which add compound literal syntax so that the user's passed lists
    // become compile-time-sized temporary arrays. C99 does not permit inline-initialised
    // VLAs, so we cannot have the macro expand to add (qcomp[matr.numRows][matr.numRows])
    // in order to preclude passing 'numQb'. We ergo accept and validate 'numQb' macro param.
    // We define private inner-functions of a macro, in lieu of writing multiline macros
    // using do-while, just to better emulate a function call for users - e.g. they
    // can wrap the macro invocations with another function call, etc.


    // the C validators check 'numQb' is consistent with the struct, but cannot check the user's passed literal sizes
    /// @private
    extern void _validateParamsToSetInlineCompMatr(CompMatr matr, int numQb);
    /// @private
    extern void _validateParamsToSetInlineDiagMatr(DiagMatr matr, int numQb);
    /// @private
    extern void _validateParamsToSetInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems);


    /// @private
    static inline void _setInlineCompMatr(CompMatr matr, int numQb, qcomp elems[1<<numQb][1<<numQb]) {
        _validateParamsToSetInlineCompMatr(matr, numQb);
        _setCompMatrFromArr(matr, elems); // validation gauranteed to pass
    }

    /// @private
    static inline void _setInlineDiagMatr(DiagMatr matr, int numQb, qcomp elems[1<<numQb]) {
        _validateParamsToSetInlineDiagMatr(matr, numQb);
        setDiagMatr(matr, elems); // 1D array decays into pointer, validation gauranteed to pass
    }

    /// @private
    static inline void _setInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems, qcomp elems[numElems]) {
        _validateParamsToSetInlineFullStateDiagMatr(matr, startInd, numElems);
        setFullStateDiagMatr(matr, startInd, elems, numElems); // 1D array decays into pointer, validation gauranteed to pass
    }


    // happily, macro arg 'numQb' must be a compile-time constant, so there is no risk of
    // unexpectedly re-evaluating user expressions due to its repetition in the macro


    /// @neverdoced
    #define setInlineCompMatr(matr, numQb, ...) \
        _setInlineCompMatr((matr), (numQb), (qcomp[1<<(numQb)][1<<(numQb)]) __VA_ARGS__)

    /// @neverdoced
    #define setInlineDiagMatr(matr, numQb, ...) \
        _setInlineDiagMatr((matr), (numQb), (qcomp[1<<(numQb)]) __VA_ARGS__)

    /// @neverdoced
    #define setInlineFullStateDiagMatr(matr, startInd, numElems, ...) \
        _setInlineFullStateDiagMatr((matr), (startInd), (numElems), (qcomp[(numElems)]) __VA_ARGS__)

    // spoofing above macros as functions to doc
    #if 0


        /** @ingroup matrices_setters
         * @notyetdoced
         * @macrodoc
         * 
         * @see
         * - reportCompMatr()
         * - createInlineCompMatr()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) examples
         */
        void setInlineCompMatr(CompMatr matr, int numQb, {{ matrix }});


        /** @ingroup matrices_setters
         * @notyetdoced
         * @macrodoc
         * 
         * @see
         * - reportDiagMatr()
         * - createInlineDiagMatr()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) examples
         */
        void setInlineDiagMatr(DiagMatr matr, int numQb, { list });


        /** @ingroup matrices_setters
         * @notyettested
         * @notyetdoced
         * @macrodoc
         * 
         * @see
         * - reportFullStateDiagMatr()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) examples
         */
        void setInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems, { list });


    #endif


#else

    // MSVC C11 does not support C99 VLAs, so the inner functions above are illegal.
    // As such, we must choose to either forego the internal validation (which 
    // checks that the passed matrix object has been prior created with e.g.
    // createDiagMatr), or expand the macro into a do-while which users cannot ergo
    // place inside another function call. We opt to preclude the latter, since it
    // seems an unlikely use-case (because the function returns void) and will give
    // a compile-time error, whereas removing validation could cause silent seg-faults
    // when users incorrectly initialise an un-created matrix.

    // Note however that because MSVC does not support C99 VLA in C11, such that
    // _setCompMatrFromArr() was not defined, so we cannot define setInlineCompMatr();
    // MSVC C users simply miss out on this convenience function. Take it up with Bill!

    /// @private
    extern void _validateParamsToSetInlineDiagMatr(DiagMatr matr, int numQb);
    /// @private
    extern void _validateParamsToSetInlineFullStateDiagMatr(FullStateDiagMatr matr, qindex startInd, qindex numElems);


    /// @neverdoced
    #define setInlineDiagMatr(matr, numQb, ...) \
        do { \
            _validateParamsToSetInlineDiagMatr((matr), (numQb)); \
            setDiagMatr((matr), (numQb), (qcomp[1<<(numQb)]) __VA_ARGS__); \
        } while (0)


    /// @neverdoced
    #define setInlineFullStateDiagMatr(matr, startInd, numElems, ...) \
        do { \
            _validateParamsToSetInlineFullStateDiagMatr((matr), (startInd), (numElems)); \
            setFullStateDiagMatr((matr), (startInd), (elems), (numElems)); \
        } while (0)

    
    // the above macros are documented in the previous #if branch

#endif



/*
 * VARIABLE-SIZE MATRIX CREATORS VIA LITERALS
 *
 * which simply combine the create*() and setInline*() functions, for
 * user convenience, and to reduce their risk of passing inconsistent params.
 * We do not define inline creators for FullStateDiagMatr, since the
 * creator automatically decides whether or not to distribute the matrix;
 * ergo the user cannot know how many elements to pass in their literal
 * (nor should they ever distribute data which fits into a single literal!)
 * 
 * These empower C and C++ users to call e.g.
 *   - CompMatr m = createInlineCompMatr(1, {{1,2},{3,4}})
 */


#if defined(__cplusplus)

    // C++ accepts vector initialiser lists


    /** @ingroup matrices_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportCompMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    CompMatr createInlineCompMatr(int numQb, std::vector<std::vector<qcomp>> elems);


    /** @ingroup matrices_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - reportDiagMatr()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.cpp) examples
     */
    DiagMatr createInlineDiagMatr(int numQb, std::vector<qcomp> elems);


#elif !defined(_MSC_VER)

    // C defines macros which add compound literal syntax so that the user's passed lists
    // become compile-time-sized temporary arrays. We use bespoke validation so that the
    // error messages reflect the name of the macro, rather than the inner called functions.
    // We define a private inner function per macro, in lieu of writing multiline macros
    // using do-while, just to better emulate a function call for users - e.g. they
    // can wrap the macro invocation with another function call.


    /// @private
    extern void _validateParamsToCreateInlineCompMatr(int numQb);
    /// @private
    extern void _validateParamsToCreateInlineDiagMatr(int numQb);


    /// @private
    static inline CompMatr _createInlineCompMatr(int numQb, qcomp elems[1<<numQb][1<<numQb]) {
        _validateParamsToCreateInlineCompMatr(numQb);
        CompMatr out = createCompMatr(numQb); // malloc failures will report 'createCompMatr', rather than 'inline' version. Alas!
        _setCompMatrFromArr(out, elems);
        return out;
    }

    /// @private
    static inline DiagMatr _createInlineDiagMatr(int numQb, qcomp elems[1<<numQb]) {
        _validateParamsToCreateInlineDiagMatr(numQb);
        DiagMatr out = createDiagMatr(numQb); // malloc failures will report 'createCompMatr', rather than 'inline' version. Alas!
        setDiagMatr(out, elems); // 1D array decays to ptr
        return out;
    }


    /// @neverdoced
    #define createInlineCompMatr(numQb, ...) \
        _createInlineCompMatr((numQb), (qcomp[1<<(numQb)][1<<(numQb)]) __VA_ARGS__)

    /// @neverdoced
    #define createInlineDiagMatr(numQb, ...) \
        _createInlineDiagMatr((numQb), (qcomp[1<<(numQb)]) __VA_ARGS__)

    // spoofing above macros as functions to doc
    #if 0


        /** @ingroup matrices_create
         * @notyetdoced
         * @macrodoc
         * 
         * @see
         * - reportCompMatr()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) examples
         */
        CompMatr createInlineCompMatr(int numQb, {{ matrix }});


        /** @ingroup matrices_create
         * @notyetdoced
         * @macrodoc
         * 
         * @see
         * - reportDiagMatr()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_matrices.c) examples
         */
        DiagMatr createInlineDiagMatr(int numQb, { list });


    #endif

#else

    // MSVC's C11 does not support C99 VLA, so we cannot use the above inner functions.
    // The nuisance of trying to create, modify then return a matrix instance using
    // MSVC-specific preprocessors is too annoying, so Windows C users miss out! :(

#endif



/*
 * SPECIAL CREATORS AND SETTERS
 */


#ifdef __cplusplus
extern "C" {
#endif


    /// @todo
    /// add std::vector<int> overloads for C++ users for the
    /// below functions (missed during original overload work)


    /// @ingroup matrices_setters
    /// @notyetdoced
    /// @notyettested
    void setDiagMatrFromMultiVarFunc(DiagMatr out, qcomp (*func)(qindex*), int* numQubitsPerVar, int numVars, int areSigned);


    /// @ingroup matrices_setters
    /// @notyetdoced
    /// @notyettested
    void setDiagMatrFromMultiDimLists(DiagMatr out, void* lists, int* numQubitsPerDim, int numDims);


    /// @ingroup matrices_create
    /// @notyetdoced
    /// @notyettested
    FullStateDiagMatr createFullStateDiagMatrFromPauliStrSum(PauliStrSum in);


    /// @ingroup matrices_setters
    /// @notyetdoced
    /// @notyettested
    void setFullStateDiagMatrFromPauliStrSum(FullStateDiagMatr out, PauliStrSum in);


    /// @ingroup matrices_setters
    /// @notyetdoced
    /// @notyettested
    void setFullStateDiagMatrFromMultiVarFunc(FullStateDiagMatr out, qcomp (*func)(qindex*), int* numQubitsPerVar, int numVars, int areSigned);


    /// @ingroup matrices_setters
    /// @notyetdoced
    /// @notyettested
    void setFullStateDiagMatrFromMultiDimLists(FullStateDiagMatr out, void* lists, int* numQubitsPerDim, int numDims);


#ifdef __cplusplus
}
#endif



/*
 * MATRIX REPORTERS
 */


// de-mangle so below are directly callable by C binary
#ifdef __cplusplus
extern "C" {
#endif


    /** @ingroup matrices_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.cpp) examples
     */
    void reportCompMatr1(CompMatr1 matrix);


    /** @ingroup matrices_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.cpp) examples
     */
    void reportCompMatr2(CompMatr2 matrix);


    /** @ingroup matrices_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.cpp) examples
     */
    void reportCompMatr(CompMatr matrix);


    /** @ingroup matrices_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.cpp) examples
     */
    void reportDiagMatr1(DiagMatr1 matrix);


    /** @ingroup matrices_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.cpp) examples
     */
    void reportDiagMatr2(DiagMatr2 matrix);


    /// @ingroup matrices_reporters
    /// @notyetdoced
    /// @notyettested
    void reportDiagMatr(DiagMatr matrix);


    /** @ingroup matrices_reporters
     * @notyetdoced
     * @notyettested
     * 
     * @see
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.c) and
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_matrices.cpp) examples
     */
    void reportFullStateDiagMatr(FullStateDiagMatr matr);


#ifdef __cplusplus
}
#endif



#endif // MATRICES_H

/** @} */ // (end file-wide doxygen defgroup)
