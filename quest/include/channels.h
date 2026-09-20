/** @file
 * Data structures for representing arbitrary channels as
 * superoperators and Kraus maps, including their constructors, 
 * getters, setters and reporters. Note the functions to
 * actually simulate these channels are exposed in decoherence.h
 * 
 * Like matrices.h, this file makes extensive use of macros to
 * overload struct initialisers for user convenience. All macros
 * herein expand to single-line definitions for safety. Some 
 * intendedly private functions are necessarily exposed here to 
 * the user, and are prefixed with an underscore.
 * 
 * Design nuances:
 * - SuperOp is a separate, independent data-structure from KrausMap 
 *   which is never assumed/validated to be CPTP. This is because the 
 *   runtime assessment of CPTP of an arbitrary superoperator is expensive, 
 *   requiring diagonalisation.
 * - KrausMap contains an internal SuperOp instance which it uses to 
 *   simulate the channel, which is re-populated whenever the constituent 
 *   Kraus operators of the map are changed by the user.
 * - KrausMap maintains an explicit list of Kraus operators, even though 
 *   only the single resulting superoperator is used for simulation. This is 
 *   so that CPTP validation can be efficiently performed at any time, and so 
 *   that KrausMap reporting can display the individual quadratically-smaller 
 *   Kraus operators, for user clarity. This is an insignificant memory waste.
 * - KrausMap must know the number of constituent Kraus operators upfront, in 
 *   order to allocate their memory. It is not possible to specify fewer Kraus 
 *   operators later when initialising the KrausMap, because tracking a "number
 *   of Kraus operators" independently of the "maximum number of Kraus operators" 
 *   is smelly and over-engineered. Changing operators requires a new KrausMap.
 * - There are no fixed-size stack-memory versions of KrausMap and SuperOp, 
 *   unlike their matrix counterparts which have (e.g.) CompMatr1. This is 
 *   because fixed-size KrausMap creation involves populating the superoperator 
 *   (and in GPU settings, copying to GPU memory) which may be an astonishingly
 *   large overhead, and expensive to copy between function stacks.
 * 
 * @author Tyson Jones
 * @author Richard Meister (aided in design)
 * @author Erich Essmann (aided in design)
 * 
 * @defgroup channels Channels
 * @ingroup api
 * @brief Data structures for representing arbitrary channels as Kraus maps and superoperators.
 * @{
 */

#ifndef CHANNELS_H
#define CHANNELS_H

#include "quest/include/types.h"

// C++ gets vector initialiser overloads, whereas C gets a macro
#ifdef __cplusplus
    #include <vector>
#endif



/*
 * unlike some other headers, we here intermix the C and C++-only
 * signatures, grouping them semantically & by their doc groups
 */



/** 
 * @defgroup channels_structs Structs
 * @brief Data structures for representing decoherence channels.
 * @{
 */


/** A superoperator which acts upon both the ket and bra space of a
 * linearised density matrix, which can represent more transformations
 * than a KrausMap.
 * 
 * An @f$n@f$-qubit superoperator is instantiated as a @f$2^{2n}\times 2^{2n}@f$
 * complex matrix, and can only be applied upon density matrices of @f$n@f$ or
 * more density matrices.
 * 
 * Like all QuEST structs, a SuperOp is safe to copy, and ergo to pass to
 * functions by-value, or be returned from them. However, its destructor
 * must only ever be called upon one such copy.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_superoperators.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_superoperators.cpp) 
 * examples of initialising SuperOp.
 * 
 * @see
 * - createSuperOp()
 * - [createInlineSuperOp()](https://quest-kit.github.io/QuEST/group__channels__create.html#ga0ee76da0f63c68a2bf26c6dda973436d)
 * - setSuperOp()
 * - syncSuperOp()
 * - reportSuperOp()
 * - mixSuperOp()
 * - destroySuperOp()
 */
typedef struct {

    /** The number of qubits of the superoperator. 
     * 
     * This is _half_ the number of qubit substates of the linearised density
     * matrix upon which the superoperator acts, but is consistent with the space
     * acted upon by an equivalent a channel or unitary.
     * 
     * The total memory costs of the superoperator scale exponentially with the
     * number of qubits; an @f$n@f$-qubit superoperator contains @f$16^n@f$
     * complex elements.
     */
    int numQubits;

    /** The dimension of the superoperator, which is a square matrix, and so
     * equals both the number of rows and columns.
     * 
     * Letting @f$n=@f$ #numQubits, then #numRows @f$=4^n@f$.
     */
    qindex numRows;
    
    /** The 2D matrix elements of the operator, stored in CPU host memory.
     * 
     * It is safest to modify this matrix through setSuperOp(), but direct modification
     * is possible; the matrix element of the `r`-th row and `c`-th colum is stored
     * at `cpuElems[r][c]`.
     * 
     * > [!IMPORTANT]
     * > It is _critical_ to call syncSuperOp() after direct modification of
     * > #cpuElems in order to update persistent superoperator properties,
     * > such as its data in GPU device memory (even when not running in
     * > GPU-accelerated mode).
     * 
     * The field #cpuElems merely aliases the 1D #cpuElemsFlat field, such that
     * modifications of #cpuElems also updates #cpuElemsFlat.
     * 
     * @see
     * - syncSuperOp()
     */
    qcomp** cpuElems;

    /** A 1D row-major form of #cpuElems.
     * 
     * Modification of the superoperator matrix should be done through #cpuElems
     * for mathematical clarity, though this 1D contiguous form may be convenient
     * when performing copying.
     * 
     * > [!IMPORTANT]
     * > It is _critical_ to call syncSuperOp() after direct modification of
     * > #cpuElemsFlat in order to update persistent superoperator properties,
     * > such as its data in GPU device memory (even when not running in
     * > GPU-accelerated mode).
     * 
     * @see
     * - syncSuperOp()
     */
    qcomp* cpuElemsFlat;

    /** The elements of the superoperator matrix, stored in GPU device memory,
     * in a 1D row-major form.
     * 
     * This is a copy of #cpuElemsFlat consulted by QuEST's GPU backend and should
     * _never_ be modified directly. Instead, it is updated by calling syncSuperOp()
     * after modifying #cpuElems or #cpuElemsFlat, which is performed automatically by
     * API functions like setSuperOp(). In this way, #gpuElemsFlat and #cpuElemsFlat
     * should never be out of sync with one another.
     * 
     * Within a GPU-enabled QuEST environment, every SuperOp allocates #gpuElemsFlat
     * in GPU memory, even if never ultimately consulted from the GPU backend.
     */
    qcomp* gpuElemsFlat;

    /** Whether the superoperator matrix elements were ever synchronised.
     * 
     * This is a heap pointer to a persistent flag which is initially @c 0 at SuperOp creation,
     * but which is permanently overwritten to @c 1 when synchronisation is performed, such as
     * via syncSuperOp() or setSuperOp(). The flag indicates whether the superoperator matrix
     * elements have been initialised (and when QuEST is GPU-accelerated, whether they have been
     * copied to GPU device memory), and ergo whether it is valid to pass the SuperOp to a
     * simulation function like mixSuperOp().
     * 
     * Note this flag can only indicate whether the matrix has _ever_ been synced; it cannot be
     * used to detect whether manual modification of #cpuElems made after an initial sync have been
     * re-synced, as required for correct behaviour in GPU mode.
     * 
     * @see
     * - syncSuperOp()
     */
    int* wasGpuSynced;

} SuperOp;


/** A Kraus map which can act upon density matrices, and can describe any
 * physical operation or channel.
 * 
 * A @f$t@f$-operator @f$n@f$-qubit KrausMap is described by @f$t@f$ complex matrices,
 * each of dimension @f$2^n\times 2^n@f$; the equivalent size of an @f$n@f$-qubit unitary.
 * A KrausMap can be applied upon density matrices of @f$n@f$ or more qubits.
 * 
 * Due to its internal representation as a SuperOp with matrix dimension
 * @f$2^{2n}\times 2^{2n}@f$, the total memory costs of a KrausMap scale exponentially 
 * with the number of qubits as @f$16^n@f$.
 * 
 * Like all QuEST structs, a KrausMap is safe to copy, and ergo to pass to
 * functions by-value, or be returned from them. However, its destructor
 * must only ever be called upon one such copy.
 * 
 * @myexample
 * 
 * See [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_krausmaps.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/initialising_krausmaps.cpp) 
 * examples of initialising KrausMap.
 * 
 * @see
 * - createKrausMap()
 * - [createInlineKrausMap()](https://quest-kit.github.io/QuEST/group__channels__create.html#gae9c49a6443896ef590ff1e4cfaa4912b)
 * - setKrausMap()
 * - syncKrausMap()
 * - reportKrausMap()
 * - mixKrausMap()
 * - destroyKrausMap()
 */
typedef struct {

    /** The number of qubits of the Kraus map. 
     * 
     * Letting @f$n=@f$ #numQubits, each Kraus operator in the map is a @f$2^n\times 2^n@f$
     * complex matrix, and the full map is described by a @f$2^{2n}\times 2^{2n}@f$ 
     * superoperator; a total of @f$16^n@f$ complex elements.
     */
    int numQubits;

    /** The number of Kraus operators in the map.
     * 
     * For example, a KrausMap form of an amplitude damping channel would contain
     * @c numMatrices=2, and a two-qubit depolarising channel would contain
     * @c numMatrices=16.
     */
    int numMatrices;

    /** The dimension of each Kraus operator, which is a square matrix, and so
     * equivalent to the number of rows and columns.
     * 
     * Letting @f$n=@f$ #numQubits, then #numRows @f$=2^n@f$.
     */
    qindex numRows;

    /** A list of each Kraus oeprator's 2D matrix, stored in CPU host memory.
     * 
     * It is safest to modify this matrix through setKrausMap(), but direct modification
     * is possible; the matrix element of the `r`-th row and `c`-th colum of the `i`-th
     * operator is stored at `matrices[i][r][c]`.
     * 
     * > [!IMPORTANT]
     * > It is _critical_ to call syncKrausMap() after direct modification of
     * > #matrices in order to update persistent KrausMap properties,
     * > such as its SuperOp data in GPU device memory (even when not running in
     * > GPU-accelerated mode).
     * 
     * Unlike other data structures (such as CompMatr), KrausMap has no GPU device memory
     * copy of #matrices, since only the superoperator (stored within #superop) is
     * needed on the device. In fact, #matrices is only used for convenient
     * initialisation of a KrausMap, for validation of CPTP, and for printing by
     * reportKrausMap().
     */
    qcomp*** matrices;

    /** A superoperator representation of the KrausMap.
     * 
     * This is computed from #matrices during syncKrausMap() (as automatically invoked
     * by setKrausMap()), and is used by QuEST's simulation backend to effect the Kraus
     * map upon a density matrix.
     * 
     * Since this is the only field of KrausMap needing synchronisation, KrausMap
     * itself lacks an explicit @c wasGpuSynced field, and instead uses that attached
     * to #superop.
     */
    SuperOp superop;

    /** Whether the Kraus map is known to be (within validation epsilon tolerance) 
     * completely positive and trace preserving (@c =1), or known to be non-CPTP (@c =0),
     * or whether it is unknown (@c =-1).
     * 
     * This is a heap pointer to a persistent flag which is initially @c =-1 at KrausMap
     * creation, and is only ever consulted and/or updated by input validation within 
     * mixKrausMap(), when such validation is enabled. Calling setKrausMap() or syncKrausMap()
     * restores #isApproxCPTP to @c =-1.
     * 
     * The property of being CPTP is measured approximately, with reference to the validation
     * epsilon as modified with setQuESTValidationEpsilon(). The flag will never be updated from
     * @c =-1 when validation is disabled. Like all epsilon-dependent fields, it is restored
     * to @c =-1 automatically whenever setValidationEpsilon() or setValidationEpsilonToDefault()
     * are called.
     * 
     * To skip CPTP validation in mixKrausMap(), users can directly mutate this field to
     * @c =1, although it is not advised.
     * 
     * @see
     * - setQuESTValidationEpsilon()
     */
    int* isApproxCPTP;

} KrausMap;


// we define no fixed-size versions (e.g. KrausMap1/2), unlike we did for CompMatr1/2
// and DiagMatr1/2. This is because the 2-qubit superoperator is 256 elements big, and
// seems inadvisably large to be passing-by-copy through the QuEST backend layers, and
// would need explicit GPU memory allocation at each invocation of mixKrausMap2() (it
// exceeds the max number of CUDA kernel args). Furthermore, repeatedly calling
// createKrausMap2() would repeatedly invoke ~256*16 flops to compute te superoperator,
// which may be an user-astonishing overhead (more astonishing than the API asymmetry).
// Finally, computing the fixed-size superoperators must be in the header (to avoid
// the issues of qcmop interoperability, just like for getCompMatr1) and could not call
// an inner function which wouldn't be user-exposed; so we would end up redefining the
// superoperator calculation THREE times!


/** @} */



// we define the remaining doc groups in advance, since their signatures are
// more naturally grouped in an implementation-specific way below. Note the
// above structs were not doc'd this way (which would be more consistent)
// because it inexplicably causes Doxygen to duplicate their section at the
// top-level under Channels (rather than under Structs). Bizarre! The order
// of declaration below will match the order shown in the html doc.
/** 
 * @defgroup channels_create Constructors
 * @brief Functions for creating channel data structures.
 * 
 * @defgroup channels_destroy Destructors
 * @brief Functions for destroying existing channel data structures.
 * 
 * @defgroup channels_reporters Reporters
 * @brief Functions for printing channels.
 * 
 * @defgroup channels_setters Setters
 * @brief Functions for overwriting the elements of channels.
 * 
 * @defgroup channels_sync Synchronisation
 * @brief Functions for overwriting a channel's GPU (VRAM) memory with its CPU (RAM) contents.
 * @details These functions are only necessary when the user wishes to manually modify the
 *          elements of a channel (in lieu of using the @ref channels_setters "Setters"), to
 *          thereafter synchronise the changes to the GPU copy of the channel. These functions 
 *          have no effect when running without GPU-acceleration, but remain legal and harmless 
 *          to call (to achieve platform agnosticism).
 */



/*
 * BASIC FUNCTIONS
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif


    /** @ingroup channels_create
     * 
     * Creates an uninitialised Kraus map.
     *
     * The returned KrausMap contains @p numOperators Kraus operators, each of which
     * spans @p numQubits many qubits. Before being passed to functions like
     * reportKrausMap() and mixKrausMap(), its elements must be populated with
     * setKrausMap() or setInlineKrausMap(), or directly modified through KrausMap::matrices,
     * though any direct modification must be followed by a call to syncKrausMap().
     *
     * The returned KrausMap should be later destroyed with destroyKrausMap().
     * 
     * > See [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.c)
     * > or [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) 
     * > examples of initialising a KrausMap.
     *
     * @param[in] numQubits     the number of qubits acted upon by the Kraus map.
     * @param[in] numOperators  the number of Kraus operators in the map.
     * @returns A new KrausMap instance.
     * @throws @validationerror
     * - if the QuEST environment is not initialised.
     * - if @p numQubits or @p numOperators are invalid.
     * - if the dimensions or memory requirements overflow.
     * - if any memory allocation fails.
     * @see
     * - [createInlineKrausMap()](https://quest-kit.github.io/QuEST/group__channels__create.html#gae9c49a6443896ef590ff1e4cfaa4912b)
     * - createSuperOp()
     * - setKrausMap()
     * - [setInlineKrausMap()](https://quest-kit.github.io/QuEST/group__channels__setters.html#ga3c60440fa9503c235e46d964bc58d3ec)
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
     * @author Tyson Jones
     */
    KrausMap createKrausMap(int numQubits, int numOperators);


    /** @ingroup channels_sync
     * 
     * Updates the internal state of @p map, necessary after manually modifying KrausMap::matrices.
     *
     * @param[in,out] map  the KrausMap to synchronise.
     * @throws @validationerror
     * - if @p map is uninitialised.
     * @see
     * - setKrausMap()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
     * @author Tyson Jones
     */
    void syncKrausMap(KrausMap map);


    /** @ingroup channels_destroy
     * Destroys a KrausMap, freeing its Kraus operators and internal SuperOp.
     *
     * Since @p map is passed by value, this function cannot nullify the caller's
     * copy of the KrausMap fields. The caller must not use @p map after destruction.
     *
     * @param[in] map  the KrausMap to destroy.
     * @throws @validationerror
     * - if @p map is uninitialised.
     * @author Tyson Jones
     */
    void destroyKrausMap(KrausMap map);


    /** @ingroup channels_reporters
     * Prints a KrausMap.
     * 
     * @myexample
     * 
     * ```cpp
        KrausMap map = createInlineKrausMap(1, 3, {
            {{1,2},{3,4}},
            {{5,5},{6,6}},
            {{1i,2i},{-3i,-4i}}
        });
        reportKrausMap(map);
     * ```
     * ```text
        KrausMap (1 qubit, 3 2x2 matrices, 1 4x4 superoperator, 528 bytes):
            [matrix 0]
                1  2  
                3  4  
            [matrix 1]
                5  5  
                6  6  
            [matrix 2]
                i    2i   
                -3i  -4i  
     * ```
     *
     * @param[in] map  the KrausMap to print.
     * @throws @validationerror
     * - if @p map is uninitialised.
     * @author Tyson Jones
     */
    void reportKrausMap(KrausMap map);


    /** @ingroup channels_create
     * 
     * Creates an uninitialised superoperator.
     *
     * The returned SuperOp represents an arbitrary linear map on vectorised density
     * matrices, spanning @p numQubits many ket-qubits and an equal number of bra-qubits.
     * Before being passed to functions
     * like reportSuperOp() and mixSuperOp(), its elements must be populated with
     * setSuperOp() or setInlineSuperOp(), or directly modified through SuperOp::cpuElems,
     * though direct modification must be followed by a call to syncSuperOp().
     *
     * The returned SuperOp should be later destroyed with destroySuperOp().
     * 
     * > See [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.c)
     * > or [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) 
     * > examples of initialising a SuperOp.
     *
     * @param[in] numQubits  the number of qubits acted upon by the superoperator.
     * @returns A new SuperOp instance.
     * @throws @validationerror
     * - if the QuEST environment is not initialised.
     * - if @p numQubits is invalid.
     * - if the dimensions or memory requirements overflow.
     * - if any memory allocation fails.
     * @see
     * - createInlineSuperOp()
     * - createKrausMap()
     * - setSuperOp()
     * - setInlineSuperOp()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) examples
     * @author Tyson Jones
     */
    SuperOp createSuperOp(int numQubits);


    /** @ingroup channels_sync
     * 
     * Updates the internal state of @p op, necessary after manually modifying SuperOp::cpuElems
     * or SuperOp::cpuElemsFlat.
     *
     * @param[in,out] op  the SuperOp to synchronise.
     * @throws @validationerror
     * - if @p op is uninitialised.
     * @see
     * - setSuperOp()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) examples
     * @author Tyson Jones
     */
    void syncSuperOp(SuperOp op);


    /** @ingroup channels_destroy
     * Destroys a SuperOp, freeing all its CPU and GPU memory.
     *
     * Since @p op is passed by value, this function cannot nullify the caller's
     * copy of the SuperOp fields. The caller must not use @p op after destruction.
     *
     * @param[in] op  the SuperOp to destroy.
     * @throws @validationerror
     * - if @p op is uninitialised.
     * @author Tyson Jones
     */
    void destroySuperOp(SuperOp op);


    /** @ingroup channels_reporters
     * Prints a SuperOp.
     * 
     * @myexample
     * 
     * ```cpp
        SuperOp op = createInlineSuperOp(1, {
            {1,2,3,4},
            {5,-(10E-2)*3.14i,7,8},
            {9,10,11,12},
            {13,14,15,16+1.23i}
        });
        reportSuperOp(op);
     * ```
     * ```text
        SuperOp (1 qubit, 4x4 qcomps, 304 bytes):
            1   2        3   4         
            5   -0.314i  7   8         
            9   10       11  12        
            13  14       15  16+1.23i 
     * ```
     *
     * @param[in] op  the SuperOp to print.
     * @throws @validationerror
     * - if @p op is uninitialised.
     * @author Tyson Jones
     */
    void reportSuperOp(SuperOp op);


#ifdef __cplusplus
}
#endif



/*
 * POINTER INITIALISERS
 *
 * which permit users to pass heap and stack pointers in both C and C++, e.g.
 *   - qcomp** ptr = malloc(...); setSuperOp(m, ptr);
 *   - qcomp* ptrs[16]; setSuperOp(m, ptrs);
 *   - qcomp*** ptr = malloc(...); setKrausMap(m, ptr);
 */


// de-mangle so below are directly callable by C and C++ binary
#ifdef __cplusplus
extern "C" {
#endif


    /** @ingroup channels_setters
     * 
     * Overwrites the Kraus operators of @p map.
     * 
     * Argument @p matrices must be a list of KrausMap::numMatrices matrices, each of 
     * dimension KrausMap::numRows by KrausMap::numRows.
     * 
     * This updates KrausMap::matrices and other internal properties. 
     *
     * @param[in,out] map       the KrausMap to overwrite.
     * @param[in]     matrices  a 3D nested list of the above dimensions.
     * @throws @validationerror
     * - if @p map is uninitialised.
     * @throws seg-fault
     * - if @p matrices is not of the expected dimensions.
     * @see
     * - [setInlineKrausMap()](https://quest-kit.github.io/QuEST/group__channels__setters.html#ga3c60440fa9503c235e46d964bc58d3ec)
     * - syncKrausMap()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
     * @author Tyson Jones
     */
    void setKrausMap(KrausMap map, qcomp*** matrices);


    /** @ingroup channels_setters
     * 
     * Overwrites the elements of @p op.
     *
     * This copies @p matrix into SuperOp::cpuElems and synchronises @p op to GPU memory
     * when relevant.
     *
     * @param[in,out] op      the SuperOp to overwrite.
     * @param[in]     matrix  a SuperOp::numRows by SuperOp::numRows matrix of new elements.
     * @throws @validationerror
     * - if @p op is uninitialised.
     * - if @p matrix is a null or invalid pointer.
     * @see
     * - setInlineSuperOp()
     * - syncSuperOp()
     * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.c) or 
     *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) examples
     * @author Tyson Jones
     */
    void setSuperOp(SuperOp op, qcomp** matrix);


#ifdef __cplusplus
}
#endif



/*
 * ARRAY, VECTOR, MATRIX INITIALISERS
 *
 * which define additional overloads for arrays, VLAs, vectors and vector initialisation lists.
 * They permit C users to additionally call e.g.
 *   - qcomp arr[16][16]; setSuperOp(m, arr);
 *   - int n=16; qcomp arr[n][n]; setSuperOp(m, arr);
 *   - setKrausMap(m, (qcomp[5][16][16]) {{{...}}});
 *   - inline temporary VLA remains impossible even in C99, however
 * and C++ users gain overloads:
 *   - int n=8; std::vector vec(n); setSuperOp(m, vec);
 *   - setKrausMap(m, {{{...}}} );
 * An unintended but harmless side-effect is the exposure of functions setKrausMapFromArr(),
 * setSuperOpFromArr(), validate_setCompMatrFromArr() and validate_setSuperOpFromArr to the user.
 */


#if defined(__cplusplus)

    // C++ overloads to accept vectors, which also enables vector initialiser literals


    /** @ingroup channels_setters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - [setInlineKrausMap()](https://quest-kit.github.io/QuEST/group__channels__setters.html#ga3c60440fa9503c235e46d964bc58d3ec)
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
     */
    void setKrausMap(KrausMap map, std::vector<std::vector<std::vector<qcomp>>> matrices);


    /** @ingroup channels_setters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - setInlineSuperOp()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) examples
     */
    void setSuperOp(SuperOp op, std::vector<std::vector<qcomp>> matrix);
    

    // C++ cannot accept VLAs so does not define 2D array overloads

#elif !defined(_MSC_VER)

    // C first defines bespoke functions which accept C99 VLAs, which we have to define here in
    // the header becauses the C++ source cannot use VLA, nor should we pass a 2D qcomp array
    // directly between C and C++ binaries (due to limited interoperability)


    // C must validate the struct fields before accessing passed 2D arrays to avoid seg-faults
    /// @private
    extern void _validateParamsToSetKrausMapFromArr(KrausMap map);
    /// @private
    extern void _validateParamsToSetSuperOpFromArr(SuperOp op);


    /// @private
    static inline void _setKrausMapFromArr(KrausMap map, qcomp matrices[map.numMatrices][map.numRows][map.numRows]) {
        _validateParamsToSetKrausMapFromArr(map);

        // create stack space for 2D collection of pointers, one to each input row
        qcomp* rows[map.numMatrices][map.numRows];
        qcomp** ptrs[map.numMatrices];

        // copy decayed array pointers into stack
        for (int n=0; n<map.numMatrices; n++) {
            for (qindex r=0; r<map.numRows; r++)
                rows[n][r] = matrices[n][r];
            ptrs[n] = rows[n];
        }

        setKrausMap(map, ptrs); // validation gauranteed to pass
    }

    /// @private
    static inline void _setSuperOpFromArr(SuperOp op, qcomp matrix[op.numRows][op.numRows]) {
        _validateParamsToSetSuperOpFromArr(op);

        // create stack space for pointers, one for each input row
        qcomp* ptrs[op.numRows];

        // copy decayed array pointers into stack
        for (qindex r=0; r<op.numRows; r++)
            ptrs[r] = matrix[r];

        setSuperOp(op, ptrs); // validation gauranteed to pass
    }


    // C then overloads setKrausMap() to call the above VLA when given arrays, using C11 Generics.
    // See the doc of getCompMatr1() in matrices.h for an explanation of Generic, and its nuances.

    /// @neverdoced
    #define setKrausMap(map, ...) \
        _Generic((__VA_ARGS__), \
            qcomp*** : setKrausMap, \
            default  : _setKrausMapFromArr \
        )((map), (__VA_ARGS__))

    /// @neverdoced
    #define setSuperOp(op, ...) \
        _Generic((__VA_ARGS__), \
            qcomp** : setSuperOp, \
            default : _setSuperOpFromArr \
        )((op), (__VA_ARGS__))

    // spoofing macros as functions
    #if 0


        /** @ingroup channels_setters
         * 
         * Overwrites the Kraus operators of @p map from a @c C array.
         * 
         * @macrodoc
         * 
         * @conly
         *
         * This is a @c C convenience macro equivalent to setKrausMap().
         * 
         * @myexample
         * 
         * ```c
            qcomp arr[2][4][4] = {
                {
                    {1,2,3,4},
                    {5,6,7,8},
                    {9,8,7,6},
                    {5,4,3,2},
                }, {
                    {1i,2i,3i},
                    {5i}
                }
            };
            KrausMap map = createKrausMap(2, 2);
            setKrausMap(map, arr);
         * ```
         *
         * @param[in,out] map       the KrausMap to overwrite.
         * @param[in]     matrices  a 3D array.
         * @throws @validationerror
         * - if @p map is uninitialised.
         * @see
         * - [setInlineKrausMap()](https://quest-kit.github.io/QuEST/group__channels__setters.html#ga3c60440fa9503c235e46d964bc58d3ec)
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.c) examples
         * @author Tyson Jones
         */
        void setKrausMap(KrausMap map, qcomp matrices[map.numMatrices][map.numRows][map.numRows]);


        /** @ingroup channels_setters
         * 
         * Overwrites the elements of @p op from a @c C array.
         * 
         * @macrodoc
         * 
         * @conly
         *
         * This is a @c C convenience macro equivalent to setSuperOp().
         * 
         * @myexample
         * 
         * ```c
            qcomp arr[4][4] = {
                {1,2,3,4},
                {5,6,7,8},
                {9,8,7,6},
                {5,4,3,2}
            };
            SuperOp a = createSuperOp(1);
            setSuperOp(a, arr);
         * ```
         *
         * @param[in,out] op      the SuperOp to overwrite.
         * @param[in]     matrix  a SuperOp::numRows by SuperOp::numRows array.
         * @throws @validationerror
         * - if @p op is uninitialised.
         * @see
         * - setSuperOp()
         * - setInlineSuperOp()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.c) examples
         * @author Tyson Jones
         */
        void setSuperOp(SuperOp op, qcomp matrix[op.numRows][op.numRows]);


    #endif


#else

    // MSVC's C11 does not support C99 VLAs, so there is no way to support _setKrausMapFromArr(),
    // and ergo no need for setKrausMap() or setSuperOp() wrappers. This sadly means MSVC C users
    // can only use the existing functions which accept qcomp*** and qcomp** respectively.

#endif



/*
 * LITERAL INITIALISERS
 *
 * which enable C users to give inline 2D and 3D array literals without having to use the
 * VLA compound literal syntax. We expose these macros to C++ too for API consistency,
 * although C++'s vector overloads achieve the same thing.
 * 
 * These empower C and C++ users to call
 *   - setInlineSuperOp(m, 1, {{...}});
 *   - setInlineKrausMap(m, 2, 16, {{{...}}});
 */


#if defined(__cplusplus)

    // C++ redirects to vector overloads, passing initialiser lists.  The args like 'numQb'
    // and 'numOps' are superfluous, but needed for consistency with the C API, so we additionally
    // validate that they match the struct dimensions (which requires validating the structs).


    /** @ingroup channels_setters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - setKrausMap()
     * - syncKrausMap()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
     */
    void setInlineKrausMap(KrausMap map, int numQb, int numOps, std::vector<std::vector<std::vector<qcomp>>> matrices);


    /** @ingroup channels_setters
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - setSuperOp()
     * - syncSuperOp()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) examples
     */
    void setInlineSuperOp(SuperOp op, int numQb, std::vector<std::vector<qcomp>> matrix);


#elif !defined(_MSC_VER)

    // C defines macros which add compound literal syntax so that the user's passed lists
    // become compile-time-sized temporary arrays. C99 does not permit inline-initialised
    // VLAs, so we cannot have the macro expand to add (qcomp[matr.numRows][matr.numRows])
    // in order to preclude passing 'numQb'. We ergo accept and validate 'numQb' macro param.
    // We define private inner-functions of a macro, in lieu of writing multiline macros
    // using do-while, just to better emulate a function call for users - e.g. they
    // can wrap the macro invocations with another function call, etc.


    // C validators check 'numQb' and 'numOps' are consistent with the struct, but cannot check the user's passed literal sizes
    /// @private
    extern void _validateParamsToSetInlineKrausMap(KrausMap map, int numQb, int numOps);
    /// @private
    extern void _validateParamsToSetInlineSuperOp(SuperOp op, int numQb);


    /// @private
    static inline void _setInlineKrausMap(KrausMap map, int numQb, int numOps, qcomp elems[numOps][1<<numQb][1<<numQb]) {
        _validateParamsToSetInlineKrausMap(map, numQb, numOps);
        _setKrausMapFromArr(map, elems);
    }

    /// @private
    static inline void _setInlineSuperOp(SuperOp op, int numQb, qcomp elems[1<<(2*numQb)][1<<(2*numQb)] ) {
        _validateParamsToSetInlineSuperOp(op, numQb);
        _setSuperOpFromArr(op, elems);
    }


    /// @neverdoced
    #define setInlineKrausMap(map, numQb, numOps, ...) \
        _setInlineKrausMap((map), (numQb), (numOps), (qcomp[(numOps)][1<<(numQb)][1<<(numQb)]) __VA_ARGS__)

    /// @neverdoced
    #define setInlineSuperOp(matr, numQb, ...) \
        _setInlineSuperOp((matr), (numQb), (qcomp[1<<(2*(numQb))][1<<(2*(numQb))]) __VA_ARGS__)

    // spoofing macros as functions
    #if 0


        /** @ingroup channels_setters
         * 
         * Overwrites the Kraus operators of @p map from an inline literal.
         * 
         * The @c {{{matrices}}} argument is a 3D array literal, of dimensions
         * @c numOps by @c 1<<numQb by @c 1<<numQb.
         * 
         * - In @c C, this is a macro, where @p numQb and @p numOps must be 
         *   compile-time literals, and turn @c {{{matrices}}} into a compound literal.
         * - In @c C++, this is a function which accepts @c {{{matrices}}} as a nested @c std::vector literal.
         *
         * @myexample
         * 
         * ```c
            KrausMap map = createKrausMap(1, 3);
            setInlineKrausMap(map, 1, 3, {
                {{1,2},{3,4}},
                {{5,5},{6,6}},
                {{1i,2i},{-3i,-4i}}
            });
         * ```
         * 
         * @param[in,out] map     the KrausMap to overwrite.
         * @param[in]     numQb   the number of qubits of the literal @c {{{matrices}}}.
         * @param[in]     numOps  the number of Kraus operators of the literal @c {{{matrices}}}.
         * @throws @validationerror
         * - if @p map is uninitialised.
         * - if @p numQb differs from the number of qubits in @p map.
         * - if @p numOps differs from the number of operators in @p map.
         * @see
         * - setKrausMap()
         * - syncKrausMap()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.c) and
         *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
         * @author Tyson Jones
         */
        void setInlineKrausMap(KrausMap map, int numQb, int numOps, {{{ matrices }}});


        /** @ingroup channels_setters
         * 
         * Overwrites the elements of @p op from an inline literal.
         * 
         * The @c {{matrix}} argument is a 2D array literal, of dimensions
         * `1<<(2*numQb)` by `1<<(2*numQb)`.
         * 
         * - In @c C, this is a macro, where @p numQb must be a compile-time literal, 
         *   amd turns @c {{matrix}} into a compound literal.
         * - In @c C++, this is a function which accepts @c {{matrix}} as a nested @c std::vector literal.
         * 
         * @myexample
         * 
         * ```c
            SuperOp a = createSuperOp(1);
            setInlineSuperOp(a, 1, {
                {1,2,3,4},
                {5,3.14i,7,8},
                {9,10,11,12},
                {13,14,15,16+1.23i}
            });
         * ```
         *
         * @param[in,out] op     the SuperOp to overwrite.
         * @param[in]     numQb  the number of qubits of the @c {{matrix}} literal.
         * @throws @validationerror
         * - if @p op is uninitialised.
         * - if @p numQb does not match the number of qubits in @p op.
         * @see
         * - setSuperOp()
         * - syncSuperOp()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.c) examples
         * @author Tyson Jones
         */
        void setInlineSuperOp(SuperOp op, int numQb, {{ matrix }});


    #endif

#else

    // MSVC's C11 does not support C99 VLA, so the inner *FromArr() functions have not
    // been defined, and ergo we cannot define setInlineKrausMap() nor setInlineSuperOp()

#endif



/*
 * LITERAL CREATORS
 *
 * which combine creators and the inline initialisation functions, so that
 * both C and C++ users can call e.g.
 *   - SuperOp op = createInlineSuperOp(2, {{...}});
 */


#if defined(__cplusplus)

    // C++ accepts vector initialiser lists


    /** @ingroup channels_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - createKrausMap()
     * - setKrausMap()
     * - syncKrausMap()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
     */
    KrausMap createInlineKrausMap(int numQubits, int numOperators, std::vector<std::vector<std::vector<qcomp>>> matrices);


    /** @ingroup channels_create
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - createSuperOp()
     * - setSuperOp()
     * - syncSuperOp()
     * - [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) examples
     */
    SuperOp createInlineSuperOp(int numQubits, std::vector<std::vector<qcomp>> matrix);


#elif !defined(_MSC_VER)

    // C defines macros which add compound literal syntax so that the user's passed lists
    // become compile-time-sized temporary arrays. We use bespoke validation so that the
    // error messages reflect the name of the macro, rather than the inner called functions.
    // We define a private inner function per macro, in lieu of writing multiline macros
    // using do-while, just to better emulate a function call for users - e.g. they
    // can wrap the macro invocation with another function call.


    /// @private
    extern void _validateParamsToCreateInlineKrausMap(int numQb, int numOps);
    /// @private
    extern void _validateParamsToCreateInlineSuperOp(int numQb);


    /// @private
    static inline KrausMap _createInlineKrausMap(int numQb, int numOps, qcomp matrices[numOps][1<<numQb][1<<numQb]) {
        _validateParamsToCreateInlineKrausMap(numQb, numOps);
        KrausMap out = createKrausMap(numQb, numOps); // malloc failures will report 'createKrausMap', rather than 'inline' version. Alas!
        _setKrausMapFromArr(out, matrices);
        return out;
    }

    /// @private
    static inline SuperOp _createInlineSuperOp(int numQb, qcomp matrix[1<<numQb][1<<numQb]) {
        _validateParamsToCreateInlineSuperOp(numQb);
        SuperOp out = createSuperOp(numQb); // malloc failures will report 'createSuperOp', rather than 'inline' version. Alas!
        _setSuperOpFromArr(out, matrix);
        return out;
    }


    /// @neverdoced
    #define createInlineKrausMap(numQb, numOps, ...) \
        _createInlineKrausMap((numQb), (numOps), (qcomp[(numOps)][1<<(numQb)][1<<(numQb)]) __VA_ARGS__)

    /// @neverdoced
    #define createInlineSuperOp(numQb, ...) \
        _createInlineSuperOp((numQb), (qcomp[1<<(2*(numQb))][1<<(2*(numQb))]) __VA_ARGS__)

    // spoofing macros as functions
    #if 0


        /** @ingroup channels_create
         * 
         * Creates and initialises a KrausMap from an inline literal.
         * 
         * This is a convenience macro which combines createKausMap() and setInlineKrausMap().
         * 
         * The @c {{{matrices}}} argument is a 3D array literal, of dimensions
         * @c numOps by @c 1<<numQb by @c 1<<numQb.
         * 
         * - In @c C, this is a macro, where @p numQb and @p numOps must be 
         *   compile-time literals, and turn @c {{{matrices}}} into a compound literal.
         * - In @c C++, this is a function which accepts @c {{{matrices}}} as a nested @c std::vector literal.
         *
         * The returned KrausMap should be later destroyed with destroyKrausMap().
         * 
         * @myexample
         * 
         * ```cpp
            KrausMap map = createInlineKrausMap(1, 3, {
                {{1,2},{3,4}},
                {{5,5},{6,6}},
                {{1i,2i},{-3i,-4i}}
            });
         * ```
         *
         * @param[in] numQb   the number of qubits acted upon by the Kraus map.
         * @param[in] numOps  the number of Kraus operators.
         * @returns A new KrausMap initialised with @p matrices.
         * @throws @validationerror
         * - if the QuEST environment is not initialised.
         * - if @p numQb or @p numOps are invalid.
         * - if dimensions or memory requirements overflow.
         * - if any memory allocation fails.
         * @see
         * - createKrausMap()
         * - setKrausMap()
         * - syncKrausMap()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.c) and
         *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_krausmaps.cpp) examples
         * @author Tyson Jones
         */
        KrausMap createInlineKrausMap(int numQb, int numOps, {{{ matrices }}});


        /** @ingroup channels_create
         * 
         * Creates and initialises a SuperOp from an inline literal.
         * 
         * This is a convenience macro which combines createSuperOp() and setInlineSuperOp().
         * 
         * The @c {{matrix}} argument is a 2D array literal, of dimensions
         * `1<<(2*numQb)` by `1<<(2*numQb)`.
         * 
         * - In @c C, this is a macro, where @p numQb must be a compile-time literal, 
         *   amd turns @c {{matrix}} into a compound literal.
         * - In @c C++, this is a function which accepts @c {{matrix}} as a nested @c std::vector literal.
         *
         * The returned SuperOp should be later destroyed with destroySuperOp().
         * 
         * @myexample
         * 
         * ```cpp
            SuperOp op = createInlineSuperOp(1, {
                {1,2,3,4},
                {5,6*3.14i,7,8},
                {9,10,11,12},
                {13,14,15,16+1.23i}
            });
         * ```
         *
         * @param[in] numQb  the number of qubits acted upon by the superoperator.
         * @returns A new SuperOp initialised with @p matrix.
         * @throws @validationerror
         * - if the QuEST environment is not initialised.
         * - if @p numQb is invalid.
         * - if dimensions or memory requirements overflow.
         * - if any memory allocation fails.
         * @see
         * - createSuperOp()
         * - setSuperOp()
         * - syncSuperOp()
         * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.c) and
         *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/initialising_superoperators.cpp) examples
         * @author Tyson Jones
         */
        SuperOp createInlineSuperOp(int numQb, {{ matrix }});


    #endif

#else

    // MSVC's C11 does not support C99 VLA, so none of the necessary inner functions are defined,
    // and ergo Windows C users cannot use createInlineKrausMap() nor createInlineSuperOp(). Tragic!

#endif



#endif // CHANNELS_H

/** @} */ // (end file-wide doxygen defgroup)
