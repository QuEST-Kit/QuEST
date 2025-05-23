/** @file
 * API signatures for creating and managing Quregs.
 * 
 * @author Tyson Jones
 * 
 * @defgroup qureg Qureg
 * @ingroup api
 * @brief Data structures for representing quantum states.
 * @{
 */

#ifndef QUREG_H
#define QUREG_H

#include "quest/include/types.h"



/*
 * These signatures are divided into three partitions; those which are
 * natively C and C++ compatible (first partition), then those which are
 * only exposed to C++ (second partition) because they return 'qcomp' 
 * which cannot cross the C++-to-C ABI, then C++only convenience functions.
 * The first partition defines the doc groups, and the latter partition 
 * functions are added into them.
 */



/*
 * C and C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/** 
 * @defgroup qureg_structs Structs
 * @brief Data structures for representing quantum registers.
 * @{
 */


/// @notyetdoced
typedef struct {

    // deployment configuration
    int isMultithreaded;
    int isGpuAccelerated;
    int isDistributed;

    // distributed configuration
    int rank;
    int numNodes;
    int logNumNodes;

    // dimension
    int isDensityMatrix;
    int numQubits;
    qindex numAmps;
    qindex logNumAmps;

    // distributed load
    qindex numAmpsPerNode;
    qindex logNumAmpsPerNode;
    qindex logNumColsPerNode;

    // amplitudes in CPU and GPU memory
    qcomp* cpuAmps;
    qcomp* gpuAmps;

    // communication buffer in CPU and GPU memory
    qcomp* cpuCommBuffer;
    qcomp* gpuCommBuffer;

} Qureg;

/** @} */



/** 
 * @defgroup qureg_create Constructors
 * @brief Functions for creating statevectors and density matrices.
 * @{
 */


/** Creates a statevector containing @p numQubits qubits, with automatically chosen deployments,
 * initialised in the zero state.
 * 
 * The chosen deployments (multithreading, GPU-acceleration and distribution) are informed by
 * which facilities are compiled, available at runtime, beneficial for the specified Qureg size,
 * and whether the necessary additional memory structures can fit in accelerators and buffers.
 * 
 * @par State
 * Let @f$ N = @f$ @p numQubits. The returned Qureg contains @f$ 2^N @f$ amplitudes, each represented
 * by a @c qcomp, initialised to state
 *   @f[ 
    \ket{0}^{\otimes N}
    \;\; = \;\;
    \{ 1, \, 0, \, 0, \, \dots, \, 0 \}.
 *   @f]
 *
 * @par Memory
 * The total allocated memory will depend upon the automatically chosen deployments, since
 * use of GPU-acceleration requires persistent device memory and distribution necessitates
 * allocating communication buffers. See createCustomQureg() for more information.
 * 
 * @equivalences
 * - This function is equivalent to calling createCustomQureg(), passing @c isDensMatr=0 and @c -1
 *   for all deployments to automate them.
 *   ```
    Qureg qureg = createCustomQureg(numQubits, 0, -1, -1, -1);
 *   ```
 * @myexample
 * ```
    Qureg qureg = createQureg(30);
    reportQuregParams(qureg);
 * ```
 * @param[in] numQubits the number of qubits in the output Qureg.
 * @returns A new Qureg instance.
 * @throws @validationerror
 * - if @p numQubits < 1
 * - if the Qureg dimensions would overflow the @c qindex type.
 * - if the total Qureg memory would overflow the @c size_t type.
 * - if the system contains insufficient RAM (or VRAM) to store the Qureg in any deployment.
 * - if any memory allocation unexpectedly fails.
 * @notyetvalidated
 * @see
 * - createDensityQureg() to create a density matrix which can additionally undergo decoherence.
 * - createForcedQureg() to create a statevector which is forced to make use of all available deployments.
 * - createCustomQureg() to explicitly set the used deployments.
 * @author Tyson Jones
 */
Qureg createQureg(int numQubits);


/** Creates a density matrix containing @p numQubits qubits, with automatically chosen deployments,
 * initialised in the zero state.
 * 
 * The chosen deployments (multithreading, GPU-acceleration and distribution) are informed by
 * which facilities are compiled, available at runtime, beneficial for the specified Qureg size,
 * and whether the necessary additional memory structures can fit in accelerators and buffers.
 * 
 * @par State
 * Let @f$ N = @f$ @p numQubits. The returned Qureg contains @f$ 2^N \times 2^N @f$ amplitudes, each 
 * represented by a @c qcomp, initialised to state
 *   @f[ 
    \ket{0}\bra{0}^{\otimes N}
    \;\; = \;\;
    \begin{pmatrix}
    1 & 0 & \dots \\
    0 & 0 &  \\
    \vdots & & \ddots 
    \end{pmatrix}.
 *   @f]
 *
 * @par Memory
 * A density matrix contains _square_ as many amplitudes as the equal-dimension statevector.
 * The total allocated memory will depend upon the automatically chosen deployments, since
 * use of GPU-acceleration requires persistent device memory and distribution necessitates
 * allocating communication buffers. See createCustomQureg() for more information.
 * 
 * @equivalences
 * - This function is equivalent to calling createCustomQureg(), passing @c isDensMatr=1 and @c -1
 *   for all deployments to automate them.
 *   ```
    Qureg qureg = createCustomQureg(numQubits, 1, -1, -1, -1);
 *   ```
 * @myexample
 * ```
    Qureg qureg = createDensityQureg(15);
    reportQuregParams(qureg);
 * ```
 * @param[in] numQubits the number of qubits in the output Qureg.
 * @returns A new Qureg instance.
 * @throws @validationerror
 * - if @p numQubits < 1
 * - if the Qureg dimensions would overflow the @c qindex type.
 * - if the total Qureg memory would overflow the @c size_t type.
 * - if the system contains insufficient RAM (or VRAM) to store the Qureg in any deployment.
 * - if any memory allocation unexpectedly fails.
 * @notyetvalidated
 * @see
 * - createQureg() to create a quadratically-smaller statevector Qureg which cannot undergo decoherence.
 * - createForcedDensityQureg() to create a density matrix which is forced to make use of all available deployments.
 * - createCustomQureg() to explicitly set the used deployments.
 * @author Tyson Jones
 */
Qureg createDensityQureg(int numQubits);


/** @notyetdoced
 * 
 * @equivalences
 * - This function is equivalent to calling createCustomQureg(), passing @c isDensMatr=0 and all
 *   deployments enabled by the QuEST environment.
 *   ```
    QuESTEnv env = getQuESTEnv();
    Qureg qureg = createCustomQureg(
        numQubits, 0, 
        env.isDistributed, 
        env.isGpuAccelerated, 
        env.isMultithreaded);
 *   ```
 */
Qureg createForcedQureg(int numQubits);


/** @notyetdoced
 * 
 * @equivalences
 * - This function is equivalent to calling createCustomQureg(), passing @c isDensMatr=1 and all
 *   deployments enabled by the QuEST environment.
 *   ```
    QuESTEnv env = getQuESTEnv();
    Qureg qureg = createCustomQureg(
        numQubits, 1, 
        env.isDistributed, 
        env.isGpuAccelerated, 
        env.isMultithreaded);
 *   ```
 */
Qureg createForcedDensityQureg(int numQubits);


/** Creates a statevector or density matrix with the specified deployments, initialised
 * in the zero state. This function is an alternative to createQureg() and createDensityQureg()
 * which permits explicitly forcing, disabling, or automating particular deployments.
 * 
 * @par State
 * Parameters @p numQubits and @p isDensMatr respectively inform the dimension of the
 * Qureg, and whether the Qureg is a density matrix or statevector. 
 * 
 * Let @f$ N = @f$ @p numQubits.
 * - When @p isDensMatr=0, the returned statevector contains @f$ 2^N @f$ amplitudes,
 *   initialised to state 
 *   @f[ 
    \ket{0}^{\otimes N}
    \;\; = \;\;
    \{ 1, \, 0, \, 0, \, \dots, \, 0 \}.
 *   @f]
 * - When @p isDensMatr=1, the returned density matrix contains @f$ 2^{N}\times 2^{N} @f$ amplitudes,
 *   initialised to state 
 *   @f[ 
    \ket{0}\bra{0}^{\otimes N}
    \;\; = \;\;
    \begin{pmatrix}
    1 & 0 & \dots \\
    0 & 0 &  \\
    \vdots & & \ddots 
    \end{pmatrix}.
 *   @f]
 * 
 * @par Deployments
 * The remaining parameters decide the deployments used to accelerate the Qureg in subsequent 
 * simulation, and the associated additional memory allocations.
 * - @p useDistrib indicates whether (@c =1) or not (@c =0) to distribute the Qureg's amplitudes
 *   across all available MPI nodes. This is suitable for Qureg which are too large to fit into a
 *   single node, and requires allocating an additional communication buffer per-node. When
 *   @c useDistrib=0 but the QuEST executable is launched in distributed settings, the Qureg 
 *   amplitudes will be duplicated upon every node.
 * - @p useGpuAccel indicates whether (@c =1) or not (@c =0) to GPU-accelerate the Qureg, and
 *   requires allocating additional persistent memory in the GPU VRAM. When combined with 
 *   @c useDistrib=1, every node will allocate both communication buffers _and_ persistent GPU 
 *   memory, and an additional persistent GPU-memory communication buffer.
 * - @p useMultithread indicates whether (@c =1) or not (@c =0) to use multithreading when
 *   subsequently modifying the Qureg with a CPU routine. This requires no additional allocations,
 *   and typically has no effect when GPU acceleration is also enabled.
 * 
 * The deployment parameters can also be @c -1 to let QuEST choose that parameter, taking into
 * account the other forced deployments. While it is always safe to _disable_ a deployment,
 * forcing a deployment which is invalid (e.g. because the device has insufficient free memory)
 * will throw a validation error.
 * 
 * @par Memory
 * The total allocated memory depends on all parameters (_except_ 
 * @p useMultithread), and the size of the variable-precision @c qcomp used to represent each
 * amplitude. This is determined by preprocessor @c FLOAT_PRECISION via
 * 
 * <center>
 * | @c FLOAT_PRECISION | @c qcomp size (bytes) | 
 * | --- | --- |
 * | 1   | 8   |
 * | 2   | 16  | 
 * | 4   | 16, 20, 32  | 
 * </center>
 * where the quad-precision size is platform specific, and is often the size of _two_ 
 * `long double` primitives.
 * 
 * Let:
 * - @f$ N = @f$ @p numQubits
 * - @f$ D=2^N @f$ or @f$ =2^{2N} @f$ (the total number of amplitudes in the state)
 * - @f$ B = @f$ @c sizeof(qcomp) (the size in bytes)
 * - @f$ W @f$ be the total number of distributed nodes (the "world size").
 * 
 * The allocated CPU memory (RAM) and GPU memory (VRAM) is
 * 
 * <center>
 * | @p useDistrib | @p useGpuAccel | RAM per node | RAM total | VRAM per node | VRAM total | memory total |
 * |---|---|---|---|---|---|---|
 * | 0 | 0 | @f$ B \, D @f$ | @f$ W  B \, D @f$ | 0 | 0 | @f$ W  B \, D @f$ |
 * | 0 | 1 | @f$ B \, D @f$ | @f$ W  B \, D @f$ | @f$ B \, D @f$ | @f$ W  B \, D @f$ | @f$ 2 \, W  B \, D @f$ |
 * | 1 | 0 | @f$ 2 \, B \, D \, / \, W @f$ | @f$ 2 \, B \, D @f$ | 0 | 0 | @f$ 2 \, B \, D @f$ |
 * | 1 | 1 | @f$ 2 \, B \, D \, / \, W @f$ | @f$ 2 \, B \, D @f$ | @f$ 2 \, B \, D \, / \, W @f$ | @f$ 2 \, B \, D @f$ | @f$ 4 \, B \, D @f$ |
 * </center>
 *
 * For illustration, using the default @c FLOAT_PRECISION=2 whereby @f$ B = 16 @f$ bytes, the <b>RAM _per node_</b>
 * over varying distributions is:
 * 
 * <center>
 * | @p isDensMatr | @p numQubits | @f$ W=1 @f$ | @f$ W=2 @f$ | @f$ W=4 @f$ | @f$ W=8 @f$ | @f$ W=16 @f$ | @f$ W=1024 @f$ |
 * | ------------- | ------------ | ----------- | ----------- | ----------- | ----------- | ------------ | ------------ |
 * | 0             | 20           | 16 MiB      | 16 MiB      | 8 MiB       | 4 MiB       | 2 MiB        | 32 KiB       |
 * | 0             | 30           | 16 GiB      | 16 GiB      | 8 GiB       | 4 GiB       | 2 GiB        | 32 MiB       |
 * | 0             | 35           | 512 GiB     | 512 GiB     | 256 GiB     | 128 GiB     | 64 GiB       | 1 GiB        |
 * | 0             | 40           | 16 TiB      | 16 TiB      | 8 TiB       | 4 TiB       | 2 TiB        | 32 GiB       |
 * | 0             | 45           | 512 TiB     | 512 TiB     | 256 TiB     | 128 TiB     | 64 TiB       | 1 TiB        |
 * | 0             | 50           | 16 PiB      | 16 PiB      | 8 PiB       | 4 PiB       | 2 PiB        | 32 TiB       |
 * | 1             | 10           | 16 MiB      | 16 MiB      | 8 MiB       | 4 MiB       | 2 MiB        | 32 KiB       |
 * | 1             | 15           | 16 GiB      | 16 GiB      | 8 GiB       | 4 GiB       | 2 GiB        | 32 MiB       |
 * | 1             | 20           | 16 TiB      | 16 TiB      | 8 TiB       | 4 TiB       | 2 TiB        | 32 GiB       |
 * | 1             | 25           | 16 PiB      | 16 PiB      | 8 PiB       | 4 PiB       | 2 PiB        | 32 TiB       |
 * </center>
 * 
 * @constraints
 * - Cannot use any deployment which has not been prior enabled during compilation, or disabled by createCustomQuESTEnv().
 * - Cannot distribute @f$ N @f$ qubits over more than @f$ 2^N @f$ nodes (regardless of @p isDensMatr).
 * - Cannot distribute when the executable was not launched using MPI (e.g. via @c mpirun).
 * - Cannot GPU-accelerate when a GPU is not available at runtime, or has insufficient memory.
 * 
 * @myexample
 * ```
    int numQubits = 30;
    int isDensMatr = 0;
    
    int useDistrib     = 1;  // use distribution
    int useMultithread = 0;  // don't use multithreading
    int useGpuAccel    = -1; // automate whether to GPU-accelerate

    Qureg qureg = createCustomQureg(
        numQubits, isDensMatr, 
        useDistrib, useGpuAccel, useMultithread);

    reportQuregParams(qureg);
 * ```
 *
 * @param[in] numQubits the number of qubits in the output Qureg.
 * @param[in] isDensMatr whether the Qureg is a density matrix (@c =1) or statevector (@c =0).
 * @param[in] useDistrib whether to force (@c =1), disable (@c =0) or automate (@c =-1) distribution.
 * @param[in] useGpuAccel whether to force (@c =1), disable (@c =0) or automate (@c =-1) GPU acceleration.
 * @param[in] useMultithread whether to force (@c =1), disable (@c =0) or automate (@c =-1) multithreading.
 * @returns A new Qureg instance of the specified dimension and deployments.
 * @throws @validationerror
 * - if @p numQubits < 1
 * - if @p isDensMatr is not @c 0 or @c 1
 * - if any of @p useDistrib, @p useGpuAccel, @p useMultithread is not @c 0, @c 1 or @c -1.
 * - if any of @p useDistrib, @p useGpuAccel, @p useMultithread is forced (@c =1) but is unsupported by the
 *   active QuESTEnv. This can happen because:
 *   - the particular deployment was disabled by initCustomQuESTEnv().
 *   - the deployment was not enabled during compilation.
 *   - @p useDistrib=1 but QuEST was not launched by MPI (e.g. via @c mpirun).
 *   - @p useGpuAccel=1 but a GPU is not accessible at runtime.
 * - if @p useDistrib=1 but the Qureg is too small to distribute over the running nodes.
 * - if the Qureg dimensions would overflow the @c qindex type.
 * - if the total Qureg memory would overflow the @c size_t type.
 * - if the system contains insufficient RAM (or VRAM) to store the Qureg.
 * - if any memory allocation unexpectedly fails.
 * @notyetvalidated
 * @see
 * - createQureg() to automate deployments (equivalent to passing @c -1).
 * - createForcedQureg() to use all available deployments.
 * @author Tyson Jones
 */
Qureg createCustomQureg(int numQubits, int isDensMatr, int useDistrib, int useGpuAccel, int useMultithread);


/// @notyetdoced
Qureg createCloneQureg(Qureg qureg);


/** @} */



/** 
 * @defgroup qureg_destroy Destructors
 * @brief Functions for destroying existing Qureg.
 * @{
 */


/// @notyetdoced
void destroyQureg(Qureg qureg);


/** @} */



/** 
 * @defgroup qureg_report Reporters
 * @brief Functions for printing Qureg states or reporting their configuration.
 * @{
 */


/** @notyetdoced
 * @notyettested
 * 
 * @see
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_quregs.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_quregs.cpp) examples
 */
void reportQuregParams(Qureg qureg);


/** @notyetdoced
 * @notyettested
 * 
 * @see
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_quregs.c) and
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/devel/examples/isolated/reporting_quregs.cpp) examples
 */
void reportQureg(Qureg qureg);


/** @} */



/** 
 * @defgroup qureg_sync Synchronisation
 * @brief Functions for copying memory between a Qureg's CPU (RAM) and GPU (VRAM) memory. 
 * @details These functions are only necessary when the user wishes to manually probe or
 *          modify the Qureg amplitudes (rather than use functions like getQuregAmps() and
 *          setQuregAmps()), to ensure that the CPU and GPU copies of the Qureg are identical.
 *          These functions have no effect when running without GPU-acceleration, but remain 
 *          legal and harmless to call, to achieve platform agnosticism.
 * @{
 */


/// @notyetdoced
/// @notyettested
void syncQuregToGpu(Qureg qureg);


/// @notyetdoced
/// @notyettested
void syncQuregFromGpu(Qureg qureg);


/// @notyetdoced
/// @notyettested
void syncSubQuregToGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps);


/// @notyetdoced
/// @notyettested
void syncSubQuregFromGpu(Qureg qureg, qindex localStartInd, qindex numLocalAmps);


/** @} */



/** 
 * @defgroup qureg_get Getters
 * @brief Functions for obtaining amplitudes from statevectors or density matrices.
 * @{
 */


/// @notyetdoced
void getQuregAmps(qcomp* outAmps, Qureg qureg, qindex startInd, qindex numAmps);


/// @notyetdoced
void getDensityQuregAmps(qcomp** outAmps, Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);


/** @} */


// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ ONLY FUNCTIONS
 *
 * which are not directly C-compatible because they pass or
 * return qcomp primitives by-value (rather than by pointer).
 * This is prohibited because the C and C++ ABI does not agree
 * on a complex type, though C's _Complex has the same memory
 * layout as C++'s std::complex<>. To work around this, the 
 * below functions have a C-compatible wrapper defined in
 * wrappers.h which passes/receives the primitives by pointer;
 * a qcomp ptr can be safely passed from the C++ source binary
 * the user's C binary. These functions use the existing doxygen
 * doc groups defined above
 */


/// @ingroup qureg_get
/// @notyetdoced
qcomp getQuregAmp(Qureg qureg, qindex index);


/// @ingroup qureg_get
/// @notyetdoced
qcomp getDensityQuregAmp(Qureg qureg, qindex row, qindex column);



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

#include <vector>


/// @ingroup qureg_get
/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cpponly
/// @see getQuregAmps()
std::vector<qcomp> getQuregAmps(Qureg qureg, qindex startInd, qindex numAmps);


/// @ingroup qureg_get
/// @notyettested
/// @notyetvalidated
/// @notyetdoced
/// @cpponly
/// @see getDensityQuregAmps()
std::vector<std::vector<qcomp>> getDensityQuregAmps(Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols);


#endif // __cplusplus



/*
 * ORPHAN DOC GROUP
 */

/** 
 * @defgroup qureg_setters Setters
 * @brief See @ref init_amps "Amplitude initialisations".
 */



#endif // QUREG_H

/** @} */ // (end file-wide doxygen defgroup)
