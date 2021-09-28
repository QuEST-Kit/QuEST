// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * The QuEST API.
 * This file contains the comments used by doxygen for generating API doc.
 *
 * @defgroup unitary Unitaries
 *      Unitary gates
 * @defgroup normgate Gates
 *      Non-unitary but norm-preserving gates, such as measurements
 * @defgroup operator Operators
 *      Non-physical operators which may be non-unitary, non-norm-preserving, even non-Hermitian
 * @defgroup decoherence Decoherence 
 *      Decoherence channels which act on density matrices to induce mixing
 * @defgroup type Data structures
 *      Precision-agnostic types and data structures for representing numbers, quantum states, and operators
 * @defgroup calc Calculations
 *      Calculations and property-getters which do not modify the studied quantum state
 * @defgroup init State initialisations
 *      Functions for preparing quantum states
 * @defgroup qasm QASM Logging
 *      Functions for recording performed gates to <a href="https://en.wikipedia.org/wiki/OpenQASM">QASM</a>
 * @defgroup debug Debugging
 *      Utilities for seeding and debugging, such as state-logging
 *
 * @author Ania Brown
 * @author Tyson Jones
 * @author Balint Koczor
 * @author Nicolas Vogt of HQS (one-qubit damping)
 */

# ifndef QUEST_H
# define QUEST_H

# include "QuEST_precision.h"

// prevent C++ name mangling
#ifdef __cplusplus
extern "C" {
#endif


/*
 * private structures
 */
    
// hide these from doxygen
/// \cond HIDDEN_SYMBOLS    

/** Codes for Z-axis phase gate variations
 *
 * @ingroup type
 * @author Ania Brown
 */
enum phaseGateType {SIGMA_Z=0, S_GATE=1, T_GATE=2};    
    
/** A logger of QASM instructions 
 *
 * @ingroup type
 * @author Tyson Jones
 */
typedef struct {
    
    char* buffer;       // generated QASM string
    int bufferSize;     // maximum number of chars before overflow
    int bufferFill;     // number of chars currently in buffer
    int isLogging;      // whether gates are being added to buffer
    
} QASMLogger;

/** Represents an array of complex numbers grouped into an array of 
 * real components and an array of coressponding complex components.
 *
 * @ingroup type
 * @author Ania Brown
 */
typedef struct ComplexArray
{
    qreal *real; 
    qreal *imag;
} ComplexArray;

/// \endcond

    

/*
 * public structures
 */
 
 /** Codes for specifying Pauli operators
  *
  * @ingroup type
  * @author Tyson Jones
  */
 enum pauliOpType {PAULI_I=0, PAULI_X=1, PAULI_Y=2, PAULI_Z=3};

/** Represents one complex number.
 *
 * @ingroup type
 * @author Ania Brown
 */
typedef struct Complex
{
    qreal real;
    qreal imag;
} Complex;

/** Represents a 2x2 matrix of complex numbers.
 *
 * In C, a ::ComplexMatrix2 can be initialised by separately specifying 
 * the real and imaginary components as nested arrays. \n
 * For example,
 * ```
 * ComplexMatrix2 m = {
 *     .real = {{1,2}, 
 *              {3,4}},
 *     .imag = {{5,6}, 
 *              {7, 8}}};
 * ```
 * specifies matrix
 * \f[
 *   m = \begin{pmatrix}
 *      1 + 5\,i & 2+6\,i \\
 *      3 + 7\,i & 4+ 8\,i
 *   \end{pmatrix}
 * \f]
 * 
 * @see 
 * - ::ComplexMatrix4
 * - createComplexMatrixN()
 *
 * @ingroup type
 * @author Balint Koczor
 * @author Tyson Jones (doc)
 */
typedef struct ComplexMatrix2
{
    qreal real[2][2];
    qreal imag[2][2];
} ComplexMatrix2;

/** Represents a 4x4 matrix of complex numbers
 *
 * In C, a ::ComplexMatrix4 can be initialised by separately specifying 
 * the real and imaginary components as nested arrays. Note that in C99, a short row 
 * that ends with a 0 with be padded with 0. \n
 * For example, 
 * ```
 * ComplexMatrix4 m = {
 *      .real = {{1,2, 3, 4},
 *               {0},
 *               {5,6,7,8},
 *               {0}},
 *      .imag = {{0},{0},{0},{1,1,1,1}}};
 * ```
 * specifies matrix
 * \f[
 *   m = \begin{pmatrix}
 *      1 & 2 & 3 & 4 \\
 *      0 & 0 & 0 & 0 \\
 *      5 & 6 & 7 & 8 \\
 *      i & i & i & i
 *   \end{pmatrix}
 * \f]
 *
 * @see 
 * - ::ComplexMatrix2
 * - createComplexMatrixN()
 *
 * @ingroup type
 * @author Balint Koczor
 * @author Tyson Jones (doc)
 */
typedef struct ComplexMatrix4
{
    qreal real[4][4];
    qreal imag[4][4];
} ComplexMatrix4;
 
 /** Represents a general 2^N by 2^N matrix of complex numbers.
  *
  * @ingroup type
  * @author Tyson Jones
  */
typedef struct ComplexMatrixN
{
  int numQubits;
  qreal **real;
  qreal **imag;
} ComplexMatrixN;

/** Represents a 3-vector of real numbers
 *
 * @ingroup type
 * @author Ania Brown
 */
typedef struct Vector
{
    qreal x, y, z;
} Vector;

/** Flags for specifying named phase functions.
 * These can be passed to functions applyNamedPhaseFunc(), applyNamedPhaseFuncOverrides(), 
 * applyParamNamedPhaseFunc(), and applyParamNamedPhaseFuncOverrides().
 *
 * Norm based phase functions:
 *    - \p NORM maps state \f$|x\rangle|y\rangle\dots\f$ to \f$\sqrt{x^2 + y^2 + \dots}\f$
 *    - \p SCALED_NORM maps state \f$|x\rangle|y\rangle\dots\f$ to \f$\text{coeff} \sqrt{x^2 + y^2 + \dots}\f$
 *    - \p INVERSE_NORM maps state \f$|x\rangle|y\rangle\dots\f$ to \f$1/\sqrt{x^2 + y^2 + \dots}\f$
 *    - \p SCALED_INVERSE_NORM maps state \f$|x\rangle|y\rangle\dots\f$ to \f$\text{coeff}/\sqrt{x^2 + y^2 + \dots}\f$
 *    - \p SCALED_INVERSE_SHIFTED_NORM maps state \f$|x\rangle|y\rangle\dots\f$ to \f$\text{coeff}/\sqrt{(x-\Delta_x)^2 + (y-\Delta_y)^2 + \dots}\f$
 *
 * Product based phase functions:
 *    - \p PRODUCT maps state \f$|x\rangle|y\rangle|z\rangle\dots\f$ to \f$x \; y \; z \dots\f$
 *    - \p SCALED_PRODUCT maps state \f$|x\rangle|y\rangle|z\rangle\dots\f$ to \f$\text{coeff} \; x \; y \; z \dots\f$
 *    - \p INVERSE_PRODUCT maps state \f$|x\rangle|y\rangle|z\rangle\dots\f$ to \f$1/(x \; y \; z \dots)\f$
 *    - \p SCALED_INVERSE_PRODUCT maps state \f$|x\rangle|y\rangle|z\rangle\dots\f$ to \f$\text{coeff}/(x \; y \; z \dots)\f$
 *
 * Euclidean distance based phase functions:
 *    - \p DISTANCE maps state \f$|x_1\rangle|x_2\rangle|y_1\rangle|y_2\rangle\dots\f$ to \f$\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + \dots}\f$
 *    - \p SCALED_DISTANCE maps state \f$|x_1\rangle|x_2\rangle|y_1\rangle|y_2\rangle\dots\f$ to \f$\text{coeff}\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + \dots}\f$
 *    - \p INVERSE_DISTANCE maps state \f$|x_1\rangle|x_2\rangle|y_1\rangle|y_2\rangle\dots\f$ to \f$1/\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + \dots}\f$
 *    - \p SCALED_INVERSE_DISTANCE maps state \f$|x_1\rangle|x_2\rangle|y_1\rangle|y_2\rangle\dots\f$ to \f$\text{coeff}/\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + \dots}\f$
 *    - \p SCALED_INVERSE_SHIFTED_DISTANCE maps state \f$|x_1\rangle|x_2\rangle|y_1\rangle|y_2\rangle\dots\f$ to \f$\text{coeff}/\sqrt{(x_1-x_2-\Delta_x)^2 + (y_1-y_2-\Delta_y)^2 + \dots}\f$
 *
 * @ingroup type 
 * @author Tyson Jones
 * @author Richard Meister (shifted functions)
 */
enum phaseFunc {
    NORM=0,     SCALED_NORM=1,      INVERSE_NORM=2,      SCALED_INVERSE_NORM=3,      SCALED_INVERSE_SHIFTED_NORM=4,
    PRODUCT=5,  SCALED_PRODUCT=6,   INVERSE_PRODUCT=7,   SCALED_INVERSE_PRODUCT=8,
    DISTANCE=9, SCALED_DISTANCE=10, INVERSE_DISTANCE=11, SCALED_INVERSE_DISTANCE=12, SCALED_INVERSE_SHIFTED_DISTANCE=13};
    
/** Flags for specifying how the bits in sub-register computational basis states 
 * are mapped to indices in functions like applyPhaseFunc().
 *
 *    - \p UNSIGNED means the bits encode an unsigned integer, hence
 * \f[ 
 * \begin{aligned}
 *     |00\rangle & \rightarrow \, 0 \\
 *     |01\rangle & \rightarrow \, 1 \\
 *     |10\rangle & \rightarrow \, 2 \\
 *     |11\rangle & \rightarrow \, 3
 * \end{aligned}
 * \f]
 *    - \p TWOS_COMPLEMENT means the bits encode a signed integer through 
 *      [two's complement](https://en.wikipedia.org/wiki/Two%27s_complement), such that
 * \f[ 
 * \begin{aligned}
 *     |000\rangle & \rightarrow \, 0 \\
 *     |001\rangle & \rightarrow \, 1 \\
 *     |010\rangle & \rightarrow \, 2 \\
 *     |011\rangle & \rightarrow \, 3 \\
 *     |100\rangle & \rightarrow \,-4 \\
 *     |101\rangle & \rightarrow \,-3 \\
 *     |110\rangle & \rightarrow \,-2 \\
 *     |111\rangle & \rightarrow \,-1
 * \end{aligned}
 * \f]
 * > Remember that the qubits specified within a sub-register, and their ordering (least to most 
 * > significant) determine the bits of a computational basis state, before intrepretation 
 * > as an encoding of an integer.
 * 
 * @ingroup type 
 * @author Tyson Jones
 */
enum bitEncoding {UNSIGNED=0, TWOS_COMPLEMENT=1};

/** A Pauli Hamiltonian, expressed as a real-weighted sum of pauli products,
 * and which can hence represent any Hermitian operator.
 *
 * @ingroup type 
 * @author Tyson Jones
 */
typedef struct PauliHamil 
{
    //! The Pauli operators acting on each qubit, flattened over every operator.
    //! This is a flat array of length PauliHamil.numSumTerms * PauliHamil.numQubits.
    enum pauliOpType* pauliCodes;
    //! The real coefficient of each Pauli product. This is an array of length PauliHamil.numSumTerms;
    qreal* termCoeffs;
    //! The number of terms in the weighted sum, or the number of Pauli products.
    int numSumTerms;
    //! The number of qubits informing the Hilbert dimension of the Hamiltonian.
    int numQubits;
} PauliHamil;

/** Represents a diagonal complex operator on the full Hilbert state of a \p Qureg.
 * The operator need not be unitary nor Hermitian (which would constrain it to
 * real values)
 *
 * @ingroup type
 * @author Tyson Jones 
 */
typedef struct DiagonalOp
{
    //! The number of qubits this operator can act on (informing its size)
    int numQubits;
    //! The number of the 2^numQubits amplitudes stored on each distributed node
    long long int numElemsPerChunk;
    //! The number of nodes between which the elements of this operator are split
    int numChunks;
    //! The position of the chunk of the operator held by this process in the full operator
    int chunkId;
    //! The real values of the 2^numQubits complex elements
    qreal *real;
    //! The imaginary values of the 2^numQubits complex elements
    qreal *imag;
    //! A copy of the elements stored persistently on the GPU
    ComplexArray deviceOperator;
} DiagonalOp;

/** Represents a system of qubits.
 * Qubits are zero-based
 *
 * @ingroup type
 * @author Ania Brown
 * @author Tyson Jones (density matrix)
 */
typedef struct Qureg
{
    //! Whether this instance is a density-state representation
    int isDensityMatrix;
    //! The number of qubits represented in either the state-vector or density matrix
    int numQubitsRepresented;
    //! Number of qubits in the state-vector - this is double the number represented for mixed states
    int numQubitsInStateVec;
    //! Number of probability amplitudes held in stateVec by this process
    //! In the non-MPI version, this is the total number of amplitudes
    long long int numAmpsPerChunk;
    //! Total number of amplitudes, which are possibly distributed among machines
    long long int numAmpsTotal;
    //! The position of the chunk of the state vector held by this process in the full state vector
    int chunkId;
    //! Number of chunks the state vector is broken up into -- the number of MPI processes used
    int numChunks;
    
    //! Computational state amplitudes - a subset thereof in the MPI version
    ComplexArray stateVec; 
    //! Temporary storage for a chunk of the state vector received from another process in the MPI version
    ComplexArray pairStateVec;
    
    //! Storage for wavefunction amplitudes in the GPU version
    ComplexArray deviceStateVec;
    //! Storage for reduction of probabilities on GPU
    qreal *firstLevelReduction, *secondLevelReduction;

    //! Storage for generated QASM output
    QASMLogger* qasmLog;
    
} Qureg;

/** Information about the environment the program is running in.
 * In practice, this holds info about MPI ranks and helps to hide MPI initialization code
 *
 * @ingroup type
 * @author Ania Brown
 * @author Tyson Jones (seeding)
 */
typedef struct QuESTEnv
{
    int rank;
    int numRanks;
    unsigned long int* seeds;
    int numSeeds;
} QuESTEnv;



/*
 * public functions
 */

/** Creates a state-vector Qureg object representing a set of qubits which will remain in a pure state.
 * 
 * Allocates space for a state-vector of complex amplitudes, which assuming a single 
 * ::qreal floating-point number requires <b>qrealBytes</b>, requires memory
 * \f[ 
 *      \text{qrealBytes} \times 2 \times 2^\text{numQubits}\;\;\text{(bytes)},
 * \f]
 * though there are additional memory costs in GPU and distributed modes.
 *
 * The returned ::Qureg begins in the zero state, as produced by initZeroState().
 *
 * Once created, the following ::Qureg fields are relevant in all backends:
 * - Qureg.numQubitsRepresented
 * - Qureg.isDensityMatrix
 *
 * > ::QuESTEnv \p env must be prior created with createQuESTEnv().
 * 
 * ### Serial
 * In serial and local (non-distributed) multithreaded modes, a state-vector \p Qureg 
 * costs only the memory above. For example, at double precision 
 * (#QuEST_PREC <b>= 2</b>, <b>qrealBytes = 8</b>), the memory costs are:
 * \p numQubits  | memory
 * ------------- | -------------
 * 10            | 16 KiB
 * 16            | 1 MiB
 * 20            | 16 MiB
 * 26            | 1 GiB
 * 30            | 16 GiB
 *
 * Individual amplitudes should be fetched and modified with functions like getAmp() and setAmps().
 * However, it is sometimes useful to access the state-vector directly, for example to 
 * create your own low-level (high performance) multithreaded functions.
 * In those instants, Qureg.stateVec can be accessed directly, storing the real 
  * and imaginary components of the state-vector amplitudes in:
 * - `Qureg.stateVec.real`
 * - `Qureg.stateVec.imag`
 *
 * The total number of amplitudes in the state-vector is
 * - Qureg.numAmpsTotal
 *
 * For example,
 * ```
 * Qureg qureg = createQureg(10, env);
 *
 * // ruin a perfectly good state-vector
 * for (long long int i=0; i<qureg.numAmpsTotal; i++) {
 *     qureg.stateVec.real[i] = rand();
 *     qureg.stateVec.imag[i] = rand(); 
 * }
 * ```
 * \n
 *
 *
 * ### GPU
 *
 * In GPU-accelerated mode, an <em>additional</em> state-vector is created in GPU memory.
 * Therefore both RAM and VRAM must be of sufficient memory to store the state-vector,
 * each of the size indicated in the Serial table above.
 *
 * > Note that many GPUs do not support quad precision ::qreal.
 *
 * Individual amplitudes of the created ::Qureg should be fetched and modified with functions like getAmp() and setAmps().
 * This is especially important since the GPU state-vector can be accessed directly, 
 * and changes to Qureg.stateVec will be ignored and overwritten.
 * To modify the state-vector "directly", one must use copyStateFromGPU() and 
 * copyStateToGPU() before and after.
 *
 * For example,
 * ```
 * Qureg qureg = createQureg(10, env);
 *
 * // ruin a perfectly good state-vector
 * for (long long int i=0; i<qureg.numAmpsTotal; i++) {
 *     qureg.stateVec.real[i] = rand();
 *     qureg.stateVec.imag[i] = rand(); 
 * }
 * copyStateToGPU();
 * ```
 * \n
 *
 *
 * ### Distributed
 *
 * In distributed mode, the state-vector is uniformly partitioned between the <b>N</b>
 * distributed nodes.
 *
 * > Only a power-of-2 number of nodes <b>N</b> may be used (e.g. <b>N = 1, 2, 4, 8, </b>...).
 * > There must additionally be at least 1 amplitude of a state-vector stored on each node.
 * > This means one cannot create a state-vector ::Qureg with fewer than \f$\log_2(\text{N})\f$ qubits.
 * 
 * In addition to Qureg.stateVec, additional memory is allocated on each node for 
 * communication buffers, of size equal to the state-vector partition.
 * Hence the total memory <em>per-node</em> required is:
 * \f[ 
 *      2 \times \text{qrealBytes} \times 2 \times 2^\text{numQubits}/N  \;\;\text{(bytes)},
 * \f]
 *
 * For example, at double precision 
 * (#QuEST_PREC <b>= 2</b>, <b>qrealBytes = 8</b>), the memory costs are:
 *
 * | \p numQubits | memory per node  |||||
 * |-------------|-------|-------|-------|--------| ------ |
 * |            | <b>N = 2</b>  | <b>N = 4</b> | <b>N = 8</b> | <b>N = 16</b> | <b>N = 32</b> |
 * | 10            | 16 KiB | 8 KiB | 4 KiB | 2 KiB  | 1 KiB  |
 * | 20           | 16 MiB | 8 MiB | 4 MiB | 2 MiB  | 1 MiB  |
 * | 30           | 16 GiB | 8 GiB | 4 GiB | 2 GiB  | 1 GiB  |
 * | 40           | 16 TiB | 8 TiB | 4 TiB | 2 TiB  | 1 TiB  |
 *
 * State-vector amplitudes should be set and modified using getAmp() and setAmps().
 * Direct modification is possible, but should be done extremely carefully, since 
 * each node only stores a <em>partition</em> of the full state-vector, which itself 
 * mightn't fit on any single node. Furthermore, an asynchronous MPI process may
 * may unexpectedly modify local amplitudes; avoid this with syncQuESTEnv().
 *
 * The fields relevant to distribution are:
 *
 * - Qureg.numAmpsPerChunk: the length of Qureg.stateVec (`.real` and `.imag`) on each node.
 * - Qureg.chunkId: the id of the node, from <b>0</b> to <b>N-1</b>.
 *
 * Therefore, this code is valid
 * ```
   syncQuESTEnv(env);
 * // set state |i> to have amplitude i
 * for (long long int i=0; i<qureg.numAmpsPerChunk; i++)
 *     qureg.stateVec.real[i] = i + qureg.chunkId * qureg.numAmpsPerChunk;
 * ```
 * while the following erroneous code would cause a segmentation fault:
 * ```
 * // incorrectly attempt to set state |i> to have amplitude i
 * for (long long int i=0; i<qureg.numAmpsTotal; i++)
 *     qureg.stateVec.real[i] = i;
 * ``` 
 * \n
 *
 *
 * @see 
 * - createDensityQureg() to create a density matrix of the equivalent number of qubits, which can enter noisy states.
 * - createCloneQureg() to create a new qureg of the size and state of an existing qureg.
 * - destroyQureg() to free the allocated ::Qureg memory.
 * - reportQuregParams() to print information about a ::Qureg.
 * - copyStateFromGPU() and copyStateToGPU() for directly modifying state-vectors in GPU mode.
 * - syncQuESTEnv() for directly modifying state-vectors in distributed mode.
 *
 * @ingroup type
 * @returns an object representing the set of qubits
 * @param[in] numQubits number of qubits in the system
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @throws invalidQuESTInputError()
 * - if \p numQubits <= 0
 * - if \p numQubits is so large that the number of amplitudes cannot fit in a long long int type, 
 * - if in distributed mode, there are more nodes than elements in the would-be state-vector
 * @throws exit 
 * - if in GPU mode, but GPU memory cannot be allocated.
 * @author Ania Brown
 * @author Tyson Jones (validation, doc)
 */
Qureg createQureg(int numQubits, QuESTEnv env);

/** Creates a density matrix Qureg object representing a set of qubits which 
 * can enter noisy and mixed states.
 * 
 * Allocates space for a matrix of complex amplitudes, which assuming a single 
 * ::qreal floating-point number requires <b>qrealBytes</b>, requires memory
 * \f[ 
 *      \text{qrealBytes} \times 2 \times 2^{2 \times\text{numQubits}}\;\;\text{(bytes)},
 * \f]
 * though there are additional memory costs in GPU and distributed modes.
 * Notice this is the memory cost of a state-vector created with createQureg()
 * of twice as many qubits.
 *
 * The returned ::Qureg begins in the zero state, as produced by initZeroState().
 *
 * Once created, the following ::Qureg fields are relevant in all backends:
 * - Qureg.numQubitsRepresented
 * - Qureg.isDensityMatrix
 *
 * Behind the scenes, density matrice are stored as state-vectors, flattened column-wise.
 * As such, individual amplitudes should be fetched with getDensityAmp(), in lieu 
 * of direct access.
 * \n
 *
 * > ::QuESTEnv \p env must be prior created with createQuESTEnv().
 * 
 * ### Serial
 * In serial and local (non-distributed) multithreaded modes, a density matrix \p Qureg 
 * costs only the memory above. For example, at double precision 
 * (#QuEST_PREC <b>= 2</b>, <b>qrealBytes = 8</b>), the memory costs are:
 * \p numQubits  | memory
 * ------------- | -------------
 * 10            | 16 MiB
 * 12            | 256 MiB
 * 14            | 4 GiB
 * 16            | 64 GiB
 * 18            | 1 TiB
 * 20            | 16 TiB
 *
 *
 * ### GPU
 *
 * In GPU-accelerated mode, an <em>additional</em> density matrix is created in GPU memory.
 * Therefore both RAM and VRAM must be of sufficient memory to store the state-vector,
 * each of the size indicated in the Serial table above.
 *
 * > Note that many GPUs do not support quad precision ::qreal.
 *
 * ### Distributed
 *
 * In distributed mode, the density matrix is uniformly partitioned between the <b>N</b>
 * distributed nodes (column-wise).
 *
 * > Only a power-of-2 number of nodes <b>N</b> may be used (e.g. <b>N = 1, 2, 4, 8, </b>...).
 * > There must additionally be at least 1 amplitude of a density matrix stored on each node.
 * > This means one cannot create a density matrix ::Qureg with fewer than \f$\log_2(\text{N})/2\f$ qubits.
 * 
 * Additional memory is allocated on each node for communication buffers, of size 
 * equal to the density matrix partition.
 * Hence the total memory <em>per-node</em> required is:
 * \f[ 
 *      2 \times \text{qrealBytes} \times 2 \times 2^{2\times\text{numQubits}}/N  \;\;\text{(bytes)},
 * \f]
 *
 * For example, at double precision 
 * (#QuEST_PREC <b>= 2</b>, <b>qrealBytes = 8</b>), the memory costs are:
 *
 * | \p numQubits | memory per node  |||||
 * |-------------|-------|-------|-------|--------| ------ |
 * |            | <b>N = 2</b>  | <b>N = 4</b> | <b>N = 8</b> | <b>N = 16</b> | <b>N = 32</b> |
 * | 10           | 16 MiB | 8 MiB | 4 MiB | 2 MiB  | 1 MiB  |
 * | 15           | 16 GiB | 8 GiB | 4 GiB | 2 GiB  | 1 GiB  |
 * | 20           | 16 TiB | 8 TiB | 4 TiB | 2 TiB  | 1 TiB  |
 *
 *
 * @see 
 * - createQureg() to create a state-vector of the equivalent number of qubits, with a square-root memory cost
 * - createCloneQureg() to create a new qureg of the size and state of an existing qureg.
 * - destroyQureg() to free the allocated \p Qureg memory.
 * - reportQuregParams() to print information about a ::Qureg.
 *
 * @ingroup type
 * @returns an object representing the set of qubits
 * @param[in] numQubits number of qubits in the system
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @throws invalidQuESTInputError()
 * - if \p numQubits <= 0
 * - if \p numQubits is so large that the number of amplitudes cannot fit in a long long int type, 
 * - if in distributed mode, there are more nodes than elements in the would-be state-vector
 * @throws exit 
 * - if in GPU mode, but GPU memory cannot be allocated.
 * @author Tyson Jones
 */
Qureg createDensityQureg(int numQubits, QuESTEnv env);

/** Create a new ::Qureg which is an exact clone of the passed qureg, which can be
 * either a state-vector or a density matrix.
 *
 * The returned \ref Qureg will have the same 
 * dimensions as the passed \p qureg and begin in an identical quantum state.
 * This must be destroyed by the user later with destroyQureg().
 *
 * @see
 * - destroyQureg() 
 * - cloneQureg()
 * - createQureg()
 * - createDensityQureg()
 * 
 * @ingroup type
 * @returns an object representing the set of qubits
 * @param[in] qureg an existing \ref Qureg to be cloned
 * @param[in] env the ::QuESTEnv
 * @author Tyson Jones
 */
Qureg createCloneQureg(Qureg qureg, QuESTEnv env);

/** Deallocate a ::Qureg, freeing its memory.
 *
 * This frees all memory bound to \p qureg, including its state-vector or 
 * density matrix in RAM, in VRAM (in GPU mode), and communication buffers 
 * (in distributed mode).
 *
 * The \p qureg must have been previously created with createQureg(), 
 * createDensityQureg() or createCloneQureg().
 *
 * @see
 * - createQureg()
 * - createDensityQureg()
 * - createCloneQureg()
 *
 * @ingroup type
 * @param[in,out] qureg the ::Qureg to be destroyed
 * @param[in] env the ::QuESTEnv
 * @author Ania Brown
 * @author Tyson Jones (improved doc)
 */
void destroyQureg(Qureg qureg, QuESTEnv env);

/** Allocate dynamic memory for a square complex matrix of any size, 
 * which can be passed to functions like multiQubitUnitary() and applyMatrixN().
 *
 * The returned matrix will have dimensions
 * \f[
 *      2^{\text{numQubits}} \times 2^{\text{numQubits}},
 * \f]
 * stored as nested arrays ComplexMatrixN.real and ComplexMatrixN.imag,
 * initialised to zero.
 *
 * Unlike a ::Qureg, the memory of a ::ComplexMatrixN is always stored in RAM, 
 * and non-distributed. Hence, elements can be directly accessed and modified:
 * ```
 * int numQubits = 5;
 * int dim = (1 << numQubits);
 * ComplexMatrixN m = createComplexMatrixN(numQubits);
 *
 * for (int r=0; r<dim; r++) {
 *     for (int c=0; c<dim; c++) {
 *         m.real[r][c] = rand();
 *         m.imag[r][c] = rand();
 *     }
 * }
 * ```
 * \n
 * A ::ComplexMatrixN can be initialised in bulk using initComplexMatrixN(),
 * though this is not C++ compatible.
 * 
 * Like ::ComplexMatrix2 and ::ComplexMatrix4 (which are incidentally stored in the stack), 
 * the returned ::ComplexMatrixN is safe to return from functions.
 *
 * > The ::ComplexMatrixN must eventually be freed using destroyComplexMatrixN(),
 * > since it is created in the dynamic heap. One can instead use getStaticComplexMatrixN()
 * > to create a ComplexMatrixN struct in the stack (which doesn't need to be later destroyed),
 * > though this may cause a stack overflow if the matrix is too large (approx 10+ qubits).
 *
 * @see
 * - destroyComplexMatrixN()
 * - getStaticComplexMatrixN()
 * - initComplexMatrixN()
 * - applyMatrixN()
 * - multiQubitUnitary()
 * - mixMultiQubitKrausMap()
 *
 * @ingroup type
 * @param[in] numQubits the number of qubits of which the returned ComplexMatrixN will correspond
 * @returns a dynamic ComplexMatrixN struct, that is one where the .real and .imag
 *  fields are arrays kept in the heap and must be later destroyed.
 * @throws invalidQuESTInputError()
 * - if \p numQubits <= 0
 * - if the memory was not allocated successfully
 * @author Tyson Jones
 */
ComplexMatrixN createComplexMatrixN(int numQubits);

/** Destroy a ComplexMatrixN instance created with createComplexMatrixN()
 *
 * It is invalid to attempt to destroy a matrix created with getStaticComplexMatrixN().
 *
 * @see
 * - getStaticComplexMatrixN()
 * - createComplexMatrixN()
 *
 * @ingroup type
 * @param[in] matr the dynamic matrix (created with createComplexMatrixN()) to deallocate
 * @throws invalidQuESTInputError()
 * - if \p matr was not yet allocated.
 * @throws malloc_error
 * -  if \p matr was static (created with getStaticComplexMatrixN())
 * @author Tyson Jones
 */
void destroyComplexMatrixN(ComplexMatrixN matr);

#ifndef __cplusplus
#ifndef _WIN32
/** Initialises a ComplexMatrixN instance to have the passed
 * \p real and \p imag values. This allows succint population of any-sized
 * ComplexMatrixN, e.g. through 2D arrays:
 *
 *     ComplexMatrixN m = createComplexMatrixN(3);
 *     initComplexMatrixN(m, 
 *         (qreal[8][8]) {{1,2,3,4,5,6,7,8}, {0}},
 *         (qreal[8][8]) {{0}});
 *
 * \p m can be created by either createComplexMatrixN() or getStaticComplexMatrixN().
 *
 * This function is only callable in C, since C++ signatures cannot
 * contain variable-length 2D arrays 
 *
 * @ingroup type
 * @param[in] m the matrix to initialise
 * @param[in] real matrix of real values; can be 2D array of array of pointers
 * @param[in] imag matrix of imaginary values; can be 2D array of array of pointers
 * @throws invalidQuESTInputError()
 * - if \p m has not been allocated (e.g. with createComplexMatrixN())
 * @author Tyson Jones
 */
void initComplexMatrixN(ComplexMatrixN m, qreal real[][1<<m.numQubits], qreal imag[][1<<m.numQubits]);
#endif
#endif

/** Dynamically allocates a Hamiltonian expressed as a real-weighted sum of products of Pauli operators.
 *
 * A ::PauliHamil is merely an encapsulation of the multiple parameters of functions 
 * like applyPauliSum().
 *
 * The Pauli operators (PauliHamil.pauliCodes) are all initialised to identity 
 * (::PAULI_I), but the coefficients (PauliHamil.termCoeffs) are not initialised.
 *
 * The Hamiltonian can be used in functions like applyPauliHamil() and applyTrotterCircuit(),
 * with \p Qureg instances of the same number of qubits.
 *
 * A ::PauliHamil can be modified directly (see ::PauliHamil), or through initPauliHamil().
 * It can furthermore be created and initialised from a file description directly with 
 * createPauliHamilFromFile().
 * 
 * > The returned dynamic \p PauliHamil instance must later be freed via destroyPauliHamil().
 *
 * @see
 * - createPauliHamilFromFile()
 * - createDiagonalOpFromPauliHamilFile()
 * - initPauliHamil()
 * - destroyPauliHamil()
 * - applyPauliSum()
 * - applyTrotterCircuit()
 * - calcExpecPauliHamil()
 *
 * @ingroup type
 * @param[in] numQubits the number of qubits on which this Hamiltonian acts 
 * @param[in] numSumTerms the number of weighted terms in the sum, or the number of Pauli products
 * @returns a dynamic \p PauliHamil struct, with fields \p pauliCodes and \p termCoeffs stored in the heap
 * @throws invalidQuESTInputError()
 * - if \p numQubits <= 0
 * - if \p numSumTerms <= 0
 * @author Tyson Jones
 */
PauliHamil createPauliHamil(int numQubits, int numSumTerms);

/** Destroy a ::PauliHamil instance, created with either createPauliHamil() or createPauliHamilFromFile().
 *
 * @ingroup type 
 * @param[in] hamil a dynamic \p PauliHamil instantiation
 * @author Tyson Jones
 */
void destroyPauliHamil(PauliHamil hamil);

/** Creates a \p PauliHamil instance, a real-weighted sum of products of Pauli operators,
 * populated with the data in filename \p fn.
 * 
 * Each line in the plaintext file is interpreted as a separate product of Pauli operators 
 * in the total sum.
 * Each line must be a space-separated list with format
 *
 *     c p1 p2 p3 ... pN
 *
 * where \p c is the real coefficient of the term, and \p p1 ... \p pN are 
 * numbers in <b>{0,1,2,3}</b> to indicate ::PAULI_I, ::PAULI_X, ::PAULI_Y, ::PAULI_Z
 * operators respectively, which act on qubits \p 0 through \p N-1 (all qubits).
 *
 * For example, the file containing
 *
 *     0.31 1 0 1 2
 *     -0.2 3 2 0 0
 *
 * encodes a two-term four-qubit Hamiltonian
 * \f[ 
 *      0.31 \, X_0 \, X_2 \, Y_3 -0.2 \, Z_0 \, Y_1 \,.
 * \f]
 *
 * The initialised ::PauliHamil can be previewed with reportPauliHamil().
 *
 * The number of qubits and terms are inferred from the file.
 * The created Hamiltonian can be used just like one created via createPauliHamil().
 * It can be modified directly (see ::PauliHamil), or through initPauliHamil().
 * 
 * > The returned dynamic \p PauliHamil instance must later be freed via destroyPauliHamil().
 *
 * @see 
 * - reportPauliHamil()
 * - destroyPauliHamil()
 * - createPauliHamil()
 * - initPauliHamil()
 * - createDiagonalOpFromPauliHamilFile()
 *
 * @ingroup type
 * @param[in] fn filename of the plaintext file specifying the pauli operators and coefficients
 * @returns a dynamic ::PauliHamil struct
 * @throws invalidQuESTInputError()
 * - if the file with name \p fn cannot be read
 * - if the file is not correctly formatted as described above
 * @author Tyson Jones
 */
PauliHamil createPauliHamilFromFile(char* fn);

/** Initialise ::PauliHamil instance \p hamil with the given term coefficients and 
 * Pauli codes (one for every qubit in every term).
 *
 * Arguments \p coeffs and \p codes encode a weighted sum of Pauli operators, with the same 
 * format as other QuEST functions (like calcExpecPauliSum()).
 *
 * This is useful to set the elements of the ::PauliHamil in batch. \n
 * For example
 * ```
 * int numQubits = 3;
 * int numTerms = 2;
 * PauliHamil hamil = createPauliHamil(numQubits, numTerms);
 * 
 * // hamil = 0.5 X0 Y1 - 0.5 Z1 X3
 * initPauliHamil(hamil, 
 *     (qreal[]) {0.5, -0.5}, 
 *     (enum pauliOpType[]) {PAULI_X,PAULI_Y,PAULI_I, PAULI_I, PAULI_Z, PAULI_X});
 * ```
 *
 * The initialised ::PauliHamil can be previewed with reportPauliHamil().
 * 
 * > \p hamil must be already created with createPauliHamil(), or createPauliHamilFromFile().
 *
 * @see
 * - reportPauliHamil()
 * - createPauliHamil()
 * - createPauliHamilFromFile()
 * 
 * @ingroup type
 * @param[in, out] hamil an existing ::PauliHamil instance to be modified
 * @param[in] coeffs an array of sum term coefficients, which must have length `hamil.numSumTerms`
 * @param[in] codes a flat array of Pauli codes, of length `hamil.numSumTerms`*`hamil.numQubits`
 * @throws invalidQuESTInputError()
 * - if \p hamil has invalid parameters (\p numQubits <= 0, \p numSumTerms <= 0)
 * - if any code in \p codes is not a valid Pauli code (::pauliOpType)
 * @author Tyson Jones
 */
void initPauliHamil(PauliHamil hamil, qreal* coeffs, enum pauliOpType* codes);

/** Creates a ::DiagonalOp representing a diagonal operator on the 
 * full Hilbert space of a ::Qureg.
 * 
 * The resulting operator need not be unitary nor Hermitian, and can be 
 * applied to any ::Qureg of a compatible number of qubits.
 *
 * This function allocates space for \f$2^{\text{numQubits}}\f$ complex amplitudes,
 * which are initially zero. This is the same cost as a local state-vector of equal 
 * number of qubits; see the Serial section of createQureg(). 
 * Note that this is a <em>paralell</em> data-type, so its 
 * ultimate memory costs depend on the hardware backends, as elaborated below.
 * 
 * The operator elements should be modified with initDiagonalOp() and setDiagonalOpElems(), 
 * and must be later freed with destroyDiagonalOp().
 *
 * ### GPU 
 *
 * In GPU-accelerated mode, this function also creates additional equally-sized 
 * persistent memory on the GPU.
 * If you wish to modify the operator elements directly (in lieu of setDiagonalOpElems()),
 * you must thereafter call syncDiagonalOp() to update the operator stored in VRAM.
 *
 * For example,
 * ```
 * DiagonalOp op = createDiagonalOp(4, env);
 * for (long long int i=0; i<op.numElemsPerChunk; i++) {
 *     op.real[i] = rand();
 *     op.imag[i] = rand();
 * }
 * syncDiagonalOp(op);
 * ```
 * 
 * ### Distribution
 * 
 * In distributed mode, the memory for the diagonal operator is divided evenly 
 * between the \f$N\f$ available nodes, such that each node contains only 
 * \f$2^{\text{numQubits}}/N\f$ complex values. This is assigned to 
 * DiagonalOp.numElemsPerChunk.
 * 
 * Users must therefore exercise care in modifying DiagonalOp.real and DiagonalOp.imag 
 * directly.
 *
 * For example, the following is valid code when when distributed between <b>N = 2</b> nodes:
 * ```
 *      // create diag({1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16})
 *      int numQubits = 4;
 *      DiagonalOp op = createDiagonalOp(numQubits4, env);
 *      for (int i=0; i<8; i++) {
 *          if (env.rank == 0)
 *              op.real[i] = (i+1);
 *          if (env.rank == 1)
 *              op.real[i] = (i+1+8);
 *      }
 * ```
 * \n
 *
 *
 * @see 
 * - createDiagonalOpFromPauliHamilFile()
 * - setDiagonalOpElems()
 * - initDiagonalOp()
 * - syncDiagonalOp()
 * - applyDiagonalOp()
 * - calcExpecDiagonalOp()
 * - destroyDiagonalOp()
 *
 * @ingroup type
 * @returns a dynamic DiagonalOp instance initialised to diag(0,0,...).
 * @param[in] numQubits number of qubits which inform the Hilbert dimension of the operator.
 * @param[in] env the ::QuESTEnv
 * @throws invalidQuESTInputError() 
 * - if \p numQubits <= 0
 * - if \p numQubits is so large that the number of elements cannot fit in a long long int type, 
 * - if in distributed mode, there are more nodes than elements in the operator
 * @throws exit 
 * - if the memory could not be allocated
 * @author Tyson Jones
 */
DiagonalOp createDiagonalOp(int numQubits, QuESTEnv env);

/** Destroys a ::DiagonalOp created with createDiagonalOp(), freeing its memory.
 *
 * @see 
 * - createDiagonalOp()
 *
 * @ingroup type
 * @param[in] op the ::DiagonalOp to destroy
 * @param[in] env the ::QuESTEnv
 * @throws invalidQuESTInputError()
 * - if \p op was not previously created
 * @author Tyson Jones
 */
void destroyDiagonalOp(DiagonalOp op, QuESTEnv env);

/** Update the GPU memory with the current values in `op.real` and `op.imag`.
 *
 * This is required after making manual changes to \p op, but is not required 
 * after functions initDiagonalOp() and setDiagonalOpElems().
 *
 * This function has no effect in other modes besides GPU mode.
 *
 * @see 
 * - ::DiagonalOp
 * - initDiagonalOp()
 * - setDiagonalOpElems()
 *
 * @ingroup type
 * @param[in,out] op the ::DiagonalOp to synch to GPU
 * @throws invalidQuESTInputError()
 * - if \p op was not created
 * @author Tyson Jones
 */
void syncDiagonalOp(DiagonalOp op);

/** Overwrites the entire ::DiagonalOp \p op with the given \p real and \p imag 
 * complex elements.
 *
 * Both \p real and \p imag must have length equal to <b>pow(2, </b>`op.numQubits`<b>)</b>.
 *
 * In GPU mode, this updates both the RAM (\p op.real and \p op.imag) <em>and</em> 
 * persistent GPU memory; there is no need to call syncDiagonalOp() afterward.
 *
 * In distributed mode, this function assumes \p real and \p imag exist fully on every 
 * node. For ::DiagonalOp which are too large to fit into a single node, use 
 * setDiagonalOpElems() or syncDiagonalOp().
 *
 * @see 
 * - setDiagonalOpElems()
 * - initDiagonalOpFromPauliHamil()
 *
 * @ingroup type
 * @param[in,out] op the diagonal operator to modify
 * @param[in] real the real components of the full set of new elements
 * @param[in] imag the imaginary components of the full set of new elements
 * @throws invalidQuESTInputError()
 * - if \p op was not created
 * @throws segmentation-fault
 * - if either \p real or \p imag have length smaller than <b>pow(2, </b>`op.numQubits`<b>)</b>
 * @author Tyson Jones
 */
void initDiagonalOp(DiagonalOp op, qreal* real, qreal* imag);

/** Populates the diagonal operator \p op to be equivalent to the given Pauli 
 * Hamiltonian \p hamil, assuming \p hamil contains only `PAULI_Z` operators.
 *
 * Given a ::PauliHamil \p hamil featuring only `PAULI_Z` and `PAULI_I`, with 
 * term coefficients \f$\{\lambda_j\}\f$, which 
 * hence has form
 * \f[
 *      \begin{aligned}
 *      \text{hamil} &= \sum\limits_j^{\text{numSumTerms}} \lambda_j
 *                      \bigotimes\limits_{k_j} \hat{Z}_k \\
 *      &\equiv \begin{pmatrix}
 *              r_1 \\ & r_2 \\ & & r_3 \\ & & & \ddots \\ & & & & r_{2^{\,\text{numQubits}}}
 *          \end{pmatrix},
 *      \end{aligned}
 * \f]
 * this function modifies \p op to 
 * \f[
 *      \text{op} \; \rightarrow \; \text{diag}
 *          \big( \; r_1, \; r_2, \; r_3, \; \dots, \; r_{2^{\,\text{numQubits}}} \, \big),
 * \f]
 * where the real amplitudes have form 
 * \f[
 *       r_i = \sum\limits_j \, s_{ij} \, \lambda_j, \;\;\;\; s_{ij} = \pm  1 \,.
 * \f]
 * This is useful since calculations with ::DiagonalOp are significantly faster than 
 * the equivalent calculations with a general ::PauliHamil. For example, 
 * applyDiagonalOp() requires a factor `numSumTerms * numQubits` fewer operations
 * than applyPauliHamil().
 * 
 * > In distributed mode, each node will contain only a sub-partition of the full diagonal matrix.  
 * > In GPU mode, both the CPU and GPU memory copies of \p op will be updated, so there 
 * > is no need to call syncDiagonalOp() afterward.
 *
 * @see
 * - createDiagonalOp()
 * - createPauliHamil()
 * - createDiagonalOpFromPauliHamilFile()
 * - initDiagonalOp()
 * - setDiagonalOpElems()
 *
 * @ingroup type
 * @param[in,out] op an existing ::DiagonalOp (e.g. created with createDiagonalOp()) to modify 
 * @param[in] hamil a ::PauliHamil of equal dimension to \p op, containing only `PAULI_Z` and `PAULI_I` operators
 * @throws invalidQuESTInputError()
 * - if \p hamil has invalid parameters (\p numQubits <= 0, \p numSumTerms <= 0)
 * - if \p op and \p hamil have unequal dimensions 
 * - if \p hamil contains any operator other than `PAULI_Z` and `PAULI_I`
 * @throws segmentation-fault 
 * - if either \p op or \p hamil have not been already created
 * @author Tyson Jones 
 * @author Milos Prokop (serial prototype)
 */
void initDiagonalOpFromPauliHamil(DiagonalOp op, PauliHamil hamil);

/** Creates and initialiases a diagonal operator from the Z Pauli Hamiltonian encoded in 
 * file with filename \p fn.
 *
 * This is a convenience function to prepare a diagonal operator from a plaintext 
 * description of an all-Z Pauli Hamiltonian. The returned ::DiagonalOp is a
 * distributed data structure, and significantly faster to use (through functions 
 * like calcExpecDiagonalOp()) than ::PauliHamil functions (like calcExpecPauliHamil()).
 *
 * - See createDiagonalOp() for info about the returned operator.
 * - See initDiagonalOpFromPauliHamil() for info about the initialised state.
 * - See createPauliHamilFromFile() for info about the required file format.
 *
 * > The returned ::DiagonalOp must be later freed with destroyDiagonalOp().  
 * > Note a ::PauliHamil from \p fn is temporarily internally created.
 *
 * This function is equivalent to
 * ```
 * // produce diagonal matrix d
 * PauliHamil h = createPauliHamilFromFile(fn);
 * DiagonalOp d = createDiagonalOp(h.numQubits, env);
 * initDiagonalOpFromPauliHamil(d, h); 
 * destroyPauliHamil(h);
 * ```
 * <br>
 *
 * @see
 * - initDiagonalOpFromPauliHamil()
 * - createPauliHamilFromFile()
 * - createDiagonalOp()
 * - destroyDiagonalOp()
 *
 * @ingroup type 
 * @param[in] fn filename of a plaintext file encoding an all-Z Pauli Hamiltonian
 * @param[in] env the session ::QuESTEnv
 * @returns a created ::DiagonalOp equivalent to the Hamiltonian in \p fn 
 * @throws invalidQuESTInputError()
 * - if file \p fn cannot be read
 * - if file \p fn does not encode a valid ::PauliHamil
 * - if the encoded ::PauliHamil consists of operators other than `PAULI_Z` and `PAULI_I`
 * @author Tyson Jones
 * @author Milos Prokop (serial prototype)
 */ 
DiagonalOp createDiagonalOpFromPauliHamilFile(char* fn, QuESTEnv env);

/** Modifies a subset (starting at index \p startInd, and ending at index 
 * \p startInd <b>+</b> \p numElems) of the elements in ::DiagonalOp \p op 
 * with the given complex numbers (passed as \p real and \p imag components).
 *
 * In GPU mode, this updates both the RAM (\p op.real and \p op.imag), and the 
 * persistent GPU memory.
 *
 * In distributed mode, this function assumes the subset \p real and \p imag exist
 * (at least) on the node containing the ultimately updated elements.\n
 * For example, below is the correct way to modify the full 8 elements of \p op 
 * when split between 2 nodes.
 * ```
 *     DiagonalOp op = createDiagonalOp(3, env);
 *     
 *     int numElems = 4;
 *     qreal re[] = {1,2,3,4};
 *     qreal im[] = {1,2,3,4};
 *     setDiagonalOpElems(op, 0, re, im, numElems);
 *     
 *     // modify re and im to the next set of elements 
 *     for (int i=0; i<4; i++) {
 *       re[i] += 4;
 *       im[i] += 4;
 *     }
 *     setDiagonalOpElems(op, 4, re, im, numElems);
 * ```
 *
 * In this way, one can avoid a single node containing all new elements which might 
 * not fit. If more elements are passed than exist on an individual node, each 
 * node merely ignores the additional elements.
 *
 * @ingroup type
 * @param[in,out] op the ::DiagonalOp to modify
 * @param[in] startInd the starting index (globally) of the subset of elements to modify
 * @param[in] real  the real components of the new elements
 * @param[in] imag  the imaginary components of the new elements
 * @param[in] numElems the number of new elements (the length of \p real and \p imag)
 * @throws invalidQuESTInputError()
 * - if \p op was not created
 * - if \p startInd is an invalid index (<0 or >=<b>pow(2,</b>`op.numQubits`<b>)</b>)
 * - if \p numElems is an invalid number of elements (<=0 or ><b>pow(2,</b>`op.numQubits`<b>)</b>)
 * - if there are fewer than \p numElems elements in \p op after \p startInd
 * @throws segmentation-fault 
 * - if either \p real or \p imag have fewer elements than \p numElems
 * @author Tyson Jones
 */
void setDiagonalOpElems(DiagonalOp op, long long int startInd, qreal* real, qreal* imag, long long int numElems);

/** Apply a diagonal operator, which is possibly non-unitary and non-Hermitian,
 * to the entire \p qureg. 
 *
 * Let \f$d_j = \text{op.real}[j] + (\text{op.imag}[j])\,\text{i} \f$, and 
 * \f[
 *  \hat{D} = \begin{pmatrix}
 *  d_0 \\
 *  & d_1 \\
 *  & & \ddots \\
 *  & & & d_{2^{\text{op.numQubits}}-1}
 *  \end{pmatrix}.
 * \f]
 * If \p qureg is a state-vector \f$|\psi\rangle\f$, this function performs
 *  \f$|\psi\rangle \rightarrow \hat{D} \, |\psi\rangle\f$. \n
 * If \p qureg is a density-matrix \f$\rho\f$, this function performs
 *  \f$\rho \rightarrow \hat{D}\, \rho\f$. Notice this has not applied \f$\hat{D}\f$
 * in the fashion of a unitary.
 * 
 * > If your operator is unitary with unit amplitudes, the phases of which can be 
 * > described by an analytic expression, you should instead use applyPhaseFunc()
 * > or applyNamedPhaseFunc() for significant memory and runtime savings.
 *
 * @see
 * - createDiagonalOp()
 * - calcExpecDiagonalOp()
 * - applyPhaseFunc()
 * - applyNamedPhaseFunc()
 *
 * @ingroup operator
 * @param[in,out] qureg the state to operate the diagonal operator upon
 * @param[in] op the diagonal operator to apply
 * @throws invalidQuESTInputError()
 * - if \p op was not created
 * - if \p op acts on a different number of qubits than \p qureg represents
 * @author Tyson Jones
 */
void applyDiagonalOp(Qureg qureg, DiagonalOp op);

/** Computes the expected value of the diagonal operator \p op for state \p qureg. 
 * Since \p op is not necessarily Hermitian, the expected value may be a complex 
 * number.
 *
 * Let \f$ D \f$ be the diagonal operator \p op, with diagonal elements \f$ d_i \f$.
 * Then if \p qureg is a state-vector \f$|\psi\rangle \f$, this function computes
 * \f[
    \langle \psi | D | \psi \rangle = \sum_i |\psi_i|^2 \, d_i
 * \f]
 * If \p qureg is a density matrix \f$ \rho \f$, this function computes 
 * \f[ 
    \text{Trace}( D \rho ) = \sum_i \rho_{ii} \, d_i
 * \f]
 *
 * @see 
 * - createDiagonalOp()
 * - applyDiagonalOp()
 * - calcExpecPauliSum()
 * - calcExpecPauliProd()
 * - calcExpecPauliHamil()
 *
 * @ingroup calc
 * @param[in] qureg a state-vector or density matrix
 * @param[in] op    the diagonal operator to compute the expected value of
 * @return the expected vaulue of the operator
 * @throws invalidQuESTInputError()
 * - if \p op was not created
 * - if \p op acts on a different number of qubits than \p qureg represents
 * @author Tyson Jones
 */
Complex calcExpecDiagonalOp(Qureg qureg, DiagonalOp op);

/** Print the current state vector of probability amplitudes for a set of qubits to file.
 * File format:
 * @verbatim
real, imag
realComponent1, imagComponent1
realComponent2, imagComponent2
...
realComponentN, imagComponentN
@endverbatim
 *
 * File naming convention:
 *
 * For each node that the program runs on, a file 'state_rank_[node_rank].csv' is generated. If there is
 * more than one node, ranks after the first do not include the header
 * @verbatim
real, imag
@endverbatim
 * so that files are easier to combine.
 *
 * @ingroup debug
 * @param[in,out] qureg object representing the set of qubits
 * @author Ania Brown
 */
void reportState(Qureg qureg);

/** Print the current state vector of probability amplitudes for a set of qubits to standard out. 
 * For debugging purposes. Each rank should print output serially. 
 * Only print output for systems <= 5 qubits
 *
 * @ingroup debug
 * @author Ania Brown
 */
void reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank);

/** Report metainformation about a set of qubits: number of qubits, number of probability amplitudes.
 *
 * @ingroup debug
 * @param[in] qureg object representing the set of qubits
 * @author Ania Brown
 */
void reportQuregParams(Qureg qureg);

/** Print the \p PauliHamil to screen. 
 * The output features a new line for each term, each with format 
 *
 *     c p1 p2 p3 ... pN
 *
 * where \p c is the real coefficient of the term, and \p p1 ... \p pN are 
 * numbers \p 0, \p 1, \p 2, \p 3 to indicate identity, pauliX, pauliY and pauliZ 
 * operators respectively, acting on qubits \p 0 through \p N-1 (all qubits).
 * A tab character separates c and p1, but single spaces separate the Pauli operators.
 *
 * @see
 * - createPauliHamil()
 * - initPauliHamil()
 * - createPauliHamilFromFile()
 *
 * @ingroup debug 
 * @param[in] hamil an instantiated PauliHamil
 * @throws invalidQuESTInputError if the parameters of \p hamil are invalid, i.e. 
 *      if \p numQubits <= 0, or if \p numSumTerms <= 0, or if \p pauliCodes 
 *      contains an invalid Pauli code.
 * @author Tyson Jones
 */ 
void reportPauliHamil(PauliHamil hamil);

/** Returns the number of qubits represented by \p qureg.
 *
 * @see 
 * - getNumAmps()
 *
 * @ingroup calc
 * @param[in] qureg a state-vecor or density matrix
 * @return `qureg.numQubitsRepresented`
 * @author Tyson Jones
 */
int getNumQubits(Qureg qureg);

/** Returns the number of complex amplitudes in a state-vector \p qureg.
 *
 * In distributed mode, this returns the total number of amplitudes in the full 
 * representation of \p qureg, and so may be larger than the number stored on 
 * each node. For the latter, refer to Qureg.numAmpsPerChunk.
 *
 * @see
 * - getNumQubits()
 *
 * @ingroup calc
 * @param[in] qureg a state-vecotor
 * @return the total number of complex amplitudes representing \p qureg.
 * @throws invalidQuESTInputError()
 * - if \p qureg is a density matrix
 * @author Tyson Jones
 */
long long int getNumAmps(Qureg qureg);

/** Initialises a qureg to have all-zero-amplitudes. This is an unphysical state 
 * useful for iteratively building a state with functions like setWeightedQureg(),
 * and should not be confused with initZeroState().
 *
 * @ingroup init
 * @param[in,out] qureg a ::Qureg of which to clear all amplitudes
 * @author Tyson Jones
 */
void initBlankState(Qureg qureg);

/** Initialise \p qureg into the zero state.
 *
 * If \p qureg is a state-vector of \f$N\f$ qubits, it is modified to state 
 * \f$ {| 0 \rangle}^{\otimes N} \f$. \n
 * If \p qureg is a density matrix of \f$ N \f$ qubits, it is modified to state 
 * \f$ {| 0 \rangle\langle 0 |}^{\otimes N} \f$
 * 
 * @ingroup init
 * @param[in,out] qureg the object representing the set of all qubits to initialise
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void initZeroState(Qureg qureg);

/** Initialise \p qureg into the plus state.
 * 
 * If \p qureg is a state-vector of \f$N\f$ qubits, it is modified to state
 * \f[ 
 *   {| + \rangle}^{\otimes N} = \frac{1}{\sqrt{2^N}} (| 0 \rangle + | 1 \rangle)^{\otimes N}\,.
 * \f]
 * If \p qureg is a density matrix of \f$N\f$ qubits, it is modified to state
 * \f[ 
 *   {| + \rangle\langle+|}^{\otimes N} = \frac{1}{{2^N}} \sum_i\sum_j |i\rangle\langle j|\,.
 * \f]
 *
 *
 * @ingroup init
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void initPlusState(Qureg qureg);

/** Initialise \p qureg into the classical state (also known as a 
 * "computational basis state") with index \p stateInd. 
 *
 * If \p qureg is a state-vector, it will become 
 *      \f$ | \text{stateInd} \rangle \f$. \n
 * If \p qureg is a density matrix, it will become
 *      \f$ | \text{stateInd} \rangle \langle \text{stateInd} | \f$.
 *
 * Classical states are indexed from zero, so that \p stateInd <b>= 0</b> produces
 * \f$ | 00 \dots 00 \rangle \f$,
 * and  \p stateInd <b>= 1</b> produces \f$ | 00 \dots 01 \rangle \f$, and 
 * \p stateInd <b>=</b> \f$ 2^N - 1 \f$ produces 
 * \f$ | 11 \dots 11 \rangle \f$.
 * 
 * Subsequent calls to getProbAmp() will yield 0 for all indices except \p stateInd,
 * and the phase of \p stateInd's amplitude will be 1 (real).
 *
 * This function can be used to initialise \p qureg into a specific binary state
 * (e.g. \p 11001) using a binary literal (supported by only some compilers):
 * ```
 *      initClassicalState(qureg, 0b11001);
 * ```
 * \n
 *
 *
 * @see
 * - initPureState()
 *
 * @ingroup init
 * @param[in,out] qureg the ::Qureg to modify
 * @param[in] stateInd the index of the basis state to modify \p qureg into
 * @throws invalidQuESTInputError()
 * - if \p stateInd is outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * @author Tyson Jones
 */
void initClassicalState(Qureg qureg, long long int stateInd);

/** Initialise \p qureg into to a given pure state of an equivalent Hilbert dimension.
 * 
 * If \p qureg is a state-vector, this merely clones \p pure into \p qureg. \n
 * If \p qureg is a density matrix, this makes \p qureg 100% likely to be in the \p pure state.
 *
 * ::Qureg \p pure is not modified.
 *
 * @see
 * - initClassicalState()
 *
 * @ingroup init
 * @param[in,out] qureg the ::Qureg to modify
 * @param[in] pure a state-vector containing the pure state into which to initialise \p qureg
 * @throws invalidQuESTInputError()
 * - if \p qureg and \p pure have mismatching dimensions
 * - if \p pure is a density matrix
 * @author Tyson Jones
 */
void initPureState(Qureg qureg, Qureg pure);

/** Initialises \p qureg to be in the un-normalised, non-physical state with 
 * with \f$n\f$-th complex amplitude given by \f$2n/10 + i(2n+1)/10\f$. 
 *
 * This is used internally for debugging and testing.
 * 
 * @ingroup debug
 * @param[in,out] qureg the register to have its amplitudes overwritten
 * @author Ania Brown
 * @author Tyson Jones (doc)
 */
void initDebugState(Qureg qureg);

/** Initialise \p qureg by specifying all amplitudes.
 * For density matrices, it is assumed the amplitudes have been flattened 
 * column-wise into the given arrays.
 *
 * The real and imaginary components of the amplitudes are passed in separate arrays,
 * \p reals and \p imags,
 * each of which must have length `qureg.numAmpsTotal`.
 * There is no automatic checking that the passed arrays are L2 normalised, so this 
 * can be used to prepare \p qureg in a non-physical state.
 *
 * In distributed mode, this would require the complete state to fit in 
 * every node. To manually prepare a state for which all amplitudes cannot fit into a single node,
 * use setAmps()
 *
 * @see
 * - setAmps()
 *
 * @ingroup init
 * @param[in,out] qureg the ::Qureg to overwrite
 * @param[in] reals array of the real components of the new amplitudes
 * @param[in] imags array of the imaginary components of the new amplitudes
 * @throws segmentation-fault
 * - if either \p reals or \p imags have fewer than `qureg.numAmpsTotal` elements
 * @author Tyson Jones
 */
void initStateFromAmps(Qureg qureg, qreal* reals, qreal* imags);

/** Overwrites a subset of the amplitudes in state-vector \p qureg, with those passed in \p reals and \p imags.
 *
 * Only amplitudes with indices in <b>[</b>\p startInd<b>,</b> \p startInd <b>+</b> \p numAmps<b>]</b> 
 * will be changed. The resulting \p qureg may not necessarily be in an L2 normalised state.
 *
 * In distributed mode, this function assumes the subset \p reals and \p imags exist
 * (at least) on the node containing the ultimately updated elements.\n
 * For example, below is the correct way to modify the full 8 elements of \p qureg 
 * when split between 2 nodes.
 * ```
 *     Qureg qureg = createQureg(3, env);
 *     
 *     long long int numAmps = 4;
 *     qreal re[] = {1,2,3,4};
 *     qreal im[] = {1,2,3,4};
 *     setAmps(qureg, 0, re, im, numAmps);
 *     
 *     // modify re and im to the next set of elements 
 *     for (int i=0; i<4; i++) {
 *       re[i] += 4;
 *       im[i] += 4;
 *     }
 *     setAmps(qureg, 4, re, im, numAmps);
 * ```
 * \n
 *
 *
 * @see
 * - setWeightedQureg()
 * - initStateFromAmps()
 * - initBlankState()
 *
 * @ingroup init
 * @param[in,out] qureg the state-vector to modify
 * @param[in] startInd the index of the first amplitude in \p qureg to modify
 * @param[in] reals array of the real components of the new amplitudes
 * @param[in] imags array of the imaginary components of the new amplitudes
 * @param[in] numAmps the length of each of the reals and imags arrays.
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a state-vector (i.e. is a density matrix)
 * - if \p startInd is outside [0, `qureg.numAmpsTotal`]
 * - if \p numAmps is outside [0, `qureg.numAmpsTotal`]
 * - if \p numAmps + \p startInd >= `qureg.numAmpsTotal`
 * @author Tyson Jones
 */
void setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps);

/** Overwrite the amplitudes of \p targetQureg with those from \p copyQureg. 
 * 
 * Registers must either both be state-vectors, or both be density matrices, and 
 * of equal dimensions.
 * Only the quantum state is cloned, while auxilary info (like recorded QASM) is unchanged.
 * copyQureg is unaffected.
 *
 * @see
 * - createCloneQureg()
 * - initPureState()
 * - setWeightedQureg()
 *
 * @ingroup init
 * @param[in, out] targetQureg the qureg to have its quantum state overwritten
 * @param[in] copyQureg the qureg to have its quantum state cloned into targetQureg.
 * @throws invalidQuESTInputError()
 * - if \p targetQureg is a state-vector while \p copyQureg is a density matrix (and vice versa)
 * - if \p targetQureg and \p copyQureg have different dimensions
 * @author Tyson Jones
 */
void cloneQureg(Qureg targetQureg, Qureg copyQureg);

/** Shift the phase between \f$ |0\rangle \f$ and \f$ |1\rangle \f$ of a single qubit by a given angle.
 * 
 * > This is equivalent to a Z-axis rotation of the Bloch-sphere up to a global phase factor.
 * For angle \f$\theta\f$, this effects single-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 & 0 \\
 * 0 & \exp(i \theta)
 * \end{pmatrix}
 * \f] 
 * with circuit diagram 
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-4, 0) {targetQubit};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_\theta$};
                \end{tikzpicture}
    \f]
 * 
 * @see
 * - controlledPhaseShift()
 * - multiControlledPhaseShift()
 * 
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to undergo a phase shift
 * @param[in] angle amount by which to shift the phase in radians
 * @throws invalidQuESTInputError()
 * - \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Tyson Jones
 */
void phaseShift(Qureg qureg, int targetQubit, qreal angle);

/** Introduce a phase factor \f$ \exp(i \theta) \f$ on state \f$ |11\rangle \f$ of qubits
 * \p idQubit1 and \p idQubit2.
 * For angle \f$\theta\f$, this effects the unitary
 * \f[
 * \begin{pmatrix}
 * 1 & & & \\
 * & 1 & & \\
 * & & 1 & \\
 * & & & \exp(i \theta)
 * \end{pmatrix}
 * \f] 
 * on \p idQubit1 and \p idQubit2.
 *     
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {qubit1};
                \node[draw=none] at (-3.5, 0) {qubit2};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_\theta$};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - phaseShift()
 * - multiControlledPhaseShift()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] idQubit1 first qubit in the state to phase shift
 * @param[in] idQubit2 second qubit in the state to phase shift
 * @param[in] angle amount by which to shift the phase in radians
 * @throws invalidQuESTInputError()
 * - if \p idQubit1 or \p idQubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if \p idQubit1 and \p idQubit2 are equal
 * @author Tyson Jones
 */
void controlledPhaseShift(Qureg qureg, int idQubit1, int idQubit2, qreal angle);

/** Introduce a phase factor \f$ \exp(i \theta) \f$ on state \f$ |1 \dots 1 \rangle \f$
 * of the passed qubits.
 *     
   \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {controls};
                \node[draw=none] at (1, .7) {$\theta$};
                
                \node[draw=none] at (0, 6) {$\vdots$};
                \draw (0, 5) -- (0, 4);
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw (0, 4) -- (0, 2); 
                
                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 0);
                
                \draw (-2,0) -- (2, 0);
                \draw[fill=black] (0, 0) circle (.2);
                \end{tikzpicture}
   \f]
 *
 * @see
 * - phaseShift()
 * - controlledPhaseShift()
 * - controlledPhaseFlip()
 * - multiControlledPhaseFlip()
 * 
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits array of qubits to phase shift
 * @param[in] numControlQubits the length of array \p controlQubits
 * @param[in] angle amount by which to shift the phase in radians
 * @throws invalidQuESTInputError()
 * - if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented])
 * - if any qubit index in \p controlQubits is outside [0, \p qureg.numQubitsRepresented])
 * - if the qubits in \p controlQubits are not unique
 * @author Tyson Jones
 */
void multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle);

/** Apply the (two-qubit) controlled phase flip gate, also known as the controlled pauliZ gate.
 * For each state, if both input qubits have value one, multiply the amplitude of that state by -1. This applies the two-qubit unitary:
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & 1 \\
 * & & & -1 
 * \end{pmatrix}
 * \f]
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {idQubit1};
                \node[draw=none] at (-3.5, 0) {idQubit2};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 0);
                
                \draw (-2,0) -- (2, 0);
                \draw[fill=black] (0, 0) circle (.2);
                \end{tikzpicture}
    \f]
 *
 * @see
 * - pauliZ()
 * - phaseShift()
 * - multiControlledPhaseFlip()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] idQubit1, idQubit2 qubits to operate upon
 * @throws invalidQuESTInputError()
 * - if \p idQubit1 or \p idQubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if \p idQubit1 and \p idQubit2 are equal
 * @author Tyson Jones
 */
void controlledPhaseFlip (Qureg qureg, int idQubit1, int idQubit2);

/** Apply the multiple-qubit controlled phase flip gate, also known as the multiple-qubit controlled pauliZ gate.
 * For each state, if all control qubits have value one, multiply the amplitude of that state by -1. This applies the many-qubit unitary:
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & \ddots \\
 * & & & 1 \\
 * & & & & -1 
 * \end{pmatrix}
 * \f]
 * on the control qubits.
 *
 * \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {controls};
                
                \node[draw=none] at (0, 6) {$\vdots$};
                \draw (0, 5) -- (0, 4);
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw (0, 4) -- (0, 2); 
                
                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 0);
                
                \draw (-2,0) -- (2, 0);
                \draw[fill=black] (0, 0) circle (.2);
                \end{tikzpicture}
   \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits array of input qubits
 * @param[in] numControlQubits number of input qubits
 * @throws invalidQuESTInputError()
 * - if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented)
 * - if any qubit in \p controlQubits is outside [0, \p qureg.numQubitsRepresented)
 * - if any qubit in \p qubits is repeated
 * @author Tyson Jones
 */
void multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits);

/** Apply the single-qubit S gate.
 * This is a rotation of \f$\pi/2\f$ around the Z-axis on the Bloch sphere, or the unitary:
 * \f[
 * \begin{pmatrix}
 * 1 & 0 \\
 * 0 & i
 * \end{pmatrix}
 * \f]
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {S};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - tGate()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate upon
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void sGate(Qureg qureg, int targetQubit);

/** Apply the single-qubit T gate.
 * This is a rotation of \f$\pi/4\f$ around the Z-axis on the Bloch sphere, or the unitary:
 * \f[
 * \begin{pmatrix}
 * 1 & 0 \\
 * 0 & \exp\left(i \frac{\pi}{4}\right)
 * \end{pmatrix}
 * \f]
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {T};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - sGate()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate upon
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void tGate(Qureg qureg, int targetQubit);

/** Create the QuEST execution environment.
 * This should be called only once, and the environment should be freed with destroyQuESTEnv at the end
 * of the user's code.
* If something needs to be done to set up the execution environment, such as 
 * initializing MPI when running in distributed mode, it is handled here.
 *
 * @see
 * - reportQuESTEnv()
 * - destroyQuESTEnv()
 * - syncQuESTEnv()
 *
 * @ingroup type
 * @return object representing the execution environment. A single instance is used for each program
 * @author Ania Brown
 */
QuESTEnv createQuESTEnv(void);

/** Destroy the QuEST environment. 
 * If something needs to be done to clean up the execution environment, such as 
 * finalizing MPI when running in distributed mode, it is handled here
 *
 * @see
 * - createQuESTEnv()
 *
 * @ingroup type
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * @author Ania Brown
 */
void destroyQuESTEnv(QuESTEnv env);

/** Guarantees that all code up to the given point has been executed on all nodes (if running in distributed mode)
 *
 * @see
 * - syncQuESTSuccess()
 *
 * @ingroup debug
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * @author Ania Brown
 */
void syncQuESTEnv(QuESTEnv env);

/** Performs a logical AND on all successCodes held by all processes. If any one process has a zero successCode
 * all processes will return a zero success code.
 *
 * @ingroup debug
 * @param[in] successCode 1 if process task succeeded, 0 if process task failed
 * @returns 1 if all processes succeeded, 0 if any one process failed
 * @author Ania Brown
 */ 
int syncQuESTSuccess(int successCode);

/** Report information about the QuEST environment
 *
 * @ingroup debug
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * @author Ania Brown
 */
void reportQuESTEnv(QuESTEnv env);

/** Sets \p str to a string containing information about the runtime environment, 
 * including whether simulation is using CUDA (for GPU), OpenMP (for multithreading)
 * and/or MPI (for distribution). The number of CPU threads and distributed ranks is
 * also reported. Note there is currently no reporting of the number of GPU cores used.
 *
 * The string format is:
 * ```
 * "CUDA=b OpenMP=b MPI=b threads=n ranks=n"
 * ```
 * where <b>b</b> is 0 or 1, and <b>n</b> are integers.
 *
 * @ingroup debug 
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * @param[out] str to be populated with the output string
 * @author Ania Brown
 * @author Tyson Jones
 */
void getEnvironmentString(QuESTEnv env, char str[200]);

/** In GPU mode, this copies the state-vector (or density matrix) from RAM 
 * (qureg.stateVec) to VRAM / GPU-memory (qureg.deviceStateVec), which is the version 
 * operated upon by other calls to the API. 
 * In CPU mode, this function has no effect.
 * In conjunction with copyStateFromGPU() (which should be called first), this allows 
 * a user to directly modify the state-vector in a harware agnostic way.
 * Note though that users should instead use setAmps() if possible.
 *
 * For example, to set the first real element to 1, one could do:
 * ```
 *     copyStateFromGPU(qureg);
 *     qureg.stateVec.real[0] = 1;
 *     copyStateToGPU(qureg);
 * ```
 *
 * Note users should never access qureg.deviceStateVec directly.
 *
 * @see
 * - copyStateFromGPU()
 *
 * @ingroup debug
 * @param[in, out] qureg the qureg of which to copy `.stateVec` to `.deviceStateVec` in GPU mode
 * @author Ania Brown 
 * @author Tyson Jones (doc)
 */
void copyStateToGPU(Qureg qureg);

/** In GPU mode, this copies the state-vector (or density matrix) from GPU memory 
 * (qureg.deviceStateVec) to RAM (qureg.stateVec), where it can be accessed/modified 
 * by the user.
 * In CPU mode, this function has no effect.
 * In conjunction with copyStateToGPU(), this allows a user to directly modify the 
 * state-vector in a harware agnostic way.
 * Note though that users should instead use setAmps() if possible.
 *
 * For example, to set the first real element to 1, one could do:
 * ```
 *     copyStateFromGPU(qureg);
 *     qureg.stateVec.real[0] = 1;
 *     copyStateToGPU(qureg);
 * ```
 *
 * Note users should never access qureg.deviceStateVec directly.
 *
 * @see
 * - copyStateToGPU()
 *
 * @ingroup debug
 * @param[in, out] qureg the qureg of which to copy `.deviceStateVec` to `.stateVec` in GPU mode
 * @author Ania Brown 
 * @author Tyson Jones (doc)
 */
void copyStateFromGPU(Qureg qureg);

/** Get the complex amplitude at a given index in the state vector.
 *
 * @see
 * - getDensityAmp()
 * - getRealAmp()
 * - getImagAmp()
 * - getProbAmp()
 * - getNumAmps()
 * - getNumQubits()
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return amplitude at index, returned as a Complex struct (with .real and .imag attributes)
 * @throws invalidQuESTInputError()
 * - if \p qureg is a density matrix
 * - if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Tyson Jones
 */
Complex getAmp(Qureg qureg, long long int index);

/** Get the real component of the complex probability amplitude at an index in the state vector.
 *
 * @see
 * - getAmp()
 * - getDensityAmp()
 * - getImagAmp()
 * - getProbAmp()
 * - getNumAmps()
 * - getNumQubits()
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return real component at that index
 * @throws invalidQuESTInputError()
 * - if \p qureg is a density matrix
 * - if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Ania Brown
 */
qreal getRealAmp(Qureg qureg, long long int index);

/** Get the imaginary component of the complex probability amplitude at an index in the state vector.
 *
 * @see
 * - getAmp()
 * - getDensityAmp()
 * - getRealAmp()
 * - getProbAmp()
 * - getNumAmps()
 * - getNumQubits()
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return imaginary component at that index
 * @throws invalidQuESTInputError()
 * - if \p qureg is a density matrix
 * - if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Ania Brown
 */
qreal getImagAmp(Qureg qureg, long long int index);

/** Get the probability of a state-vector at an index in the full state vector.
 *
 * @see
 * - getAmp()
 * - getDensityAmp()
 * - getRealAmp()
 * - getImagAmp()
 * - getNumAmps()
 * - getNumQubits()
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return realEl*realEl + imagEl*imagEl
 * @throws invalidQuESTInputError()
 * - if \p qureg is a density matrix
 * - if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Ania Brown
 */
qreal getProbAmp(Qureg qureg, long long int index);

/** Get an amplitude from a density matrix at a given row and column.
 *
 * @see
 * - getAmp()
 * - getRealAmp()
 * - getImagAmp()
 * - getProbAmp()
 * - getNumAmps()
 * - getNumQubits()
 *
 * @ingroup calc
 * @param[in] qureg object representing a density matrix
 * @param[in] row row of the desired amplitude in the density matrix
 * @param[in] col column of the desired amplitude in the density matrix
 * @return a Complex scalar representing the desired amplitude
 * @throws invalidQuESTInputError()
 * - if \p qureg is a state-vector, 
 * - if \p row or \p col are outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Tyson Jones
 */
Complex getDensityAmp(Qureg qureg, long long int row, long long int col);

/** A debugging function which calculates the probability of the qubits in \p qureg 
 * being in any state, which should always be 1 for correctly normalised states 
 * (hence returning a real number).
 * For state-vectors \f$ \psi \f$, this is the norm of the entire state-vector 
 * (the sum of the absolute-value-squared of every amplitude):
 * \f[
 *      \sum\limits_i |\psi_i|^2
 * \f]
 * and for density matrices \f$ \rho \f$, it is the trace:
 * \f[
 *      \text{Trace}(\rho) = \sum\limits_i \rho_{i,i} \;
 * \f]
 *
 * For un-normalised density matrices (those directly modified or initialised by the user), 
 * this function returns the real component of the trace.
 *
 * Note this calculation utilises Kahan summation for greater accuracy, and hence is
 * not parallelised and so will be slower than other functions.
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @return the total probability of the qubits in \p qureg being in any state
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
qreal calcTotalProb(Qureg qureg);

/** Apply a single-qubit unitary parameterised by two given complex scalars.
 * Given valid complex numbers \f$\alpha\f$ and \f$\beta\f$, applies the unitary
 * \f[
 * U =
 * \begin{pmatrix}
 * \alpha & -\beta^* \\
 * \beta & \alpha^*
 * \end{pmatrix}
 * \f]
 * which is general up to a global phase factor.               
 * Valid \f$\alpha\f$, \f$\beta\f$ satisfy \f$|\alpha|^2 + |\beta|^2 = 1\f$. 
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - controlledCompactUnitary()
 * - unitary()
 * - twoQubitUnitary()
 * - multiQubitUnitary()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] alpha complex unitary parameter (row 1, column 1)
 * @param[in] beta complex unitary parameter (row 2, column 1)
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if \p alpha, \p beta don't satisfy |`alpha`|^2 + |`beta`|^2 = 1
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void compactUnitary(Qureg qureg, int targetQubit, Complex alpha, Complex beta);

/** Apply a general single-qubit unitary (including a global phase factor).
 * The passed 2x2 ComplexMatrix must be unitary, otherwise an error is thrown.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    \f]
 * 
 * If \p qureg is a state-vector, then the resulting state is \f$ u \, |\text{qureg}\rangle \f$.\n
 * If \p qureg is a density-matrix \f$ \rho \f$, then the resulting state is \f$ u \, \rho \, u^\dagger \f$.
 *
 * > Use applyMatrix2() to left-multiply a non-unitary ::ComplexMatrix2
 *
 * @see
 * - ::ComplexMatrix2
 * - controlledUnitary()
 * - multiControlledUnitary()
 * - multiStateControlledUnitary()
 * - twoQubitUnitary()
 * - multiQubitUnitary()
 * - applyMatrix2()
 * - compactUnitary()
 *
 * @ingroup unitary                                                              
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] u unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if matrix \p u is not unitary
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void unitary(Qureg qureg, int targetQubit, ComplexMatrix2 u);

/** Rotate a single qubit by a given angle around the X-axis of the Bloch-sphere. For angle \f$\theta\f$, applies
 * \f[
 * \begin{pmatrix}
 * \cos\theta/2 & -i \sin \theta/2\\
 * -i \sin \theta/2 & \cos \theta/2
 * \end{pmatrix}
 * \f]
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {rot};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_x(\theta)$};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - controlledRotateX()
 * - rotateY()
 * - rotateZ()
 * - rotateAroundAxis()
 * - multiRotateZ()
 * - multiRotatePauli()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws invalidQuESTInputError()
 * - if \p rotQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void rotateX(Qureg qureg, int rotQubit, qreal angle);

/** Rotate a single qubit by a given angle around the Y-axis of the Bloch-sphere. 
 * For angle \f$\theta\f$, applies
 * \f[
 * \begin{pmatrix}
 * \cos\theta/2 & - \sin \theta/2\\
 * \sin \theta/2 & \cos \theta/2
 * \end{pmatrix}
 * \f]            
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {rot};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_y(\theta)$};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - controlledRotateY()
 * - rotateX()
 * - rotateZ()
 * - rotateAroundAxis()
 * - multiRotateZ()
 * - multiRotatePauli()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws invalidQuESTInputError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc, debug)
 */
void rotateY(Qureg qureg, int rotQubit, qreal angle);

/** Rotate a single qubit by a given angle around the Z-axis of the Bloch-sphere (also known as a phase shift gate).   
 * For angle \f$\theta\f$, applies
 * \f[
 * \begin{pmatrix}
 * \exp(-i \theta/2) & 0 \\
 * 0 & \exp(i \theta/2)
 * \end{pmatrix}
 * \f] 
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {rot};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_z(\theta)$};
                \end{tikzpicture}
    \f]
 * 
 * @see
 * - multiRotateZ()
 * - controlledRotateZ()
 * - rotateY()
 * - rotateX()
 * - rotateAroundAxis()
 * - multiRotatePauli()
 * - phaseShift()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws invalidQuESTInputError()
 * - if \p rotQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void rotateZ(Qureg qureg, int rotQubit, qreal angle);

/** Rotate a single qubit by a given angle around a given \ref Vector on the Bloch-sphere.      
 * The vector must not be zero (else an error is thrown), but needn't be unit magnitude, since 
 * it will be normalised by QuEST.
 *
 * For angle \f$\theta\f$ and axis vector \f$\vec{n}\f$, applies \f$R_{\hat{n}} = \exp \left(- i \frac{\theta}{2} \hat{n} \cdot \vec{\sigma} \right) \f$
 * where \f$\vec{\sigma}\f$ is the vector of Pauli matrices.
 *
 * @see
 * - controlledRotateAroundAxis()
 * - rotateX()
 * - rotateY()
 * - rotateZ()
 * - multiRotateZ()
 * - multiRotatePauli()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @param[in] axis vector around which to rotate (can be non-unit; will be normalised)
 * @throws invalidQuESTInputError()
 * - if \p rotQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if \p axis is the zero vector
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void rotateAroundAxis(Qureg qureg, int rotQubit, qreal angle, Vector axis);

/** Applies a controlled rotation by a given angle around the X-axis of the Bloch-sphere. 
 * The target qubit is rotated in states where the control qubit has value 1.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_x(\theta)$};
                \end{tikzpicture}
    \f] 
 *
 * @see
 * - rotateX()
 * - controlledRotateY()
 * - controlledRotateZ()
 * - controlledRotateAroundAxis()
 * - controlledPhaseShift()
 * - multiRotateZ()
 * - multiRotatePauli()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubit and \p targetQubit are equal
 * @author Tyson Jones
 */
void controlledRotateX(Qureg qureg, int controlQubit, int targetQubit, qreal angle);

/** Applies a controlled rotation by a given angle around the Y-axis of the Bloch-sphere. 
 * The target qubit is rotated in states where the control qubit has value 1.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_y(\theta)$};
                \end{tikzpicture}
    \f] 
 *
 * - rotateY()
 * - controlledRotateX()
 * - controlledRotateZ()
 * - controlledRotateAroundAxis()
 * - controlledPhaseShift()
 * - multiRotateZ()
 * - multiRotatePauli()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubit and \p targetQubit are equal
 * @author Tyson Jones
 */
void controlledRotateY(Qureg qureg, int controlQubit, int targetQubit, qreal angle);

/** Applies a controlled rotation by a given angle around the Z-axis of the Bloch-sphere. 
 * The target qubit is rotated in states where the control qubit has value 1.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_z(\theta)$};
                \end{tikzpicture}
    \f] 
 *
 * @see
 * - rotateZ()
 * - controlledRotateX()
 * - controlledRotateY()
 * - controlledRotateAroundAxis()
 * - controlledPhaseShift()
 * - multiRotateZ()
 * - multiRotatePauli()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubit and \p targetQubit are equal
 * @author Tyson Jones
 */
void controlledRotateZ(Qureg qureg, int controlQubit, int targetQubit, qreal angle);

/** Applies a controlled rotation by a given angle around a given vector on the Bloch-sphere.      
 * The vector must not be zero (else an error is thrown), but needn't be unit magnitude.
 *
 * For angle \f$\theta\f$ and axis vector \f$\vec{n}\f$, applies \f$R_{\hat{n}} = \exp \left(- i \frac{\theta}{2} \hat{n} \cdot \vec{\sigma} \right) \f$ to states where the target qubit is 1 
 * (\f$\vec{\sigma}\f$ is the vector of Pauli matrices).
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_{\hat{n}}(\theta)$};
                \end{tikzpicture}
    \f]
 *
 * @see 
 * - rotateAroundAxis()
 * - multiRotatePauli()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit with value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @param[in] axis vector around which to rotate (can be non-unit; will be normalised)
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubit and \p targetQubit are equal
 * - if \p axis is the zero vector
 * @author Tyson Jones
 */
void controlledRotateAroundAxis(Qureg qureg, int controlQubit, int targetQubit, qreal angle, Vector axis);

/** Apply a controlled unitary (single control, single target) parameterised by two given complex scalars.
 * Given valid complex numbers \f$\alpha\f$ and \f$\beta\f$, applies the two-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\
 * & & \alpha & -\beta^* \\
 * & & \beta & \alpha^*
 * \end{pmatrix}
 * \f]
 * to the control and target qubits.
 * Valid \f$\alpha\f$, \f$\beta\f$ satisfy \f$|\alpha|^2 + |\beta|^2 = 1\f$. 
 * The target unitary is general up to a global phase factor.         
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$U_{\alpha, \beta}$};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - compactUnitary()
 * - controlledUnitary()
 * - multiControlledUnitary()
 * - multiStateControlledUnitary()
 *
 * @ingroup unitary                                                             
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit apply the target unitary if this qubit has value 1
 * @param[in] targetQubit qubit on which to apply the target unitary
 * @param[in] alpha complex unitary parameter (row 1, column 1)
 * @param[in] beta complex unitary parameter (row 2, column 1)
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubit and \p targetQubit are equal
 * - if \p alpha, \p beta don't satisfy |`alpha`|^2 + |`beta`|^2 = 1
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void controlledCompactUnitary(Qureg qureg, int controlQubit, int targetQubit, Complex alpha, Complex beta);

/** Apply a general controlled unitary (single control, single target), which can include a global phase factor.
 * The given unitary is applied to the target qubit if the control qubit has value 1,
 * effecting the two-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\
 * & & u_{00} & u_{01}\\
 * & & u_{10} & u_{11}
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *      
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - ::ComplexMatrix2
 * - unitary()
 * - multiControlledUnitary()
 * - multiStateControlledUnitary()
 * - twoQubitUnitary()
 * - multiQubitUnitary()
 *
 * @ingroup unitary                                                           
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit apply unitary if this qubit is 1
 * @param[in] targetQubit qubit to operate on
 * @param[in] u single-qubit unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubit and \p targetQubit are equal
 * - if \p u is not unitary
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void controlledUnitary(Qureg qureg, int controlQubit, int targetQubit, ComplexMatrix2 u);

/** Apply a general multiple-control single-target unitary, which can include
 * a global phase factor. Any number of control qubits can be specified,
 * and if all have value 1, the given unitary is applied to the target qubit.
 * This effects the many-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & \ddots \\
 * & & & u_{00} & u_{01}\\
 * & & & u_{10} & u_{11}
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 * The given 2x2 ComplexMatrix must be unitary, otherwise an error is thrown.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 3) {controls};
                \node[draw=none] at (-3.5, 0) {target};

                \node[draw=none] at (0, 6) {$\vdots$};
                \draw (0, 5) -- (0, 4);
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw (0, 4) -- (0, 2);         
                
                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    \f]
 *
 * @see 
 * - ::ComplexMatrix2
 * - unitary()
 * - controlledUnitary()
 * - multiStateControlledUnitary()
 * - twoQubitUnitary()
 * - multiQubitUnitary()
 *
 * @ingroup unitary                                                            
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits applies unitary if all qubits in this array equal 1
 * @param[in] numControlQubits number of control qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] u single-qubit unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented])
 * - if any qubit index (\p targetQubit or one in \p controlQubits) is outside [0, \p qureg.numQubitsRepresented])
 * - if any qubit in \p controlQubits is repeated
 * - if \p controlQubits contains \p targetQubit
 * - if \p u is not unitary
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void multiControlledUnitary(Qureg qureg, int* controlQubits, int numControlQubits, int targetQubit, ComplexMatrix2 u);

/** Apply the single-qubit Pauli-X (also known as the X, sigma-X, NOT or bit-flip) gate.
 * This is a rotation of \f$\pi\f$ around the x-axis on the Bloch sphere. I.e. 
 * \f[
 * \begin{pmatrix}
 * 0 & 1 \\
 * 1 & 0
 * \end{pmatrix}
 * \f]   
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (2, 0);
                \draw (0, 0) circle (.5);
                \draw (0, .5) -- (0, -.5);
                \end{tikzpicture}
    \f]   
 *
 * @see
 * - rotateX()
 * - pauliY()
 * - pauliZ()
 * - controlledNot()
 * - multiQubitNot()
 * - multiControlledMultiQubitNot()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void pauliX(Qureg qureg, int targetQubit);

/** Apply the single-qubit Pauli-Y (also known as the Y or sigma-Y) gate.
 * This is a rotation of \f$\pi\f$ around the Y-axis on the Bloch sphere. I.e. 
 * \f[
 * \begin{pmatrix}
 * 0 & -i \\
 * i & 0
 * \end{pmatrix}
 * \f]  
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$\sigma_y$};
                \end{tikzpicture}
    \f]      
 *
 * @see
 * - rotateY()
 * - pauliX()
 * - pauliZ()
 * 
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void pauliY(Qureg qureg, int targetQubit);

/** Apply the single-qubit Pauli-Z (also known as the Z, sigma-Z or phase-flip) gate.
 * This is a rotation of \f$\pi\f$ around the Z-axis (a phase shift) on the Bloch sphere. I.e. 
 * \f[
 * \begin{pmatrix}
 * 1 & 0 \\
 * 0 & -1
 * \end{pmatrix}
 * \f]   
 * with circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$\sigma_z$};
                \end{tikzpicture}
    \f]     
 *
 * @see
 * - phaseShift()
 * - rotateZ()
 * - pauliX()
 * - pauliY()
 * - controlledPhaseFlip()
 * - multiControlledPhaseFlip()
 * 
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void pauliZ(Qureg qureg, int targetQubit);

/** Apply the single-qubit Hadamard gate.
 * This takes \f$|0\rangle\f$ to \f$|+\rangle\f$ and \f$|1\rangle\f$ to \f$|-\rangle\f$, and is equivalent to a rotation of
 * \f$\pi\f$ around the x-axis then \f$\pi/2\f$ about the y-axis on the Bloch-sphere. I.e.
 * \f[
 * \frac{1}{\sqrt{2}}
 * \begin{pmatrix}
 * 1 & 1 \\
 * 1 & -1
 * \end{pmatrix}
 * \f]  
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {H};
                \end{tikzpicture}
    \f]  
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void hadamard(Qureg qureg, int targetQubit);

/** Apply the controlled not (single control, single target) gate, also
 * known as the c-X, c-sigma-X, c-Pauli-X and c-bit-flip gate.
 * This applies pauliX to the target qubit if the control qubit has value 1.
 * This effects the two-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & & 1 \\
 * & & 1
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, -.5);
                
                \draw (-2,0) -- (2, 0);
                \draw (0, 0) circle (.5);
                \end{tikzpicture}
    \f]  
 *
 * @see
 * - multiControlledMultiQubitNot()
 * - pauliX()
 *
 * @ingroup unitary
 * @param[in,out] qureg the state-vector or density matrix to modify
 * @param[in] controlQubit nots the target if this qubit is 1
 * @param[in] targetQubit qubit to not
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * - if \p controlQubit and \p targetQubit are equal
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void controlledNot(Qureg qureg, int controlQubit, int targetQubit);

/** Apply a NOT (or Pauli X) gate with multiple control and target qubits. 
 * This applies pauliX to qubits \p targs on every basis state for which the 
 * control qubits \p ctrls are all in the \f$|1\rangle\f$ state. The ordering within 
 * each of \p ctrls and \p targs has no effect on the operation.
 * > This function is equivalent, but significantly faster (approximately \p numTargs times)
 * > than applying controlled NOTs on each qubit in \p targs in turn, since:
 * > \f[
 * >     C_{a, \,b, \,\dots}( X_c \otimes X_d \otimes \dots ) \equiv
 * >     C_{a, \,b, \,\dots}( X_c) \; \otimes \; C_{a, \,b, \,\dots}(X_d) \; \otimes \; \dots
 * > \f]
 *
 * The effected unitary, if \p targs and \p ctrls happened to be contiguous, has matrix:
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & \ddots \\
 * & & & &   &    & {{\scriptstyle\cdot}^{{\scriptstyle\cdot}^{{\scriptstyle\cdot}}}} \\
 * & & & &   & 1  &   \\
 * & & & & 1 &    &  \\
 * & & & {{\scriptstyle\cdot}^{{\scriptstyle\cdot}^{{\scriptstyle\cdot}}}} & & &
 * \end{pmatrix}
 * \f]
 * and circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 1) {targets};
                \node[draw=none] at (-3.5, 5) {controls};
                
                \node[draw=none] at (0, 8) {$\vdots$};
                \draw (0, 7) -- (0, 6);
                
                \draw (-2, 6) -- (2, 6);
                \draw[fill=black] (0, 6) circle (.2);
                \draw (0, 6) -- (0, 4);         
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw(0, 4) -- (0, -1);
                
                \draw (-2,2) -- (2, 2);
                \draw (0, 2) circle (.4);

                \draw (-2,0) -- (2, 0);
                \draw (0, 0) circle (.4);
                
                \node[draw=none] at (0, -1.5) {$\vdots$};
                \end{tikzpicture}
    \f]
 * > In distributed mode, this operation requires at most a single round of pair-wise 
 * > communication between nodes, and hence is as efficient as pauliX().
 *
 * @see
 * - multiQubitNot()
 * - controlledNot()
 * - pauliX()
 *
 * @ingroup unitary
 * @param[in,out] qureg a state-vector or density matrix to modify
 * @param[in] ctrls a list of the control qubit indices
 * @param[in] numCtrls the length of list \p ctrls 
 * @param[in] targs a list of the qubits to be targeted by the X gates
 * @param[in] numTargs the length of list \p targs
 * @throws invalidQuESTInputError()
 * - if any qubit in \p ctrls and \p targs is invalid, i.e. outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * - if \p ctrls or \p targs contain any repetitions
 * - if any qubit in \p ctrls is also in \p targs (and vice versa)
 * - if \p numTargs <b>< 1</b>
 * - if \p numCtrls <b>< 1</b> (use multiQubitNot() for no controls)
 * @throws segmentation-fault
 * - if \p ctrls contains fewer elements than \p numCtrls
 * - if \p targs contains fewer elements than \p numTargs
 * @author Tyson Jones
 */
void multiControlledMultiQubitNot(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs);

/** Apply a NOT (or Pauli X) gate with multiple target qubits, which has the same 
 * effect as (but is much faster than) applying each single-qubit NOT gate in turn.
 *
 * The ordering within \p targs has no effect on the operation.
 * > This function is equivalent, but significantly faster (approximately \p numTargs times)
 * > than applying pauliX() on each qubit in \p targs in turn.
 * > \f[
 * >     X_a \otimes X_b \otimes \dots 
 * > \f]
 *
 * The effected unitary, if \p targs happen to be contiguous, has matrix:
 * \f[
 * \begin{pmatrix}
 *  &   &    & {{\scriptstyle\cdot}^{{\scriptstyle\cdot}^{{\scriptstyle\cdot}}}} \\
 *  &   & 1  &   \\
 *  & 1 &    &  \\
 *  {{\scriptstyle\cdot}^{{\scriptstyle\cdot}^{{\scriptstyle\cdot}}}} & & &
 * \end{pmatrix}
 * \f]
 * and circuit diagram:
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 1) {targets};
                \draw (0, -1) -- (0, 2.4);
                
                \draw (-2,2) -- (2, 2);
                \draw (0, 2) circle (.4);

                \draw (-2,0) -- (2, 0);
                \draw (0, 0) circle (.4);
                
                \node[draw=none] at (0, -1.5) {$\vdots$};
                \end{tikzpicture}
    \f]
 * > In distributed mode, this operation requires at most a single round of pair-wise 
 * > communication between nodes, and hence is as efficient as pauliX().
 *
 * @see
 * - multiControlledMultiQubitNot()
 * - controlledNot()
 * - pauliX()
 *
 * @ingroup unitary
 * @param[in,out] qureg a state-vector or density matrix to modify
 * @param[in] targs a list of the qubits to be targeted by the X gates
 * @param[in] numTargs the length of list \p targs
 * @throws invalidQuESTInputError()
 * - if any qubit in \p targs is invalid, i.e. outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * - if \p targs contain any repetitions
 * - if \p numTargs <b>< 1</b> 
 * - if \p numTargs <b>></b>`qureg.numQubitsRepresented`
 * @throws segmentation-fault
 * - if \p targs contains fewer elements than \p numTargs
 * @author Tyson Jones
 */
void multiQubitNot(Qureg qureg, int* targs, int numTargs);

/** Apply the controlled pauliY (single control, single target) gate, also
 * known as the c-Y and c-sigma-Y gate.
 * This applies pauliY to the target qubit if the control qubit has value 1.
 * This effects the two-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & & -i \\
 * & & i
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {Y};
                \end{tikzpicture}
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit applies pauliY to the target if this qubit is 1
 * @param[in] targetQubit qubit to not
 * @throws invalidQuESTInputError()
 * - if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubit and \p targetQubit are equal
 * @author Tyson Jones
 * @author Ania Brown (debug)
 */
void controlledPauliY(Qureg qureg, int controlQubit, int targetQubit);

/** Gives the probability of a specified qubit being measured in the given outcome (0 or 1).
 * This performs no actual measurement and does not change the state of the qubits.
 * 
 * For state-vectors, this function works by summing the absolute-value-squared of every 
 * amplitude in the state-vector for which \p measureQubit \p = \p 0. 
 * If \p outcome \p = \p 1, it returns \p 1 minus this value. Hence for unnormalised 
 * state-vectors, this result will differ from the absolute-value-squared of every 
 * amplitude where \p measureQubit \p = \p outcome.
 *
 * For density matrices, this function sums the diagonal values (should be real)
 * corresponding to \p measureQubit \p = \p 0 (returning 1 minus this if \p outcome \p = \p 1).
 *
 * @see 
 * - calcProbOfAllOutcomes()
 * 
 * @see
 * - measure()
 * - measureWithStats()
 * - collapseToOutcome()
 * - calcTotalProb()
 *
 * @ingroup calc
 * @param[in] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to study
 * @param[in] outcome for which to find the probability of the qubit being measured in
 * @return probability of qubit measureQubit being measured in the given outcome
 * @throws invalidQuESTInputError()
 * - if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if \p outcome is not in {0, 1}
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
qreal calcProbOfOutcome(Qureg qureg, int measureQubit, int outcome);

/** Populates \p outcomeProbs with the probabilities of every outcome of the sub-register
 * contained in \p qubits. 
 *
 * > This performs no actual measurement and does not modify \p qureg.
 *
 * - For \p qubits <b>= {a, b, c, ...}</b>, \p outcomeProbs
 *   is populated to
 *   \f[
 *     \text{outcomeProbs} = \{ \; |\alpha_0|^2, \; |\alpha_1|^2, \; |\alpha_2|^2, \; |\alpha_3|^2, \; ... \; \},
 *   \f]
 *   where \f$|\alpha_j|^2\f$ are the probabilities of the respective outcome states (interpreting
 *   \p qubits as ordered least to most significant)
 *   \f[
 *      |\dots\textbf{c\,b\,a}\rangle_i \; \; = \;\; |000\rangle,  \;\; |001\rangle \;\; |010\rangle \;\; |011\rangle, \;\; \dots
 *   \f]
 *   understood in a state-vector \p qureg \f$|\psi\rangle\f$ as
 *   \f[
 *      |\psi\rangle = \sum\limits_i^{\text{numQubits}} \alpha_i \; |\dots\textbf{c\,b\,a}\rangle_i 
 *        \; \otimes \; |\phi\rangle_i,
 *   \f]
 *   or in a density matrix \p qureg \f$\rho\f$ as
 *   \f[
 *      \begin{aligned}
 *      \rho &= \sum\limits_{i,j}^{\text{numQubits}} \; \beta_{ij} \; |\dots\textbf{c\,b\,a}\rangle_i\,\langle\dots\textbf{c\,b\,a}|_j 
 *            \; \otimes \; \mu_{ij} \\
 *           &= \sum\limits_i^{\text{numQubits}} \; |\alpha_i|^2 \;  |\dots\textbf{c\,b\,a}\rangle\langle\dots\textbf{c\,b\,a}|_i  \;\; + \, 
 *            \sum\limits_{i \ne j}^{\text{numQubits}} \; \beta_{ij} \; |\dots\textbf{c\,b\,a}\rangle_i\,\langle\dots\textbf{c\,b\,a}|_j 
 *            \; \otimes \; \mu_{ij},
 *       \end{aligned}
 *   \f]
 *   where \f$|\phi\rangle_i\f$ and \f$\mu_{ij}\f$ are understood to be the separated states of the remaining qubits.
 *   
 * > \p outcomeProbs must be a pre-allocated array of length \f$2^{\text{numQubits}}\f$,
 * > which is equivalent to <b>1<<</b>\p numQubits. In distributed mode, every node receives 
 * > the full list of outcome probabilities.
 *
 * - Note that the indices in \p qubits need not be adjacent nor ordered. The order of
 *   \p qubits determines the order of \p outcomeProbs, whereby \p qubits are treated
 *   as <em>increasing</em> significance.
 *   \n\n
 *
 * - For example, given a \p qureg initialised into state-vector
 *   \f[
 *      |\psi\rangle = 
 *          \alpha_0 |000\rangle \;+\; \alpha_1 |001\rangle \;+\; 
 *          \alpha_2 |010\rangle \;+\; \alpha_3 |011\rangle \;+\;
 *          \alpha_4 |100\rangle \;+\; \alpha_5 |101\rangle \;+\; 
 *          \alpha_6 |110\rangle \;+\; \alpha_7 |111\rangle,
 *   \f]
 *   then executing
 *   ```
 *   int qubits[] = {2, 0};
 *   int numQubits = 2;
 *
 *   qreal outcomeProbs[1<<numQubits];
 *   calcProbOfAllOutcomes(outcomeProbs, qureg, qubits, numQubits);
 *   ```
 *   would populate \p outcomeProbs with
 *   \f[
 *      \text{outcomeProbs} = \{ \;\; |\alpha_0|^2+|\alpha_2|^2, \;\; |\alpha_4|^2+|\alpha_6|^2, \;\;
 *                                 |\alpha_1|^2+|\alpha_3|^2, \;\; |\alpha_5|^2+|\alpha_7|^2  \;\; \}.
 *   \f]
 *
 * > Since all probability amplitudes of a state-vector are ultimately involved in 
 * > the output probabilities, calcProbOfAllOutcomes() works as expected for 
 * > un-normalised states. This is similarly true for density matrices, where all
 * > diagonal elements are involved, although only the real values of the diagonal elements 
 * > will be consulted.
 *
 * @see 
 * - calcProbOfOutcome()
 *
 * @ingroup calc
 * @param[out] outcomeProbs a pre-allocated array of length <b>1<<</b>\p numQubits, 
 *      which will be modified to contain all outcome probabilities
 * @param[in] qureg a state-vector or density matrix to study
 * @param[in] qubits a list of qubits to study
 * @param[in] numQubits the length of list \p qubits
 * @throws invalidQuESTInputError()
 * - if \p numQubits <= 0
 * - if any index in \p qubits is invalid, i.e. outside <b>[0,</b> \p qureg.numQubitsRepresented <b>)</b>
 * - if \p qubits contains any repetitions
 * @throws segmentation-fault
 * - if \p outcomeProbs is not pre-allocated
 * - if \p outcomeProbs contains space for fewer than <b>1<<</b>\p numQubits elements
 * @author Tyson Jones
 */
void calcProbOfAllOutcomes(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits);

/** Updates \p qureg to be consistent with measuring \p measureQubit in the given 
 * \p outcome (0 or 1), and returns the probability of such a measurement outcome. 
 * This is effectively performing a renormalising projection, or a measurement with a forced outcome.
 * This is an irreversible change to the state, whereby computational states
 * inconsistant with the outcome are given zero amplitude and the \p qureg is renormalised.
 * The given outcome must not have a near zero probability, else it cannot be
 * collapsed into.
 *
 * Note that the collapse probably used for renormalisation is calculated for 
 * \p outcome \p = \p 0, and assumed 1 minus this probability if \p outcome \p = \p 1.
 * Hence this routine will not correctly project un-normalised quregs onto 
 * \p outcome \p = \p 1.
 *
 * To avoid renormalisation after projection, or force projection into non-physical 
 * states with very small probability, use applyProjector().
 * 
 * @see
 * - measure()
 * - measureWithStats()
 * - applyProjector()
 *
 * @ingroup normgate
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @param[in] outcome to force the measure qubit to enter
 * @return probability of the (forced) measurement outcome
 * @throws invalidQuESTInputError()
 * - if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if \p outcome is not in {0, 1}
 * - if the probability of \p outcome is zero (within machine epsilon)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
qreal collapseToOutcome(Qureg qureg, int measureQubit, int outcome);

/** Measures a single qubit, collapsing it randomly to 0 or 1.
 *
 * Outcome probabilities are weighted by the state vector, which is irreversibly
 * changed after collapse to be consistent with the outcome.
 *
 * > The random outcome generator is seeded by seedQuESTDefault() within 
 * > createQuESTEnv(), unless later overridden by seedQuEST().
 * 
 * @see
 * - measureWithStats()
 * - collapseToOutcome()
 * - seedQuEST()
 * - seedQuESTDefault()
 * 
 * @ingroup normgate
 * @param[in, out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @return the measurement outcome, 0 or 1
 * @throws invalidQuESTInputError()
 * - if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
int measure(Qureg qureg, int measureQubit);

/** Measures a single qubit, collapsing it randomly to 0 or 1, and
 * additionally gives the probability of that outcome.
 * Outcome probabilities are weighted by the state vector, which is irreversibly
 * changed after collapse to be consistent with the outcome.
 *
 * > The random outcome generator is seeded by seedQuESTDefault() within 
 * > createQuESTEnv(), unless later overridden by seedQuEST().
 *
 * @see 
 * - measure()
 * - collapseToOutcome()
 * - seedQuESTDefault()
 * - seedQuEST()
 *
 * @ingroup normgate
 * @param[in, out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @param[out] outcomeProb a pointer to a qreal which is set to the probability of the occurred outcome
 * @return the measurement outcome, 0 or 1
 * @throws invalidQuESTInputError()
 * - if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
int measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb);

/** Computes the inner product \f$ \langle \text{bra} | \text{ket} \rangle \f$ of two 
 * equal-size state vectors, given by 
 * \f[ 
     \langle \text{bra} | \text{ket} \rangle = \sum_i {\text{bra}_i}^* \; \times \; \text{ket}_i 
 * \f]
 * The same \p qureg may be passed as both \p bra and \p ket, 
 * though we recommend users check state-vector normalisation with \p calcTotalProb which 
 * employs Kahan summation for greater accuracy.
 * Neither state-vector is modified.
 *
 * This function returns the correct inner product even if \p bra and \p ket are 
 * not correctly normalised states.
 *
 * @see
 * - calcDensityInnerProduct()
 *
 * @ingroup calc
 * @param[in] bra qureg to be the 'bra' (i.e. have its values conjugate transposed) in the inner product 
 * @param[in] ket qureg to be the 'ket' in the inner product 
 * @return the complex inner product of \p bra and \p ket 
 * @throws invalidQuESTInputError()
 * - if either \p bra and \p ket are not both state-vectors
 * - if \p bra and \p ket do not have equal dimensions
 * @author Tyson Jones
 */
Complex calcInnerProduct(Qureg bra, Qureg ket);

/** Computes the Hilbert-Schmidt scalar product
 * (which is equivalent to the Frobenius inner product of matrices) 
 * of two density matrices \p rho1 and \p rho2 of equivalent size.
 * That is, we define the Hilbert-Schmidt scalar product
 * \f[
    ((\rho_1, \rho_2))_{HS} := \text{Tr}[ \rho_1^\dagger \rho_2 ],
 * \f]
 * which is equivalent to the sum of products of matrix elemets, i.e.,
 * \f[
    ((\rho_1, \rho_2))_{HS} = \sum\limits_i \sum\limits_j  (\rho_1)_{ij}^* (\rho_2)_{ij}
 * \f]
 * Assuming that both density matrices are Hermitian,
 * the resulting scalar product is real and invariant under
 * reordering its arguments as 
 * \f[
    ((\rho_1, \rho_2))_{HS} = ((\rho_2, \rho_1))_{HS} = \text{Tr}[\rho_1 \rho_2]
 * \f]
 * If both \p rho1 and \p rho2 are density matrices of pure states
 * \p bra and \p ket, then the equality holds
 * \f[
    ((\rho_1, \rho_2))_{HS} = |\langle \text{bra} | \text{ket} \rangle|^2.
 * \f]
 * If either or both of \p rho1 and \p rho2 are non Hermitian (i.e. invalid density 
 * matrices), then this function returns the real component of the scalar product, 
 * and discards the imaginary component. That is, it returns 
 * \f[
     \text{Re}\{ \text{Tr}[ \rho_1^\dagger \rho_2 ] \} = \text{Re}\{ \text{Tr}[ \rho_2^\dagger \rho_1 ] \}.
 * \f]
 * This is still sometimes useful, e.g. in calculating the inner product with an 
 * anti-commutator, e.g. (for Hermitian \f$ \sigma \f$, \f$ \rho \f$, \f$ H \f$)
 * \f[
 *      ((\sigma, H \rho + \rho H))_{HS} = 2 \; \text{Re} \{ ((\sigma, H \rho))_{HS} \} 
 * \f]
 * where \f$ H \rho \f$ could be a weighted sum of Pauli products applied to \f$ \rho \f$ 
 * through applyPauliSum().
 * 
 * @see
 * - calcInnerProduct()
 * - calcHilbertSchmidtDistance()
 *
 * @ingroup calc
 * @param[in] rho1 qureg as a density matrix (to have its values conjugate transposed)
 * @param[in] rho2 qureg as a density matrix
 * @returns the real Hilbert-Schmidt scalar product of density matrices
            \p rho1 and \p rho2 (assuming Hermiticity)
 * @throws invalidQuESTInputError()
 * - if \p rho1 and \p rho2 are not both density matrices
 * - if \p rho1 and \p rho2 have mismatching dimensions
 * @author Balint Koczor (CPU)
 * @author Tyson Jones (GPU)
 */
qreal calcDensityInnerProduct(Qureg rho1, Qureg rho2);

/** Seeds the random number generator with the (master node) current time and process ID.
 * 
 * This is the default seeding used by createQuESTEnv(), and determines the 
 * outcomes in functions like measure() and measureWithStats().
 * 
 * In distributed mode, every node agrees on the seed (nominated by the master node)
 * such that every node generates the same sequence of pseudorandom numbers.
 *
 * > QuEST uses the 
 * > <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html">Mersenne Twister</a>
 * > for random number generation. 
 *
 * @see
 * - Use seedQuEST() to provide a custom seed, overriding the default.
 * - Use getQuESTSeeds() to obtain the seeds currently being used for RNG.
 *
 * @ingroup debug
 * @param[in] env a pointer to the ::QuESTEnv runtime environment
 * @author Ania Brown
 * @author Balint Koczor (Windows compatibility)
 * @author Tyson Jones (doc)
 **/
void seedQuESTDefault(QuESTEnv *env);

/** Seeds the random number generator with a custom array of key(s), overriding the
 * default keys.
 *
 * This determines the sequence of outcomes in functions like measure() and measureWithStats().
 *
 * In distributed mode, the key(s) passed to the master node will be broadcast to all 
 * other nodes, such that every node generates the same sequence of pseudorandom numbers.
 *
 * This function will copy the contents of \p seedArray into a permanent array 
 * `env.seeds`, so \p seedArray is afterward safe to free.
 *
 * > QuEST uses the 
 * > <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html">Mersenne Twister</a>
 * > for random number generation. 
 *
 * @see
 * - Use seedQuESTDefault() to seed via the current timestamp and process id.
 * - Use getQuESTSeeds() to obtain the seeds currently being used for RNG.
 *
 * @ingroup debug
 * @param[in] env a pointer to the ::QuESTEnv runtime environment
 * @param[in] seedArray Array of integers to use as seed. 
 *  This allows the MT to be initialised with more than a 32-bit integer if required
 * @param[in] numSeeds Length of seedArray
 * @author Ania Brown
 * @author Tyson Jones (doc)
 **/
void seedQuEST(QuESTEnv *env, unsigned long int *seedArray, int numSeeds);

/** Obtain the seeds presently used in random number generation.
 *
 * This function sets argument \p seeds to the address of the array of keys
 * which have seeded QuEST's
 * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html">Mersenne Twister</a>
 * random number generator. \p numSeeds is set to the length of \p seeds.
 * These are the seeds which inform the outcomes of random functions like 
 * measure(), and are set using seedQuEST() and seedQuESTDefault().
 *
 * > The output \p seeds array <b>must not</b> be freed, and should not be modified.
 *
 * Obtaining QuEST's seeds is useful for seeding your own random number generators,
 * so that a simulation (with random QuEST measurements, and your own random decisions)
 * can be precisely repeated later, just by calling seedQuEST().
 *
 * Note this function merely sets the arguments to the attributes for \p env. 
 * I.e.
 * ```
 *     unsigned long int* seeds;
 *     int numSeeds;
 *     getQuESTSeeds(env, &seeds, &numSeeds);
 *     
 *     func(seeds, numSeeds);
 * ```
 * is equivalent to
 * ```
 *     func(env.seeds, env.numSeeds);
 * ```
 * However, one should not rely upon their local pointer from getQuESTSeeds() to be 
 * automatically updated after a subsequent call to seedQuEST() or seedQuESTDefault().
 * Instead, getQuESTSeeds() should be recalled.
 *
 * @see
 * - seedQuEST()
 * - seedQuESTDefault()
 *
 * @ingroup debug
 * @param[in] env the ::QuESTEnv runtime environment
 * @param[in] seeds a pointer to an unitialised array to be modified
 * @param[in] numSeeds a pointer to an integer to be modified
 * @author Tyson Jones
 **/
void getQuESTSeeds(QuESTEnv env, unsigned long int** seeds, int* numSeeds);

/** Enable QASM recording. Gates applied to qureg will here-after be added to a
 * growing log of QASM instructions, progressively consuming more memory until 
 * disabled with stopRecordingQASM(). The QASM log is bound to this qureg instance.
 *
 * @ingroup qasm
 * @param[in,out] qureg The qureg to begin recording subsequent operations upon
 * @author Tyson Jones
 */
void startRecordingQASM(Qureg qureg);

/** Disable QASM recording. The recorded QASM will be maintained in qureg
 * and continue to be appended to if startRecordingQASM is recalled.
 *
 * Has no effect if \p qureg was not already recording operations.
 *
 * @ingroup qasm
 * @param[in,out] qureg The qureg to halt recording subsequent operations upon
 * @author Tyson Jones
 */
void stopRecordingQASM(Qureg qureg);

/** Clear all QASM so far recorded. This does not start or stop recording.
 *
 * @ingroup qasm
 * @param[in,out] qureg The qureg of which to clear the QASM log
 * @author Tyson Jones
 */
void clearRecordedQASM(Qureg qureg);

/** Print recorded QASM to stdout. This does not clear the QASM log, nor does it 
 * start or stop QASM recording.
 * 
 * @ingroup qasm
 * @param[in] qureg Prints the QASM recorded for this qureg.
 * @author Tyson Jones
 */
void printRecordedQASM(Qureg qureg);

/** Writes recorded QASM to a file, throwing an error if inaccessible.
 *
 * @ingroup qasm
 * @param[in] qureg Writes the QASM recorded for this qureg to file
 * @param[in] filename The filename of the file to contain the recorded QASM
 * @throws invalidQuESTInputError()
 * - if \p filename cannot be written to
 * @author Tyson Jones
 */
void writeRecordedQASMToFile(Qureg qureg, char* filename);

/** Mixes a density matrix \p qureg to induce single-qubit dephasing noise.
 * With probability \p prob, applies Pauli Z to \p targetQubit.
 *
 * This transforms \p qureg = \f$\rho\f$ into the mixed state
 * \f[
 * (1 - \text{prob}) \, \rho + \text{prob} \; Z_q \, \rho \, Z_q
 * \f]
 * where q = \p targetQubit.
 * \p prob cannot exceed 1/2, which maximally mixes \p targetQubit.
 *
 * @see
 * - mixTwoQubitDephasing()
 * - mixDamping()
 * - mixDepolarising()
 * - mixKrausMap()
 * - mixPauli()
 * - mixDensityMatrix()
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce dephasing noise
 * @param[in] prob the probability of the phase error occuring
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if \p prob is not in [0, 1/2]
 * @author Tyson Jones (GPU, doc)
 * @author Ania Brown (CPU, distributed)
 */
void mixDephasing(Qureg qureg, int targetQubit, qreal prob);

/** Mixes a density matrix \p qureg to induce two-qubit dephasing noise.
 * With probability \p prob, applies Pauli Z to either or both qubits.
 *
 * This transforms \p qureg = \f$\rho\f$ into the mixed state
 * \f[
 * (1 - \text{prob}) \, \rho + \frac{\text{prob}}{3} \; \left( 
 *      Z_a \, \rho \, Z_a + 
 *      Z_b \, \rho \, Z_b + 
 *      Z_a Z_b \, \rho \, Z_a Z_b
 * \right)
 * \f]
 * where a = \p qubit1, b = \p qubit2.
 * \p prob cannot exceed 3/4, at which maximal mixing occurs.
 *
 * @see
 * - mixDephasing()
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] qubit1 qubit upon which to induce dephasing noise
 * @param[in] qubit2 qubit upon which to induce dephasing noise
 * @param[in] prob the probability of the phase error occuring
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if either \p qubit1 or \p qubit2 is outside [0, \p qureg.numQubitsRepresented)
 * - if \p qubit1 = \p qubit2
 * - if \p prob is not in [0, 3/4]
 * @author Tyson Jones (GPU, doc)
 * @author Ania Brown (CPU, distributed)
 */
void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob);

/** Mixes a density matrix \p qureg to induce single-qubit homogeneous depolarising noise.
 * This is equivalent to, with probability \p prob, uniformly randomly applying 
 * either Pauli X, Y, or Z to \p targetQubit.
 *
 * This transforms \p qureg = \f$\rho\f$ into the mixed state
 * \f[
 * (1 - \text{prob}) \, \rho + \frac{\text{prob}}{3} \; \left( 
 *      X_q \, \rho \, X_q + 
 *      Y_q \, \rho \, Y_q + 
 *      Z_q \, \rho \, Z_q
 * \right)
 * \f]
 * where q = \p targetQubit.
 * \p prob cannot exceed 3/4, at which maximal mixing occurs.
 * The produced state is equivalently expressed as
 * \f[
 *      \left( 1 - \frac{4}{3} \text{prob} \right) \rho +
 *      \left( \frac{4}{3} \text{prob} \right) \frac{\vec{\bf{1}}}{2}
 * \f]
 * where \f$ \frac{\vec{\bf{1}}}{2} \f$ is the maximally mixed state of the target 
 * qubit.
 *
 * @see
 * - mixTwoQubitDepolarising()
 * - mixDephasing()
 * - mixDamping()
 * - mixKrausMap()
 * - mixPauli()
 * - mixDensityMatrix()
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if \p prob is not in [0, 3/4]
 * @author Tyson Jones (GPU, doc)
 * @author Ania Brown (CPU, distributed)
 */
void mixDepolarising(Qureg qureg, int targetQubit, qreal prob);

/** Mixes a density matrix \p qureg to induce single-qubit amplitude damping (decay to 0 state).
 * With probability \p prob, applies damping (transition from 1 to 0 state).
 *
 * This transforms \p qureg = \f$\rho\f$ into the mixed state
 * \f[
    K_0 \rho K_0^\dagger + K_1 \rho K_1^\dagger
 * \f]
 * where q = \p targetQubit and \f$K_0\f$ and \f$K_1\f$ are Kraus operators
 * \f[
        K_0 = \begin{pmatrix} 1 & 0 \\ 0 & \sqrt{1-\text{prob}} \end{pmatrix}, \;\;
        K_1 = \begin{pmatrix} 0 & \sqrt{\text{prob}} \\ 0 & 0 \end{pmatrix}.
 * \f]
 * \p prob cannot exceed 1, at which total damping/decay occurs. Note that unlike
 * mixDephasing() and mixDepolarising(), this function can increase the purity of a 
 * mixed state (by, as \p prob becomes 1, gaining certainty that the qubit is in
 * the 0 state).
 *
 * @see
 * - mixDephasing()
 * - mixDepolarising()
 * - mixKrausMap()
 * - mixPauli()
 * - mixDensityMatrix()
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce amplitude damping
 * @param[in] prob the probability of the damping
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if \p prob is not in [0, 1]
 * @author Nicolas Vogt of HQS (local CPU)
 * @author Ania Brown (GPU, patched local CPU)
 * @author Tyson Jones (distributed, doc)
 */
void mixDamping(Qureg qureg, int targetQubit, qreal prob);

/** Mixes a density matrix \p qureg to induce two-qubit homogeneous depolarising noise.
 * With probability \p prob, applies to \p qubit1 and \p qubit2 any operator of the set
 * \f$\{ IX, IY, IZ, XI, YI, ZI, XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ \}\f$.
 * Note this is the set of all two-qubit Pauli gates excluding \f$II\f$.
 *
 * This transforms \p qureg = \f$\rho\f$ into the mixed state
 * \f[
 * (1 - \text{prob}) \, \rho \; + \; \frac{\text{prob}}{15} \; \left( 
 *      \sum \limits_{\sigma_a \in \{X_a,Y_a,Z_a,I_a\}}
 *      \sum \limits_{\sigma_b \in \{X_b,Y_b,Z_b,I_b\}}
 *      \sigma_a \sigma_b \; \rho \; \sigma_a \sigma_b
 * \right)
 * - \frac{\text{prob}}{15} I_a I_b \; \rho \; I_a I_b
 * \f]
 * or verbosely
 * \f[
 * (1 - \text{prob}) \, \rho + \frac{\text{prob}}{15} \; \left( 
 * \begin{aligned}
 *      &X_a \, \rho \, X_a + 
 *      X_b \, \rho \, X_b + 
 *      Y_a \, \rho \, Y_a + 
 *      Y_b \, \rho \, Y_b + 
 *      Z_a \, \rho \, Z_a + 
 *      Z_b \, \rho \, Z_b 
 *   \\
 *    + &X_a X_b \, \rho \, X_a X_b +
 *      X_a Y_b \, \rho \, X_a Y_b +
 *      X_a Z_b \, \rho \, X_a Z_b +
 *      Y_a X_b \, \rho \, Y_a X_b
 * \\
 *   + &Y_a Y_b \, \rho \, Y_a Y_b +
 *      Y_a Z_b \, \rho \, Y_a Z_b +
 *      Z_a X_b \, \rho \, Z_a X_b + 
 *      Z_a Y_b \, \rho \, Z_a Y_b + 
 *      Z_a Z_b \, \rho \, Z_a Z_b
 * \end{aligned}
 * \right)
 * \f]
 * where a = \p qubit1, b = \p qubit2.
 *
 * \p prob cannot exceed 15/16, at which maximal mixing occurs.
 *
 * The produced state is equivalently expressed as
 * \f[ 
 *      \left( 1 - \frac{16}{15} \text{prob} \right) \rho + \left( \frac{16}{15} \text{prob} \right) \frac{\vec{\bf{1}}}{2}
 * \f]
 * where \f$ \frac{\vec{\bf{1}}}{2} \f$ is the maximally mixed state of the two
 * target qubits.
 *
 * @see
 * - mixDepolarising()
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] qubit1 qubit upon which to induce depolarising noise
 * @param[in] qubit2 qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if either \p qubit1 or \p qubit2 is outside [0, \p qureg.numQubitsRepresented)
 * - if \p qubit1 = \p qubit2
 * - if \p prob is not in [0, 15/16]
 * @author Tyson Jones (GPU, doc)
 * @author Ania Brown (CPU, distributed)
 */
void mixTwoQubitDepolarising(Qureg qureg, int qubit1, int qubit2, qreal prob);

/** Mixes a density matrix \p qureg to induce general single-qubit Pauli noise.
 * With probabilities \p probX, \p probY and \p probZ, applies Pauli X, Y, and Z
 * respectively to \p targetQubit.
 *
 * This transforms \p qureg = \f$\rho\f$ into the mixed state
 * \f[
 * (1 - \text{probX} - \text{probY} - \text{probZ}) \, \rho + \;\;\;
 *      (\text{probX})\; X_q \, \rho \, X_q + \;\;\;
 *      (\text{probY})\; Y_q \, \rho \, Y_q + \;\;\;
 *      (\text{probZ})\; Z_q \, \rho \, Z_q
 * \f]
 * where q = \p targetQubit.
 * Each of \p probX, \p probY and \p probZ cannot exceed the chance of no error: 
 * 1 - \p probX - \p probY - \p probZ
 *
 * This function operates by first converting the given Pauli probabilities into 
 * a single-qubit Kraus map (four 2x2 operators).
 *
 * @see
 * - mixDephasing()
 * - mixDepolarising()
 * - mixDamping()
 * - mixKrausMap()
 * - mixDensityMatrix()
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit to decohere
 * @param[in] probX the probability of inducing an X error
 * @param[in] probY the probability of inducing an Y error
 * @param[in] probZ the probability of inducing an Z error
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * - if any of \p probX, \p probY or \p probZ are not in [0, 1]
 * - if any of p in {\p probX, \p probY or \p probZ} don't satisfy p <= (1 - \p probX - \p probY - \p probZ)
 * @author Balint Koczor
 * @author Tyson Jones (refactored, doc)
 */
void mixPauli(Qureg qureg, int targetQubit, qreal probX, qreal probY, qreal probZ);

/** Modifies combineQureg to become (1-\p prob)\p combineProb + \p prob \p otherQureg.
 * Both registers must be equal-dimension density matrices, and prob must be in [0, 1].
 *
 * @see
 * - mixDephasing()
 * - mixDepolarising()
 * - mixDamping()
 * - mixKrausMap()
 * - mixPauli()
 *
 * @ingroup decoherence
 * @param[in,out] combineQureg a density matrix to be modified
 * @param[in] prob the probability of \p otherQureg in the modified \p combineQureg
 * @param[in] otherQureg a density matrix to be mixed into \p combineQureg
 * @throws invalidQuESTInputError()
 * - if either \p combineQureg or \p otherQureg are not density matrices
 * - if the dimensions of \p combineQureg and \p otherQureg do not match
 * - if \p prob is not in [0, 1]
 * @author Tyson Jones
 */
void mixDensityMatrix(Qureg combineQureg, qreal prob, Qureg otherQureg);

/** Calculates the purity of a density matrix, by the trace of the density matrix squared.
 * Returns \f$\text{Tr}(\rho^2)\f$.
 * For a pure state, this =1.
 * For a mixed state, the purity is less than 1 and is lower bounded by 1/2^n, where
 * n is the number of qubits. The minimum purity is achieved for the maximally mixed state identity/2^n.
 *
 * This function does not accept state-vectors, which clearly have purity 1.
 * 
 * Note this function will give incorrect results for non-Hermitian Quregs (i.e. 
 * invalid density matrices), which will disagree with \f$\text{Tr}(\rho^2)\f$.
 * Instead, this function returns \f$\sum_{ij} |\rho_{ij}|^2 \f$.
 *
 * @see
 * - calcFidelity()
 * - calcHilbertSchmidtDistance()
 * - calcTotalProb()
 *
 * @ingroup calc
 * @param[in] qureg a density matrix of which to measure the purity
 * @return the purity
 * @throws invalidQuESTInputError()
 * - if either \p combineQureg or \p otherQureg are not density matrices
 * - if the dimensions of \p combineQureg and \p otherQureg do not match
 * - if \p prob is not in [0, 1]
 * @author Tyson Jones
 */
qreal calcPurity(Qureg qureg);

/** Calculates the fidelity of \p qureg (a state-vector or density matrix) against 
 * a reference pure state (necessarily a state-vector).
 * If \p qureg is a state-vector, this function computes 
 * \f[ 
    |\langle \text{qureg} | \text{pureState} \rangle|^2
 * \f]
 * If \p qureg is a density matrix, this function computes 
 * \f[ 
    \langle \text{pureState} | \text{qureg} | \text{pureState} \rangle
 * \f]
 * In either case, the returned fidelity lies in [0, 1] (assuming both input 
 * states have valid normalisation). If any of the input \p Quregs are not 
 * normalised, this function will return the real component of the correct 
 * linear algebra calculation.
 *
 * The number of qubits represented in \p qureg and \p pureState must match.
 * 
 * @see
 * - calcHilbertSchmidtDistance()
 * - calcPurity()
 *
 * @ingroup calc
 * @param[in] qureg a density matrix or state vector
 * @param[in] pureState a state vector
 * @return the fidelity between the input registers
 * @throws invalidQuESTInputError()
 * - if the second argument (\p pureState) is not a state-vector
 * - if the number of qubits in \p qureg and \p pureState do not match
 * @author Tyson Jones
 */
qreal calcFidelity(Qureg qureg, Qureg pureState);

/** Performs a SWAP gate between \p qubit1 and \p qubit2.
 * This effects
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & & 1 \\\
 * &  1 \\
 * & & & 1
 * \end{pmatrix}
 * \f]
 * on the designated qubits, though is performed internally by three CNOT gates.
 *
   \f[
               \begin{tikzpicture}[scale=.5]
               \node[draw=none] at (-3.5, 2) {qubit1};
               \node[draw=none] at (-3.5, 0) {qubit2};

               \draw (-2, 2) -- (2, 2);
               \draw (0, 2) -- (0, 0);
               \draw (-2,0) -- (2, 0);

               \draw (-.35,-.35) -- (.35,.35);
               \draw (-.35,.35) -- (.35,-.35);

               \draw (-.35,-.35 + 2) -- (.35,.35 + 2);
               \draw (-.35,.35 + 2) -- (.35,-.35 + 2);

               \end{tikzpicture}
   \f]
 *
 * @see
 * - sqrtSwapGate()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] qubit1 qubit to swap
 * @param[in] qubit2 other qubit to swap
 * @throws invalidQuESTInputError()
 * - if either \p qubit1 or \p qubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if \p qubit1 and \p qubit2 are equal
 * @author Tyson Jones
 */
void swapGate(Qureg qureg, int qubit1, int qubit2);


/** Performs a sqrt SWAP gate between \p qubit1 and \p qubit2.
 * This effects
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & \frac{1}{2}(1+i) & \frac{1}{2}(1-i) \\\
 * & \frac{1}{2}(1-i) & \frac{1}{2}(1+i) \\
 * & & & 1
 * \end{pmatrix}
 * \f]
 * on the designated qubits, though is performed internally by three CNOT gates.
 *
   \f[
               \begin{tikzpicture}[scale=.5]
               \node[draw=none] at (-3.5, 2) {qubit1};
               \node[draw=none] at (-3.5, 0) {qubit2};

               \draw (-2, 2) -- (2, 2);
               \draw (0, 2) -- (0, 0);
               \draw (-2,0) -- (2, 0);

               \draw (-.35,-.35) -- (.35,.35);
               \draw (-.35,.35) -- (.35,-.35);

               \draw (-.35,-.35 + 2) -- (.35,.35 + 2);
               \draw (-.35,.35 + 2) -- (.35,-.35 + 2);
               
               \draw[fill=white] (0, 1) circle (.5);
               \node[draw=none] at (0, 1) {1/2};

               \end{tikzpicture}
   \f]
 *
 * @see
 * - swapGate()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] qb1 qubit to sqrt swap
 * @param[in] qb2 other qubit to sqrt swap
 * @throws invalidQuESTInputError()
 * - if either \p qubit1 or \p qubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if \p qubit1 and \p qubit2 are equal
 * @author Tyson Jones
 */
void sqrtSwapGate(Qureg qureg, int qb1, int qb2);

/** Apply a general single-qubit unitary with multiple control qubits, conditioned 
 * upon a specific bit sequence.
 *
 * Any number of control qubits can be specified, along with their classical
 * state (`0` or `1`) to condition upon. Only amplitudes of computational basis states 
 * for which `controlQubits` have corresponding bit values `controlState` are modified 
 * by `u`.
 * 
 * > This function is equivalent (albeit faster) to applying pauliX() on each of the control qubits
 * > which are conditioned on outcome `0`, calling multiControlledUnitary(), then 
 * > re-appplying pauliX() on the same qubits.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 3) {controls};
                \node[draw=none] at (-3.5, 0) {target};

                \node[draw=none] at (0, 6) {$\vdots$};
                \draw (0, 5) -- (0, 4);
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw (0, 4) -- (0, 2);         
                
                \draw (-2, 2) -- (2, 2);
                \draw[fill=white] (0, 2) circle (.2);
                \draw (0, 2-.2) -- (0, 1);
                
                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    \f]
 *
 * @see
 * - ::ComplexMatrix2
 * - unitary()
 * - controlledUnitary()
 * - multiControlledUnitary()
 * - twoQubitUnitary()
 * - multiQubitUnitary()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits the indices of the control qubits 
 * @param[in] controlState the bit values (0 or 1) of each control qubit, upon which to condition
 * @param[in] numControlQubits number of control qubits
 * @param[in] targetQubit qubit to operate the unitary upon
 * @param[in] u single-qubit unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented])
 * - if any qubit index (\p targetQubit or one in \p controlQubits) is outside [0, \p qureg.numQubitsRepresented]),
 * - if any qubit in \p controlQubits is repeated
 * - if \p controlQubits contains \p targetQubit
 * - if any element of \p controlState is not a bit (`0` or `1`)
 * - if \p u is not unitary
 * @author Tyson Jones
 */
void multiStateControlledUnitary(Qureg qureg, int* controlQubits, int* controlState, int numControlQubits, int targetQubit, ComplexMatrix2 u);

/** Apply a multi-qubit Z rotation, also known as a phase gadget, on a selected number of qubits. 
 * This is the unitary 
 * \f[ 
 *    \exp \left( - i \, \frac{\theta}{2} \; \bigotimes_{j}^{\text{numQubits}} Z_j\right)
 * \f]
 * where the Pauli Z gates operate the qubits listed in `qubits`, and cause 
 * rotations of \f$\theta =\f$ \p angle.
 *
 * > All qubits not appearing in \p qubits are assumed to receive the identity operator.
 *
 * This has the effect of premultiplying every amplitude with 
 * \f$\exp(\pm i \theta/2)\f$ where the sign is determined by the parity of
 * the target qubits for that amplitude.
 *
 * @see
 * - multiControlledMultiRotateZ()
 * - multiRotatePauli()
 * - controlledRotateZ()
 * - rotateZ()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] qubits a list of the indices of the target qubits 
 * @param[in] numQubits number of target qubits
 * @param[in] angle the angle by which the multi-qubit state is rotated around the Z axis
 * @throws invalidQuESTInputError()
 * - if \p numQubits is outside [1, \p qureg.numQubitsRepresented])
 * - if any qubit in \p qubits is outside [0, \p qureg.numQubitsRepresented])
 * - if any qubit in \p qubits is repeated
 * @throws segmentation-fault
 * - if \p qubits contains fewer elements than \p numQubits
 * @author Tyson Jones
 */
void multiRotateZ(Qureg qureg, int* qubits, int numQubits, qreal angle);

/** Apply a multi-qubit multi-Pauli rotation, also known as a Pauli gadget,
 * on a selected number of qubits. 
 * This is the unitary 
 * \f[ 
 *    \exp \left( - i \, \frac{\theta}{2} \; \bigotimes_{j}^{\text{numTargets}} \hat{\sigma}_j\right)
 * \f]
 * where  \f$\theta = \f$\p angle and \f$\hat{\sigma}_j \in \{X, Y, Z\}\f$ is a Pauli operator
 * ::pauliOpType operating upon the corresponding qubit `targetQubits`.
 *
 * For example:
 * ```
 *     multiRotatePauli(qureg, (int[]) {4,5,8,9}, (int[]) {0,1,2,3}, 4, .1)
 * ```
 * effects 
 * \f[
 *    \exp \left( - i \, (0.1/2) \; X_5 \, Y_8 \, Z_9 \right) 
 * \f] 
 * on \p qureg, where unspecified qubits (along with those targeted by `PAULI_I`) are 
 * assumed to receive the identity operator (excluded from exponentiation). 
 * > This means specifying `PAULI_I` does *not* induce a global phase factor \f$\exp(-i \theta/2)\f$.
 * > Hence, if all \p targetPaulis are identity, then this function does nothing to \p qureg.
 * > Specifying `PAULI_I` on a qubit is superfluous but allowed for convenience.
 *
 * This function effects the Pauli gadget by first rotating the qubits which are 
 * nominated to receive `X` or `Y` Paulis into alternate basis, performing 
 * multiRotateZ() on all target qubits, then restoring 
 * the original basis. 
 *
 * @see
 * - multiControlledMultiRotatePauli()
 * - multiRotateZ()
 * - rotateX()
 * - rotateY()
 * - rotateZ()
 * - rotateAroundAxis()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubits a list of the indices of the target qubits 
 * @param[in] targetPaulis a list of the Pauli operators (::pauliOpType)
 *      to apply to the corresponding qubits in \p targetQubits
 * @param[in] numTargets number of target qubits, i.e. the length of \p targetQubits and \p targetPaulis
 * @param[in] angle the angle by which the multi-qubit state is rotated
 * @throws invalidQuESTInputError()
 * - if \p numTargets is outside [1, \p qureg.numQubitsRepresented)
 * - if any qubit in \p targetQubits is outside [0, \p qureg.numQubitsRepresented)
 * - if any qubit in \p targetQubits is repeated
 * - if any element of \p targetPaulis is not one of `PAULI_I`, `PAULI_X`, `PAULI_Y`, `PAULI_Z`
 * @throws segmentation-fault
 * - if \p targetQubits contains fewer elements than \p numTargets
 * - if \p targetPaulis contains fewer elements than \p numTargets
 * @author Tyson Jones
 */
void multiRotatePauli(Qureg qureg, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle);

/** Apply a multi-controlled multi-target Z rotation, also known as a controlled phase gadget.
 * This is the unitary 
 * \f[ 
 *    |1\rangle\langle 1|^{\otimes\, \text{numControls}} \; \otimes \,
 *     \exp \left( - i \, \frac{\theta}{2} \; \bigotimes_{j}^{\text{numTargets}} Z_j\right)
 *     \;\;+\;\; \sum\limits_{k=0}^{2^{\,\text{numControls}} - 2} |k\rangle\langle k| \otimes \text{I}
 * \f]
 * where the Pauli Z gates operate upon the qubits in `targetQubits`, and cause 
 * rotations of \f$\theta =\f$ \p angle.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-4, 1) {targets};
                \node[draw=none] at (-4, 5) {controls};
                
                \node[draw=none] at (0, 8) {$\vdots$};
                \draw (0, 7) -- (0, 6);
                
                \draw (-2.5, 6) -- (2.5, 6);
                \draw[fill=black] (0, 6) circle (.2);
                \draw (0, 6) -- (0, 4);         
                
                \draw (-2.5, 4) -- (2.5, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw(0, 4) -- (0, 3);

                \draw (-2.5,0) -- (-1.5, 0);
                \draw (1.5, 0) -- (2.5, 0);
                \draw (-2.5,2) -- (-1.5, 2);
                \draw (1.5, 2) -- (2.5, 2);
                \draw (-1.5,-1)--(-1.5,3)--(1.5,3)--(1.5,-1);
                \node[draw=none] at (0, 1) {$e^{-i\frac{\theta}{2}Z^{\otimes}}$};
                \node[draw=none] at (0, -1) {$\vdots$};
                
                \end{tikzpicture}
    \f]
 * 
 * > All qubits not appearing in \p targetQubits and \p controlQubits are assumed to receive the identity operator.
 *
 * This has the effect of premultiplying all amplitudes (for which the control qubits are `1`) 
 * with \f$\exp(\pm i \theta/2)\f$, where the sign is determined by the parity of
 * the target qubits for that amplitude.
 *
 * @see
 * - multiControlledMultiRotatePauli()
 * - multiRotatePauli()
 * - multiRotateZ()
 * - controlledRotateZ()
 * - rotateZ()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits list of the indices of qubits to control upon
 * @param[in] numControls length of length `controlQubits`
 * @param[in] targetQubits a list of the indices of the target qubits 
 * @param[in] numTargets length of list `targetQubits`
 * @param[in] angle the angle by which the multi-qubit state is rotated around the Z axis
 * @throws invalidQuESTInputError()
 * - if any qubit in \p controlQubits and \p targetQubits is invalid, i.e. outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * - if \p controlQubits or \p targetQubits contain any repetitions
 * - if any qubit in \p controlQubits is also in \p targetQubits (and vice versa)
 * - if \p numTargets <b>< 1</b>
 * - if \p numControls <b>< 1</b> (use multiRotateZ() for no controls)
 * @throws segmentation-fault
 * - if \p controlQubits contains fewer elements than \p numControls
 * - if \p targetQubits contains fewer elements than \p numTargets
 * @author Tyson Jones
 */
void multiControlledMultiRotateZ(Qureg qureg, int* controlQubits, int numControls, int* targetQubits, int numTargets, qreal angle);

/** Apply a multi-controlled multi-target multi-Pauli rotation, also known as a 
 * controlled Pauli gadget.
 * This is the unitary 
 * \f[ 
 *    |1\rangle\langle 1|^{\otimes\, \text{numControls}} \; \otimes \,
 *     \exp \left( - i \, \frac{\theta}{2} \; \bigotimes_{j}^{\text{numTargets}} \hat{\sigma}_j\right)
 *     \;\;+\;\; \sum\limits_{k=0}^{2^{\,\text{numControls}} - 2} |k\rangle\langle k| \otimes \text{I}
 * \f]
 * where \f$\hat{\sigma}_j\f$ are the Pauli operators (::pauliOpType) in `targetPaulis`, which operate 
 * upon the corresponding qubits in `targetQubits`.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-4, 1) {targets};
                \node[draw=none] at (-4, 5) {controls};
                
                \node[draw=none] at (0, 8) {$\vdots$};
                \draw (0, 7) -- (0, 6);
                
                \draw (-2.5, 6) -- (2.5, 6);
                \draw[fill=black] (0, 6) circle (.2);
                \draw (0, 6) -- (0, 4);         
                
                \draw (-2.5, 4) -- (2.5, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw(0, 4) -- (0, 3);

                \draw (-2.5,0) -- (-1.5, 0);
                \draw (1.5, 0) -- (2.5, 0);
                \draw (-2.5,2) -- (-1.5, 2);
                \draw (1.5, 2) -- (2.5, 2);
                \draw (-1.5,-1)--(-1.5,3)--(1.5,3)--(1.5,-1);
                \node[draw=none] at (0, 1) {$e^{-i\frac{\theta}{2} \bigotimes\limits_j \hat{\sigma}_j }$};
                \node[draw=none] at (0, -1) {$\vdots$};
                
                \end{tikzpicture}
    \f]
 * 
 * > All qubits not appearing in \p targetQubits and \p controlQubits are assumed to receive the identity operator.
 *
 * For example:
 * ```
 *     int numCtrls = 1;
 *     int numTargs = 4;
 *     int ctrls[] = {4};
 *     int targs[] = {0,1,2,3};
 *     
 *     pauliOpType paulis[] = {PAULI_X, PAULI_Y, PAULI_Z, PAULI_I};
 *     
 *     multiControlledMultiRotatePauli(
 *         qureg, ctrls, numCtrls, targs, paulis, numTargs, 0.1);
 * ```
 * effects 
 * \f[
 *    |1\rangle\langle 1 | \otimes \exp\left( -i \, (0.1/2) \, X_0 \, Y_1 \, Z_2 \right) \, \text{I}_3
 *    \;\; + \;\; |0\rangle\langle 0| \otimes \text{I}^{\otimes 4}
 * \f] 
 * on \p qureg, where unspecified qubits (along with those targeted by `PAULI_I`) are 
 * assumed to receive the identity operator (excluded from exponentiation). 
 *
 * > This means specifying `PAULI_I` does *not* induce a global phase factor \f$\exp(-i \theta/2)\f$.
 * > Hence, if all \p targetPaulis are identity, then this function does nothing to \p qureg.
 * > Specifying `PAULI_I` on a qubit is superfluous but allowed for convenience.
 *
 * This function effects the controlled Pauli gadget by first (controlled) 
 * rotating the qubits which are targeted with either `X` or `Y` into alternate basis, 
 * performing multiControlledMultiRotateZ() on all target qubits, then restoring 
 * the original basis. 
 *
 * @see
 * - multiControlledMultiRotateZ()
 * - multiRotatePauli()
 * - multiRotateZ()
 * - rotateX()
 * - rotateY()
 * - rotateZ()
 * - rotateAroundAxis()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits list of the indices of qubits to control upon
 * @param[in] numControls length of length `controlQubits`
 * @param[in] targetQubits a list of the indices of the target qubits 
 * @param[in] targetPaulis a list of the Pauli operators around which to rotate the target qubits
 * @param[in] numTargets length of list `targetQubits`
 * @param[in] angle the angle by which the multi-qubit state is rotated around the Z axis
 * @throws invalidQuESTInputError()
 * - if any qubit in \p controlQubits and \p targetQubits is invalid, i.e. outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * - if \p controlQubits or \p targetQubits contain any repetitions
 * - if any qubit in \p controlQubits is also in \p targetQubits (and vice versa)
 * - if \p numTargets <b>< 1</b>
 * - if \p numControls <b>< 1</b> (use multiRotateZ() for no controls)
 * - if any element of \p targetPaulis is not one of `PAULI_I`, `PAULI_X`, `PAULI_Y`, `PAULI_Z`
 * @throws segmentation-fault
 * - if \p controlQubits contains fewer elements than \p numControls
 * - if \p targetQubits contains fewer elements than \p numTargets
 * - if \p targetPaulis contains fewer elements than \p numTargets
 * @author Tyson Jones
 */
void multiControlledMultiRotatePauli(Qureg qureg, int* controlQubits, int numControls, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle);

/** Computes the expected value of a product of Pauli operators.
 * Letting \f$ \sigma = \otimes_j \hat{\sigma}_j \f$ be the operators indicated by \p pauliCodes 
 * and acting on qubits \p targetQubits, this function computes \f$ \langle \psi | \sigma | \psi \rangle \f$ 
 * if \p qureg = \f$ \psi \f$ is a state-vector, and computes \f$ \text{Trace}(\sigma \rho) \f$ 
 * if \p qureg = \f$ \rho \f$ is a density matrix.
 * 
 * \p pauliCodes is an array of length \p numTargets which specifies which Pauli operators to 
 * enact on the corresponding qubits in \p targetQubits, where 0 = \p PAULI_I, 1 = \p PAULI_X, 
 * 2 = \p PAULI_Y, 3 = \p PAULI_Z. The target qubits must be unique, and at most \p qureg.numQubitsRepresented
 * may be specified. For example, on a 7-qubit state-vector,
 * ```
 *     calcExpecPauliProd(qureg, {4,5,6}, {PAULI_X, PAULI_I, PAULI_Z}, 3, workspace);
 * ```
 * will compute \f$ \langle \psi | I I I I X I Z | \psi \rangle \f$ (where in this notation, the left-most operator
 * applies to the least-significant qubit, i.e. that with index 0).
 *
 * \p workspace must be a register with the same type (state-vector vs density matrix) and dimensions 
 * (number of represented qubits) as \p qureg, and is used as working space. When this function returns, \p qureg 
 * will be unchanged and \p workspace will be set to \f$ \sigma | \psi \rangle \f$ (if \p qureg is a state-vector)
 * or \f$ \sigma \rho \f$ (if \p qureg is a density matrix). NOTE that this last quantity is NOT the result of applying 
 * the paulis as unitaries, \f$ \sigma^\dagger \rho \sigma \f$, but is instead the result of their direct 
 * multiplication with the density matrix. It is therefore itself not a valid density matrix.
 *
 * This function works by cloning the \p qureg state into \p workspace, applying the specified 
 * Pauli operators to \p workspace then computing its inner product with \p qureg (for state-vectors)
 * or its trace (for density matrices). It therefore should scale linearly in time with the number of 
 * specified non-identity Pauli operators, which is bounded by the number of represented qubits.
 *
 * @see
 * - calcExpecDiagonalOp()
 * - calcExpecPauliSum()
 * - calcExpecPauliHamil()
 *
 * @ingroup calc
 * @param[in] qureg the register of which to find the expected value, which is unchanged by this function
 * @param[in] targetQubits a list of the indices of the target qubits 
 * @param[in] pauliCodes a list of the Pauli codes (0=PAULI_I, 1=PAULI_X, 2=PAULI_Y, 3=PAULI_Z) 
 *      to apply to the corresponding qubits in \p targetQubits
 * @param[in] numTargets number of target qubits, i.e. the length of \p targetQubits and \p pauliCodes
 * @param[in,out] workspace a working-space qureg with the same dimensions as \p qureg, which is modified 
 *      to be the result of multiplying the state with the pauli operators
 * @throws invalidQuESTInputError()
 * - if \p numTargets is outside [1, \p qureg.numQubitsRepresented])
 * - if any qubit in \p targetQubits is outside [0, \p qureg.numQubitsRepresented))
 * - if any qubit in \p targetQubits is repeated
 * - if any code in \p pauliCodes is not in {0,1,2,3}
 * - if \p workspace is not of the same type and dimensions as \p qureg
 * @author Tyson Jones
 */
qreal calcExpecPauliProd(Qureg qureg, int* targetQubits, enum pauliOpType* pauliCodes, int numTargets, Qureg workspace);

/** Computes the expected value of a sum of products of Pauli operators.
 * Let \f$ H = \sum_i c_i \otimes_j^{N} \hat{\sigma}_{i,j} \f$ be 
 * the operators indicated by \p allPauliCodes (where \f$ c_i \in \f$ \p termCoeffs 
 * and \f$ N = \f$ \p qureg.numQubitsRepresented).
 * This function computes \f$ \langle \psi | H | \psi \rangle \f$ 
 * if \p qureg = \f$ \psi \f$ is a state-vector, and computes \f$ \text{Trace}(H \rho) =\text{Trace}(\rho H) \f$ 
 * if \p qureg = \f$ \rho \f$ is a density matrix.
 *
 * \p allPauliCodes is an array of length \p numSumTerms*\p qureg.numQubitsRepresented
 * which specifies which Pauli operators to apply, where 0 = \p PAULI_I, 1 = \p PAULI_X, 
 * 2 = \p PAULI_Y, 3 = \p PAULI_Z. For each sum term, a Pauli operator must be specified for 
 * EVERY qubit in \p qureg; each set of \p numSumTerms operators will be grouped into a product.
 * \p termCoeffs is an arrray of length \p numSumTerms containing the term coefficients.
 * For example, on a 3-qubit state-vector,
 * ```
 *     int paulis[6] = {PAULI_X, PAULI_I, PAULI_I,  PAULI_X, PAULI_Y, PAULI_Z};
 *     qreal coeffs[2] = {1.5, -3.6};
 *     calcExpecPauliSum(qureg, paulis, coeffs, 2, workspace);
 * ```
 * will compute \f$ \langle \psi | (1.5 X I I - 3.6 X Y Z) | \psi \rangle \f$ (where in this notation, the left-most operator
 * applies to the least-significant qubit, i.e. that with index 0).
 * 
 * \p workspace must be a register with the same type (state-vector vs density matrix) and dimensions 
 * (number of represented qubits) as \p qureg, and is used as working space. When this function returns, \p qureg 
 * will be unchanged and \p workspace will be set to \p qureg pre-multiplied with the final Pauli product.
 * NOTE that if \p qureg is a density matrix, \p workspace will become \f$ \hat{\sigma} \rho \f$ 
 * which is itself not a density matrix (it is distinct from \f$ \hat{\sigma} \rho \hat{\sigma}^\dagger \f$).
 *
 * This function works by cloning the \p qureg state into \p workspace, applying each of the specified
 * Pauli products to \p workspace (one Pauli operation at a time), then computing its inner product with \p qureg (for state-vectors)
 * or its trace (for density matrices) multiplied with the corresponding coefficient, and summing these contributions. 
 * It therefore should scale linearly in time with the total number of non-identity specified Pauli operators.
 *
 * @see
 * - calcExpecDiagonalOp()
 * - calcExpecPauliProd()
 * - calcExpecPauliHamil()
 *
 * @ingroup calc
 * @param[in] qureg the register of which to find the expected value, which is unchanged by this function
 * @param[in] allPauliCodes a list of the Pauli codes (0=PAULI_I, 1=PAULI_X, 2=PAULI_Y, 3=PAULI_Z) 
 *      of all Paulis involved in the products of terms. A Pauli must be specified for each qubit 
 *      in the register, in every term of the sum.
 * @param[in] termCoeffs The coefficients of each term in the sum of Pauli products
 * @param[in] numSumTerms The total number of Pauli products specified
 * @param[in,out] workspace a working-space qureg with the same dimensions as \p qureg, which is modified 
 *      to be the result of multiplying the state with the final specified Pauli product
 * @throws invalidQuESTInputError()
 * - if any code in \p allPauliCodes is not in {0,1,2,3}
 * - if \p numSumTerms <= 0,
 * - if \p workspace is not of the same type and dimensions as \p qureg
 * @author Tyson Jones
 */
qreal calcExpecPauliSum(Qureg qureg, enum pauliOpType* allPauliCodes, qreal* termCoeffs, int numSumTerms, Qureg workspace);

/** Computes the expected value of \p qureg under Hermitian operator \p hamil.
 * Represent \p hamil as \f$ H = \sum_i c_i \otimes_j^{N} \hat{\sigma}_{i,j} \f$
 *  (where \f$ c_i \in \f$ \p hamil.termCoeffs and \f$ N = \f$ \p hamil.numQubits).
 * This function computes \f$ \langle \psi | H | \psi \rangle \f$ 
 * if \p qureg = \f$ \psi \f$ is a state-vector, and computes \f$ \text{Trace}(H \rho) =\text{Trace}(\rho H) \f$ 
 * if \p qureg = \f$ \rho \f$ is a density matrix.
 *
 * This function is merely an encapsulation of calcExpecPauliSum() - refer to the doc 
 * there for an elaboration.
 * 
 * \p workspace must be a register with the same type (state-vector vs density matrix) and dimensions 
 * (number of represented qubits) as \p qureg and \p hamil, and is used as working space. 
 * When this function returns, \p qureg  will be unchanged and \p workspace will be set to
 * \p qureg pre-multiplied with the final Pauli product in \p hamil.
 * NOTE that if \p qureg is a density matrix, \p workspace will become \f$ \hat{\sigma} \rho \f$ 
 * which is itself not a density matrix (it is distinct from \f$ \hat{\sigma} \rho \hat{\sigma}^\dagger \f$).
 *
 * This function works by cloning the \p qureg state into \p workspace, applying each of the specified
 * Pauli products in \p hamil to \p workspace (one Pauli operation at a time), then computing its inner product with \p qureg (for state-vectors)
 * or its trace (for density matrices) multiplied with the corresponding coefficient, and summing these contributions. 
 * It therefore should scale linearly in time with the total number of non-identity specified Pauli operators.
 *
 * @see 
 * - createPauliHamil()
 * - calcExpecDiagonalOp()
 * - calcExpecPauliSum()
 * - calcExpecPauliProd()
 *
 * @ingroup calc
 * @param[in] qureg the register of which to find the expected value, which is unchanged by this function
 * @param[in] hamil a \p PauliHamil created with createPauliHamil() or createPauliHamilFromFile()
 * @param[in,out] workspace a working-space qureg with the same dimensions as \p qureg, which is modified 
 *      to be the result of multiplying the state with the final specified Pauli product
 * @throws invalidQuESTInputError()
 * - if any code in \p hamil.pauliCodes is not a valid Pauli code
 * - if \p hamil.numSumTerms <= 0
 * - if \p workspace is not of the same type and dimensions as \p qureg and \p hamil
 * @author Tyson Jones
 */
qreal calcExpecPauliHamil(Qureg qureg, PauliHamil hamil, Qureg workspace);

/** Apply a general two-qubit unitary (including a global phase factor).
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target2};
                \node[draw=none] at (-3.5, 2) {target1};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-2,2) -- (-1, 2);
                \draw (1, 2) -- (2, 2);
                \draw (-1,-1)--(-1,3)--(1,3)--(1,-1)--cycle;
                \node[draw=none] at (0, 1) {U};
                \end{tikzpicture}
    \f]
 *
 * \p targetQubit1 is treated as the \p least significant qubit in \p u, such that 
 * a row in \p u is dotted with the vector
 * \f$ |\text{targetQubit2} \;\; \text{targetQubit1}\rangle : \{ |00\rangle, |01\rangle, |10\rangle, |11\rangle \} \f$
 *
 * For example, 
 * ```
 *     twoQubitUnitary(qureg, a, b, u);
 * ```
 * will invoke multiplication
 * \f[
 * \begin{pmatrix}
 * u_{00} & u_{01} & u_{02} & u_{03} \\
 * u_{10} & u_{11} & u_{12} & u_{13} \\
 * u_{20} & u_{21} & u_{22} & u_{23} \\
 * u_{30} & u_{31} & u_{32} & u_{33}
 * \end{pmatrix}
 * \begin{pmatrix}
 * |ba\rangle = |00\rangle \\
 * |ba\rangle = |01\rangle \\
 * |ba\rangle = |10\rangle \\
 * |ba\rangle = |11\rangle 
 * \end{pmatrix}
 * \f]
 *
 * The passed ::ComplexMatrix4 must be unitary, otherwise an error is thrown.
 * > Use applyMatrix4() to left-multiply a non-unitary ::ComplexMatrix4.
 *                 
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q/4 nodes.
 * 
 * @see 
 * - ::ComplexMatrix4
 * - controlledTwoQubitUnitary()
 * - multiControlledTwoQubitUnitary()
 * - multiQubitUnitary()
 * - applyMatrix4()
 *
 * @ingroup unitary                                            
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit1 first qubit to operate on, treated as least significant in \p u
 * @param[in] targetQubit2 second qubit to operate on, treated as most significant in \p u
 * @param[in] u unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p targetQubit1 or \p targetQubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if \p targetQubit1 equals \p targetQubit2
 * - if matrix \p u is not unitary
 * - if each node cannot fit 4 amplitudes in distributed mode
 * @author Tyson Jones
 */
void twoQubitUnitary(Qureg qureg, int targetQubit1, int targetQubit2, ComplexMatrix4 u);

/** Apply a general controlled two-qubit unitary (including a global phase factor).
 * The given unitary is applied to the target amplitudes where the control qubit has value 1.
 * This effects the many-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\
 * & & 1 \\
 * & & & 1 \\
 * & & & & u_{00} & u_{01} & u_{02} & u_{03} \\
 * & & & & u_{10} & u_{11} & u_{12} & u_{13} \\
 * & & & & u_{20} & u_{21} & u_{22} & u_{23} \\
 * & & & & u_{30} & u_{31} & u_{32} & u_{33}
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *
 * \p targetQubit1 is treated as the \p least significant qubit in \p u, such that 
 * a row in \p u is dotted with the vector
 * \f$ |\text{targetQubit2} \;\; \text{targetQubit1}\rangle : \{ |00\rangle, |01\rangle, |10\rangle, |11\rangle \} \f$
 *
 * The passed 4x4 ComplexMatrix must be unitary, otherwise an error is thrown.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target1};
                \node[draw=none] at (-3.5, 2) {target2};
                \node[draw=none] at (-3.5, 4) {control};      
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw(0, 4) -- (0, 3);

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-2,2) -- (-1, 2);
                \draw (1, 2) -- (2, 2);
                \draw (-1,-1)--(-1,3)--(1,3)--(1,-1)--cycle;
                \node[draw=none] at (0, 1) {U};
                \end{tikzpicture}
    \f]
 *
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q/4 nodes.
 *
 * @see
 * - ::ComplexMatrix4
 * - multiControlledTwoQubitUnitary() 
 * - multiQubitUnitary()
 * - unitary()
 *
 * @ingroup unitary                                                          
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit the control qubit which must be in state 1 to effect the given unitary
 * @param[in] targetQubit1 first qubit to operate on, treated as least significant in \p u
 * @param[in] targetQubit2 second qubit to operate on, treated as most significant in \p u
 * @param[in] u unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p controlQubit, \p targetQubit1 or \p targetQubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if any of \p controlQubit, \p targetQubit1 and \p targetQubit2 are equal
 * - if matrix \p u is not unitary
 *  - if each node cannot fit 4 amplitudes in distributed mode.
 * @author Tyson Jones
 */
void controlledTwoQubitUnitary(Qureg qureg, int controlQubit, int targetQubit1, int targetQubit2, ComplexMatrix4 u);

/** Apply a general multi-controlled two-qubit unitary (including a global phase factor).
 * Any number of control qubits can be specified, and if all have value 1, 
  * the given unitary is applied to the target qubit.
 * This effects the many-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & \ddots \\
 * & & & u_{00} & u_{01} & u_{02} & u_{03} \\
 * & & & u_{10} & u_{11} & u_{12} & u_{13} \\
 * & & & u_{20} & u_{21} & u_{22} & u_{23} \\
 * & & & u_{30} & u_{31} & u_{32} & u_{33}
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 
 * \p targetQubit1 is treated as the \p least significant qubit in \p u, such that 
 * a row in \p u is dotted with the vector
 * \f$ |\text{targetQubit2} \;\; \text{targetQubit1}\rangle : \{ |00\rangle, |01\rangle, |10\rangle, |11\rangle \} \f$
 * 
 * The passed 4x4 ComplexMatrix must be unitary, otherwise an error is thrown.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target1};
                \node[draw=none] at (-3.5, 2) {target2};
                \node[draw=none] at (-3.5, 5) {controls};
                
                \node[draw=none] at (0, 8) {$\vdots$};
                \draw (0, 7) -- (0, 6);
                
                \draw (-2, 6) -- (2, 6);
                \draw[fill=black] (0, 6) circle (.2);
                \draw (0, 6) -- (0, 4);         
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw(0, 4) -- (0, 3);

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-2,2) -- (-1, 2);
                \draw (1, 2) -- (2, 2);
                \draw (-1,-1)--(-1,3)--(1,3)--(1,-1)--cycle;
                \node[draw=none] at (0, 1) {U};
                \end{tikzpicture}
    \f]
 *
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q/4 nodes.
 *
 * @see
 * - ::ComplexMatrix4
 * - twoQubitUnitary()
 * - controlledTwoQubitUnitary()
 * - multiQubitUnitary()
 * - unitary()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits the control qubits which all must be in state 1 to effect the given unitary
 * @param[in] numControlQubits the number of control qubits
 * @param[in] targetQubit1 first target qubit, treated as least significant in \p u
 * @param[in] targetQubit2 second target qubit, treated as most significant in \p u
 * @param[in] u unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p targetQubit1 or \p targetQubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if \p targetQubit1 equals \p targetQubit2
 * - if any qubit in \p controlQubits is outside [0, \p qureg.numQubitsRepresented)
 * - if \p controlQubits are not unique
 * - if either \p targetQubit1 and \p targetQubit2 are in \p controlQubits
 * - if matrix \p u is not unitary
 * - if each node cannot fit 4 amplitudes in distributed mode
 * @author Tyson Jones
 */
void multiControlledTwoQubitUnitary(Qureg qureg, int* controlQubits, int numControlQubits, int targetQubit1, int targetQubit2, ComplexMatrix4 u);

/** Apply a general multi-qubit unitary (including a global phase factor) with any number of target qubits.
 *
 * The first target qubit in \p targs is treated as \b least significant in \p u.
 * For example, 
 * ```
 *     multiQubitUnitary(qureg, (int []) {a, b, c}, 3, u);
 * ```
 * will invoke multiplication
 * \f[
 * \begin{pmatrix}
 * u_{00} & u_{01} & u_{02} & u_{03} & u_{04} & u_{05} & u_{06} & u_{07} \\
 * u_{10} & u_{11} & u_{12} & u_{13} & u_{14} & u_{15} & u_{16} & u_{17} \\
 * u_{20} & u_{21} & u_{22} & u_{23} & u_{24} & u_{25} & u_{26} & u_{27} \\
 * u_{30} & u_{31} & u_{32} & u_{33} & u_{34} & u_{35} & u_{36} & u_{37} \\
 * u_{40} & u_{41} & u_{42} & u_{43} & u_{44} & u_{45} & u_{46} & u_{47} \\
 * u_{50} & u_{51} & u_{52} & u_{53} & u_{54} & u_{55} & u_{56} & u_{57} \\
 * u_{60} & u_{61} & u_{62} & u_{63} & u_{64} & u_{65} & u_{66} & u_{67} \\
 * u_{70} & u_{71} & u_{72} & u_{73} & u_{74} & u_{75} & u_{76} & u_{77} \\
 * \end{pmatrix}
 * \begin{pmatrix}
 * |cba\rangle = |000\rangle \\
 * |cba\rangle = |001\rangle \\
 * |cba\rangle = |010\rangle \\
 * |cba\rangle = |011\rangle \\
 * |cba\rangle = |100\rangle \\
 * |cba\rangle = |101\rangle \\
 * |cba\rangle = |110\rangle \\
 * |cba\rangle = |111\rangle 
 * \end{pmatrix}
 * \f]
 * 
 * The passed ComplexMatrix must be unitary and be a compatible size with the specified number of
 * target qubits, otherwise an error is thrown.
 * > To left-multiply a non-unitary ::ComplexMatrixN, use applyMatrixN().
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 1) {targets};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-2,2) -- (-1, 2);
                \draw (1, 2) -- (2, 2);
                \draw (-1,-1)--(-1,3)--(1,3)--(1,-1);
                \node[draw=none] at (0, 1) {U};
                \node[draw=none] at (0, -1) {$\vdots$};
                
                \end{tikzpicture}
    \f]
 *
 * Note that in multithreaded mode, each thread will clone 2^\p numTargs amplitudes,
 * and store these in the runtime stack.
 * Using t threads, the total memory overhead of this function is t*2^\p numTargs.
 * For many targets (e.g. 16 qubits), this may cause a stack-overflow / seg-fault 
 * (e.g. on a 1 MiB stack).
 * 
 * Note too that in distributed mode, this routine requires that each node contains 
 * at least 2^\p numTargs amplitudes in the register. This means an q-qubit register (state vector or density matrix) 
 * can be distributed by at most 2^q / 2^\p numTargs nodes.
 *
 * @see
 * - createComplexMatrixN()
 * - controlledMultiQubitUnitary()
 * - multiControlledMultiQubitUnitary()
 * - applyMatrixN()
 * - twoQubitUnitary()
 * - unitary()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targs a list of the target qubits, ordered least significant to most in \p u
 * @param[in] numTargs the number of target qubits
 * @param[in] u unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if any index in \p targs is outside of [0, \p qureg.numQubitsRepresented)
 * - if \p targs are not unique
 * - if matrix \p u is not unitary
 * - if \p u is not of a compatible size with \p numTargs
 * - if a node cannot fit the required number of target amplitudes in distributed mode
 * @author Tyson Jones
 */
void multiQubitUnitary(Qureg qureg, int* targs, int numTargs, ComplexMatrixN u);

/** Apply a general controlled multi-qubit unitary (including a global phase factor).
 * One control and any number of target qubits can be specified.
 * This effects the many-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & 1 \\
 * & & & 1 \\
 * & & & & u_{00} & u_{01} & \dots  \\
 * & & & & u_{10} & u_{11} & \dots \\
 * & & & & \vdots & \vdots & \ddots
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *
 * The target qubits in \p targs are treated as ordered least significant 
 * to most significant in \p u.
 *
 * The passed ComplexMatrix must be unitary and be a compatible size with the specified number of
 * target qubits, otherwise an error is thrown.
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 1) {targets};
                \node[draw=none] at (-3.5, 4) {control};      
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw(0, 4) -- (0, 3);

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-2,2) -- (-1, 2);
                \draw (1, 2) -- (2, 2);
                \draw (-1,-1)--(-1,3)--(1,3)--(1,-1);
                \node[draw=none] at (0, 1) {U};
                \node[draw=none] at (0, -1) {$\vdots$};
                
                \end{tikzpicture}
    \f]
 *
 * Note that in multithreaded mode, each thread will clone 2^\p numTargs amplitudes,
 * and store these in the runtime stack.
 * Using t threads, the total memory overhead of this function is t*2^\p numTargs.
 * For many targets (e.g. 16 qubits), this may cause a stack-overflow / seg-fault 
 * (e.g. on a 1 MiB stack).
 *
 * Note too that in distributed mode, this routine requires that each node contains at least 2^\p numTargs amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q / 2^\p numTargs nodes.
 *
 * @see
 * - createComplexMatrixN()
 * - multiQubitUnitary()
 * - multiControlledMultiQubitUnitary()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] ctrl the control qubit
 * @param[in] targs a list of the target qubits, ordered least to most significant
 * @param[in] numTargs the number of target qubits
 * @param[in] u unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p ctrl or any index in \p targs is outside of [0, \p qureg.numQubitsRepresented)
 * - if \p targs are not unique
 * - if \p targs contains \p ctrl
 * - if matrix \p u is not unitary
 * - if a node cannot fit the required number of target amplitudes in distributed mode
 * @author Tyson Jones
 */
void controlledMultiQubitUnitary(Qureg qureg, int ctrl, int* targs, int numTargs, ComplexMatrixN u);

/** Apply a general multi-controlled multi-qubit unitary (including a global phase factor).
 * Any number of control and target qubits can be specified.
 * This effects the many-qubit unitary
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & \ddots \\
 * & & & u_{00} & u_{01} & \dots  \\
 * & & & u_{10} & u_{11} & \dots \\
 * & & & \vdots & \vdots & \ddots
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *
 * The target qubits in \p targs are treated as ordered least significant 
 * to most significant in \p u.
 *
 * The passed ::ComplexMatrixN must be unitary and be a compatible size with the specified number of
 * target qubits, otherwise an error is thrown.
 * > To left-multiply a non-unitary ::ComplexMatrixN, including control qubits,
 * > use applyMultiControlledMatrixN()
 *
    \f[
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 1) {targets};
                \node[draw=none] at (-3.5, 5) {controls};
                
                \node[draw=none] at (0, 8) {$\vdots$};
                \draw (0, 7) -- (0, 6);
                
                \draw (-2, 6) -- (2, 6);
                \draw[fill=black] (0, 6) circle (.2);
                \draw (0, 6) -- (0, 4);         
                
                \draw (-2, 4) -- (2, 4);
                \draw[fill=black] (0, 4) circle (.2);
                \draw(0, 4) -- (0, 3);

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-2,2) -- (-1, 2);
                \draw (1, 2) -- (2, 2);
                \draw (-1,-1)--(-1,3)--(1,3)--(1,-1);
                \node[draw=none] at (0, 1) {U};
                \node[draw=none] at (0, -1) {$\vdots$};
                
                \end{tikzpicture}
    \f]
 *
 * Note that in multithreaded mode, each thread will clone 2^\p numTargs amplitudes,
 * and store these in the runtime stack.
 * Using t threads, the total memory overhead of this function is t*2^\p numTargs.
 * For many targets (e.g. 16 qubits), this may cause a stack-overflow / seg-fault 
 * (e.g. on a 1 MiB stack).
 *
 * Note that in distributed mode, this routine requires that each node contains at least 2^\p numTargs amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q / 2^\p numTargs nodes.
 *
 * @see
 * - createComplexMatrixN()
 * - applyMultiControlledMatrixN()
 * - multiControlledMultiQubitNot()
 * - controlledMultiQubitUnitary()
 * - multiQubitUnitary()
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] ctrls a list of the control qubits
 * @param[in] numCtrls the number of control qubits
 * @param[in] targs a list of the target qubits, ordered least to most significant
 * @param[in] numTargs the number of target qubits
 * @param[in] u unitary matrix to apply
 * @throws invalidQuESTInputError()
 * - if any qubit in \p ctrls and \p targs is invalid, i.e. outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * - if \p ctrls or \p targs contain any repetitions
 * - if any qubit in \p ctrls is also in \p targs (and vice versa)
 * - if \p numTargs <b>< 1</b>
 * - if \p numCtrls <b>< 1</b> (use multiQubitUnitary() for no controls)
 * - if matrix \p u is not unitary
 * - if a node cannot fit the required number of target amplitudes in distributed mode
 * @throws segmentation-fault
 * - if \p ctrls contains fewer elements than \p numCtrls
 * - if \p targs contains fewer elements than \p numTargs
 * @author Tyson Jones
 */
void multiControlledMultiQubitUnitary(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, ComplexMatrixN u);

/** Apply a general single-qubit Kraus map to a density matrix, as specified by at most 
 * four Kraus operators, \f$K_i\f$ (\p ops). A Kraus map is also referred to as 
 * a "operator-sum representation" of a quantum channel, and enables the simulation of 
 * general single-qubit noise process,
 * by effecting 
 * \f[
    \rho \to \sum\limits_i^{\text{numOps}} K_i \rho K_i^\dagger
 * \f]
 *
 * The Kraus map must be completely positive and trace preserving, which constrains each 
 * \f$ K_i \f$ in \p ops by
 * \f[
    \sum \limits_i^{\text{numOps}} K_i^\dagger K_i = I
 * \f]
 * where \f$ I \f$ is the identity matrix.
 *
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register can be distributed by at most 2^(q-2) numTargs nodes.
 *
 * @see
 * - ::ComplexMatrix2
 * - mixTwoQubitKrausMap()
 * - mixMultiQubitKrausMap()
 * - mixDephasing()
 * - mixDepolarising()
 * - mixDamping()
 * - mixPauli()
 * - mixDensityMatrix()
 *
 * @ingroup decoherence
 * @param[in,out] qureg the density matrix to which to apply the map
 * @param[in] target the target qubit of the map
 * @param[in] ops an array of at most 4 Kraus operators
 * @param[in] numOps the number of operators in \p ops which must be >0 and <= 4.
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if \p target is outside of [0, \p qureg.numQubitsRepresented)
 * - if \p numOps is outside [1, 4]
 * - if \p ops do not create a completely positive, trace preserving map
 * - if a node cannot fit 4 amplitudes in distributed mode
 * @author Balint Koczor
 * @author Tyson Jones (refactored, doc)
 */
void mixKrausMap(Qureg qureg, int target, ComplexMatrix2 *ops, int numOps);

/** Apply a general two-qubit Kraus map to a density matrix, as specified by at most 
 * sixteen Kraus operators. A Kraus map is also referred to as a "operator-sum representation"
 * of a quantum channel. This allows one to simulate a general two-qubit noise process.
 *
 * The Kraus map must be completely positive and trace preserving, which constrains each 
 * \f$ K_i \f$ in \p ops by
 * \f[
    \sum \limits_i^{\text{numOps}} K_i^\dagger K_i = I
 * \f]
 * where \f$ I \f$ is the identity matrix.
 *
 * \p targetQubit1 is treated as the \p least significant qubit in each op in \p ops.
 *
 * Note that in distributed mode, this routine requires that each node contains at least 16 amplitudes.
 * This means an q-qubit register can be distributed by at most 2^(q-4) numTargs nodes.
 *
 * @see
 * - ::ComplexMatrix4
 * - mixMultiQubitKrausMap()
 * - mixKrausMap()
 *
 * @ingroup decoherence
 * @param[in,out] qureg the density matrix to which to apply the map
 * @param[in] target1 the least significant target qubit in \p ops
 * @param[in] target2 the most significant target qubit in \p ops
 * @param[in] ops an array of at most 16 Kraus operators
 * @param[in] numOps the number of operators in \p ops which must be >0 and <= 16.
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if either \p target1 or \p target2 is outside of [0, \p qureg.numQubitsRepresented)
 * - if \p target1 = \p target2
 * - if \p numOps is outside [1, 16]
 * - if \p ops do not create a completely positive, trace preserving map
 * - if a node cannot fit 16 amplitudes in distributed mode
 * @author Balint Koczor
 * @author Tyson Jones (refactored, doc)
 */
void mixTwoQubitKrausMap(Qureg qureg, int target1, int target2, ComplexMatrix4 *ops, int numOps);

/** Apply a general N-qubit Kraus map to a density matrix, as specified by at most (2N)^2
 * Kraus operators. This allows one to simulate a general noise process.
 *
 * The Kraus map must be completely positive and trace preserving, which constrains each 
 * \f$ K_i \f$ in \p ops by
 * \f[
    \sum \limits_i^{\text{numOps}} K_i^\dagger K_i = I
 * \f]
 * where \f$ I \f$ is the identity matrix.
 *
 * The first qubit in \p targets is treated as the \p least significant qubit in each op in \p ops.
 *
 * Note that in distributed mode, this routine requires that each node contains at least (2N)^2 amplitudes.
 * This means an q-qubit register can be distributed by at most 2^(q-2)/N^2 nodes.
 *
 * Note too that this routine internally creates a 'superoperator'; a complex matrix of dimensions
 * 2^(2*numTargets) by 2^(2*numTargets). Therefore, invoking this function incurs, 
 * for numTargs={1,2,3,4,5, ...}, an additional memory overhead of (at double-precision)
 * {0.25 KiB, 4 KiB, 64 KiB, 1 MiB, 16 MiB, ...} (respectively).
 * At quad precision (usually 10 B per number, but possibly 16 B due to alignment),
 * this costs at most double the amount of memory. 
 * For numTargets < 4, this superoperator will be created in the runtime 
 * stack. For numTargs >= 4, the superoperator will be allocated in the heap and 
 * therefore this routine may suffer an anomalous slowdown.
 *
 * @see
 * - createComplexMatrixN()
 * - initComplexMatrixN()
 * - mixKrausMap()
 * - mixTwoQubitKrausMap()
 *
 * @ingroup decoherence
 * @param[in,out] qureg the density matrix to which to apply the map
 * @param[in] targets a list of target qubit indices, the first of which is treated as least significant in each op in \p ops
 * @param[in] numTargets the length of \p targets
 * @param[in] ops an array of at most (2N)^2 Kraus operators
 * @param[in] numOps the number of operators in \p ops which must be >0 and <= (2N)^2.
 * @throws invalidQuESTInputError()
 * - if \p qureg is not a density matrix
 * - if any target in \p targets is outside of [0, \p qureg.numQubitsRepresented)
 * - if any qubit in \p targets is repeated
 * - if \p numOps is outside [1, (2 \p numTargets)^2]
 * - if any ComplexMatrixN in \p ops does not have op.numQubits == \p numTargets
 * - if \p ops do not create a completely positive, trace preserving map
 * - if a node cannot fit (2N)^2 amplitudes in distributed mode
 * @author Tyson Jones
 * @author Balint Koczor
 */
void mixMultiQubitKrausMap(Qureg qureg, int* targets, int numTargets, ComplexMatrixN* ops, int numOps);

/** Computes the Hilbert Schmidt distance between two density matrices \p a and \p b, 
 * defined as the Frobenius norm of the difference between them.
 * That is, we define the Hilbert Schmidt distance
 * \f[
    D(a, b) = \| a - b \|_F = \sqrt{  \text{Tr}[ (a-b)(a-b)^\dagger ]   }
 * \f]
 * This is equivalent to the square-root of the sum of the absolute value squared of the 
 * element-differences of the matrices, i.e.
 * \f[
    D(a, b) = \sqrt{ \sum\limits_i \sum\limits_j | a_{ij} - b_{ij} |^2 }
 * \f]
 * We caution this may differ by some definitions of the Hilbert Schmidt distance 
 * by a square-root.
 *
 * This function correctly returns the result of the above formulations even when 
 * \p a and \p b are incorrectly normalised (i.e. are general matrices).
 *
 * @see
 * - calcDensityInnerProduct()
 * - calcFidelity()
 * - calcPurity()
 *
 * @ingroup calc
 * @param[in] a a density matrix
 * @param[in] b an equally-sized density matrix
 * @throws invalidQuESTInputError()
 * - if either \p a or \p b are not density matrices
 * - if \p a and \p have mismatching dimensions
 * @author Balint Koczor
 * @author Tyson Jones (refactored, doc)
 */
qreal calcHilbertSchmidtDistance(Qureg a, Qureg b);

/** Modifies qureg \p out to the result of (\p facOut \p out + \p fac1 \p qureg1 + \p fac2 \p qureg2), 
 * imposing no constraints on normalisation. Works for both state-vectors and density matrices.
 * Note that afterward, \p out may not longer be normalised and ergo no longer a valid 
 * state-vector or density matrix. Users must therefore be careful passing \p out to
 * other QuEST functions which assume normalisation in order to function correctly.
 *
 * \p qureg1, \p qureg2 and \p out must be all state-vectors, or all density matrices,
 * of equal dimensions. \p out can be one (or both) of \p qureg1 and \p qureg2.
 *
 * @ingroup init
 * @param[in] fac1 the complex number by which to scale \p qureg1 in the output state 
 * @param[in] qureg1 the first qureg to add to \p out, which is itself unmodified
 * @param[in] fac2 the complex number by which to scale \p qureg2 in the output state 
 * @param[in] qureg2 the second qureg to add to \p out, which is itself unmodified
 * @param[in] facOut the complex factor by which to multiply the current elements of \p out.
 *      \p out is completely overwritten if \p facOut is set to (Complex) {.real=0,.imag=0}
 * @param[in,out] out the qureg to be modified, to be scaled by \p facOut then have \p fac1 \p qureg1 and
 *      \p fac2 \p qureg2 added to it.
 * @throws invalidQuESTInputError()
 * - if \p qureg1, \p qureg2 and \p aren't all state-vectors or all density matrices
 * - if the dimensions of \p qureg1, \p qureg2 and \p aren't equal
 * @author Tyson Jones
 */
void setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out);

/** Modifies \p outQureg to be the result of applying the weighted sum of Pauli products (a Hermitian but not 
 * necessarily unitary operator) to \p inQureg. Note that afterward, \p outQureg may no longer be normalised and ergo not a
 * state-vector or density matrix. Users must therefore be careful passing \p outQureg to
 * other QuEST functions which assume normalisation in order to function correctly.
 *
 * Letting \f$ \alpha = \sum_i c_i \otimes_j^{N} \hat{\sigma}_{i,j} \f$ be 
 * the operators indicated by \p allPauliCodes (where \f$ c_i \in \f$ \p termCoeffs and \f$ N = \f$ \p qureg.numQubitsRepresented), 
 * this function effects \f$ \alpha | \psi \rangle \f$ on state-vector \f$ |\psi\rangle \f$
 * and \f$\alpha \rho\f$ (left matrix multiplication) on density matrix \f$ \rho \f$.
 *
 * \p allPauliCodes is an array of length \p numSumTerms*\p qureg.numQubitsRepresented
 *  which specifies which Pauli operators to apply, where 0 = \p PAULI_I, 1 = \p PAULI_X, 
 * 2 = \p PAULI_Y, 3 = \p PAULI_Z. For each sum term, a Pauli operator must be specified for 
 * EVERY qubit in \p qureg; each set of \p numSumTerms operators will be grouped into a product.
 * \p termCoeffs is an arrray of length \p numSumTerms containing the term coefficients.
 * For example, on a 3-qubit state-vector,
 * ```
 *     int paulis[6] = {PAULI_X, PAULI_I, PAULI_I,  PAULI_X, PAULI_Y, PAULI_Z};
 *     qreal coeffs[2] = {1.5, -3.6};
 *     applyPauliSum(inQureg, paulis, coeffs, 2, outQureg);
 * ```
 * will apply Hermitian operation \f$ (1.5 X I I - 3.6 X Y Z) \f$ 
 * (where in this notation, the left-most operator applies to the least-significant qubit, i.e. that with index 0).
 *
 * In theory, \p inQureg is unchanged though its state is temporarily 
 * modified and is reverted by re-applying Paulis (XX=YY=ZZ=I), so may see a change by small numerical errors.
 * The initial state in \p outQureg is not used.
 *
 * \p inQureg and \p outQureg must both be state-vectors, or both density matrices,
 * of equal dimensions. \p inQureg cannot be \p outQureg.
 *
 * This function works by applying each Pauli product to \p inQureg in turn, 
 * and adding the resulting state (weighted by a coefficient in \p termCoeffs)
 * to the initially-blanked \p outQureg. Ergo it should scale with the total number 
 * of Pauli operators specified (excluding identities), and the qureg dimension. 
 *
 * @see
 * - calcExpecPauliSum()
 * - applyPauliHamil()
 *
 * @ingroup operator
 * @param[in] inQureg the register containing the state which \p outQureg will be set to, under
 *      the action of the Hermitiain operator specified by the Pauli codes. \p inQureg should be 
 *      unchanged, though may vary slightly due to numerical error.
 * @param[in] allPauliCodes a list of the Pauli codes (0=PAULI_I, 1=PAULI_X, 2=PAULI_Y, 3=PAULI_Z) 
 *      of all Paulis involved in the products of terms. A Pauli must be specified for each qubit 
 *      in the register, in every term of the sum.
 * @param[in] termCoeffs The coefficients of each term in the sum of Pauli products
 * @param[in] numSumTerms The total number of Pauli products specified
 * @param[out] outQureg the qureg to modify to be the result of applyling the weighted Pauli sum operator 
 *      to the state in \p inQureg
 * @throws invalidQuESTInputError()
 * - if any code in \p allPauliCodes is not in {0,1,2,3}
 * - if numSumTerms <= 0
 * - if \p inQureg is not of the same type and dimensions as \p outQureg
 * @author Tyson Jones
 */
void applyPauliSum(Qureg inQureg, enum pauliOpType* allPauliCodes, qreal* termCoeffs, int numSumTerms, Qureg outQureg);

/** Modifies \p outQureg to be the result of applying \p PauliHamil (a Hermitian but not 
 * necessarily unitary operator) to \p inQureg. Note that afterward, \p outQureg may no longer be normalised and ergo not a
 * state-vector or density matrix. Users must therefore be careful passing \p outQureg to
 * other QuEST functions which assume normalisation in order to function correctly.
 *
 * This is merely an encapsulation of applyPauliSum(), which can refer to for elaborated doc.
 *
 * Letting \p hamil be expressed as \f$ \alpha = \sum_i c_i \otimes_j^{N} \hat{\sigma}_{i,j} \f$ 
 * (where \f$ c_i \in \f$ \p hamil.termCoeffs and \f$ N = \f$ \p hamil.numQubits), 
 * this function effects \f$ \alpha | \psi \rangle \f$ on state-vector \f$ |\psi\rangle \f$
 * and \f$\alpha \rho\f$ (left matrix multiplication) on density matrix \f$ \rho \f$.
 *
 * In theory, \p inQureg is unchanged though its state is temporarily 
 * modified and is reverted by re-applying Paulis (XX=YY=ZZ=I), so may see a change by small numerical errors.
 * The initial state in \p outQureg is not used.
 *
 * \p inQureg and \p outQureg must both be state-vectors, or both density matrices,
 * of equal dimensions to \p hamil.
 * \p inQureg cannot be \p outQureg.
 *
 * This function works by applying each Pauli product in \p hamil to \p inQureg in turn, 
 * and adding the resulting state (weighted by a coefficient in \p termCoeffs)
 * to the initially-blanked \p outQureg. Ergo it should scale with the total number 
 * of Pauli operators specified (excluding identities), and the qureg dimension. 
 *
 * @see
 * - createPauliHamil()
 * - createPauliHamilFromFile()
 * - calcExpecPauliHamil()
 * - applyTrotterCircuit()
 *
 * @ingroup operator
 * @param[in] inQureg the register containing the state which \p outQureg will be set to, under
 *      the action of \p hamil. \p inQureg should be unchanged, though may vary slightly due to numerical error.
 * @param[in] hamil a weighted sum of products of pauli operators
 * @param[out] outQureg the qureg to modify to be the result of applyling \p hamil to the state in \p inQureg
 * @throws invalidQuESTInputError()
 * - if any code in \p hamil.pauliCodes is not a valid Pauli code
 * - if \p numSumTerms <= 0
 * - if \p inQureg is not of the same type and dimensions as \p outQureg and \p hamil
 * @author Tyson Jones
 */
void applyPauliHamil(Qureg inQureg, PauliHamil hamil, Qureg outQureg);

/** Applies a trotterisation of unitary evolution \f$ \exp(-i \, \text{hamil} \, \text{time}) \f$
 * to \p qureg. This is a sequence of unitary operators, effected by multiRotatePauli(),
 * which together approximate the action of full unitary-time evolution under the given Hamiltonian.
 *
 * Notate \f$ \text{hamil} = \sum_j^N c_j \, \hat \sigma_j \f$ where \f$c_j\f$ is a real 
 * coefficient in \p hamil, \f$\hat \sigma_j\f$ is the corresponding product of Pauli operators,
 * of which there are a total \f$N\f$.
 * Then, \p order=1 performs first-order Trotterisation, whereby
 * \f[
 *   \exp(-i \, \text{hamil} \, \text{time})
 *      \approx 
 *    \prod\limits^{\text{reps}} \prod\limits_{j=1}^{N} \exp(-i \, c_j \, \text{time} \, \hat\sigma_j / \text{reps})
 * \f]
 * \p order=2 performs the lowest order "symmetrized" Suzuki decomposition, whereby 
 * \f[
 *   \exp(-i \, \text{hamil} \, \text{time})
 *      \approx 
 *    \prod\limits^{\text{reps}} \left[
 *         \prod\limits_{j=1}^{N} \exp(-i \, c_j \, \text{time} \, \hat\sigma_j / (2 \, \text{reps}))
 *          \prod\limits_{j=N}^{1} \exp(-i \, c_j \, \text{time} \, \hat\sigma_j / (2 \, \text{reps}))
 *     \right]
 * \f]
 * Greater even values of \p order specify higher-order symmetrized decompositions 
 * \f$ S[\text{time}, \text{order}, \text{reps}] \f$ which satisfy 
 * \f[
 *      S[\text{time}, \text{order}, 1] = 
 *          \left( \prod\limits^2 S[p \, \text{time}, \text{order}-2, 1] \right)
 *          S[ (1-4p)\,\text{time}, \text{order}-2, 1]
 *          \left( \prod\limits^2 S[p \, \text{time}, \text{order}-2, 1] \right)
 * \f]
 * and 
 * \f[
 *      S[\text{time}, \text{order}, \text{reps}] = 
 *          \prod\limits^{\text{reps}} S[\text{time}/\text{reps}, \text{order}, 1]
 * \f]
 * where \f$ p = \left( 4 - 4^{1/(\text{order}-1)} \right)^{-1} \f$.
 * 
 * These formulations are taken from 'Finding Exponential Product Formulas
 * of Higher Orders', Naomichi Hatano and Masuo Suzuki (2005) (<a href="https://arxiv.org/abs/math-ph/0506007">arXiv</a>).
 *
 * Note that the applied Trotter circuit is captured by QASM, if QASM logging is enabled
 * on \p qureg. \n
 * For example:
 * ```
 * startRecordingQASM(qureg);
 * applyTrotterCircuit(qureg, hamil, 1, 2, 1);
 * printRecordedQASM(qureg);
 * ```
 * may show 
 * ```
 * // Beginning of Trotter circuit (time 1, order 2, 1 repetitions).
 * // Here, a multiRotatePauli with angle 0.5 and paulis X Y I  was applied.
 * // Here, a multiRotatePauli with angle -0.5 and paulis I Z X  was applied.
 * // Here, a multiRotatePauli with angle -0.5 and paulis I Z X  was applied.
 * // Here, a multiRotatePauli with angle 0.5 and paulis X Y I  was applied.
 * // End of Trotter circuit
 * ``` 
 * \n
 *
 *
 * @see
 * - createPauliHamil()
 *
 * @ingroup operator
 * @param[in,out] qureg the register to modify under the approximate unitary-time evolution
 * @param[in] hamil the hamiltonian under which to approxiamte unitary-time evolution
 * @param[in] time the target evolution time, which is permitted to be both positive and negative.
 * @param[in] order the order of Trotter-Suzuki decomposition to use. Higher orders (necessarily even)
 *      are more accurate but prescribe an exponentially increasing number of gates.
 * @param[in] reps the number of repetitions of the decomposition of the given order. This 
 *      improves the accuracy but prescribes a linearly increasing number of gates.
 * @throws invalidQuESTInputError()
 * - if \p qureg.numQubitsRepresented != \p hamil.numQubits
 * - if \p hamil contains invalid parameters or Pauli codes,
 * - if \p order is not in {1, 2, 4, 6, ...}
 * - or if \p reps <= 0
 * @author Tyson Jones
 */
void applyTrotterCircuit(Qureg qureg, PauliHamil hamil, qreal time, int order, int reps);

/** Apply a general 2-by-2 matrix, which may be non-unitary. The matrix is 
 * left-multiplied onto the state, for both state-vectors and density matrices.
 * 
 * Note this differs from the action of unitary() on a density matrix.
 * 
 * This function may leave \p qureg is an unnormalised state.
 *
 * @see
 * - ::ComplexMatrix2 
 * - unitary()
 *
 * @ingroup operator                                                              
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate \p u upon
 * @param[in] u matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Tyson Jones
 */
void applyMatrix2(Qureg qureg, int targetQubit, ComplexMatrix2 u);

/** Apply a general 4-by-4 matrix, which may be non-unitary. The matrix is 
 * left-multiplied onto the state, for both state-vectors and density matrices.
 * 
 * Note this differs from the action of twoQubitUnitary() on a density matrix.
 *
 * \p targetQubit1 is treated as the \p least significant qubit in \p u, such that 
 * a row in \p u is dotted with the vector
 * \f$ |\text{targetQubit2} \;\; \text{targetQubit1}\rangle : \{ |00\rangle, |01\rangle, |10\rangle, |11\rangle \} \f$
 *
 * For example, 
 * ```
 *     applyMatrix4(qureg, a, b, u);
 * ```
 * will invoke multiplication
 * \f[
 * \begin{pmatrix}
 * u_{00} & u_{01} & u_{02} & u_{03} \\
 * u_{10} & u_{11} & u_{12} & u_{13} \\
 * u_{20} & u_{21} & u_{22} & u_{23} \\
 * u_{30} & u_{31} & u_{32} & u_{33}
 * \end{pmatrix}
 * \begin{pmatrix}
 * |ba\rangle = |00\rangle \\
 * |ba\rangle = |01\rangle \\
 * |ba\rangle = |10\rangle \\
 * |ba\rangle = |11\rangle 
 * \end{pmatrix}
 * \f]
 *
 * This function may leave \p qureg is an unnormalised state.
 *                 
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q/4 nodes.
 *
 * @see
 * - ::ComplexMatrix4
 * - twoQubitUnitary()
 * 
 * @ingroup operator                                            
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit1 first qubit to operate on, treated as least significant in \p u
 * @param[in] targetQubit2 second qubit to operate on, treated as most significant in \p u
 * @param[in] u matrix to apply
 * @throws invalidQuESTInputError()
 * - if \p targetQubit1 or \p targetQubit2 are outside [0, \p qureg.numQubitsRepresented)
 * - if \p targetQubit1 equals \p targetQubit2
 * - if each node cannot fit 4 amplitudes in distributed mode
 * @author Tyson Jones
 */
void applyMatrix4(Qureg qureg, int targetQubit1, int targetQubit2, ComplexMatrix4 u);

/** Apply a general N-by-N matrix, which may be non-unitary, on any number of target qubits.
 * The matrix is left-multiplied onto the state, for both state-vectors and density matrices.
 * Note this differs from the action of multiQubitUnitary() on a density matrix.
 *
 * The first target qubit in \p targs is treated as \b least significant in \p u.
 * For example, 
 * ```
 *     applyMatrixN(qureg, (int []) {a, b, c}, 3, u);
 * ```
 * will invoke multiplication
 * \f[
 * \begin{pmatrix}
 * u_{00} & u_{01} & u_{02} & u_{03} & u_{04} & u_{05} & u_{06} & u_{07} \\
 * u_{10} & u_{11} & u_{12} & u_{13} & u_{14} & u_{15} & u_{16} & u_{17} \\
 * u_{20} & u_{21} & u_{22} & u_{23} & u_{24} & u_{25} & u_{26} & u_{27} \\
 * u_{30} & u_{31} & u_{32} & u_{33} & u_{34} & u_{35} & u_{36} & u_{37} \\
 * u_{40} & u_{41} & u_{42} & u_{43} & u_{44} & u_{45} & u_{46} & u_{47} \\
 * u_{50} & u_{51} & u_{52} & u_{53} & u_{54} & u_{55} & u_{56} & u_{57} \\
 * u_{60} & u_{61} & u_{62} & u_{63} & u_{64} & u_{65} & u_{66} & u_{67} \\
 * u_{70} & u_{71} & u_{72} & u_{73} & u_{74} & u_{75} & u_{76} & u_{77} \\
 * \end{pmatrix}
 * \begin{pmatrix}
 * |cba\rangle = |000\rangle \\
 * |cba\rangle = |001\rangle \\
 * |cba\rangle = |010\rangle \\
 * |cba\rangle = |011\rangle \\
 * |cba\rangle = |100\rangle \\
 * |cba\rangle = |101\rangle \\
 * |cba\rangle = |110\rangle \\
 * |cba\rangle = |111\rangle 
 * \end{pmatrix}
 * \f]
 *
 * This function may leave \p qureg is an unnormalised state.
 * 
 * The passed ComplexMatrix must be a compatible size with the specified number of
 * target qubits, otherwise an error is thrown.
 *
 * Note that in multithreaded mode, each thread will clone 2^\p numTargs amplitudes,
 * and store these in the runtime stack.
 * Using t threads, the total memory overhead of this function is t*2^\p numTargs.
 * For many targets (e.g. 16 qubits), this may cause a stack-overflow / seg-fault 
 * (e.g. on a 1 MiB stack).
 * 
 * Note too that in distributed mode, this routine requires that each node contains 
 * at least 2^\p numTargs amplitudes in the register. This means an q-qubit register (state vector or density matrix) 
 * can be distributed by at most 2^q / 2^\p numTargs nodes.
 *
 * @see
 * - createComplexMatrixN()
 * - getStaticComplexMatrixN()
 * - applyMultiControlledMatrixN()
 * - multiQubitUnitary()
 *
 * @ingroup operator
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targs a list of the target qubits, ordered least significant to most in \p u
 * @param[in] numTargs the number of target qubits
 * @param[in] u matrix to apply
 * @throws invalidQuESTInputError()
 * - if any index in \p targs is outside of [0, \p qureg.numQubitsRepresented)
 * - if \p targs are not unique
 * - if \p u is not of a compatible size with \p numTargs
 * - if a node cannot fit the required number of target amplitudes in distributed mode
 * @author Tyson Jones
 */
void applyMatrixN(Qureg qureg, int* targs, int numTargs, ComplexMatrixN u);

/** Apply a general N-by-N matrix, which may be non-unitary, with additional controlled qubits.
 * The matrix is left-multiplied onto the state, for both state-vectors and density matrices.
 * Hence, this function differs from multiControlledMultiQubitUnitary() by more than just permitting a non-unitary 
 * matrix.
 *
 * This function may leave \p qureg is an unnormalised state.
*
 * Any number of control and target qubits can be specified.
 * This effects the many-qubit matrix
 * \f[
 * \begin{pmatrix}
 * 1 \\
 * & 1 \\\
 * & & \ddots \\
 * & & & u_{00} & u_{01} & \dots  \\
 * & & & u_{10} & u_{11} & \dots \\
 * & & & \vdots & \vdots & \ddots
 * \end{pmatrix}
 * \f]
 * on the control and target qubits.
 *
 * The target qubits in \p targs are treated as ordered least significant 
 * to most significant in \p u.
 *
 * The passed ComplexMatrix must be a compatible size with the specified number of
 * target qubits, otherwise an error is thrown.
 *
 * Note that in multithreaded mode, each thread will clone 2^\p numTargs amplitudes,
 * and store these in the runtime stack.
 * Using t threads, the total memory overhead of this function is t*2^\p numTargs.
 * For many targets (e.g. 16 qubits), this may cause a stack-overflow / seg-fault 
 * (e.g. on a 1 MiB stack).
 *
 * Note that in distributed mode, this routine requires that each node contains at least 2^\p numTargs amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q / 2^\p numTargs nodes.
 *
 * @ingroup operator
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] ctrls a list of the control qubits
 * @param[in] numCtrls the number of control qubits
 * @param[in] targs a list of the target qubits, ordered least to most significant
 * @param[in] numTargs the number of target qubits
 * @param[in] u matrix to apply
 * @throws invalidQuESTInputError()
 * - if any index in \p ctrls and \p targs is outside of [0, \p qureg.numQubitsRepresented)
 * - if \p ctrls and \p targs are not unique
 * - if matrix \p u is not a compatible size with \p numTargs
 * - if a node cannot fit the required number of target amplitudes in distributed mode
 * @author Tyson Jones
 */
void applyMultiControlledMatrixN(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, ComplexMatrixN u);

/** An internal function called when invalid arguments are passed to a QuEST API
 * call, which the user can optionally override by redefining. This function is 
 * a weak symbol, so that users can choose how input errors are handled, by 
 * redefining it in their own code. Users must ensure that the triggered API 
 * call does not continue (e.g. the user exits or throws an exception), else 
 * QuEST will continue with the valid input and likely trigger a seg-fault.
 * This function is triggered before any internal state-change, hence it is 
 * safe to interrupt with exceptions.
 *
 * E.g. in C
 * @code
void invalidQuESTInputError(const char* errMsg, const char* errFunc) {
 
     // log to file
     printf("ERROR! Writing to file...\n");
     FILE *fp = fopen("errorlog.txt", "w");
     fprintf(fp, "incorrect usage of function '%s': %s", errFunc, errMsg);
     fclose(fp);
     
     // exit
     exit(1);
} @endcode
 * 
 * This function is compatible with C++ exceptions, though note the user of 
 * extern "C":
 * @code
extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {

    string err = "in function " + string(errFunc) + ": " + string(errMsg);
    throw std::invalid_argument(err);
} @endcode
 *
 * @ingroup debug
 * @param[in] errMsg a string describing the nature of the argument error
 * @param[in] errFunc the name of the invalidly-called API function
 * @throws invalidQuESTInputError()
 * - unless overriden by the user
 * @author Tyson Jones
 */
void invalidQuESTInputError(const char* errMsg, const char* errFunc);
 
#ifndef __cplusplus
#ifndef _WIN32
 // hide this function from doxygen
 /// \cond HIDDEN_SYMBOLS
/** Creates a ComplexMatrixN struct with .real and .imag arrays kept entirely 
 * in the stack. 
 * This function should not be directly called by the user; instead, users should 
 * call the macro getStaticComplexMatrixN.
 *
 * The passed 'reStorage' and 'imStorage' must be arrays with 
 * length (1<<numQubits) and are populated with pointers to rows of the 2D 
 * arrays 're' and 'im', and then attached to the returned ComplexMatrixN instance.
 * For example:
 * ```
 *     ComplexMatrixN m = bindArraysToStackComplexMatrixN(
 *         1, 
 *         (qreal[][2]) {{1,0},{0,1}}, (qreal[][2]) {{0}}, 
 *         (qreal*[2]) {0}, (qreal*[2]) {0}
 *     );
 * ```
 * Note that this ComplexMatrixN instance, since kept in the stack, cannot be *returned*
 * beyond the calling scope which would result in a dangling pointer.
 * This is unlike a ComplexMatrixN instance created with createComplexMatrixN, which 
 * is dynamic (lives in heap) and can be returned, through needs explicit freeing 
 * with destroyComplexMatrixN.
 *
 * This function is only callable in C, since C++ signatures cannot contain 
 * variable-length 2D arrays.
 *
 * @ingroup type
 * @param[in] numQubits the number of qubits that the ComplexMatrixN corresponds to.
 *  note the .real and .imag arrays of the returned ComplexMatrixN will have 
 *  2^numQubits rows and columns.
 * @param[in] re a 2D array (2^numQubits by 2^numQubits) of the real components
 * @param[in] im a 2D array (2^numQubits by 2^numQubits) of the imag components
 * @param[in] reStorage a 1D array of length 2^numQubits
 * @param[in] imStorage a 1D array of length 2^numQubits
 * @returns a ComplexMatrixN struct with stack-stored .real and .imag arrays,
 *  which are actually the passed \p reStorage and \p imStorage arrays, populated 
 *  with pointers to the rows of \p re and \p im
 * @author Tyson Jones
 */
ComplexMatrixN bindArraysToStackComplexMatrixN(
    int numQubits, qreal re[][1<<numQubits], qreal im[][1<<numQubits], 
    qreal** reStorage, qreal** imStorage);
#endif
#endif
/// \endcond

// hide this function from doxygen
/// \cond HIDDEN_SYMBOLS
#define UNPACK_ARR(...) __VA_ARGS__
/// \endcond

#ifndef __cplusplus
/** Creates a ComplexMatrixN struct which lives in the stack and so does not 
 * need freeing, but cannot be returned beyond the calling scope. That is, 
 * the .real and .imag arrays of the returned ComplexMatrixN live in the stack
 * as opposed to that returned by createComplexMatrixN() (which live in the heap).
 * Note the real and imag components must be wrapped in paranthesis, e.g.
 * ```
 *     ComplexMatrixN u = getStaticComplexMatrixN(1, ({{1,2},{3,4}}), ({{0}}));
 * ```
 *
 * Here is an example of an incorrect usage, since a 'local' ComplexMatrixN cannot
 * leave the calling scope (otherwise inducing dangling pointers):
 * ```
 *     ComplexMatrixN getMyMatrix(void) {
 *         return getStaticComplexMatrixN(1, ({{1,2},{3,4}}), ({{0}}));
 *     }
 * ```
 *
 * This function is actually a single-line anonymous macro, so can be safely 
 * invoked within arguments to other functions, e.g.
 * ```
 *      multiQubitUnitary(
 *          qureg, (int[]) {0}, 1, 
 *          getStaticComplexMatrixN(1, ({{1,0},{0,1}}), ({{0}}))
 *      );
 * ```
 *
 * The returned ComplexMatrixN can be accessed and modified in the same way as
 * that returned by createComplexMatrixN(), e.g.
 * ```
 *      ComplexMatrixN u = getStaticComplexMatrixN(3, ({{0}}), ({{0}}));
 *      for (int i=0; i<8; i++)
 *          for (int j=0; j<8; j++)
 *              u.real[i][j] = .1;
 * ```
 * \n
 *
 * > Note that the first argument \p numQubits must be a literal.
 *
 * > This macro is only callable in C, since it invokes the function 
 * > bindArraysToStackComplexMatrixN() which is only callable in C.
 *
 * @see 
 * - createComplexMatrixN()
 *
 * @ingroup type
 * @author Tyson Jones
 */
#define getStaticComplexMatrixN(numQubits, re, im) \
    bindArraysToStackComplexMatrixN( \
        numQubits, \
        (qreal[1<<numQubits][1<<numQubits]) UNPACK_ARR re, \
        (qreal[1<<numQubits][1<<numQubits]) UNPACK_ARR im, \
        (double*[1<<numQubits]) {NULL}, (double*[1<<numQubits]) {NULL} \
    )
#endif

/** Induces a phase change upon each amplitude of \p qureg, determined by the passed 
 * exponential polynomial "phase function". This effects a diagonal unitary of unit complex scalars,
 * targeting the nominated \p qubits.
 *
 * - Arguments \p coeffs and \p exponents together specify a real exponential polynomial \f$f(r)\f$ 
 *   with \p numTerms terms, of the form
 *   \f[ 
 *    f(r) = \sum\limits_{i}^{\text{numTerms}} \text{coeffs}[i] \; r^{\, \text{exponents}[i]}\,,
 *   \f] 
 *   where both \p coeffs and \p exponents can be negative, positive and fractional. 
 *   For example,
 *   ```
 *      qreal coeffs[] = {1, -3.14};
 *      qreal exponents[] = {2, -5.5};
 *      int numTerms = 2;
 *   ```
 *   constitutes the function 
 *   \f[
 *       f(r) =  1 \, r^2 - 3.14 \, r^{-5.5}.
 *   \f]
 *   Note you cannot use fractional exponents with \p encoding <b>=</b> ::TWOS_COMPLEMENT, 
 *   since the negative indices would generate (illegal) complex phases, and must 
 *   be overriden with applyPhaseFuncOverrides(). \n
 * > If your function \f$f(r)\f$ diverges at one or more \f$r\f$ values, you must instead 
 * > use applyPhaseFuncOverrides() and specify explicit phase changes for these values.
 * > Otherwise, the corresponding amplitudes of the state-vector will become indeterminate (like `NaN`).
 * > Note that use of any negative exponent will result in divergences at \f$r=0\f$.
 *
 * - The function \f$f(r)\f$ specifies the phase change to induce upon amplitude \f$\alpha\f$ 
 *   of computational basis state with index \f$r\f$, such that
 *   \f[
 *    \alpha \, |r\rangle \rightarrow \, \exp(i f(r)) \; \alpha \, |r\rangle.
 *   \f]
 *   The index \f$r\f$ associated with each computational basis state is determined by 
 *   the binary value of the specified \p qubits (ordered least to most significant), 
 *   interpreted under the given ::bitEncoding \p encoding. \n\n
 *   For example, under \p encoding <b>=</b> \p UNSIGNED and \p qubits <b>= {0,1}</b>,
 *   \f[ 
 *   \begin{aligned}
 *     |0\mathbf{00}\rangle & \rightarrow \, e^{i f(0)}\,|0\mathbf{00}\rangle \\
 *     |0\mathbf{01}\rangle & \rightarrow \, e^{i f(1)}\,|0\mathbf{01}\rangle \\
 *     |0\mathbf{10}\rangle & \rightarrow \, e^{i f(2)}\,|0\mathbf{10}\rangle \\
 *     |0\mathbf{11}\rangle & \rightarrow \, e^{i f(3)}\,|0\mathbf{11}\rangle \\
 *     |1\mathbf{00}\rangle & \rightarrow \, e^{i f(0)}\,|1\mathbf{00}\rangle \\
 *     |1\mathbf{01}\rangle & \rightarrow \, e^{i f(1)}\,|1\mathbf{01}\rangle \\
 *     |1\mathbf{10}\rangle & \rightarrow \, e^{i f(2)}\,|1\mathbf{10}\rangle \\
 *     |1\mathbf{11}\rangle & \rightarrow \, e^{i f(3)}\,|1\mathbf{11}\rangle
 *   \end{aligned}
 *   \f]
 *
 * - If \p qureg is a density matrix \f$\rho\f$, this function modifies \p qureg to
 *   \f[
 *      \rho \rightarrow \hat{D} \, \rho \, \hat{D}^\dagger
 *   \f]
 *   where \f$\hat{D}\f$ is the diagonal unitary operator 
 *   \f[
 *      \hat{D} = \text{diag} \, \{ \; e^{i f(r_0)}, \; e^{i f(r_1)}, \;  \dots \; \}.
 *   \f]
 *   This means element \f$\rho_{jk}\f$ is modified to
 *   \f[
 *      \alpha \, |j\rangle\langle k| \; \rightarrow \; e^{i (f(r_j) - f(r_k))} \; \alpha \, |j\rangle\langle k|
 *   \f]\n
 *
 * - The interpreted phase function can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyPhaseFunc(qureg, ...);
 *   printRecordedQASM(qureg);
 *   ```
 *   may show 
 *   ```
 *   // Here, applyPhaseFunc() multiplied a complex scalar of the form
 *   //     exp(i (1 x^3))
 *   //   upon every substate |x>, informed by qubits (under an unsigned binary encoding)
 *   //     {4, 1, 2, 0}
 *   ``` 
 *
 * > This function may become numerically imprecise for quickly growing phase functions 
 * > which admit very large phases, for example of 10^10. 
 *
 * @see 
 * - applyPhaseFuncOverrides() to override the phase function for specific states.
 * - applyMultiVarPhaseFunc() for multi-variable exponential polynomial phase functions.
 * - applyNamedPhaseFunc() for a set of specific phase functions.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 *
 * @ingroup operator
 * @param[in,out] qureg the state-vector or density matrix to be modified
 * @param[in] qubits a list of the indices of the qubits which will inform \f$r\f$ for each amplitude in \p qureg
 * @param[in] numQubits the length of list \p qubits
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r\f$ from the bits of \p qubits in each basis state of \p qureg
 * @param[in] coeffs the coefficients of the exponential polynomial phase function \f$f(r)\f$
 * @param[in] exponents the exponents of the exponential polynomial phase function \f$f(r)\f$
 * @param[in] numTerms the length of list \p coeffs, which must be the same as that of \p exponents
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique
 * - if \p numQubits < 0 or \p numQubits >= `qureg.numQubitsRepresented` 
 * - if \p encoding is not a valid ::bitEncoding
 * - if \p encoding is not compatible with \p numQubits (e.g. \p TWOS_COMPLEMENT with only 1 qubit)
 * - if \p exponents contains a fractional number despite \p encoding <b>=</b> ::TWOS_COMPLEMENT (you must instead use applyPhaseFuncOverrides() and override all negative indices)
 * - if \p exponents contains a negative power (you must instead use applyPhaseFuncOverrides() and override the zero index)
 * - if \p numTerms <= 0
 * @author Tyson Jones
 */
void applyPhaseFunc(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms);

/** Induces a phase change upon each amplitude of \p qureg, determined by the passed 
 * exponential polynomial "phase function", and an explicit set of 'overriding' values at specific 
 * state indices.
 *
 * See applyPhaseFunc() first for a full description.
 *
 * - As in applyPhaseFunc(), the arguments \p coeffs and \p exponents specify a phase 
 *   function \f$f(r)\f$, where \f$r\f$ is determined by \p qubits and \p encoding for 
 *   each basis state of \p qureg.\n\n
 * - Additionally, \p overrideInds is a list of length \p numOverrides which specifies 
 *   the values of \f$r\f$ for which to explicitly set the induced phase change.\n
 *   The overriding phase changes are specified in the corresponding elements of \p overridePhases.\n\n
 *   For example, 
 *   ```
 *      int qubits[] = {0,1};
 *      enum bitEncoding encoding = UNSIGNED;
 *
 *      long long int overrideInds[] = {2};
 *      qreal overridePhases[] = {M_PI};
 *
 *      applyPhaseFuncOverrides(...);
 *   ```
 *   would effect the same diagonal unitary of applyPhaseFunc(), <em>except</em> that all 
 *   instance of \f$f(r=2)\f$ are overriden with phase \f$\pi\f$. \n I.e.
 *   \f[ 
 *   \begin{aligned}
 *     |0\mathbf{00}\rangle & \rightarrow \, e^{i f(0)}\,|0\mathbf{00}\rangle \\
 *     |0\mathbf{01}\rangle & \rightarrow \, e^{i f(1)}\,|0\mathbf{01}\rangle \\
 *     |0\mathbf{10}\rangle & \rightarrow \, e^{i \pi} \hspace{12pt} |0\mathbf{10}\rangle \\
 *     |0\mathbf{11}\rangle & \rightarrow \, e^{i f(3)}\,|0\mathbf{11}\rangle \\
 *     |1\mathbf{00}\rangle & \rightarrow \, e^{i f(0)}\,|1\mathbf{00}\rangle \\
 *     |1\mathbf{01}\rangle & \rightarrow \, e^{i f(1)}\,|1\mathbf{01}\rangle \\
 *     |1\mathbf{10}\rangle & \rightarrow \, e^{i \pi} \hspace{12pt} |1\mathbf{10}\rangle \\
 *     |1\mathbf{11}\rangle & \rightarrow \, e^{i f(3)}\,|1\mathbf{11}\rangle
 *   \end{aligned}
 *   \f]
 *   Note that if \p encoding <b>=</b> ::TWOS_COMPLEMENT, \a and \f$f(r)\f$ features a 
 *   fractional exponent, then every negative phase index must be overriden. This 
 *   is checked and enforced by QuEST's validation, \a unless there are more than 
 *   16 targeted qubits, in which case valid input is assumed (due to an otherwise 
 *   prohibitive performance overhead).
 *   \n
 * > Overriding phases are checked at each computational basis state of \p qureg <em>before</em>
 * > evaluating the phase function \f$f(r)\f$, and hence are useful for avoiding 
 * > singularities or errors at diverging values of \f$r\f$.
 *
 * - If \p qureg is a density matrix \f$\rho\f$, the overrides determine the diagonal unitary matrix 
 *   \f$\hat{D}\f$, which is then applied to \p qureg as
 *   \f[
 *      \rho \; \rightarrow \; \hat{D} \, \rho \hat{D}^\dagger.
 *   \f]
 *   This means that with overrides \f$f(r_j) \rightarrow \theta\f$ and \f$f(r_k) \rightarrow \phi\f$,
 *   element \f$\rho_{jk}\f$ is modified to
 *   \f[
 *      \alpha \, |j\rangle\langle k| \; \rightarrow \; 
 *          \exp(\, i \, (\theta - \phi) \, ) \; \alpha \, |j\rangle\langle k|.
 *   \f]\n
 *
 * - The interpreted phase function and list of overrides can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyPhaseFunc(qureg, ...);
 *   printRecordedQASM(qureg);
 *   ```
 *   may show 
 *   ```
 *   // Here, applyPhaseFunc() multiplied a complex scalar of the form
 *   //     exp(i (0.3 x^(-5) + 4 x^1 + 1 x^3))
 *   //   upon every substate |x>, informed by qubits (under a two's complement binary encoding)
 *   //     {4, 1, 2, 0}
 *   //   though with overrides
 *   //     |0> -> exp(i 3.14159)
 *   //     |1> -> exp(i (-3.14159))
 *   //     |2> -> exp(i 0)
 *   ```
 * \n
 *
 *
 * @see 
 * - applyPhaseFunc() for full doc on how \f$f(r)\f$ is evaluated.
 * - applyMultiVarPhaseFunc() for multi-variable exponential polynomial phase functions.
 * - applyNamedPhaseFunc() for a set of specific phase functions.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 * 
 * @ingroup operator
 * @param[in,out] qureg the state-vector or density matrix to be modified
 * @param[in] qubits a list of the indices of the qubits which will inform \f$r\f$ for each amplitude in \p qureg
 * @param[in] numQubits the length of list \p qubits
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r\f$ from the bits of \p qubits in each basis state of \p qureg
 * @param[in] coeffs the coefficients of the exponential polynomial phase function \f$f(r)\f$
 * @param[in] exponents the exponents of the exponential polynomial phase function \f$f(r)\f$
 * @param[in] numTerms the length of list \p coeffs, which must be the same as that of \p exponents
 * @param[in] overrideInds a list of sub-state indices (values of \f$r\f$) of which to explicit set the phase change
 * @param[in] overridePhases a list of replacement phase changes, for the corresponding \f$r\f$ values in \p overrideInds (one to one)
 * @param[in] numOverrides the lengths of lists \p overrideInds and \p overridePhases
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique
 * - if \p numQubits < 0 or \p numQubits >= `qureg.numQubitsRepresented` 
 * - if \p encoding is not a valid ::bitEncoding
 * - if \p encoding is not compatible with \p numQubits (i.e. \p TWOS_COMPLEMENT with 1 qubit)
 * - if \p numTerms <= 0
 * - if any value in \p overrideInds is not producible by \p qubits under the given \p encoding (e.g. 2 unsigned qubits cannot represent index 9)
 * - if \p numOverrides < 0
 * - if \p exponents contains a negative power and the (consequently diverging) zero index is not contained in \p overrideInds 
 * - if \p encoding is ::TWOS_COMPLEMENT, and \p exponents contains a fractional number, but \p overrideInds does not contain every possible negative index (checked only up to 16 targeted qubits)
 * @author Tyson Jones
 */
void applyPhaseFuncOverrides(Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int numTerms, long long int* overrideInds, qreal* overridePhases, int numOverrides);

/** Induces a phase change upon each amplitude of \p qureg, determined by a
 * multi-variable exponential polynomial "phase function". 
 *
 * This is a multi-variable extension of applyPhaseFunc(), whereby multiple sub-registers inform 
 * separate variables in the exponential polynomial function, and effects a diagonal unitary 
 * operator.
 *
 * - Arguments \p coeffs, \p exponents and \p numTermsPerReg together specify a real 
 *   exponential polynomial \f$f(\vec{r})\f$ of the form
 *   \f[ 
 *    f(r_1, \; \dots, \; r_{\text{numRegs}}) = \sum\limits_j^{\text{numRegs}} \; \sum\limits_{i}^{\text{numTermsPerReg}[j]} \; c_{i,j} \; {r_j}^{\; p_{i,j}}\,,
 *   \f] 
 *   where both coefficients \f$c_{i,j}\f$ and exponents \f$p_{i,j}\f$ can be any real number, subject to constraints described below.
 *   \n\n
 *   While \p coeffs and \p exponents are flat lists, they should be considered grouped into 
 *   #`numRegs` sublists with lengths given by \p numTermsPerReg (which itself has length \p numRegs). \n\n
 *   For example,
 *   ```
 *      int numRegs = 3;
 *      qreal coeffs[] =        {1,  2, 4,  -3.14};
 *      qreal exponents[] =     {2,  1, 5,   0.5 };
 *      int numTermsPerReg[] =  {1,  2,      1   };
 *   ```
 *   constitutes the function
 *   \f[
 *      f(\vec{r}) =  1 \, {r_1}^2 + 2 \, {r_2} + 4 \, {r_2}^{5} - 3.14 \, {r_3}^{0.5}.
 *   \f] \n
 *   > This means lists \p coeffs and \p exponents should both be of length equal to the sum of \p numTermsPerReg.
 *   Unlike applyPhaseFunc(), applyMultiVarPhaseFunc() places additional constraints on the 
 *   exponents in \f$f(\vec{r})\f$, due to the exponentially growing costs of overriding 
 *   diverging indices. Namely:\n
 *   -# \p exponents must not contain a negative number, since this would result in a divergence 
 *         when that register is zero, which would need to be overriden for every other register 
 *         basis state. If \f$f(\vec{r})\f$ must contain a negative exponent, you should instead 
 *         call applyPhaseFuncOverrides() once for each register/variable, and override the 
 *         zero index for the relevant variable. This works, because 
 *         \f[  \exp( i \sum_j f_j(r_j) ) = \prod_j \exp(i f_j(r_j) ). \f]
 *   -# \p exponents must not contain a fractional number if \p endoding <b>=</b> ::TWOS_COMPLEMENT, 
 *         because such a term would produce illegal complex values at negative register indices.
 *         Similar to the problem above, each negative register index would require overriding at 
 *         every index of the other registers, and hence require an exponential number of overrides.
 *         Therefore, if \f$f(\vec{r})\f$ must contain a negative exponent, you should instead 
 *         call applyPhaseFuncOverrides() once for each register/variable, and override every 
 *         negative index of each register in turn.
 * \n\n
 * - Lists \p qubits and \p numQubitsPerReg together describe #`numRegs` sub-registers of \p qureg,
 *   which can each contain a different number of qubits. \n
 *   Although \p qubits is a flat list of unique qubit indices, it should be imagined grouped into #`numRegs` sub-lists, 
 *   of lengths given by \p numQubitsPerReg. \n\n
 *   For example,
 *   ```
 *      int qubits[] =          {0,1,  3,4,5,  7}
 *      int numQubitsPerReg[] = {2,    3,      1};
 *      int numRegs = 3;
 *   ```
 *   describes three sub-registers, which are bolded below in an eight-qubit zero-state.
 *   \f[
 *      |r_3\rangle \; |0\rangle \; |r_2\rangle \; |0\rangle \; |r_1\rangle = 
 *      |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{00}\rangle
 *   \f]
 *   Note that the qubits need not be ordered increasing, and qubits within each sub-register 
 *   are assumed ordered least to most significant in that sub-register.\n\n
 *   > List \p qubits should have length equal to the sum of elements in \p numQubitsPerReg.
 *
 * - Each sub-register is associated with a variable \f$r_j\f$ in phase function \f$f(\vec{r})\f$. \n
 *   For a given computational basis state of \p qureg, the value of each variable is determined 
 *   by the binary value in the corresponding sub-register, when intepreted with ::bitEncoding \p encoding. \n
 *   See ::bitEncoding for more information.\n\n
 *
 * - The function \f$f(\vec{r})\f$ specifies the phase change to induce upon amplitude \f$\alpha\f$ 
 *   of computational basis state with the nominated sub-registers encoding values \f$r_1, \; \dots\f$.
 *   \f[
 *    \alpha \, |r_{\text{numRegs}}, \; \dots, \; r_2, \; r_1 \rangle \rightarrow \, \exp(i f(\vec{r}\,)) \; \alpha \, |r_{\text{numRegs}}, \; \dots, \; r_2, \; r_1 \rangle.
 *   \f]
 *   For example, using the sub-registers in the previous example and \p encoding <b>=</b> \p UNSIGNED, the 
 *   following states receive amplitude factors:
 *   \f[ 
 *   \begin{aligned}
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{00}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=0)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{01}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=1)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{10}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=2)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{11}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=3)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |1\rangle \; |\mathbf{00}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=0)} \\
 *   & \;\;\;\vdots \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{111}\rangle \; |0\rangle \; |\mathbf{01}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=7,r_1=1)} \\
 *   & \;\;\;\vdots \\
 *     |\mathbf{1}\rangle \; |0\rangle \; |\mathbf{111}\rangle \; |0\rangle \; |\mathbf{11}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=1,r_2=7,r_1=3)}
 *   \end{aligned}
 *   \f]
 *
 * - If \p qureg is a density matrix \f$\rho\f$, then its elements are modified as 
 *   \f[
 *      \alpha \, |j\rangle\langle k| \; \rightarrow \;
 *          \exp(i \, (f(\vec{r}_j) - f(\vec{r}_k)) \, ) \; \alpha \, |j\rangle\langle k|,
 *   \f]
 *   where \f$f(\vec{r}_j)\f$ and \f$f(\vec{r}_k)\f$ are determined as above.\n\n
 *
 * - The interpreted phase function can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyMultiVarPhaseFunc(qureg, ...);
 *   printRecordedQASM(qureg);
 *   ```
 *   would show, for the above example,
 *   ```
 *   // Here, applyMultiVarPhaseFunc() multiplied a complex scalar of the form
 *   //     exp(i (
 *   //          + 1 x^2
 *   //          + 2 y + 4 y^(-1)
 *   //          - 3.14 z^0.5 ))
 *   //   upon substates informed by qubits (under an unsigned binary encoding)
 *   //     |x> = {0, 1}
 *   //     |y> = {3, 4, 5}
 *   //     |z> = {7}
 *   ```
 * \n
 *
 *
 * @see 
 * - applyMultiVarPhaseFuncOverrides() to additionally specify explicit phases for specific sub-register values.
 * - applyNamedPhaseFunc() for a set of specific and potentially multi-variable phase functions.
 * - applyPhaseFunc() for a single-variable polynomial exponential phase function, which is approximately twice as fast.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 *
 * @ingroup operator
 * @param[in,out] qureg the state-vector or density matrix to be modified
 * @param[in] qubits a list of all the qubit indices contained in each sub-register
 * @param[in] numQubitsPerReg a list of the lengths of each sub-list in \p qubits
 * @param[in] numRegs the number of sub-registers, which is the length of both \p numQubitsPerReg and \p numTermsPerReg
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r_j\f$ from the bits of a sub-register
 * @param[in] coeffs the coefficients of all terms of the exponential polynomial phase function \f$f(\vec{r})\f$
 * @param[in] exponents the exponents of all terms of the exponential polynomial phase function \f$f(\vec{r})\f$
 * @param[in] numTermsPerReg a list of the number of \p coeff and \p exponent terms supplied for each variable/sub-register
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique (including if sub-registers overlap)
 * - if \p numRegs <= 0 or \p numRegs > 100 (constrained by `MAX_NUM_REGS_APPLY_ARBITRARY_PHASE` in QuEST_precision.h)
 * - if \p encoding is not a valid ::bitEncoding
 * - if the size of any sub-register is incompatible with \p encoding (e.g. contains fewer than two qubits in \p encoding <b>=</b> \p TWOS_COMPLEMENT)
 * - if any element of \p numTermsPerReg is < 1
 * - if \p exponents contains a negative number 
 * - if \p exponents contains a fractional number despite \p encoding <b>=</b> ::TWOS_COMPLEMENT
 * @author Tyson Jones
 */
void applyMultiVarPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg);

/** Induces a phase change upon each amplitude of \p qureg, determined by a
 * multi-variable exponential polynomial "phase function", and an explicit set of 'overriding' 
 * values at specific state indices.
 
 * See applyMultiVarPhaseFunc() first for a full description.
 *
 * - As in applyMultiVarPhaseFunc(), the arguments \p coeffs and \p exponents specify a 
 *   multi-variable phase function \f$f(\vec{r})\f$, where \f$\vec{r}\f$ is determined by 
 *   the sub-registers in \p qubits, and ::bitEncoding \p encoding for each basis state of \p qureg.\n\n
 * 
 * - Additionally, \p overrideInds is a list of length \p numOverrides which specifies 
 *   the values of \f$\vec{r}\f$ for which to explicitly set the induced phase change.\n
 *   While flat, \p overrideInds should be imagined grouped into sub-lists of length
 *   \p numRegs, which specify the full \f$\{r_1,\; \dots \;r_{\text{numRegs}} \} \f$ coordinate to override. \n
 *   Each sublist corresponds to a single element of \p overridePhases. \n
 *   For example,
 *   ```
 *   int numRegs = 3;
 *   int numOverrides = 2;
 *   long long int overrideInds[] = { 0,0,0,   1,2,3  };
 *   qreal overridePhases[]       = { M_PI,   - M_PI };
 *   ```
 *   denotes that any basis state of \p qureg with sub-register values \f$\{r_3,r_2,r_1\} = \{0, 0, 0\}\f$
 *   (or \f$\{r_3,r_2,r_1\} = \{1,2,3\}\f$) should receive phase change \f$\pi\f$ (or \f$-\pi\f$)
 *   in lieu of \f$\exp(i f(r_3=0,r_2=0,r_1=0))\f$.\n\n
 *   > Note that you cannot use applyMultiVarPhaseFuncOverrides() to override divergences 
 *   > in \f$f(\vec{r})\f$, since each diverging value \f$r_j\f$ would need to be overriden 
 *   > as an \f$\vec{r}\f$ coordinate for every basis state of the other registers; the number 
 *   > of overrides grows exponentially. Ergo, if \p exponents contains a negative number 
 *   > (diverging at \f$r_j=0\f$), or \p exponents contains a fractional number despite 
 *   > \p encoding <b>=</b> ::TWOS_COMPLEMENT (producing complex phases at negative indices),
 *   > you must instead call applyPhaseFuncOverrides() for each variable in turn and 
 *   > override the diverging \f$r_j\f$ (each independently of the other registers).
 *
 * - The interpreted overrides can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyMultiVarPhaseFuncOverrides(qureg, ...);
 *   printRecordedQASM(qureg);
 *   ```
 *   may show
 *   ```
 *   // Here, applyMultiVarPhaseFunc() multiplied ...
 *   //   though with overrides
 *   //     |x=0, y=0, z=0> -> exp(i 3.14159)
 *   //     |x=1, y=2, z=3> -> exp(i (-3.14159))
 *   ```
 * \n
 *
 *
 * @see 
 * - applyNamedPhaseFunc() for a set of specific and potentially multi-variable phase functions.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 *
 * @ingroup operator
 * @param[in,out] qureg the state-vector or density-matrix to be modified
 * @param[in] qubits a list of all the qubit indices contained in each sub-register
 * @param[in] numQubitsPerReg a list of the lengths of each sub-list in \p qubits
 * @param[in] numRegs the number of sub-registers, which is the length of both \p numQubitsPerReg and \p numTermsPerReg
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r_j\f$ from the bits of a sub-register
 * @param[in] coeffs the coefficients of all terms of the exponential polynomial phase function \f$f(\vec{r})\f$
 * @param[in] exponents the exponents of all terms of the exponential polynomial phase function \f$f(\vec{r})\f$
 * @param[in] numTermsPerReg a list of the number of \p coeff and \p exponent terms supplied for each variable/sub-register
 * @param[in] overrideInds a flattened list of sub-register coordinates (values of \f$\vec{r}\f$) of which to explicit set the phase change
 * @param[in] overridePhases a list of replacement phase changes, for the corresponding \f$\vec{r}\f$ values in \p overrideInds
 * @param[in] numOverrides the lengths of list \p overridePhases (but not necessarily of \p overrideInds)
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique (including if sub-registers overlap)
 * - if \p numRegs <= 0 or \p numRegs > 100 (constrained by `MAX_NUM_REGS_APPLY_ARBITRARY_PHASE` in QuEST_precision.h)
 * - if \p encoding is not a valid ::bitEncoding
 * - if the size of any sub-register is incompatible with \p encoding (e.g. contains fewer than two qubits in \p encoding <b>=</b> \p TWOS_COMPLEMENT)
 * - if any element of \p numTermsPerReg is < 1
 * - if \p exponents contains a negative number 
 * - if \p exponents contains a fractional number despite \p encoding <b>=</b> ::TWOS_COMPLEMENT
 * - if any value in \p overrideInds is not producible by its corresponding sub-register under the given \p encoding (e.g. 2 unsigned qubits cannot represent index 9)
 * - if \p numOverrides < 0
 * @author Tyson Jones
 */
void applyMultiVarPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, qreal* coeffs, qreal* exponents, int* numTermsPerReg, long long int* overrideInds, qreal* overridePhases, int numOverrides);

/** Induces a phase change upon each amplitude of \p qureg, determined by a
 * named (and potentially multi-variable) phase function.
 *
 * This effects a diagonal unitary operator, with a phase function \f$f(\vec{r})\f$ which may not be 
 * simply expressible as an exponential polynomial in functions applyPhaseFunc() and applyMultiVarPhaseFunc().
 *
 * Arguments \p qubits and \p numQubitsPerReg encode sub-registers of \p qureg in the same 
 * manner as in applyMultiVarPhaseFunc():
 * - Lists \p qubits and \p numQubitsPerReg together describe #`numRegs` sub-registers of \p qureg,
 *   which can each contain a different number of qubits. \n
 *   Although \p qubits is a flat list of unique qubit indices, it should be imagined grouped into #`numRegs` sub-lists, 
 *   of lengths given by \p numQubitsPerReg. \n\n
 *   For example,
 *   ```
 *      int qubits[] =          {0,1,  3,4,5,  7}
 *      int numQubitsPerReg[] = {2,    3,      1};
 *      int numRegs = 3;
 *   ```
 *   describes three sub-registers, which are bolded below in an eight-qubit zero-state.
 *   \f[
 *      |r_3\rangle \; |0\rangle \; |r_2\rangle \; |0\rangle \; |r_1\rangle = 
 *      |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{00}\rangle
 *   \f]
 *   Note that the qubits need not be ordered increasing, and qubits within each sub-register 
 *   are assumed ordered least to most significant in that sub-register.\n\n
 *   > List \p qubits should have length equal to the sum of elements in \p numQubitsPerReg.
 *
 * - Each sub-register is associated with a variable \f$r_j\f$ in phase function \f$f(\vec{r})\f$. \n
 *   For a given computational basis state of \p qureg, the value of each variable is determined 
 *   by the binary value in the corresponding sub-register, when intepreted with ::bitEncoding \p encoding. \n
 *   See ::bitEncoding for more information.\n\n
 *
 * - Argument \p functionNameCode determines the phase function \f$f(\vec{r})\f$.\n
 *   For example, 
 *   ```
 *   int numRegs = 3;
 *   enum phaseFunc functionNameCode = NORM;
 *   ```
 *   describes phase function 
 *   \f[
 *      f(\vec{r}) = \sqrt{ {r_1}^2 + {r_2}^2 + {r_3} ^2 }.
 *   \f]
 *   See ::phaseFunc for a list and description of all named phase functions. \n
 *   Some phase functions, like \p SCALED_NORM, require passing additional parameters, through 
 *   the function applyParamNamedPhaseFunc().\n\n
 *   > If the phase function \f$f(\vec{r})\f$ diverges at one or more \f$\vec{r}\f$ values, you should instead 
 *   > use applyNamedPhaseFuncOverrides() and specify explicit phase changes for these coordinates.
 *   > Otherwise, the corresponding amplitudes of \p qureg will become indeterminate (like `NaN`). \n
 *
 * - The function \f$f(\vec{r})\f$ specifies the phase change to induce upon amplitude \f$\alpha\f$ 
 *   of computational basis state with the nominated sub-registers encoding values \f$r_1, \; \dots\f$.
 *   \f[
 *    \alpha \, |r_{\text{numRegs}}, \; \dots, \; r_2, \; r_1 \rangle \rightarrow \, \exp(i f(\vec{r}\,)) \; \alpha \, |r_{\text{numRegs}}, \; \dots, \; r_2, \; r_1 \rangle.
 *   \f]
 *   For example, using the sub-registers in the above example and \p encoding <b>=</b> \p UNSIGNED, the 
 *   following states receive amplitude factors:
 *   \f[ 
 *   \begin{aligned}
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{00}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=0)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{01}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=1)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{10}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=2)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |0\rangle \; |\mathbf{11}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=3)} \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{000}\rangle \; |1\rangle \; |\mathbf{00}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=0,r_1=0)} \\
 *   & \;\;\;\vdots \\
 *     |\mathbf{0}\rangle \; |0\rangle \; |\mathbf{111}\rangle \; |0\rangle \; |\mathbf{01}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=0,r_2=7,r_1=1)} \\
 *   & \;\;\;\vdots \\
 *     |\mathbf{1}\rangle \; |0\rangle \; |\mathbf{111}\rangle \; |0\rangle \; |\mathbf{11}\rangle & 
 *        \rightarrow \, 
 *            e^{i f(r_3=1,r_2=7,r_1=3)}
 *   \end{aligned}
 *   \f]\n
 *
 * - If \p qureg is a density matrix, its elements are modified to
 *   \f[
        \alpha \, |j\rangle\langle k| \; \rightarrow \;
            \exp(i (f(\vec{r}_j) \, - \, f(\vec{r}_k))) \; \alpha \, |j\rangle\langle k|
 *   \f]
 *   where \f$f(\vec{r}_j)\f$ and \f$f(\vec{r}_k)\f$ are determined as above. This is equivalent 
 *   to modification
 *   \f[
 *          \rho \; \rightarrow \; \hat{D} \, \rho \, \hat{D}^\dagger
 *   \f]
 *   where \f$\hat{D}\f$ is the diagonal unitary 
 *   \f[
 *      \hat{D} = \text{diag}\, \{ \; e^{i f(\vec{r_0})}, \; e^{i f(\vec{r_1})}, \; \dots \; \}.
 *   \f]\n
 *
 * - The interpreted phase function can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyNamedPhaseFunc(qureg, ..., INVERSE_DISTANCE, ... );
 *   printRecordedQASM(qureg);
 *   ```
 *   may show
 *   ```
 *   // Here, applyNamedPhaseFunc() multiplied a complex scalar of form
 *   //     exp(i 1 / sqrt((x-y)^2 + (z-t)^2))
 *   ```
 * \n
 *
 *
 * @see 
 * - applyNamedPhaseFuncOverrides() to additionally specify phase values for specific sub-register indices.
 * - applyParamNamedPhaseFunc() to specify named phase functions which require additional parameters.
 * - applyPhaseFunc() to specify a general single-variable exponential polynomial phase function.
 * - applyMultiVarPhaseFunc() to specify a general multi-variable exponential polynomial phase function.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 * 
 * @ingroup operator
 * @param[in,out] qureg the state-vector or density-matrix to be modified
 * @param[in] qubits a list of all the qubit indices contained in each sub-register
 * @param[in] numQubitsPerReg a list of the lengths of each sub-list in \p qubits
 * @param[in] numRegs the number of sub-registers, which is the length of both \p numQubitsPerReg and \p numTermsPerReg
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r_j\f$ from the bits of a sub-register
 * @param[in] functionNameCode the ::phaseFunc \f$f(\vec{r})\f$
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique (including if sub-registers overlap)
 * - if \p numRegs <= 0 or \p numRegs > 100 (constrained by `MAX_NUM_REGS_APPLY_ARBITRARY_PHASE` in QuEST_precision.h)
 * - if \p encoding is not a valid ::bitEncoding
 * - if the size of any sub-register is incompatible with \p encoding (e.g. contains fewer than two qubits in \p encoding <b>=</b> \p TWOS_COMPLEMENT)
 * - if \p functionNameCode is not a valid ::phaseFunc
 * - if \p functionNameCode requires additional parameters, which must instead be passed with applyParamNamedPhaseFunc()
 * @author Tyson Jones
 */
void applyNamedPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode);

/** Induces a phase change upon each amplitude of \p qureg, determined by a
 * named (and potentially multi-variable) phase function, and an explicit set of 'overriding' 
 * values at specific state indices.
 *
 * See applyNamedPhaseFunc() first for a full description.
 *
 * - As in applyNamedPhaseFunc(), \p functionNameCode specifies a
 *   multi-variable phase function \f$f(\vec{r})\f$, where \f$\vec{r}\f$ is determined by 
 *   the sub-registers in \p qubits, and ::bitEncoding \p encoding for each basis state of \p qureg.\n\n
 * 
 * - Additionally, \p overrideInds is a list of length \p numOverrides which specifies 
 *   the values of \f$\vec{r}\f$ for which to explicitly set the induced phase change.\n
 *   While flat, \p overrideInds should be imagined grouped into sub-lists of length
 *   \p numRegs, which specify the full \f$\{r_1,\; \dots \;r_{\text{numRegs}} \} \f$ coordinate to override. \n
 *   Each sublist corresponds to a single element of \p overridePhases. \n
 *   For example,
 *   ```
 *   int numRegs = 3;
 *   int numOverrides = 2;
 *   long long int overrideInds[] = { 0,0,0,   1,2,3  };
 *   qreal overridePhases[]       = { M_PI,   - M_PI };
 *   ```
 *   denotes that any basis state of \p qureg with sub-register values \f$\{r_3,r_2,r_1\} = \{0, 0, 0\}\f$
 *   (or \f$\{r_3,r_2,r_1\} = \{1,2,3\}\f$) should receive phase change \f$\pi\f$ (or \f$-\pi\f$)
 *   in lieu of \f$\exp(i f(r_3=0,r_2=0,r_1=0))\f$.\n\n
 *
 * - The interpreted overrides can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyNamedPhaseFuncOverrides(qureg, ...);
 *   printRecordedQASM(qureg);
 *   ```
 *   may show
 *   ```
 *   // Here, applyNamedPhaseFunc() multiplied ...
 *   //   though with overrides
 *   //     |x=0, y=0, z=0> -> exp(i 3.14159)
 *   //     |x=1, y=2, z=3> -> exp(i (-3.14159))
 *   ```
 * \n
 *
 *
 * @see 
 * - applyParamNamedPhaseFuncOverrides() to specify <em>parameterised</em> named phase functions, with phase overrides.
 * - applyPhaseFunc() to specify a general single-variable exponential polynomial phase function.
 * - applyMultiVarPhaseFunc() to specify a general multi-variable exponential polynomial phase function.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 * 
 * @ingroup operator
 * @param[in,out] qureg the state-vector pr density-matrix to be modified
 * @param[in] qubits a list of all the qubit indices contained in each sub-register
 * @param[in] numQubitsPerReg a list of the lengths of each sub-list in \p qubits
 * @param[in] numRegs the number of sub-registers, which is the length of both \p numQubitsPerReg and \p numTermsPerReg
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r_j\f$ from the bits of a sub-register
 * @param[in] functionNameCode the ::phaseFunc \f$f(\vec{r})\f$
 * @param[in] overrideInds a flattened list of sub-register coordinates (values of \f$\vec{r}\f$) of which to explicit set the phase change
 * @param[in] overridePhases a list of replacement phase changes, for the corresponding \f$\vec{r}\f$ values in \p overrideInds
 * @param[in] numOverrides the lengths of list \p overridePhases (but not necessarily of \p overrideInds)
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique (including if sub-registers overlap)
 * - if \p numRegs <= 0 or \p numRegs > 100 (constrained by `MAX_NUM_REGS_APPLY_ARBITRARY_PHASE` in QuEST_precision.h)
 * - if \p encoding is not a valid ::bitEncoding
 * - if the size of any sub-register is incompatible with \p encoding (e.g. contains fewer than two qubits in \p encoding <b>=</b> \p TWOS_COMPLEMENT)
 * - if \p functionNameCode is not a valid ::phaseFunc
 * - if \p functionNameCode requires additional parameters, which must instead be passed with applyParamNamedPhaseFunc()
 * - if any value in \p overrideInds is not producible by its corresponding sub-register under the given \p encoding (e.g. 2 unsigned qubits cannot represent index 9)
 * - if \p numOverrides < 0
 * @author Tyson Jones
 */
void applyNamedPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, long long int* overrideInds, qreal* overridePhases, int numOverrides);

/** Induces a phase change upon each amplitude of \p qureg, determined by a
 * named, paramaterized (and potentially multi-variable) phase function.
 *
 * See applyNamedPhaseFunc() for full documentation. \n
 * This function merely accepts additional ::phaseFunc names which accept one (or more) parameters.
 *
 * - Argument \p functionNameCode, which determines the phase function \f$f(\vec{r}, \vec{\theta})\f$,
 *   can include parameterised ::phaseFunc names like \p SCALED_NORM, which require additional 
 *   parameters \f$\vec{\theta}\f$ passed via list \p params.\n
 *   For example, 
 *   ```
 *   enum phaseFunc functionNameCode = SCALED_PRODUCT;
 *   qreal params[] = {0.5};
 *   int numParams = 1;
 *   applyParamNamedPhaseFunc(..., functionNameCode, params, numParams);
 *   ```
 *   invokes phase function 
 *   \f[
 *      f(\vec{r}, \theta)|_{\theta=0.5} \; = \; 0.5 \prod_j^{\text{numRegs}} \; r_j\,.
 *   \f] 
 *   See ::phaseFunc for all named phased functions.
 *
 * - Functions with divergences, like \p INVERSE_NORM and \p SCALED_INVERSE_DISTANCE, must accompany 
 *   an extra parameter to specify an overriding phase at the divergence. For example,
 *   ```
 *   enum phaseFunc functionNameCode = SCALED_INVERSE_NORM;
 *   qreal params[] = {0.5, M_PI};
 *   int numParams = 2;
 *   applyParamNamedPhaseFunc(..., functionNameCode, params, numParams);
 *   ```
 *   invokes phase function 
 *   \f[
 *      f(\vec{r}, \theta)|_{\theta=0.5} \; = \; \begin{cases} \pi & \;\;\; \vec{r}=\vec{0} \\ \displaystyle 0.5 \left[ \sum_j^{\text{numRegs}} {r_j}^2 \right]^{-1/2} & \;\;\;\text{otherwise} \end{cases}.
 *   \f] 
 *   Notice the order of the parameters matches the order of the words in the \p phaseFunc.
 *   > Functions \p SCALED_INVERSE_SHIFTED_NORM and \p SCALED_INVERSE_SHIFTED_DISTANCE,
 *   > which can have denominators arbitrarily close to zero, will invoke the 
 *   > divergence parameter whenever the denominator is smaller than (or equal to)
 *   > machine precision `REAL_EPS`.
 *
 * - Functions allowing the shifting of sub-register values, which are \p SCALED_INVERSE_SHIFTED_NORM
 *   and \p SCALED_INVERSE_SHIFTED_DISTANCE, need these shift values to be passed in the \p params
 *   argument _after_ the scaling and divergence override parameters listed above. The function
 *   \p SCALED_INVERSE_SHIFTED_NORM needs as many extra parameters, as there are sub-registers;
 *   \p SCALED_INVERSE_SHIFTED_DISTANCE needs one extra parameter for each pair of sub-registers.
 *   For example,
 *   ```
 *   enum phaseFunc functionNameCode = SCALED_INVERSE_SHIFTED_NORM;
 *   int qubits[] = {0,1,2,3, 4,5,6,7};
 *   int qubitsPerReg[] = {4, 4};
 *   qreal params[] = {0.5, M_PI, 0.8, -0.3};
 *   int numParams = 4;
 *   applyParamNamedPhaseFunc(..., qubits, qubitsPerReg, 2, ..., functionNameCode, params, numParams);
 *   ```
 *   invokes phase function 
 *   \f[
 *      f(\vec{r}) \; = \; \begin{cases} \pi & \;\;\; \vec{r}=\vec{0} \\ \displaystyle 0.5 \left[(r_1-0.8)^2 + (r_2+0.3)^2\right]^{-1/2} & \;\;\;\text{otherwise} \end{cases}.
 *   \f] 
 *   and
 *   ```
 *   enum phaseFunc functionNameCode = SCALED_INVERSE_SHIFTED_DISTANCE;
 *   int qubits[] = {0,1, 2,3, 4,5, 6,7};
 *   int qubitsPerReg[] = {2, 2, 2, 2};
 *   qreal params[] = {0.5, M_PI, 0.8, -0.3};
 *   int numParams = 4;
 *   applyParamNamedPhaseFunc(..., qubits, qubitsPerReg, 4, ..., functionNameCode, params, numParams);
 *   ```
 *   invokes phase function 
 *   \f[
 *      f(\vec{r}) \; = \; \begin{cases} \pi & \;\;\; \vec{r}=\vec{0} \\ \displaystyle 0.5 \left[(r_1-r_2-0.8)^2 + (r_3-r_4+0.3)^2\right]^{-1/2} & \;\;\;\text{otherwise} \end{cases}.
 *   \f] 
 * 
 *   > You can further override \f$f(\vec{r}, \vec{\theta})\f$ at one or more \f$\vec{r}\f$ values
 *   > via applyParamNamedPhaseFuncOverrides().
 *
 * - The interpreted parameterised phase function can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyParamNamedPhaseFunc(...);
 *   printRecordedQASM(qureg);
 *   ```
 *   may show
 *   ```
 *   // Here, applyNamedPhaseFunc() multiplied a complex scalar of form
 *   //     exp(i (-0.5) / (x y z))
 *   ```
 * \n
 *
 *
 * @see 
 * - applyParamNamedPhaseFuncOverrides() to additionally specify phase values for specific sub-register indices.
 * - applyPhaseFunc() to specify a general single-variable exponential polynomial phase function.
 * - applyMultiVarPhaseFunc() to specify a general multi-variable exponential polynomial phase function.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 * 
 * @ingroup operator
 * @param[in,out] qureg the state-vector or density-matrix to be modified
 * @param[in] qubits a list of all the qubit indices contained in each sub-register
 * @param[in] numQubitsPerReg a list of the lengths of each sub-list in \p qubits
 * @param[in] numRegs the number of sub-registers, which is the length of both \p numQubitsPerReg and \p numTermsPerReg
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r_j\f$ from the bits of a sub-register
 * @param[in] functionNameCode the ::phaseFunc \f$f(\vec{r}, \vec{\theta})\f$
 * @param[in] params a list of any additional parameters needed by the ::phaseFunc \p functionNameCode
 * @param[in] numParams the length of list \p params
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique (including if sub-registers overlap)
 * - if \p numRegs <= 0 or \p numRegs > 100 (constrained by `MAX_NUM_REGS_APPLY_ARBITRARY_PHASE` in QuEST_precision.h)
 * - if \p encoding is not a valid ::bitEncoding
 * - if the size of any sub-register is incompatible with \p encoding (e.g. contains fewer than two qubits in \p encoding <b>=</b> \p TWOS_COMPLEMENT)
 * - if \p functionNameCode is not a valid ::phaseFunc
 * - if \p numParams is incompatible with \p functionNameCode (for example, no parameters were passed to \p SCALED_PRODUCT)
 * @author Tyson Jones
 * @author Richard Meister (shifted functions)
 */
void applyParamNamedPhaseFunc(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, qreal* params, int numParams);

/** Induces a phase change upon each amplitude of \p qureg, determined by a
 * named, parameterised (and potentially multi-variable) phase function, and an explicit set of 'overriding' 
 * values at specific state indices.
 *
 * See applyParamNamedPhaseFunc() and applyNamedPhaseFunc() first for a full description.
 *
 * - As in applyParamNamedPhaseFunc(), \p functionNameCode specifies a parameterised
 *   multi-variable phase function \f$f(\vec{r}, \vec{\theta})\f$, where \f$\vec{\theta}\f$ is
 *   passed in list \p params, and \f$\vec{r}\f$ is determined both by 
 *   the sub-registers in \p qubits, and ::bitEncoding \p encoding for each basis state of \p qureg.\n\n
 * 
 * - Additionally, \p overrideInds is a list of length \p numOverrides which specifies 
 *   the values of \f$\vec{r}\f$ for which to explicitly set the induced phase change.\n
 *   While flat, \p overrideInds should be imagined grouped into sub-lists of length
 *   \p numRegs, which specify the full \f$\{r_1,\; \dots \;r_{\text{numRegs}} \} \f$ coordinate to override. \n
 *   Each sublist corresponds to a single element of \p overridePhases. \n
 *   For example,
 *   ```
 *   int numRegs = 3;
 *   int numOverrides = 2;
 *   long long int overrideInds[] = { 0,0,0,   1,2,3  };
 *   qreal overridePhases[]       = { M_PI,   - M_PI };
 *   ```
 *   denotes that any basis state of \p qureg with sub-register values \f$\{r_3,r_2,r_1\} = \{0, 0, 0\}\f$
 *   (or \f$\{r_3,r_2,r_1\} = \{1,2,3\}\f$) should receive phase change \f$\pi\f$ (or \f$-\pi\f$)
 *   in lieu of \f$\exp(i f(r_3=0,r_2=0,r_1=0, \vec{\theta}))\f$.\n\n
 *
 * - The interpreted overrides can be previewed in the QASM log, as a comment. \n
 *   For example:
 *   ```
 *   startRecordingQASM(qureg);
 *   applyParamNamedPhaseFuncOverrides(qureg, ...);
 *   printRecordedQASM(qureg);
 *   ```
 *   may show
 *   ```
 *   // Here, applyParamNamedPhaseFunc() multiplied ...
 *   //   though with overrides
 *   //     |x=0, y=0, z=0> -> exp(i 3.14159)
 *   //     |x=1, y=2, z=3> -> exp(i (-3.14159))
 *   ```
 * \n
 *
 *
 * @see 
 * - applyPhaseFunc() to specify a general single-variable exponential polynomial phase function.
 * - applyMultiVarPhaseFunc() to specify a general multi-variable exponential polynomial phase function.
 * - applyDiagonalOp() to apply a non-unitary diagonal operator.
 * 
 * @ingroup operator
 * @param[in,out] qureg the state-vector or density-matrix to be modified
 * @param[in] qubits a list of all the qubit indices contained in each sub-register
 * @param[in] numQubitsPerReg a list of the lengths of each sub-list in \p qubits
 * @param[in] numRegs the number of sub-registers, which is the length of both \p numQubitsPerReg and \p numTermsPerReg
 * @param[in] encoding the ::bitEncoding under which to infer the binary value \f$r_j\f$ from the bits of a sub-register
 * @param[in] functionNameCode the ::phaseFunc \f$f(\vec{r}, \vec{\theta})\f$
 * @param[in] params a list of any additional parameters \f$\vec{\theta}\f$ needed by the ::phaseFunc \p functionNameCode
 * @param[in] numParams the length of list \p params
 * @param[in] overrideInds a flattened list of sub-register coordinates (values of \f$\vec{r}\f$) of which to explicit set the phase change
 * @param[in] overridePhases a list of replacement phase changes, for the corresponding \f$\vec{r}\f$ values in \p overrideInds
 * @param[in] numOverrides the lengths of list \p overridePhases (but not necessarily of \p overrideInds)
 * @exception invalidQuESTInputError()
 * - if any qubit in \p qubits has an invalid index (i.e. does not satisfy 0 <= qubit < `qureg.numQubitsRepresented`)
 * - if the elements of \p qubits are not unique (including if sub-registers overlap)
 * - if \p numRegs <= 0 or \p numRegs > 100 (constrained by `MAX_NUM_REGS_APPLY_ARBITRARY_PHASE` in QuEST_precision.h)
 * - if \p encoding is not a valid ::bitEncoding
 * - if the size of any sub-register is incompatible with \p encoding (e.g. contains fewer than two qubits in \p encoding <b>=</b> \p TWOS_COMPLEMENT)
 * - if \p functionNameCode is not a valid ::phaseFunc
 * - if \p numParams is incompatible with \p functionNameCode (for example, no parameters were passed to \p SCALED_PRODUCT)
 * - if any value in \p overrideInds is not producible by its corresponding sub-register under the given \p encoding (e.g. 2 unsigned qubits cannot represent index 9)
 * - if \p numOverrides < 0
 * @author Tyson Jones
 */
void applyParamNamedPhaseFuncOverrides(Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding, enum phaseFunc functionNameCode, qreal* params, int numParams, long long int* overrideInds, qreal* overridePhases, int numOverrides);

/** Applies the quantum Fourier transform (QFT) to the entirety of \p qureg. 
 * The effected unitary circuit (shown here for 4 qubits, bottom qubit is <b>0</b>) resembles
 \f[
             \begin{tikzpicture}[scale=.5]
             \draw (-2, 5) -- (23, 5);
             \draw (-2, 3) -- (23, 3);
             \draw (-2, 1) -- (23, 1);
             \draw (-2, -1) -- (23, -1);

             \draw[fill=white] (-1, 4) -- (-1, 6) -- (1, 6) -- (1,4) -- cycle;
             \node[draw=none] at (0, 5) {H};

             \draw(2, 5) -- (2, 3);
             \draw[fill=black] (2, 5) circle (.2);
             \draw[fill=black] (2, 3) circle (.2);
             \draw(4, 5) -- (4, 1);
             \draw[fill=black] (4, 5) circle (.2);
             \draw[fill=black] (4, 1) circle (.2);
             \draw(6, 5) -- (6, -1);
             \draw[fill=black] (6, 5) circle (.2);
             \draw[fill=black] (6, -1) circle (.2);

             \draw[fill=white] (-1+8, 4-2) -- (-1+8, 6-2) -- (1+8, 6-2) -- (1+8,4-2) -- cycle;
             \node[draw=none] at (8, 5-2) {H};

             \draw(10, 5-2) -- (10, 3-2);
             \draw[fill=black] (10, 5-2) circle (.2);
             \draw[fill=black] (10, 3-2) circle (.2);
             \draw(12, 5-2) -- (12, 3-4);
             \draw[fill=black] (12, 5-2) circle (.2);
             \draw[fill=black] (12, 3-4) circle (.2);

             \draw[fill=white] (-1+8+6, 4-4) -- (-1+8+6, 6-4) -- (1+8+6, 6-4) -- (1+8+6,4-4) -- cycle;
             \node[draw=none] at (8+6, 5-4) {H};

             \draw(16, 5-2-2) -- (16, 3-4);
             \draw[fill=black] (16, 5-2-2) circle (.2);
             \draw[fill=black] (16, 3-4) circle (.2);

             \draw[fill=white] (-1+8+6+4, 4-4-2) -- (-1+8+6+4, 6-4-2) -- (1+8+6+4, 6-4-2) -- (1+8+6+4,4-4-2) -- cycle;
             \node[draw=none] at (8+6+4, 5-4-2) {H};

             \draw (20, 5) -- (20, -1);
             \draw (20 - .35, 5 + .35) -- (20 + .35, 5 - .35);
             \draw (20 - .35, 5 - .35) -- (20 + .35, 5 + .35);
             \draw (20 - .35, -1 + .35) -- (20 + .35, -1 - .35);
             \draw (20 - .35, -1 - .35) -- (20 + .35, -1 + .35);
             \draw (22, 3) -- (22, 1);
             \draw (22 - .35, 3 + .35) -- (22 + .35, 3 - .35);
             \draw (22 - .35, 3 - .35) -- (22 + .35, 3 + .35);
             \draw (22 - .35, 1 + .35) -- (22 + .35, 1 - .35);
             \draw (22 - .35, 1 - .35) -- (22 + .35, 1 + .35);
             \end{tikzpicture}
 \f]
 * though is performed more efficiently.
 * 
 * - If \p qureg is a state-vector, the output amplitudes are the discrete Fourier 
 *   transform (DFT) of the input amplitudes, in the exact ordering. This is true 
 *   even if \p qureg is unnormalised. \n
 *   Precisely,
 *   \f[
 *      \text{QFT} \, \left(  \sum\limits_{x=0}^{2^N-1} \alpha_x |x\rangle \right) 
 *      = 
 *      \frac{1}{\sqrt{2^N}}
 *       \sum\limits_{x=0}^{2^N-1} \left( 
 *           \sum\limits_{y=0}^{2^N-1} e^{2 \pi \, i \, x \, y / 2^N} \; \alpha_y
 *       \right) |x\rangle
 *   \f]
 *
 * - If \p qureg is a density matrix \f$\rho\f$, it will be changed under the unitary action 
 *   of the QFT. This can be imagined as each mixed state-vector undergoing the DFT 
 *   on its amplitudes. This is true even if \p qureg is unnormalised.
 *   \f[
 *      \rho \; \rightarrow \; \text{QFT} \; \rho \; \text{QFT}^{\dagger}
 *   \f]
 *
 * > This function merges contiguous controlled-phase gates into single invocations 
 * > of applyNamedPhaseFunc(), and hence is significantly faster than performing 
 * > the QFT circuit directly.
 *
 * > Furthermore, in distributed mode, this function requires only \f$\log_2(\text{\#nodes})\f$ 
 * > rounds of pair-wise
 * > communication, and hence is exponentially faster than directly performing the 
 * > DFT on the amplitudes of \p qureg.
 *
 * @see 
 * - applyQFT() to apply the QFT to a sub-register of \p qureg.
 * 
 * @ingroup operator
 * @param[in,out] qureg a state-vector or density matrix to modify
 * @author Tyson Jones
 */
void applyFullQFT(Qureg qureg);

/** Applies the quantum Fourier transform (QFT) to a specific subset of qubits 
 * of the register \p qureg. 
 *
 * The order of \p qubits affects the ultimate unitary.
 * The canonical full-state QFT (applyFullQFT()) is achieved by targeting every 
 * qubit in increasing order.
 *
 * The effected unitary circuit (shown here for \p numQubits <b>= 4</b>) resembles
 * \f[
        \begin{tikzpicture}[scale=.5]
        \draw (-2, 5) -- (23, 5);    \node[draw=none] at (-4,5) {qubits[3]};
        \draw (-2, 3) -- (23, 3);    \node[draw=none] at (-4,3) {qubits[2]};
        \draw (-2, 1) -- (23, 1);     \node[draw=none] at (-4,1) {qubits[1]};
        \draw (-2, -1) -- (23, -1);  \node[draw=none] at (-4,-1) {qubits[0]};

        \draw[fill=white] (-1, 4) -- (-1, 6) -- (1, 6) -- (1,4) -- cycle;
        \node[draw=none] at (0, 5) {H};

        \draw(2, 5) -- (2, 3);
        \draw[fill=black] (2, 5) circle (.2);
        \draw[fill=black] (2, 3) circle (.2);
        \draw(4, 5) -- (4, 1);
        \draw[fill=black] (4, 5) circle (.2);
        \draw[fill=black] (4, 1) circle (.2);
        \draw(6, 5) -- (6, -1);
        \draw[fill=black] (6, 5) circle (.2);
        \draw[fill=black] (6, -1) circle (.2);

        \draw[fill=white] (-1+8, 4-2) -- (-1+8, 6-2) -- (1+8, 6-2) -- (1+8,4-2) -- cycle;
        \node[draw=none] at (8, 5-2) {H};

        \draw(10, 5-2) -- (10, 3-2);
        \draw[fill=black] (10, 5-2) circle (.2);
        \draw[fill=black] (10, 3-2) circle (.2);
        \draw(12, 5-2) -- (12, 3-4);
        \draw[fill=black] (12, 5-2) circle (.2);
        \draw[fill=black] (12, 3-4) circle (.2);

        \draw[fill=white] (-1+8+6, 4-4) -- (-1+8+6, 6-4) -- (1+8+6, 6-4) -- (1+8+6,4-4) -- cycle;
        \node[draw=none] at (8+6, 5-4) {H};

        \draw(16, 5-2-2) -- (16, 3-4);
        \draw[fill=black] (16, 5-2-2) circle (.2);
        \draw[fill=black] (16, 3-4) circle (.2);

        \draw[fill=white] (-1+8+6+4, 4-4-2) -- (-1+8+6+4, 6-4-2) -- (1+8+6+4, 6-4-2) -- (1+8+6+4,4-4-2) -- cycle;
        \node[draw=none] at (8+6+4, 5-4-2) {H};

        \draw (20, 5) -- (20, -1);
        \draw (20 - .35, 5 + .35) -- (20 + .35, 5 - .35);
        \draw (20 - .35, 5 - .35) -- (20 + .35, 5 + .35);
        \draw (20 - .35, -1 + .35) -- (20 + .35, -1 - .35);
        \draw (20 - .35, -1 - .35) -- (20 + .35, -1 + .35);
        \draw (22, 3) -- (22, 1);
        \draw (22 - .35, 3 + .35) -- (22 + .35, 3 - .35);
        \draw (22 - .35, 3 - .35) -- (22 + .35, 3 + .35);
        \draw (22 - .35, 1 + .35) -- (22 + .35, 1 - .35);
        \draw (22 - .35, 1 - .35) -- (22 + .35, 1 + .35);
        \end{tikzpicture}
 * \f]
 * though is performed more efficiently.
 *
 * - If \p qureg is a state-vector, the output amplitudes are a kronecker product of 
 *   the discrete Fourier transform (DFT) acting upon the targeted amplitudes, and the 
 *   remaining. \n
 *   Precisely, 
 *   - let \f$|x,r\rangle\f$ represent a computational basis state where 
 *     \f$x\f$ is the binary value of the targeted qubits, and \f$r\f$ is the binary value 
 *     of the remaining qubits. 
 *   - let \f$|x_j,r_j\rangle\f$ be the \f$j\text{th}\f$ such state.
 *   - let \f$n =\f$ \p numQubits, and \f$N =\f$ `qureg.numQubitsRepresented`.\n
 * Then, this function effects
 *   \f[
 *      (\text{QFT}\otimes 1) \, \left(  \sum\limits_{j=0}^{2^N-1} \alpha_j \, |x_j,r_j\rangle \right) 
 *      = 
 *      \frac{1}{\sqrt{2^n}}
 *       \sum\limits_{j=0}^{2^N-1} \alpha_j \left( 
 *           \sum\limits_{y=0}^{2^n-1} e^{2 \pi \, i \, x_j \, y / 2^n} \; 
 *           |y,r_j \rangle
 *       \right)
 *   \f]
 *
 * - If \p qureg is a density matrix \f$\rho\f$, it will be changed under the unitary action 
 *   of the QFT. This can be imagined as each mixed state-vector undergoing the DFT 
 *   on its amplitudes. This is true even if \p qureg is unnormalised.
 *   \f[
 *      \rho \; \rightarrow \; \text{QFT} \; \rho \; \text{QFT}^{\dagger}
 *   \f]
 *
 * > This function merges contiguous controlled-phase gates into single invocations 
 * > of applyNamedPhaseFunc(), and hence is significantly faster than performing 
 * > the QFT circuit directly.
 *
 * > Furthermore, in distributed mode, this function requires only \f$\log_2(\text{\#nodes})\f$ 
 * > rounds of pair-wise
 * > communication, and hence is exponentially faster than directly performing the 
 * > DFT on the amplitudes of \p qureg.
 *
 * @see 
 * - applyFullQFT() to apply the QFT to the entirety of \p qureg.
 *
 * @ingroup operator
 * @param[in,out] qureg a state-vector or density matrix to modify
 * @param[in] qubits a list of the qubits to operate the QFT upon
 * @param[in] numQubits the length of list \p qubits
 * @throws invalidQuESTInputError()
 * - if any qubit in \p qubits is invalid, i.e. outside <b>[0, </b>`qureg.numQubitsRepresented`<b>)</b>
 * - if \p qubits contain any repetitions
 * - if \p numQubits <b>< 1</b> 
 * - if \p numQubits <b>></b>`qureg.numQubitsRepresented`
 * @throws segmentation-fault
 * - if \p qubits contains fewer elements than \p numQubits
 * @author Tyson Jones
 */
void applyQFT(Qureg qureg, int* qubits, int numQubits);

/** Force the target \p qubit of \p qureg into the given classical \p outcome, via a 
 * non-renormalising projection.
 *
 * This function zeroes all amplitudes in the state-vector or density-matrix which 
 * correspond to the opposite \p outcome given. Unlike collapseToOutcome(), it does 
 * not thereafter normalise \p qureg, and hence may leave it in a non-physical state.
 *
 * Note there is no requirement that the \p outcome state has a non-zero proability, and hence 
 * this function may leave \p qureg in a blank state, like that produced by initBlankState().
 * 
 * @see
 * - collapseToOutcome() for a norm-preserving equivalent, like a forced measurement
 *
 * @ingroup operator
 * @param[in,out] qureg a state-vector or density matrix to modify
 * @param[in] qubit the qubit to which to apply the projector 
 * @param[in] the single-qubit outcome (`0` or `1`) to project \p qubit into
 * @throws invalidQuESTInputError()
 * - if \p qubit is outside [0, `qureg.numQubitsRepresented`)
 * - if \p outcome is not in {0,1}
 * @author Tyson Jones
 */
void applyProjector(Qureg qureg, int qubit, int outcome);

// end prevention of C++ name mangling
#ifdef __cplusplus
}
#endif

#endif // QUEST_H

