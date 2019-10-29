// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * The QuEST API.
 * This file contains the comments used by doxygen for generating API doc
 *
 * The API is split into categories
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
 * @author Balint Koczor (Kraus maps)
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

/** Represents a 2x2 matrix of complex numbers
 * 
 * @ingroup type
 * @author Balint Koczor
 */
typedef struct ComplexMatrix2
{
    qreal real[2][2];
    qreal imag[2][2];
} ComplexMatrix2;

/** Represents a 4x4 matrix of complex numbers
 *
 * @ingroup type
 * @author Balint Koczor
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
 */
typedef struct QuESTEnv
{
    int rank;
    int numRanks;
} QuESTEnv;



/*
 * public functions
 */

/** Create a Qureg object representing a set of qubits which will remain in a pure state.
 * Allocate space for state vector of probability amplitudes, including space for temporary values to be copied from
 * one other chunk if running the distributed version. Define properties related to the size of the set of qubits.
 * The qubits are initialised in the zero state (i.e. initZeroState is automatically called)
 *
 * @ingroup type
 * @returns an object representing the set of qubits
 * @param[in] numQubits number of qubits in the system
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @throws exitWithError if \p numQubits <= 0
 * @author Ania Brown
 */
Qureg createQureg(int numQubits, QuESTEnv env);

/** Create a Qureg for qubits which are represented by a density matrix, and can be in mixed states.
 * Allocates space for a density matrix of probability amplitudes, including space for temporary values to be copied from
 * one other chunk if running the distributed version. Define properties related to the size of the set of qubits.
 * initZeroState is automatically called allocation, so that the density qureg begins in the zero state |0><0|.
 *
 * @ingroup type
 * @returns an object representing the set of qubits
 * @param[in] numQubits number of qubits in the system
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @throws exitWithError if \p numQubits <= 0
 * @author Tyson Jones
 */
Qureg createDensityQureg(int numQubits, QuESTEnv env);

/** Create a new Qureg which is an exact clone of the passed qureg, which can be
 * either a statevector or a density matrix. That is, it will have the same 
 * dimensions as the passed qureg and begin in an identical quantum state.
 * This must be destroyed by the user later with destroyQureg()
 *
 * @ingroup type
 * @returns an object representing the set of qubits
 * @param[in] qureg an existing qureg to be cloned
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @author Tyson Jones
 */
Qureg createCloneQureg(Qureg qureg, QuESTEnv env);

/** Deallocate a Qureg object representing a set of qubits.
 * Free memory allocated to state vector of probability amplitudes, including temporary vector for
 * values copied from another chunk if running the distributed version.
 *
 * @ingroup type
 * @param[in,out] qureg object to be deallocated
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @author Ania Brown
 */
void destroyQureg(Qureg qureg, QuESTEnv env);

/** Create (dynamically) a square complex matrix which can be passed to the multi-qubit general unitary functions.
 * The matrix will have dimensions (2^\p numQubits) by (2^\p numQubits), and all elements
 * of .real and .imag are initialised to zero.
 * The ComplexMatrixN must eventually be freed using destroyComplexMatrixN().
 * Like ComplexMatrix2 and ComplexMatrix4, the returned ComplexMatrixN is safe to 
 * return from functions.
 * 
 * One can instead use getStaticComplexMatrixN() to create a ComplexMatrixN struct 
 * in the stack (which doesn't need to be explicitly freed).
 *
 * @ingroup type
 * @param[in] numQubits the number of qubits of which the returned ComplexMatrixN will correspond
 * @returns a dynamic ComplexMatrixN struct, that is one where the .real and .imag
 *  fields are arrays kept in the heap and must be freed.
 * @author Tyson Jones
 */
ComplexMatrixN createComplexMatrixN(int numQubits);

/** Destroy a ComplexMatrixN instance created with createComplexMatrixN()
 *
 * It is invalid to attempt to destroy a matrix created with getStaticComplexMatrixN().
 *
 * @ingroup type
 * @param[in] matr the dynamic matrix (created with createComplexMatrixN()) to deallocate
 * @throws exitWithError if \p matr was not yet allocated.
 * @throws malloc_error if \p matr was static (created with getStaticComplexMatrixN())
 * @author Tyson Jones
 */
void destroyComplexMatrixN(ComplexMatrixN matr);

#ifndef __cplusplus
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
 * @throws exitWithError if \p m has not been allocated (e.g. with createComplexMatrixN())
 * @author Tyson Jones
 */
void initComplexMatrixN(ComplexMatrixN m, qreal real[][1<<m.numQubits], qreal imag[][1<<m.numQubits]);
#endif 

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

/** Get the number of qubits in a qureg object
 *
 * @ingroup calc
 * @author Tyson Jones
 */
int getNumQubits(Qureg qureg);

/** Get the number of probability amplitudes in a qureg object, given by 2^numQubits
 *
 * @ingroup calc
 * @author Tyson Jones
 */
long long int getNumAmps(Qureg qureg);

/** Initialises a qureg to have all-zero-amplitudes. This is an unphysical state 
 * useful for iteratively building a state with e.g. \p setWeightedQureg,
 * and should not be confused with the zero state |0...0>
 *
 * @ingroup init
 * @param[in,out] qureg the object representing the set of all qubits to initialise
 * @author Tyson Jones
 */
void initBlankState(Qureg qureg);

/** Initialise a set of \f$ N \f$ qubits to the classical zero state 
 * \f$ {| 0 \rangle}^{\otimes N} \f$.
 *
 * @ingroup init
 * @param[in,out] qureg the object representing the set of all qubits to initialise
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void initZeroState(Qureg qureg);

/** Initialise a set of \f$ N \f$ qubits to the plus state
 * \f$ {| + \rangle}^{\otimes N} = \frac{1}{\sqrt{2^N}} (| 0 \rangle + | 1 \rangle)^{\otimes N} \f$.
 * This is the product state of \f$N\f$ qubits where every classical state is uniformly 
 * populated with real coefficient \f$\frac{1}{\sqrt{2^N}}\f$.
 * This is equivalent to applying a Hadamard to every qubit in the zero state: 
 * \f$ \hat{H}^{\otimes N} {|0\rangle}^{\otimes N} \f$
 *
 * @ingroup init
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void initPlusState(Qureg qureg);

/** Initialise a set of \f$ N \f$ qubits to the classical state with index \p stateInd.
 * Note \f$ | 00 \dots 00 \rangle \f$ has \p stateInd 0, \f$ | 00 \dots 01 \rangle \f$ has \p stateInd 1, 
 * \f$ | 11 \dots 11 \rangle \f$ has \p stateInd \f$ 2^N - 1 \f$, etc.
 * Subsequent calls to getProbAmp will yield 0 for all indices except \p stateInd.
 *
 * @ingroup init
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @param[in] stateInd the index (0 to the number of amplitudes, exclusive) of the state to give probability 1
 * @author Tyson Jones
 */
void initClassicalState(Qureg qureg, long long int stateInd);

/** Initialise a set of \f$ N \f$ qubits, which can be a state vector or density matrix, to a given pure state.
 * If \p qureg is a state-vector, this merely makes \p qureg an identical copy of \p pure.
 * If \p qureg is a density matrix, this makes \p qureg 100% likely to be in the \p pure state.
 *
 * @ingroup init
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @param[in] pure the pure state to be copied or to give probability 1 in qureg
 * @author Tyson Jones
 */
void initPureState(Qureg qureg, Qureg pure);

/** Initialises \p qureg to be in the un-normalised, non-physical state with 
 * with n-th complex amplitude (2n/10 + i(2n+1)/10). This is used internally for 
 * debugging and testing.
 * 
 * @ingroup debug
 * @param[in,out] qureg the register to have its amplitudes overwritten
 * @author Ania Brown
 * @author Tyson Jones (doc)
 */
void initDebugState(Qureg qureg);

/** Initialise qureg by specifying the complete statevector.
 * The real and imaginary components of the amplitudes are passed in separate arrays,
 * each of which must have length \p qureg.numAmpsTotal.
 * There is no automatic checking that the passed arrays are L2 normalised.
 *
 * @ingroup init
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @param[in] reals array of the real components of the new amplitudes
 * @param[in] imags array of the imaginary components of the new amplitudes
 * @throws exitWithError
 *      if \p qureg is not a statevector (i.e. is a density matrix)
 * @author Tyson Jones
 */
void initStateFromAmps(Qureg qureg, qreal* reals, qreal* imags);

/** Edit qureg to adopt the subset of amplitudes passed in \p reals and \p imags.
 * Only amplitudes with indices in [\p startInd, \p startInd + \p numAmps] will be changed, which means
 * the new state may not be L2 normalised. This allows the user to initialise a custom state by 
 * setting batches of amplitudes.
 *
 * @ingroup init
 * @param[in,out] qureg the statevector to modify
 * @param[in] startInd the index of the first amplitude in \p qureg's statevector to modify
 * @param[in] reals array of the real components of the new amplitudes
 * @param[in] imags array of the imaginary components of the new amplitudes
 * @param[in] numAmps the length of each of the reals and imags arrays.
 * @throws exitWithError
 *      if \p qureg is not a statevector (i.e. is a density matrix)
 * @author Tyson Jones
 */
void setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps);

/** Set targetQureg to be a clone of copyQureg. 
 * Registers must either both be state-vectors, or both be density matrices.
 * Only the quantum state is cloned, auxilary info (like recorded QASM) is unchanged.
 * copyQureg is unaffected.
 *
 * @ingroup init
 * @param[in, out] targetQureg the qureg to have its quantum state overwritten
 * @param[in] copyQureg the qureg to have its quantum state cloned in targetQureg.
 * @author Tyson Jones
 */
void cloneQureg(Qureg targetQureg, Qureg copyQureg);

/** Shift the phase between \f$ |0\rangle \f$ and \f$ |1\rangle \f$ of a single qubit by a given angle.
 * This is equivalent to a rotation Z-axis of the Bloch-sphere up to a global phase factor.
 * For angle \f$\theta\f$, applies
 * \f[
 * \begin{pmatrix}
 * 1 & 0 \\
 * 0 & \exp(i \theta)
 * \end{pmatrix}
 * \f] 
 *     
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {rot};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_\theta$};
                \end{tikzpicture}
    }
    \f]
 * 
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to undergo a phase shift
 * @param[in] angle amount by which to shift the phase in radians
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Tyson Jones
 */
void phaseShift(Qureg qureg, const int targetQubit, qreal angle);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] idQubit1 first qubit in the state to phase shift
 * @param[in] idQubit2 second qubit in the state to phase shift
 * @param[in] angle amount by which to shift the phase in radians
 * @throws exitWithError
 *  if \p idQubit1 or \p idQubit2 are outside [0, \p qureg.numQubitsRepresented), or are equal
 * @author Tyson Jones
 */
void controlledPhaseShift(Qureg qureg, const int idQubit1, const int idQubit2, qreal angle);

/** Introduce a phase factor \f$ \exp(i \theta) \f$ on state \f$ |1 \dots 1 \rangle \f$
 * of the passed qubits.
 *     
   \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
   \f]
 * 
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits array of qubits to phase shift
 * @param[in] numControlQubits the length of array \p controlQubits
 * @param[in] angle amount by which to shift the phase in radians
 * @throws exitWithError
 *      if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit index in \p controlQubits is outside
 *      [0, \p qureg.numQubitsRepresented]), 
 *      or if any qubit in \p controlQubits is repeated.
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
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {idQubit1};
                \node[draw=none] at (-3.5, 0) {idQubit2};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, 0);
                
                \draw (-2,0) -- (2, 0);
                \draw[fill=black] (0, 0) circle (.2);
                \end{tikzpicture}
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] idQubit1, idQubit2 qubits to operate upon
 * @throws exitWithError 
 *  if \p idQubit1 or \p idQubit2 are outside [0, \p qureg.numQubitsRepresented), or are equal
 * @author Tyson Jones
 */
void controlledPhaseFlip (Qureg qureg, const int idQubit1, const int idQubit2);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
   \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits array of input qubits
 * @param[in] numControlQubits number of input qubits
 * @throws exitWithError 
 *      if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented),
 *      or if any qubit in \p controlQubits is outside [0, \p qureg.numQubitsRepresented),
 *      or if any qubit in \p qubits is repeated.
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
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {S};
                \end{tikzpicture}
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate upon
 * @throws exitWithError if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void sGate(Qureg qureg, const int targetQubit);

/** Apply the single-qubit T gate.
 * This is a rotation of \f$\pi/4\f$ around the Z-axis on the Bloch sphere, or the unitary:
 * \f[
 * \begin{pmatrix}
 * 1 & 0 \\
 * 0 & \exp\left(i \frac{\pi}{4}\right)
 * \end{pmatrix}
 * \f]
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {T};
                \end{tikzpicture}
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate upon
 * @throws exitWithError if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void tGate(Qureg qureg, const int targetQubit);

/** Create the QuEST execution environment.
 * This should be called only once, and the environment should be freed with destroyQuESTEnv at the end
 * of the user's code.
* If something needs to be done to set up the execution environment, such as 
 * initializing MPI when running in distributed mode, it is handled here.
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
 * @ingroup type
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * @author Ania Brown
 */
void destroyQuESTEnv(QuESTEnv env);

/** Guarantees that all code up to the given point has been executed on all nodes (if running in distributed mode)
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

/** Sets \p str to a string containing the number of qubits in \p qureg, and the 
 * hardware facilities used (e.g. GPU, MPI and/or OMP).
 * 
 * @ingroup debug 
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * @param[in] qureg the qureg of which to query the simulating hardware
 * @param[out] str to be populated with the output string
 * @author Ania Brown
 */
void getEnvironmentString(QuESTEnv env, Qureg qureg, char str[200]);

/** In GPU mode, this copies the state-vector (or density matrix) from RAM 
 * (qureg.stateVec) to VRAM / GPU-memory (qureg.deviceStateVec), which is the version 
 * operated upon by other calls to the API. 
 * In CPU mode, this function has no effect.
 * In conjunction with copyStateFromGPU() (which should be called first), this allows 
 * a user to directly modify the state-vector in a harware agnostic way.
 * Note though that users should instead use setAmps() if possible.
 *
 * For example, to set the first real element to 1, one could do:
 * 
 *     copyStateFromGPU(qureg);
 *     qureg.stateVec.real[0] = 1;
 *     copyStateToGPU(qureg);
 *
 * Note users should never access qureg.deviceStateVec directly.
 *
 * @ingroup debug
 * @param[in, out] qureg the qureg of which to copy .stateVec to .deviceStateVec in GPU mode
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
 * 
 *     copyStateFromGPU(qureg);
 *     qureg.stateVec.real[0] = 1;
 *     copyStateToGPU(qureg);
 *
 * Note users should never access qureg.deviceStateVec directly.
 *
 * @ingroup debug
 * @param[in, out] qureg the qureg of which to copy .deviceStateVec to .stateVec in GPU mode
 * @author Ania Brown 
 * @author Tyson Jones (doc)
 */
void copyStateFromGPU(Qureg qureg);

/** Get the complex amplitude at a given index in the state vector.
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return amplitude at index, returned as a Complex struct (with .real and .imag attributes)
 * @throws exitWithError
 *      if \p qureg is a density matrix,
 *      or if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Tyson Jones
 */
Complex getAmp(Qureg qureg, long long int index);

/** Get the real component of the complex probability amplitude at an index in the state vector.
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return real component at that index
 * @throws exitWithError
 *      if \p qureg is a density matrix,
 *      or if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Ania Brown
 */
qreal getRealAmp(Qureg qureg, long long int index);

/** Get the imaginary component of the complex probability amplitude at an index in the state vector.
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return imaginary component at that index
 * @throws exitWithError
 *      if \p qureg is a density matrix,
 *      or if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Ania Brown
 */
qreal getImagAmp(Qureg qureg, long long int index);

/** Get the probability of a state-vector at an index in the full state vector.
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return realEl*realEl + imagEl*imagEl
 * @throws exitWithError
 *      if \p qureg is a density matrix,
 *      or if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Ania Brown
 */
qreal getProbAmp(Qureg qureg, long long int index);

/** Get an amplitude from a density matrix at a given row and column.
 *
 * @ingroup calc
 * @param[in] qureg object representing a density matrix
 * @param[in] row row of the desired amplitude in the density matrix
 * @param[in] col column of the desired amplitude in the density matrix
 * @return a Complex scalar representing the desired amplitude
 * @throws exitWithError
 *      if \p qureg is a statevector, 
 *      or if \p row or \p col are outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 * @author Tyson Jones
 */
Complex getDensityAmp(Qureg qureg, long long int row, long long int col);

/** A debugging function which calculates the probability of being in any state, which should always be 1.
 * For pure states, this is the norm of the entire state vector and for mixed states, is the trace of
 * the density matrix.
 * Note this calculation utilises Kahan summation for greaster accuracy, but is
 * not parallelised and so will be slower than other functions.
 *
 * @ingroup calc
 * @param[in] qureg object representing a set of qubits
 * @return total probability
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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] alpha complex unitary parameter (row 1, column 1)
 * @param[in] beta complex unitary parameter (row 2, column 1)
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p alpha, \p beta don't satisfy |\p alpha|^2 + |\p beta|^2 = 1.
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void compactUnitary(Qureg qureg, const int targetQubit, Complex alpha, Complex beta);

/** Apply a general single-qubit unitary (including a global phase factor).
 * The passed 2x2 ComplexMatrix must be unitary, otherwise an error is thrown.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {U};
                \end{tikzpicture}
    }
    \f]
 *
 * @ingroup unitary                                                              
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or matrix \p u is not unitary.
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void unitary(Qureg qureg, const int targetQubit, ComplexMatrix2 u);

/** Rotate a single qubit by a given angle around the X-axis of the Bloch-sphere. For angle \f$\theta\f$, applies
 * \f[
 * \begin{pmatrix}
 * \cos\theta/2 & -i \sin \theta/2\\
 * -i \sin \theta/2 & \cos \theta/2
 * \end{pmatrix}
 * \f]
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {rot};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_x(\theta)$};
                \end{tikzpicture}
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void rotateX(Qureg qureg, const int rotQubit, qreal angle);

/** Rotate a single qubit by a given angle around the Y-axis of the Bloch-sphere. 
 * For angle \f$\theta\f$, applies
 * \f[
 * \begin{pmatrix}
 * \cos\theta/2 & - \sin \theta/2\\
 * \sin \theta/2 & \cos \theta/2
 * \end{pmatrix}
 * \f]            
 * 
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {rot};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_y(\theta)$};
                \end{tikzpicture}
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc, debug)
 */
void rotateY(Qureg qureg, const int rotQubit, qreal angle);

/** Rotate a single qubit by a given angle around the Z-axis of the Bloch-sphere (also known as a phase shift gate).   
 * For angle \f$\theta\f$, applies
 * \f[
 * \begin{pmatrix}
 * \exp(-i \theta/2) & 0 \\
 * 0 & \exp(i \theta/2)
 * \end{pmatrix}
 * \f] 
 *     
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {rot};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$R_z(\theta)$};
                \end{tikzpicture}
    }
    \f]
 * 
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void rotateZ(Qureg qureg, const int rotQubit, qreal angle);

/** Rotate a single qubit by a given angle around a given \ref Vector on the Bloch-sphere.      
 * The vector must not be zero (else an error is thrown), but needn't be unit magnitude.
 *
 * For angle \f$\theta\f$ and axis vector \f$\vec{n}\f$, applies \f$R_{\hat{n}} = \exp \left(- i \frac{\theta}{2} \hat{n} \cdot \vec{\sigma} \right) \f$
 * where \f$\vec{\sigma}\f$ is the vector of Pauli matrices.
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @param[in] axis vector around which to rotate (can be non-unit; will be normalised)
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p axis is the zero vector
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void rotateAroundAxis(Qureg qureg, const int rotQubit, qreal angle, Vector axis);


/** Applies a controlled rotation by a given angle around the X-axis of the Bloch-sphere. 
 * The target qubit is rotated in states where the control qubit has value 1.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f] 
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal.
 * @author Tyson Jones
 */
void controlledRotateX(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle);

/** Applies a controlled rotation by a given angle around the Y-axis of the Bloch-sphere. 
 * The target qubit is rotated in states where the control qubit has value 1.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f] 
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal.
 * @author Tyson Jones
 */
void controlledRotateY(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle);

/** Applies a controlled rotation by a given angle around the Z-axis of the Bloch-sphere. 
 * The target qubit is rotated in states where the control qubit has value 1.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f] 
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal.
 * @author Tyson Jones
 */
void controlledRotateZ(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle);

/** Applies a controlled rotation by a given angle around a given vector on the Bloch-sphere.      
 * The vector must not be zero (else an error is thrown), but needn't be unit magnitude.
 *
 * For angle \f$\theta\f$ and axis vector \f$\vec{n}\f$, applies \f$R_{\hat{n}} = \exp \left(- i \frac{\theta}{2} \hat{n} \cdot \vec{\sigma} \right) \f$ to states where the target qubit is 1 
 * (\f$\vec{\sigma}\f$ is the vector of Pauli matrices).
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit with value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @param[in] axis vector around which to rotate (can be non-unit; will be normalised)
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal
 *      or if \p axis is the zero vector
 * @author Tyson Jones
 */
void controlledRotateAroundAxis(Qureg qureg, const int controlQubit, const int targetQubit, qreal angle, Vector axis);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * @ingroup unitary                                                             
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit apply the target unitary if this qubit has value 1
 * @param[in] targetQubit qubit on which to apply the target unitary
 * @param[in] alpha complex unitary parameter (row 1, column 1)
 * @param[in] beta complex unitary parameter (row 2, column 1)
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal,
 *      or if \p alpha, \p beta don't satisfy |\p alpha|^2 + |\p beta|^2 = 1.
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void controlledCompactUnitary(Qureg qureg, const int controlQubit, const int targetQubit, Complex alpha, Complex beta);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * @ingroup unitary                                                           
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit apply unitary if this qubit is 1
 * @param[in] targetQubit qubit to operate on
 * @param[in] u single-qubit unitary matrix to apply
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal,
 *      or if \p u is not unitary.
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void controlledUnitary(Qureg qureg, const int controlQubit, const int targetQubit, ComplexMatrix2 u);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * @ingroup unitary                                                            
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits applies unitary if all qubits in this array equal 1
 * @param[in] numControlQubits number of control qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] u single-qubit unitary matrix to apply
 * @throws exitWithError
 *      if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit index (\p targetQubit or one in \p controlQubits) is outside
 *      [0, \p qureg.numQubitsRepresented]), 
 *      or if any qubit in \p controlQubits is repeated,
 *      or if \p controlQubits contains \p targetQubit,
 *      or if \p u is not unitary.
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void multiControlledUnitary(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit, ComplexMatrix2 u);

/** Apply the single-qubit Pauli-X (also known as the X, sigma-X, NOT or bit-flip) gate.
 * This is a rotation of \f$\pi\f$ around the x-axis on the Bloch sphere. I.e. 
 * \f[
 * \begin{pmatrix}
 * 0 & 1 \\
 * 1 & 0
 * \end{pmatrix}
 * \f]   
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (2, 0);
                \draw (0, 0) circle (.5);
                \draw (0, .5) -- (0, -.5);
                \end{tikzpicture}
    }
    \f]   
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void pauliX(Qureg qureg, const int targetQubit);

/** Apply the single-qubit Pauli-Y (also known as the Y or sigma-Y) gate.
 * This is a rotation of \f$\pi\f$ around the Y-axis on the Bloch sphere. I.e. 
 * \f[
 * \begin{pmatrix}
 * 0 & -i \\
 * i & 0
 * \end{pmatrix}
 * \f]  
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$\sigma_y$};
                \end{tikzpicture}
    }
    \f]      
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void pauliY(Qureg qureg, const int targetQubit);

/** Apply the single-qubit Pauli-Z (also known as the Z, sigma-Z or phase-flip) gate.
 * This is a rotation of \f$\pi\f$ around the Z-axis (a phase shift) on the Bloch sphere. I.e. 
 * \f[
 * \begin{pmatrix}
 * 1 & 0 \\
 * 0 & -1
 * \end{pmatrix}
 * \f]   
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {$\sigma_z$};
                \end{tikzpicture}
    }
    \f]     
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void pauliZ(Qureg qureg, const int targetQubit);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2,0) -- (-1, 0);
                \draw (1, 0) -- (2, 0);
                \draw (-1,-1)--(-1,1)--(1,1)--(1,-1)--cycle;
                \node[draw=none] at (0, 0) {H};
                \end{tikzpicture}
    }
    \f]  
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void hadamard(Qureg qureg, const int targetQubit);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
                \begin{tikzpicture}[scale=.5]
                \node[draw=none] at (-3.5, 2) {control};
                \node[draw=none] at (-3.5, 0) {target};

                \draw (-2, 2) -- (2, 2);
                \draw[fill=black] (0, 2) circle (.2);
                \draw (0, 2) -- (0, -.5);
                
                \draw (-2,0) -- (2, 0);
                \draw (0, 0) circle (.5);
                \end{tikzpicture}
    }
    \f]  
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit nots the target if this qubit is 1
 * @param[in] targetQubit qubit to not
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented), or are equal.
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix, doc)
 */
void controlledNot(Qureg qureg, const int controlQubit, const int targetQubit);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit applies pauliY to the target if this qubit is 1
 * @param[in] targetQubit qubit to not
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented), or are equal.
 * @author Tyson Jones
 * @author Ania Brown (debug)
 */
void controlledPauliY(Qureg qureg, const int controlQubit, const int targetQubit);

/** Gives the probability of a specified qubit being measured in the given outcome (0 or 1).
 * This performs no actual measurement and does not change the state of the qubits.
 * 
 * @ingroup calc
 * @param[in] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to study
 * @param[in] outcome for which to find the probability of the qubit being measured in
 * @return probability of qubit measureQubit being measured in the given outcome
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p outcome is not in {0, 1}.
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
qreal calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome);

/** Updates the state vector to be consistent with measuring the measure qubit in the given outcome (0 or 1), and returns the probability of such a measurement outcome. 
 * This is effectively performing a measurement and forcing the outcome.
 * This is an irreversible change to the state vector, whereby incompatible states
 * in the state vector are given zero amplitude and the remaining states are renormalised.
 * Exits with error if the given outcome has ~zero probability, and so cannot be
 * collapsed into.
 * 
 * @ingroup normgate
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @param[in] outcome to force the measure qubit to enter
 * @return probability of the (forced) measurement outcome
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p outcome is not in {0, 1},
 *      or if the probability of \p outcome is zero (within machine epsilon)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
qreal collapseToOutcome(Qureg qureg, const int measureQubit, int outcome);

/** Measures a single qubit, collapsing it randomly to 0 or 1.
 * Outcome probabilities are weighted by the state vector, which is irreversibly
 * changed after collapse to be consistent with the outcome.
 *
 * @ingroup normgate
 * @param[in, out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @return the measurement outcome, 0 or 1
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
int measure(Qureg qureg, int measureQubit);

/** Measures a single qubit, collapsing it randomly to 0 or 1, and
 * additionally gives the probability of that outcome.
 * Outcome probabilities are weighted by the state vector, which is irreversibly
 * changed after collapse to be consistent with the outcome.
 *
 * @ingroup normgate
 * @param[in, out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @param[out] outcomeProb a pointer to a qreal which is set to the probability of the occurred outcome
 * @return the measurement outcome, 0 or 1
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 * @author Ania Brown (state-vector)
 * @author Tyson Jones (density matrix)
 */
int measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb);

/** Computes the inner product \f$ \langle \text{bra} | \text{ket} \rangle \f$ of two 
 * equal-size state vectors. The same \p qureg may be passed as both \p bra and \p ket, 
 * though we recommend users check state-vector normalisation with \p calcTotalProb which 
 * employs Kahan summation for greater accuracy.
 * Neither state-vector is modified.
 *
 * @ingroup calc
 * @param[in] bra qureg to be the 'bra' (i.e. have its values conjugate transposed) in the inner product 
 * @param[in] ket qureg to be the 'ket' in the inner product 
 * @return the complex inner product of \p bra and \p ket 
 * @throws exitWithError
 *      if either \p bra or \p ket are not state-vectors, 
 *      or if \p bra and \p ket do not have equal dimensions.
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
    ((\rho_1, \rho_2))_{HS} = ((\rho_2, \rho_1))_{HS}
 * \f]
 * Also note that if both \p rho1 and \p rho2 are density matrices of pure states
 * \p bra and \p ket, then the equality holds
 * \f[
    ((\rho_1, \rho_2))_{HS} = |\langle \text{bra} | \text{ket} \rangle|^2.
 * \f]
 *
 * @ingroup calc
 * @param[in] rho1 qureg as a density matrix (to have its values conjugate transposed)
 * @param[in] rho2 qureg as a density matrix
 * @returns the real Hilbert-Schmidt scalar product of density matrices
            \p rho1 and \p rho2 (assuming Hermiticity)
 * @throws exitWithError
 *      if \p rho1 and \p rho2 are not density matrices or
 *      have mismatching dimensions.
 * @author Balint Koczor (CPU)
 * @author Tyson Jones (GPU)
 */
qreal calcDensityInnerProduct(Qureg rho1, Qureg rho2);

/** Seed the Mersenne Twister used for random number generation in the QuEST environment with an example
 * defualt seed.
 * This default seeding function uses the mt19937 init_by_array function with two keys -- 
 * time and pid. Subsequent calls to mt19937 genrand functions will use this seeding. 
 * For a multi process code, the same seed is given to all process, therefore this seeding is only
 * appropriate to use for functions such as measure where all processes require the same random value.
 *
 * For more information about the MT, see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
 * 
 * @ingroup debug
 * @author Ania Brown
 * @author Balint Koczor (Windows compatibility)
 **/
void seedQuESTDefault(void);

/** Seed the Mersenne Twister used for random number generation in the QuEST environment with
 * a user defined seed.
 * This function uses the mt19937 init_by_array function with numSeeds keys supplied by the user.
 * Subsequent calls to mt19937 genrand functions will use this seeding. 
 * For a multi process code, the same seed is given to all process, therefore this seeding is only
 * appropriate to use for functions such as measure where all processes require the same random value.
 *
 * For more information about the MT, see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
 *
 * @ingroup debug
 * @param[in] seedArray Array of integers to use as seed. 
 *  This allows the MT to be initialised with more than a 32-bit integer if required
 * @param[in] numSeeds Length of seedArray
 * @author Ania Brown
 **/
void seedQuEST(unsigned long int *seedArray, int numSeeds);

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
 * @throws exitWithError if \p filename cannot be written to
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
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce dephasing noise
 * @param[in] prob the probability of the phase error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p prob is not in [0, 1/2]
 * @author Tyson Jones (GPU, doc)
 * @author Ania Brown (CPU, distributed)
 */
void mixDephasing(Qureg qureg, const int targetQubit, qreal prob);

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
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] qubit1 qubit upon which to induce dephasing noise
 * @param[in] qubit2 qubit upon which to induce dephasing noise
 * @param[in] prob the probability of the phase error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if either \p qubit1 or \p qubit2 is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p qubit1 = \p qubit2,
 *      or if \p prob is not in [0, 3/4]
 * @author Tyson Jones (GPU, doc)
 * @author Ania Brown (CPU, distributed)
 */
void mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal prob);

/** Mixes a density matrix \p qureg to induce single-qubit homogeneous depolarising noise.
 * With probability \p prob, applies (uniformly) either Pauli X, Y, or Z to \p targetQubit.
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
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p prob is not in [0, 3/4]
 * @author Tyson Jones (GPU, doc)
 * @author Ania Brown (CPU, distributed)
 */
void mixDepolarising(Qureg qureg, const int targetQubit, qreal prob);

/** Mixes a density matrix \p qureg to induce single-qubit damping (decay to 0 state).
 * With probability \p prob, applies damping (transition from 1 to 0 state).
 *
 * This transforms \p qureg = \f$\rho\f$ into the mixed state
 * \f[
 * (1 - \text{prob}) \, \rho + \text{prob} \; \left( 
 *      \sigma^{-} \, \rho \, sigma^{+} 
 * \right)
 * \f]
 * where q = \p targetQubit.
 * \p prob cannot exceed 1, at which total damping/decay occurs.
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p prob is not in [0, 1]
 * @author Nicolas Vogt (HQS)
 * @author Ania Brown (patched)
 */
void mixDamping(Qureg qureg, const int targetQubit, qreal prob);

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
 *      X_a \, \rho \, X_a + 
 *      X_b \, \rho \, X_b + 
 *      Y_a \, \rho \, Y_a + 
 *      Y_b \, \rho \, Y_b + 
 *      Z_a \, \rho \, Z_a + 
 *      Z_b \, \rho \, Z_b + 
 *      X_a X_b \, \rho \, X_a X_b +
 *      X_a Y_b \, \rho \, X_a Y_b +
 *      X_a Z_b \, \rho \, X_a Z_b +
 *      Y_a X_b \, \rho \, Y_a X_b +
 *      Y_a Y_b \, \rho \, Y_a Y_b +
 *      Y_a Z_b \, \rho \, Y_a Z_b +
 *      Z_a X_b \, \rho \, Z_a X_b + 
 *      Z_a Y_b \, \rho \, Z_a Y_b + 
 *      Z_a Z_b \, \rho \, Z_a Z_b
 * \right)
 * \f]
 * where a = \p qubit1, b = \p qubit2.
 * \p prob cannot exceed 15/16, at which maximal mixing occurs.
 *
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] qubit1 qubit upon which to induce depolarising noise
 * @param[in] qubit2 qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if either \p qubit1 or \p qubit2 is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p qubit1 = \p qubit2,
 *      or if \p prob is not in [0, 15/16]
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
 * @ingroup decoherence
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit to decohere
 * @param[in] probX the probability of inducing an X error
 * @param[in] probX the probability of inducing an Y error
 * @param[in] probX the probability of inducing an Z error
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if any of \p probX, \p probY or \p probZ are not in [0, 1],
 *      or if any of p in {\p probX, \p probY or \p probZ} don't satisfy
 *      p <= (1 - \p probX - \p probY - \p probZ)
 * @author Balint Koczor
 * @author Tyson Jones (refactored, doc)
 */
void mixPauli(Qureg qureg, int targetQubit, qreal probX, qreal probY, qreal probZ);

/** Modifies combineQureg to become (1-prob)combineProb + prob otherQureg.
 * Both registers must be equal-dimension density matrices, and prob must be in [0, 1].
 *
 * @ingroup decoherence
 * @param[in,out] combineQureg a density matrix to be modified
 * @param[in] prob the probability of \p otherQureg in the modified \p combineQureg
 * @param[in] otherQureg a density matrix to be mixed into \p combineQureg
 * @throws exitWithError
 *      if either \p combineQureg or \p otherQureg are not density matrices,
 *      or if the dimensions of \p combineQureg and \p otherQureg do not match,
 *      or if \p prob is not in [0, 1]
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
 * @ingroup calc
 * @param[in] qureg a density matrix of which to measure the purity
 * @return the purity
 * @throws exitWithError
 *      if either \p combineQureg or \p otherQureg are not density matrices,
 *      or if the dimensions of \p combineQureg and \p otherQureg do not match,
 *      or if \p prob is not in [0, 1]
 * @author Tyson Jones
 */
qreal calcPurity(Qureg qureg);

/** Calculates the fidelity of qureg (a statevector or density matrix) against 
 * a reference pure state (necessarily a statevector).
 * For two pure states, this computes 
 * \f[ 
    |\langle \text{qureg} | \text{pureState} \rangle|^2
 * \f]
 * For a mixed and pure state, this computes 
 * \f[ 
    \langle \text{pureState} | \text{qureg} | \text{pureState} \rangle
 * \f]
 * In either case, the fidelity lies in [0, 1].
 * The number of qubits represented in \p qureg and \p pureState must match.
 * 
 * @ingroup calc
 * @param[in] qureg a density matrix or state vector
 * @param[in] pureState a state vector
 * @return the fidelity between the input registers
 * @throws exitWithError
 *      if the second argument (\p pureState) is not a statevector, 
 *      or if the number of qubits in \p qureg and \p pureState do not match
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
   \setlength{\fboxrule}{0.01pt}
   \fbox{
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
   }
   \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] qubit1 qubit to swap
 * @param[in] qubit2 other qubit to swap
 * @throws exitWithError
 *      if either \p qubit1 or \p qubit2 are outside [0, \p qureg.numQubitsRepresented), or are equal.
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
   \setlength{\fboxrule}{0.01pt}
   \fbox{
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
   }
   \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] qubit1 qubit to sqrt swap
 * @param[in] qubit2 other qubit to sqrt swap
 * @throws exitWithError
 *      if either \p qubit1 or \p qubit2 are outside [0, \p qureg.numQubitsRepresented), or are equal.
 * @author Tyson Jones
 */
void sqrtSwapGate(Qureg qureg, int qb1, int qb2);

/** Apply a general multiple-control, conditioned on a specific bit sequence,
 *  single-target unitary, which can include a global phase factor. 
 * Any number of control qubits can be specified, along with which of their 
 * states (0 or 1) to condition upon; when the specified controls are in the 
 * specified state, the given unitary is applied to the target qubit.
 * This is equivalent to NOTing the control bits which are conditioned on 0, 
 * calling multiControlledUnitary then NOTing the same control bits.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits the indices of the control qubits 
 * @param[in] controlState the bit values (0 or 1) of each control qubit, upon which to condition
 * @param[in] numControlQubits number of control qubits
 * @param[in] targetQubit qubit to operate the unitary upon
 * @param[in] u single-qubit unitary matrix to apply
 * @throws exitWithError
 *      if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit index (\p targetQubit or one in \p controlQubits) is outside
 *      [0, \p qureg.numQubitsRepresented]), 
 *      or if any qubit in \p controlQubits is repeated.,
 *      or if \p controlQubits contains \p targetQubit,
 *      or if any element of controlState is not a bit (0 or 1),
 *      or if \p u is not unitary.
 * @author Tyson Jones
 */
void multiStateControlledUnitary(
    Qureg qureg, int* controlQubits, int* controlState, const int numControlQubits, 
    const int targetQubit, ComplexMatrix2 u
);

/** Apply a multi-qubit Z rotation on a selected number of qubits. 
 * This is the unitary 
 * \f[ 
    \exp \left( - i \theta/2 \bigotimes_{j} Z_j\right)
 * \f]
 * where the Pauli Z gates operate upon the passed list \f$j \in\f$ \p qubits, and cause 
 * rotations of \f$\theta =\f$ \p angle.
 * All qubits not appearing in \p qubits are assumed to receive the identity operator.
 * This has the effect of premultiplying every amplitude with 
 * \f$\exp(\pm i \theta/2)\f$ where the sign is determined by the parity of
 * the target qubits for that amplitude.
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] qubits a list of the indices of the target qubits 
 * @param[in] numQubits number of target qubits
 * @param[in] angle the angle by which the multi-qubit state is rotated around the Z axis
 * @throws exitWithError
 *      if \p numQubits is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit in \p qubits is outside [0, \p qureg.numQubitsRepresented])
 *      or if any qubit in \p qubits is repeated.
 * @author Tyson Jones
 */
void multiRotateZ(Qureg qureg, int* qubits, int numQubits, qreal angle);

/** Apply a multi-qubit multi-Pauli rotation on a selected number of qubits. 
 * This is the unitary 
 * \f[ 
    \exp \left( - i \theta/2 \bigotimes_{j} \hat{\sigma}_j\right)
 * \f]
 * where \f$\hat{\sigma}_j \in \{1, X, Y, Z\}\f$ is a Pauli operator (indicated by
 * codes 0, 1, 2, 3 respectively in \p targetPaulis, or by enums PAULI_I,
 * PAULI_X, PAULI_Y and PAULI_Z) operating upon the qubit 
 * \p targetQubits[j], and \f$\theta\f$ is the passed \p angle.
 *  The operators specified in \p targetPaulis act on the corresponding qubit in \p targetQubits. 
 * For example:
 * 
 *     multiRotatePauli(qureg, (int[]) {4,5,8,9}, (int[]) {0,1,2,3}, 4, .1)
 *
 * effects \f$ \exp \left( - i .1/2 X_5 Y_8 Z_9 \right) \f$ on \p qureg, 
 * where unspecified qubits (along with those specified with Pauli code 0) are 
 * assumed to receive the identity operator. Note that specifying the identity 
 * Pauli (code=0) on a qubit is superfluous but allowed for convenience.
 *
 * This function effects this unitary by first rotating the qubits which are 
 * nominated to receive X or Y Paulis into alternate basis, performing 
 * multiRotateZ on all target qubits receiving X, Y or Z Paulis, then restoring 
 * the original basis. In the worst case, this means that 1+2*\p numTargets
 * primitive unitaries are performed on the statevector, and double this on 
 * density matrices.
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubits a list of the indices of the target qubits 
 * @param[in] targetPaulis a list of the Pauli codes (0=PAULI_I, 1=PAULI_X, 2=PAULI_Y, 3=PAULI_Z) 
 *      to apply to the corresponding qubits in \p targetQubits
 * @param[in] numTargets number of target qubits, i.e. the length of \p targetQubits and \p targetPaulis
 * @param[in] angle the angle by which the multi-qubit state is rotated
 * @throws exitWithError
 *      if \p numQubits is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit in \p qubits is outside [0, \p qureg.numQubitsRepresented))
 *      or if any qubit in \p qubits is repeated.
 * @author Tyson Jones
 */
void multiRotatePauli(Qureg qureg, int* targetQubits, enum pauliOpType* targetPaulis, int numTargets, qreal angle);

/** Computes the expected value of a product of Pauli operators.
 * Letting \f$ \sigma = \otimes_j \hat{\sigma}_j \f$ be the operators indicated by \p pauliCodes 
 * and acting on qubits \p targetQubits, this function computes \f$ \langle \psi | \sigma | \psi \rangle \f$ 
 * if \p qureg = \f$ \psi \f$ is a statevector, and computes \f$ \text{Trace}(\sigma \rho) \f$ 
 * if \p qureg = \f$ \rho \f$ is a density matrix.
 * 
 * \p pauliCodes is an array of length \p numTargets which specifies which Pauli operators to 
 * enact on the corresponding qubits in \p targetQubits, where 0 = \p PAULI_I, 1 = \p PAULI_X, 
 * 2 = \p PAULI_Y, 3 = \p PAULI_Z. The target qubits must be unique, and at most \p qureg.numQubitsRepresented
 * may be specified. For example, on a 7-qubit statevector,
 * 
 *     calcExpecPauliProd(qureg, {4,5,6}, {PAULI_X, PAULI_I, PAULI_Z}, 3, workspace);
 *
 * will compute \f$ \langle \psi | I I I I X I Z | \psi \rangle \f$ (where in this notation, the left-most operator
 * applies to the least-significant qubit, i.e. that with index 0).
 *
 * \p workspace must be a register with the same type (statevector vs density matrix) and dimensions 
 * (numQubitsRepresented) as \p qureg, and is used as working space. When this function returns, \p qureg 
 * will be unchanged and \p workspace will be set to \f$ \sigma | \psi \rangle \f$ (if \p qureg is a statevector)
 * or \f$ \sigma \rho \f$ (if \p qureg is a density matrix). NOTE that this last quantity is NOT the result of applying 
 * the paulis as unitaries, \f$ \sigma^\dagger \rho \sigma \f$, but is instead the result of their direct 
 * multiplication with the density matrix. It is therefore itself not a valid density matrix.
 *
 * This function works by cloning the \p qureg state into \p workspace, applying the specified 
 * Pauli operators to \p workspace then computing its inner product with \p qureg (for statevectors)
 * or its trace (for density matrices). It therefore should scale linearly in time with the number of 
 * specified non-identity Pauli operators, which is bounded by the number of represented qubits.
 *
 * @ingroup calc
 * @param[in] qureg the register of which to find the expected value, which is unchanged by this function
 * @param[in] targetQubits a list of the indices of the target qubits 
 * @param[in] pauliCodes a list of the Pauli codes (0=PAULI_I, 1=PAULI_X, 2=PAULI_Y, 3=PAULI_Z) 
 *      to apply to the corresponding qubits in \p targetQubits
 * @param[in] numTargets number of target qubits, i.e. the length of \p targetQubits and \p pauliCodes
 * @param[in,out] workspace a working-space qureg with the same dimensions as \p qureg, which is modified 
 *      to be the result of multiplying the state with the pauli operators
 * @throws exitWithError
 *      if \p numTargets is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit in \p targetQubits is outside [0, \p qureg.numQubitsRepresented))
 *      or if any qubit in \p targetQubits is repeated,
 *      or if any code in \p pauliCodes is not in {0,1,2,3},
 *      or if \p workspace is not of the same type and dimensions as \p qureg
 * @author Tyson Jones
 */
qreal calcExpecPauliProd(Qureg qureg, int* targetQubits, enum pauliOpType* pauliCodes, int numTargets, Qureg workspace);

/** Computes the expected value of a sum of products of Pauli operators.
 * Letting \f$ \alpha = \sum_i c_i \otimes_j^{N} \hat{\sigma}_{i,j} \f$ be 
 * the operators indicated by \p allPauliCodes (where \f$ c_i \in \f$ \p termCoeffs and \f$ N = \f$ \p qureg.numQubitsRepresented), 
 * this function computes \f$ \langle \psi | \alpha | \psi \rangle \f$ 
 * if \p qureg = \f$ \psi \f$ is a statevector, and computes \f$ \text{Trace}(\alpha \rho) \f$ 
 * if \p qureg = \f$ \rho \f$ is a density matrix.
 *
 * \p allPauliCodes is an array of length \p numSumTerms*\p qureg.numQubitsRepresented
 * which specifies which Pauli operators to apply, where 0 = \p PAULI_I, 1 = \p PAULI_X, 
 * 2 = \p PAULI_Y, 3 = \p PAULI_Z. For each sum term, a Pauli operator must be specified for 
 * EVERY qubit in \p qureg; each set of \p numSumTerms operators will be grouped into a product.
 * \p termCoeffs is an arrray of length \p numSumTerms containing the term coefficients.
 * For example, on a 3-qubit statevector,
 *
 *     int paulis[6] = {PAULI_X, PAULI_I, PAULI_I,  PAULI_X, PAULI_Y, PAULI_Z};
 *     qreal coeffs[2] = {1.5, -3.6};
 *     calcExpecPauliSum(qureg, paulis, coeffs, 2, workspace);
 *
 * will compute \f$ \langle \psi | (1.5 X I I - 3.6 X Y Z) | \psi \rangle \f$ (where in this notation, the left-most operator
 * applies to the least-significant qubit, i.e. that with index 0).
 * 
 * \p workspace must be a register with the same type (statevector vs density matrix) and dimensions 
 * (numQubitsRepresented) as \p qureg, and is used as working space. When this function returns, \p qureg 
 * will be unchanged and \p workspace will be set to \p qureg pre-multiplied with the final Pauli product.
 * NOTE that if \p qureg is a density matrix, \p workspace will become \f$ \hat{\sigma} \rho \f$ 
 * which is itself not a density matrix (it is distinct from \f$ \hat{\sigma}^\dagger \rho \hat{\sigma} \f$).
 *
 * This function works by cloning the \p qureg state into \p workspace, applying each of the specified
 * Pauli products to \p workspace (one Pauli operation at a time), then computing its inner product with \p qureg (for statevectors)
 * or its trace (for density matrices) multiplied with the corresponding coefficient, and summing these contributions. 
 * It therefore should scale linearly in time with the total number of non-identity specified Pauli operators.
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
 * @throws exitWithError
 *      if any code in \p allPauliCodes is not in {0,1,2,3},
 *      or if numSumTerms <= 0,
 *      or if \p workspace is not of the same type and dimensions as \p qureg
 * @author Tyson Jones
 */
qreal calcExpecPauliSum(Qureg qureg, enum pauliOpType* allPauliCodes, qreal* termCoeffs, int numSumTerms, Qureg workspace);

/** Apply a general two-qubit unitary (including a global phase factor).
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * \p targetQubit1 is treated as the \p least significant qubit in \p u, such that 
 * a row in \p u is dotted with the vector
 * \f$ |\text{targetQubit2} \;\; \text{targetQubit1}\rangle : \{ |00\rangle, |01\rangle, |10\rangle, |11\rangle \} \f$
 *
 * For example, 

 *     twoQubitUnitary(qureg, a, b, u);
 *
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
 * The passed 4x4 ComplexMatrix must be unitary, otherwise an error is thrown.
 *                 
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q/4 nodes.
 * 
 * @ingroup unitary                                            
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit1 first qubit to operate on, treated as least significant in \p u
 * @param[in] targetQubit2 second qubit to operate on, treated as most significant in \p u
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if \p targetQubit1 or \p targetQubit2 are outside [0, \p qureg.numQubitsRepresented),
 *      or if \p targetQubit1 equals \p targetQubit2,
 *      or matrix \p u is not unitary, 
 *      or if each node cannot fit 4 amplitudes in distributed mode.
 * @author Tyson Jones
 */
void twoQubitUnitary(Qureg qureg, const int targetQubit1, const int targetQubit2, ComplexMatrix4 u);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q/4 nodes.
 *
 * @ingroup unitary                                                          
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit the control qubit which must be in state 1 to effect the given unitary
 * @param[in] targetQubit1 first qubit to operate on, treated as least significant in \p u
 * @param[in] targetQubit2 second qubit to operate on, treated as most significant in \p u
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if \p controlQubit, \p targetQubit1 or \p targetQubit2 are outside [0, \p qureg.numQubitsRepresented),
 *      or if any of \p controlQubit, \p targetQubit1 and \p targetQubit2 are equal,
 *      or matrix \p u is not unitary, 
 *     or if each node cannot fit 4 amplitudes in distributed mode.
 * @author Tyson Jones
 */
void controlledTwoQubitUnitary(Qureg qureg, const int controlQubit, const int targetQubit1, const int targetQubit2, ComplexMatrix4 u);

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
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * Note that in distributed mode, this routine requires that each node contains at least 4 amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q/4 nodes.
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits the control qubits which all must be in state 1 to effect the given unitary
 * @param[in] numControlQubits the number of control qubits
 * @param[in] targetQubit1 first target qubit, treated as least significant in \p u
 * @param[in] targetQubit2 second target qubit, treated as most significant in \p u
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if \p targetQubit1 or \p targetQubit2 are outside [0, \p qureg.numQubitsRepresented),
 *      or if \p targetQubit1 equals \p targetQubit2,
 *      or if any qubit in \p controlQubits is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p controlQubits are not unique, or if either \p targetQubit1 and \p targetQubit2
 *      are in \p controlQubits,
 *      or if matrix \p u is not unitary,
 *      or if each node cannot fit 4 amplitudes in distributed mode.
 * @author Tyson Jones
 */
void multiControlledTwoQubitUnitary(Qureg qureg, int* controlQubits, const int numControlQubits, const int targetQubit1, const int targetQubit2, ComplexMatrix4 u);

/** Apply a general multi-qubit unitary (including a global phase factor) with any number of target qubits.
 *
 * The first target qubit in \p targs is treated as \b least sigifnicant in \p u.
 * For example, 

 *     multiQubitUnitary(qureg, (int []) {a, b, c}, 3, u);
 *
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
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * Note that in distributed mode, this routine requires that each node contains at least 2^\p numTargs amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q / 2^\p numTargs nodes.
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targs a list of the target qubits, ordered least significant to most in \p u
 * @param[in] numTargs the number of target qubits
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if any index in \p targs is outside of [0, \p qureg.numQubitsRepresented),
 *      or if \p targs are not unique,
 *      or if matrix \p u is not unitary,
 *      or if a node cannot fit the required number of target amplitudes in distributed mode.
 * @author Tyson Jones
 */
void multiQubitUnitary(Qureg qureg, int* targs, const int numTargs, ComplexMatrixN u);

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
 * The target qubits in \p targs are treated as ordered least sigifnicant 
 * to most significant in \p u.
 *
 * The passed ComplexMatrix must be unitary and be a compatible size with the specified number of
 * target qubits, otherwise an error is thrown.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * Note that in distributed mode, this routine requires that each node contains at least 2^\p numTargs amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q / 2^\p numTargs nodes.
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] ctrl the control qubit
 * @param[in] targs a list of the target qubits, ordered least to most significant
 * @param[in] numTargs the number of target qubits
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if \p ctrl or any index in \p targs is outside of [0, \p qureg.numQubitsRepresented),
 *      or if \p targs are not unique,
 *      or if \p targs contains \p ctrl,
 *      or if matrix \p u is not unitary,
 *      or if a node cannot fit the required number of target amplitudes in distributed mode.
 * @author Tyson Jones
 */
void controlledMultiQubitUnitary(Qureg qureg, int ctrl, int* targs, const int numTargs, ComplexMatrixN u);

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
 * The target qubits in \p targs are treated as ordered least sigifnicant 
 * to most significant in \p u.
 *
 * The passed ComplexMatrix must be unitary and be a compatible size with the specified number of
 * target qubits, otherwise an error is thrown.
 *
    \f[
    \setlength{\fboxrule}{0.01pt}
    \fbox{
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
    }
    \f]
 *
 * Note that in distributed mode, this routine requires that each node contains at least 2^\p numTargs amplitudes.
 * This means an q-qubit register (state vector or density matrix) can be distributed 
 * by at most 2^q / 2^\p numTargs nodes.
 *
 * @ingroup unitary
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] ctrls a list of the control qubits
 * @param[in] numCtrls the number of control qubits
 * @param[in] targs a list of the target qubits, ordered least to most significant
 * @param[in] numTargs the number of target qubits
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if any index in \p ctrls and \p targs is outside of [0, \p qureg.numQubitsRepresented),
 *      or if \p ctrls and \p targs are not unique,
 *      or if matrix \p u is not unitary,
 *      or if a node cannot fit the required number of target amplitudes in distributed mode.
 * @author Tyson Jones
 */
void multiControlledMultiQubitUnitary(Qureg qureg, int* ctrls, const int numCtrls, int* targs, const int numTargs, ComplexMatrixN u);

/** Apply a general single-qubit Kraus map to a density matrix, as specified by at most 
 * four Kraus operators. A Kraus map is also referred to as a "operator-sum representation"
 * of a quantum channel. This allows one to simulate a general single-qubit noise process.
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
 * @ingroup decoherence
 * @param[in,out] qureg the density matrix to which to apply the map
 * @param[in] target the target qubit of the map
 * @param[in] ops an array of at most 4 Kraus operators
 * @param[in] numOps the number of operators in \p ops which must be >0 and <= 4.
 * @throws exitWithError
 *      if \p qureg is not a density matrix, 
 *      or if \p target is outside of [0, \p qureg.numQubitsRepresented),
 *      or if \p numOps is outside [1, 4],
 *      or if \p ops do not create a completely positive, trace preserving map,
 *      or if a node cannot fit 4 amplitudes in distributed mode.
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
 * @ingroup decoherence
 * @param[in,out] qureg the density matrix to which to apply the map
 * @param[in] target1 the least significant target qubit in \p ops
 * @param[in] target2 the most significant target qubit in \p ops
 * @param[in] ops an array of at most 16 Kraus operators
 * @param[in] numOps the number of operators in \p ops which must be >0 and <= 16.
 * @throws exitWithError
 *      if \p qureg is not a density matrix, 
 *      or if either \p target1 or \p target2 is outside of [0, \p qureg.numQubitsRepresented),
 *      or if \p target1 = \p target2,
 *      or if \p numOps is outside [1, 16],
 *      or if \p ops do not create a completely positive, trace preserving map,
 *      or if a node cannot fit 16 amplitudes in distributed mode.
 * @author Balint Koczor
 * @author Tyson Jones (refactored, doc)
 */
void mixTwoQubitKrausMap(Qureg qureg, int target1, int target2, ComplexMatrix4 *ops, int numOps);

/** Apply a general N-qubit Kraus map to a density matrix, as specified by at most (2N)^2
 * Kraus operators. A Kraus map is also referred to as a "operator-sum representation"
 * of a quantum channel. This allows one to simulate a general two-qubit noise process.
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
 * @ingroup decoherence
 * @param[in,out] qureg the density matrix to which to apply the map
 * @param[in] targets a list of target qubit indices, the first of which is treated as least significant in each op in \p ops
 * @param[in] numTargets the length of \p targets
 * @param[in] ops an array of at most (2N)^2 Kraus operators
 * @param[in] numOps the number of operators in \p ops which must be >0 and <= (2N)^2.
 * @throws exitWithError
 *      if \p qureg is not a density matrix, 
 *      or if any target in \p targets is outside of [0, \p qureg.numQubitsRepresented),
 *      or if any qubit in \p targets is repeated,
 *      or if \p numOps is outside [1, (2 \p numTargets)^2],
 *      or if any ComplexMatrixN in \ops does not have op.numQubits == \p numTargets,
 *      or if \p ops do not create a completely positive, trace preserving map,
 *      or if a node cannot fit (2N)^2 amplitudes in distributed mode.
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
 * @ingroup calc
 * @param[in] a a density matrix
 * @param[in] b an equally-sized density matrix
 * @throws exitWithError
 *      if either \p a or \p b are not density matrices,
 *      or if \p a and \p have mismatching dimensions.
 * @author Balint Koczor
 * @author Tyson Jones (refactored, doc)
 */
qreal calcHilbertSchmidtDistance(Qureg a, Qureg b);

/** Modifies qureg \p out to the result of (\p facOut \p out + \p fac1 \p qureg1 + \p fac2 \p qureg2), 
 * imposing no constraints on normalisation. Works for both statevectors and density matrices.
 * Note that afterward, \p out may not longer be normalised and ergo no longer a valid 
 * statevector or density matrix. Users must therefore be careful passing \p out to
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
 * @throws exitWithError
 *      if \p qureg1, \p qureg2 and \p aren't all state-vectors or all density matrices,
 *      or if the dimensions of \p qureg1, \p qureg2 and \p aren't equal
 * @author Tyson Jones
 */
void setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out);

/** Modifies \p outQureg to be the result of applying the weighted sum of Pauli products (a Hermitian but not 
 * necessarily unitary operator) to \p inQureg. Note that afterward, \p outQureg may not longer be normalised and ergo not a
 * statevector or density matrix. Users must therefore be careful passing \p outQureg to
 * other QuEST functions which assume normalisation in order to function correctly.

 * Letting \f$ \alpha = \sum_i c_i \otimes_j^{N} \hat{\sigma}_{i,j} \f$ be 
 * the operators indicated by \p allPauliCodes (where \f$ c_i \in \f$ \p termCoeffs and \f$ N = \f$ \p qureg.numQubitsRepresented), 
 * this function effects \f$ \alpha | \psi \rangle \f$ on statevector \f$ |\psi\rangle \f$
 * and $\alpha \rho$ (matrix multiplication) on density matrix \f$ \rho \f$.
 *
 * \p allPauliCodes is an array of length \p numSumTerms*\p qureg.numQubitsRepresented
 * which specifies which Pauli operators to apply, where 0 = \p PAULI_I, 1 = \p PAULI_X, 
 * 2 = \p PAULI_Y, 3 = \p PAULI_Z. For each sum term, a Pauli operator must be specified for 
 * EVERY qubit in \p qureg; each set of \p numSumTerms operators will be grouped into a product.
 * \p termCoeffs is an arrray of length \p numSumTerms containing the term coefficients.
 * For example, on a 3-qubit statevector,
 *
 *     int paulis[6] = {PAULI_X, PAULI_I, PAULI_I,  PAULI_X, PAULI_Y, PAULI_Z};
 *     qreal coeffs[2] = {1.5, -3.6};
 *     applyPauliSum(inQureg, paulis, coeffs, 2, outQureg);
 *
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
 * @throws exitWithError
 *      if any code in \p allPauliCodes is not in {0,1,2,3},
 *      or if numSumTerms <= 0,
 *      or if \p inQureg is not of the same type and dimensions as \p outQureg
 * @author Tyson Jones
 */
void applyPauliSum(Qureg inQureg, enum pauliOpType* allPauliCodes, qreal* termCoeffs, int numSumTerms, Qureg outQureg);
 
#ifndef __cplusplus
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
 *
 *     ComplexMatrixN m = bindArraysToStackComplexMatrixN(
 *         1, 
 *         (qreal[][2]) {{1,0},{0,1}}, (qreal[][2]) {{0}}, 
 *         (qreal*[2]) {0}, (qreal*[2]) {0}
 *     );
 *
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
 *
 *     ComplexMatrixN u = getStaticComplexMatrixN(1, ({{1,2},{3,4}}), ({{0}}));
 *
 * Here is an example of an incorrect usage, since a 'local' ComplexMatrixN cannot
 * leave the calling scope (otherwise inducing dangling pointers):
 *
 *     ComplexMatrixN getMyMatrix(void) {
 *         return getStaticComplexMatrixN(1, ({{1,2},{3,4}}), ({{0}}));
 *     }
 * 
 * This function is actually a single-line anonymous macro, so can be safely 
 * invoked within arguments to other functions, e.g.
 *
 *      multiQubitUnitary(
 *          qureg, (int[]) {0}, 1, 
 *          getStaticComplexMatrixN(1, ({{1,0},{0,1}}), ({{0}}))
 *      );
 *
 * The returned ComplexMatrixN can be accessed and modified in the same way as
 * that returned by createComplexMatrixN(), e.g.
 *
 *      ComplexMatrixN u = getStaticComplexMatrixN(3, ({{0}}), ({{0}}));
 *      for (int i=0; i<8; i++)
 *          for (int j=0; j<8; j++)
 *              u.real[i][j] = .1;
 *
 * Note that the first argument \p numQubits must be a literal.
 *
 * This macro is only callable in C, since it invokes the function 
 * bindArraysToStackComplexMatrixN() which is only callable in C.
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


// end prevention of C++ name mangling
#ifdef __cplusplus
}
#endif

#endif // QUEST_H

