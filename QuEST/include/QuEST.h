// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * The QuEST API.
 * This file contains the comments used by doxygen for generating API doc
*/

# ifndef QUEST_H
# define QUEST_H

# include "QuEST_precision.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * private structures
 */
    
// hide these from doxygen
/// \cond HIDDEN_SYMBOLS    

// Codes for Z-axis phase gate variations
enum phaseGateType {SIGMA_Z=0, S_GATE=1, T_GATE=2};    
    
/** A logger of QASM instructions */
typedef struct {
    
    char* buffer;       // generated QASM string
    int bufferSize;     // maximum number of chars before overflow
    int bufferFill;     // number of chars currently in buffer
    int isLogging;      // whether gates are being added to buffer
    
} QASMLogger;

/** Represents an array of complex numbers grouped into an array of 
 * real components and an array of coressponding complex components.
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

/** Represents one complex number.
 */
typedef struct Complex
{
    qreal real;
    qreal imag;
} Complex;

/** Represents a 2x2 matrix of complex numbers
 */
typedef struct ComplexMatrix2
{
    Complex r0c0, r0c1;
    Complex r1c0, r1c1;
} ComplexMatrix2;

/** Represents a 3-vector of real numbers
 */
typedef struct Vector
{
    qreal x, y, z;
} Vector;

/** Represents a system of qubits.
 * Qubits are zero-based
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
In practice, this holds info about MPI ranks and helps to hide MPI initialization code
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
 * @returns an object representing the set of qubits
 * @param[in] numQubits number of qubits in the system
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @throws exitWithError if \p numQubits <= 0
 */
Qureg createQureg(int numQubits, QuESTEnv env);

/** Create a Qureg for qubits which are represented by a density matrix, and can be in mixed states.
 * Allocate space for a density matrix of probability amplitudes, including space for temporary values to be copied from
 * one other chunk if running the distributed version. Define properties related to the size of the set of qubits.
 * initZeroState should be called after this to initialise the qubits to the zero pure state.
 *
 * @returns an object representing the set of qubits
 * @param[in] numQubits number of qubits in the system
 * @param[in] env object representing the execution environment (local, multinode etc)
 * @throws exitWithError if \p numQubits <= 0
 */
Qureg createDensityQureg(int numQubits, QuESTEnv env);

/** Deallocate a Qureg object representing a set of qubits.
 * Free memory allocated to state vector of probability amplitudes, including temporary vector for
 * values copied from another chunk if running the distributed version.
 *
 * @param[in,out] qureg object to be deallocated
 * @param[in] env object representing the execution environment (local, multinode etc)
 */
void destroyQureg(Qureg qureg, QuESTEnv env);

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

 * @param[in,out] qureg object representing the set of qubits
 */
void reportState(Qureg qureg);

/** Print the current state vector of probability amplitudes for a set of qubits to standard out. 
 * For debugging purposes. Each rank should print output serially. 
 * Only print output for systems <= 5 qubits
 */
void reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank);

/** Report metainformation about a set of qubits: number of qubits, number of probability amplitudes.

 * @param[in] qureg object representing the set of qubits
 */
void reportQuregParams(Qureg qureg);

/** Get the number of qubits in a qureg object
 */
int getNumQubits(Qureg qureg);

/** Get the number of probability amplitudes in a qureg object, given by 2^numQubits
 */
int getNumAmps(Qureg qureg);

/**
 * Initialise a set of \f$ N \f$ qubits to the classical zero state 
 * \f$ {| 0 \rangle}^{\otimes N} \f$.
 *
 * @param[in,out] qureg the object representing the set of all qubits to initialise
 */
void initZeroState(Qureg qureg);

/** Initialise a set of \f$ N \f$ qubits to the plus state
 * \f$ {| + \rangle}^{\otimes N} = \frac{1}{\sqrt{2^N}} (| 0 \rangle + | 1 \rangle)^{\otimes N} \f$.
 * This is the product state of \f$N\f$ qubits where every classical state is uniformly 
 * populated with real coefficient \f$\frac{1}{\sqrt{2^N}}\f$.
 * This is equivalent to applying a Hadamard to every qubit in the zero state: 
 * \f$ \hat{H}^{\otimes N} {|0\rangle}^{\otimes N} \f$
 *
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 */
void initPlusState(Qureg qureg);

/** Initialise a set of \f$ N \f$ qubits to the classical state with index \p stateInd.
 * Note \f$ | 00 \dots 00 \rangle \f$ has \p stateInd 0, \f$ | 00 \dots 01 \rangle \f$ has \p stateInd 1, 
 * \f$ | 11 \dots 11 \rangle \f$ has \p stateInd \f$ 2^N - 1 \f$, etc.
 * Subsequent calls to getProbAmp will yield 0 for all indices except \p stateInd.
 *
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @param[in] stateInd the index (0 to the number of amplitudes, exclusive) of the state to give probability 1
 */
void initClassicalState(Qureg qureg, long long int stateInd);

/** Initialise a set of \f$ N \f$ qubits, which can be a state vector or density matrix, to a given pure state.
 * If \p qureg is a state-vector, this merely makes \p qureg an identical copy of \p pure.
 * If \p qureg is a density matrix, this makes \p qureg 100% likely to be in the \p pure state.
 *
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @param[in] pure the pure state to be copied or to give probability 1 in qureg
 */
void initPureState(Qureg qureg, Qureg pure);

/** Initialise qureg by specifying the complete statevector.
 * The real and imaginary components of the amplitudes are passed in separate arrays,
 * each of which must have length \p qureg.numAmpsTotal.
 * There is no automatic checking that the passed arrays are L2 normalised.
 *
 * @param[in,out] qureg the object representing the set of qubits to be initialised
 * @param[in] reals array of the real components of the new amplitudes
 * @param[in] imags array of the imaginary components of the new amplitudes
 * @throws exitWithError
 *      if \p qureg is not a statevector (i.e. is a density matrix)
 */
void initStateFromAmps(Qureg qureg, qreal* reals, qreal* imags);

/** Edit qureg to adopt the subset of amplitudes passed in \p reals and \p imags.
 * Only amplitudes with indices in [\p startInd, \p startInd + \p numAmps] will be changed, which means
 * the new state may not be L2 normalised. This allows the user to initialise a custom state by 
 * setting batches of amplitudes.
 *
 * @param[in,out] qureg the statevector to modify
 * @param[in] startInd the index of the first amplitude in \p qureg's statevector to modify
 * @param[in] reals array of the real components of the new amplitudes
 * @param[in] imags array of the imaginary components of the new amplitudes
 * @param[in] numAmps the length of each of the reals and imags arrays.
 * @throws exitWithError
 *      if \p qureg is not a statevector (i.e. is a density matrix)
 */
void setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps);

/** Set targetQureg to be a clone of copyQureg. 
 * Registers must either both be state-vectors, or both be density matrices.
 * Only the quantum state is cloned, auxilary info (like recorded QASM) is unchanged.
 * copyQureg is unaffected.
 *
 * @param[in, out] targetQureg the qureg to have its quantum state overwritten
 * @param[in] copyQureg the qureg to have its quantum state cloned in targetQureg.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to undergo a phase shift
 * @param[in] angle amount by which to shift the phase in radians
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] idQubit1 first qubit in the state to phase shift
 * @param[in] idQubit2 second qubit in the state to phase shift
 * @param[in] angle amount by which to shift the phase in radians
 * @throws exitWithError
 *  if \p idQubit1 or \p idQubit2 are outside [0, \p qureg.numQubitsRepresented), or are equal
 */
void controlledPhaseShift(Qureg qureg, const int idQubit1, const int idQubit2, qreal angle);

/** Introduce a phase factor \f$ \exp(i \theta) \f$ on state \f$ |1 \dots 1 \rangle \f$
 * of the passed qubits.
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
                \node[draw=none] at (0, 0) {$R_\theta$};
                \end{tikzpicture}
    }
    \f]
 * 
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits array of qubits to phase shift
 * @param[in] numControlQubits the length of array \p controlQubits
 * @param[in] angle amount by which to shift the phase in radians
 * @throws exitWithError
 *      if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit index in \p controlQubits is outside
 *      [0, \p qureg.numQubitsRepresented])
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] idQubit1, idQubit2 qubits to operate upon
 * @throws exitWithError 
 *  if \p idQubit1 or \p idQubit2 are outside [0, \p qureg.numQubitsRepresented), or are equal
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits array of input qubits
 * @param[in] numControlQubits number of input qubits
 * @throws exitWithError 
 *      if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented) 
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate upon
 * @throws exitWithError if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate upon
 * @throws exitWithError if \p targetQubit is outside [0, \p qureg.numQubitsRepresented)
 */
void tGate(Qureg qureg, const int targetQubit);

/** Create the QuEST execution environment.
 * This should be called only once, and the environment should be freed with destroyQuESTEnv at the end
 * of the user's code.
* If something needs to be done to set up the execution environment, such as 
 * initializing MPI when running in distributed mode, it is handled here.
 *
 * @return object representing the execution environment. A single instance is used for each program
 */
QuESTEnv createQuESTEnv(void);

/** Destroy the QuEST environment. 
 * If something needs to be done to clean up the execution environment, such as 
 * finalizing MPI when running in distributed mode, it is handled here
 *
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void destroyQuESTEnv(QuESTEnv env);

/** Guarantees that all code up to the given point has been executed on all nodes (if running in distributed mode)
 *
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void syncQuESTEnv(QuESTEnv env);

/** Performs a logical AND on all successCodes held by all processes. If any one process has a zero successCode
 * all processes will return a zero success code.
 *
 * @param[in] successCode 1 if process task succeeded, 0 if process task failed
 * @returns 1 if all processes succeeded, 0 if any one process failed
 */ 
int syncQuESTSuccess(int successCode);

/** Report information about the QuEST environment
 *
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void reportQuESTEnv(QuESTEnv env);

void getEnvironmentString(QuESTEnv env, Qureg qureg, char str[200]);

/** Get the complex amplitude at a given index in the state vector.
 *
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return amplitude at index, returned as a Complex struct (with .real and .imag attributes)
 * @throws exitWithError
 *      if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 */
Complex getAmp(Qureg qureg, long long int index);

/** Get the real component of the complex probability amplitude at an index in the state vector.
 * For debugging purposes.
 *
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return real component at that index
 * @throws exitWithError
 *      if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 */
qreal getRealAmp(Qureg qureg, long long int index);

/** Get the imaginary component of the complex probability amplitude at an index in the state vector.
 * For debugging purposes.
 *
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return imaginary component at that index
 * @throws exitWithError
 *      if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 */
qreal getImagAmp(Qureg qureg, long long int index);

/** Get the probability of the state at an index in the full state vector.
 *
 * @param[in] qureg object representing a set of qubits
 * @param[in] index index in state vector of probability amplitudes
 * @return realEl*realEl + imagEl*imagEl
 * @throws exitWithError
 *      if \p index is outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 */
qreal getProbAmp(Qureg qureg, long long int index);

/** Get an amplitude from a density matrix at a given row and column.
 *
 * @param[in] qureg object representing a density matrix
 * @param[in] row row of the desired amplitude in the density matrix
 * @param[in] col column of the desired amplitude in the density matrix
 * @return a Complex scalar representing the desired amplitude
 * @throws exitWithError
 *      if \p qureg is a statevector, or if \p row or \p col are outside [0, \f$2^{N}\f$) where \f$N = \f$ \p qureg.numQubitsRepresented
 */
Complex getDensityAmp(Qureg qureg, long long int row, long long int col);

/** A debugging function which calculates the probability of being in any state, which should always be 1.
 * For pure states, this is the norm of the entire state vector and for mixed states, is the trace of
 * the density matrix.
 * Note this calculation utilises Kahan summation for greaster accuracy, but is
 * not parallelised and so will be slower than other functions.
 *
 * @param[in] qureg object representing a set of qubits
 * @return total probability
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] alpha complex unitary parameter (row 1, column 1)
 * @param[in] beta complex unitary parameter (row 2, column 1)
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p alpha, \p beta don't satisfy |\p alpha|^2 + |\p beta|^2 = 1.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] u unitary matrix to apply
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or matrix \p u is not unitary.
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
 *
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented).
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented).
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented).
 */
void rotateZ(Qureg qureg, const int rotQubit, qreal angle);

/** Rotate a single qubit by a given angle around a given vector on the Bloch-sphere.      
 * The vector must not be zero (else an error is thrown), but needn't be unit magnitude.
 *
 * For angle \f$\theta\f$ and axis vector \f$\vec{n}\f$, applies \f$R_{\hat{n}} = \exp \left(- i \frac{\theta}{2} \hat{n} \cdot \vec{\sigma} \right) \f$
 * where \f$\vec{\sigma}\f$ is the vector of Pauli matrices.
 *
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] rotQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @param[in] axis vector around which to rotate (can be non-unit; will be normalised)
 * @throws exitWithError
 *      if \p rotQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p axis is the zero vector
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit which has value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate the target qubit in radians
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit qubit with value 1 in the rotated states
 * @param[in] targetQubit qubit to rotate
 * @param[in] angle angle by which to rotate in radians
 * @param[in] axis vector around which to rotate (can be non-unit; will be normalised)
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal
 *      or if \p axis is the zero vector
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit apply the target unitary if this qubit has value 1
 * @param[in] targetQubit qubit on which to apply the target unitary
 * @param[in] alpha complex unitary parameter (row 1, column 1)
 * @param[in] beta complex unitary parameter (row 2, column 1)
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal,
 *      or if \p alpha, \p beta don't satisfy |\p alpha|^2 + |\p beta|^2 = 1.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit apply unitary if this qubit is 1
 * @param[in] targetQubit qubit to operate on
 * @param[in] u single-qubit unitary matrix to apply
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented) or are equal,
 *      or if \p u is not unitary.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubits applies unitary if all qubits in this array equal 1
 * @param[in] numControlQubits number of control qubits
 * @param[in] targetQubit qubit to operate on
 * @param[in] u single-qubit unitary matrix to apply
 * @throws exitWithError
 *      if \p numControlQubits is outside [1, \p qureg.numQubitsRepresented]),
 *      or if any qubit index (\p targetQubit or one in \p controlQubits) is outside
 *      [0, \p qureg.numQubitsRepresented]), 
 *      or if \p controlQubits contains \p targetQubit,
 *      or if \p u is not unitary.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] targetQubit qubit to operate on
 * @throws exitWithError
 *      if \p targetQubit is outside [0, \p qureg.numQubitsRepresented).
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit nots the target if this qubit is 1
 * @param[in] targetQubit qubit to not
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented), or are equal.
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
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] controlQubit applies pauliY to the target if this qubit is 1
 * @param[in] targetQubit qubit to not
 * @throws exitWithError
 *      if either \p controlQubit or \p targetQubit are outside [0, \p qureg.numQubitsRepresented), or are equal.
 */
void controlledPauliY(Qureg qureg, const int controlQubit, const int targetQubit);

/** Gives the probability of a specified qubit being measured in the given outcome (0 or 1).
 * This performs no actual measurement and does not change the state of the qubits.
 * 
 * @param[in] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to study
 * @param[in] outcome for which to find the probability of the qubit being measured in
 * @return probability of qubit measureQubit being measured in the given outcome
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p outcome is not in {0, 1}.
 */
qreal calcProbOfOutcome(Qureg qureg, const int measureQubit, int outcome);

/** Updates the state vector to be consistent with measuring the measure qubit in the given outcome (0 or 1), and returns the probability of such a measurement outcome. 
 * This is effectively performing a measurement and forcing the outcome.
 * This is an irreversible change to the state vector, whereby incompatible states
 * in the state vector are given zero amplitude and the remaining states are renormalised.
 * Exits with error if the given outcome has ~zero probability, and so cannot be
 * collapsed into.
 * 
 * @param[in,out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @param[in] outcome to force the measure qubit to enter
 * @return probability of the (forced) measurement outcome
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p outcome is not in {0, 1},
 *      or if the probability of \p outcome is zero (within machine epsilon)
 */
qreal collapseToOutcome(Qureg qureg, const int measureQubit, int outcome);

/** Measures a single qubit, collapsing it randomly to 0 or 1.
 * Outcome probabilities are weighted by the state vector, which is irreversibly
 * changed after collapse to be consistent with the outcome.
 *
 * @param[in, out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @return the measurement outcome, 0 or 1
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 */
int measure(Qureg qureg, int measureQubit);

/** Measures a single qubit, collapsing it randomly to 0 or 1, and
 * additionally gives the probability of that outcome.
 * Outcome probabilities are weighted by the state vector, which is irreversibly
 * changed after collapse to be consistent with the outcome.
 *
 * @param[in, out] qureg object representing the set of all qubits
 * @param[in] measureQubit qubit to measure
 * @param[out] outcomeProb a pointer to a qreal which is set to the probability of the occurred outcome
 * @return the measurement outcome, 0 or 1
 * @throws exitWithError
 *      if \p measureQubit is outside [0, \p qureg.numQubitsRepresented)
 */
int measureWithStats(Qureg qureg, int measureQubit, qreal *outcomeProb);

/** Computes <bra|ket> */
Complex calcInnerProduct(Qureg bra, Qureg ket);

/** Seed the Mersenne Twister used for random number generation in the QuEST environment with an example
 * defualt seed.
 * This default seeding function uses the mt19937 init_by_array function with two keys -- 
 * time and pid. Subsequent calls to mt19937 genrand functions will use this seeding. 
 * For a multi process code, the same seed is given to all process, therefore this seeding is only
 * appropriate to use for functions such as measure where all processes require the same random value.
 *
 * For more information about the MT, see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
 **/
void seedQuESTDefault(void);

/** Seed the Mersenne Twister used for random number generation in the QuEST environment with
 * a user defined seed.
 * This function uses the mt19937 init_by_array function with numSeeds keys supplied by the user.
 * Subsequent calls to mt19937 genrand functions will use this seeding. 
 * For a multi process code, the same seed is given to all process, therefore this seeding is only
 * appropriate to use for functions such as measure where all processes require the same random value.
 *
 * @param[in] seedArray Array of integers to use as seed. This allows the MT to be initialised with more 
 * than a 32-bit integer if required
 * @param[in] numSeeds Length of seedArray
 *
 * For more information about the MT, see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
 **/
void seedQuEST(unsigned long int *seedArray, int numSeeds);

/** Enable QASM recording. Gates applied to qureg will here-after be added to growing QASM instructions,
 * progressively consuming more memory until stopped
 */
void startRecordingQASM(Qureg qureg);

/** Disable QASM recording. The recorded QASM will be maintained in qureg
and continue to be
 * added to if startRecordingQASM is recalled.
 */
void stopRecordingQASM(Qureg qureg);

/** Clear all QASM so far recorded. This does not start or stop recording
 */
void clearRecordedQASM(Qureg qureg);

/** Print recorded QASM to stdout */
void printRecordedQASM(Qureg qureg);

/** Writes recorded QASM to a file, throwing an error if inaccessible */
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
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce dephasing noise
 * @param[in] prob the probability of the phase error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p prob is not in [0, 1/2]
 */
void applyOneQubitDephaseError(Qureg qureg, const int targetQubit, qreal prob);

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
 * @param[in,out] qureg a density matrix
 * @param[in] qubit1 qubit upon which to induce dephasing noise
 * @param[in] qubit2 qubit upon which to induce dephasing noise
 * @param[in] prob the probability of the phase error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if either \p qubit1 or \p qubit2 is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p qubit1 = \p qubit2,
 *      or if \p prob is not in [0, 3/4]
 */
void applyTwoQubitDephaseError(Qureg qureg, const int qubit1, const int qubit2, qreal prob);

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
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p prob is not in [0, 3/4]
 */
void applyOneQubitDepolariseError(Qureg qureg, const int targetQubit, qreal prob);

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
 * @param[in,out] qureg a density matrix
 * @param[in] targetQubit qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if \p targetQubit is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p prob is not in [0, 1]
 */
void applyOneQubitDampingError(Qureg qureg, const int targetQubit, qreal prob);

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
 * @param[in,out] qureg a density matrix
 * @param[in] qubit1 qubit upon which to induce depolarising noise
 * @param[in] qubit2 qubit upon which to induce depolarising noise
 * @param[in] prob the probability of the depolarising error occuring
 * @throws exitWithError
 *      if \p qureg is not a density matrix,
 *      or if either \p qubit1 or \p qubit2 is outside [0, \p qureg.numQubitsRepresented),
 *      or if \p qubit1 = \p qubit2,
 *      or if \p prob is not in [0, 15/16]
 */
void applyTwoQubitDepolariseError(Qureg qureg, const int qubit1, const int qubit2, qreal prob);

/** Modifies combineQureg to become (1-prob)combineProb + prob otherQureg.
 * Both registers must be equal-dimension density matrices, and prob must be in [0, 1].
 *
 * @param[in,out] combineQureg a density matrix to be modified
 * @param[in] prob the probability of \p otherQureg in the modified \p combineQureg
 * @param[in] otherQureg a density matrix to be mixed into \p combineQureg
 * @throws exitWithError
 *      if either \p combineQureg or \p otherQureg are not density matrices,
 *      or if the dimensions of \p combineQureg and \p otherQureg do not match,
 *      or if \p prob is not in [0, 1]
 */
void addDensityMatrix(Qureg combineQureg, qreal prob, Qureg otherQureg);

/** Calculates the purity of a density matrix, by the trace of the density matrix squared.
 * Returns \f$\text{Tr}(\rho^2)\f$.
 * For a pure state, this =1.
 * For a mixed state, the purity is less than 1 and is lower bounded by 1/2^n, where
 * n is the number of qubits. The minimum purity is achieved for the maximally mixed state identity/2^n.
 *
 * This function does not accept state-vectors, which clearly have purity 1.
 *
 * @param[in] qureg a density matrix of which to measure the purity
 * @return the purity
 * @throws exitWithError
 *      if either \p combineQureg or \p otherQureg are not density matrices,
 *      or if the dimensions of \p combineQureg and \p otherQureg do not match,
 *      or if \p prob is not in [0, 1]
 */
qreal calcPurity(Qureg qureg);

/** Calculates the fidelity of qureg (a statevector or density matrix) against 
 * a reference pure state (necessarily a statevector).
 * For two pure states, this is |<qureg|pureState>|^2.
 * For a mixed and pure state, this is <pureState|qureg|pureState>.
 * In either case, the fidelity lies in [0, 1].
 * The number of qubits represented in \p qureg and \p pureState must match.
 * 
 * @param[in] qureg a density matrix or state vector
 * @param[in] pureState a state vector
 * @return the fidelity between the input registers
 * @throws exitWithError
 *      if the dimensions \p qureg and \p pureState do not match
 */
qreal calcFidelity(Qureg qureg, Qureg pureState);


#ifdef __cplusplus
}
#endif

#endif // QUEST_H

