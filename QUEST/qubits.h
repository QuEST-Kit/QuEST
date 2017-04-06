# ifndef QUBITS
# define QUBITS
/** @file
 * The QUEST library API and objects. 
*/

/** Represents an array of complex numbers grouped into an array of real components and an array of coressponding complex components.
*/
typedef struct ComplexArray
{
	double *real; 
	double *imag;
} ComplexArray;

/** Represents one complex number.
*/
typedef struct Complex
{
	double real;
	double imag;
} Complex;

/** Represents a system of qubits.
Qubits are zero-based and the the first qubit is the rightmost
*/
typedef struct MultiQubit
{
	//! Probablilty amplitudes for the multi qubit state
	ComplexArray stateVec; 
	//! Temporary storage for a chunk of the state vector received from another process in the MPI version
	ComplexArray pairStateVec;
	//! Number of qubits in the state
	int numQubits;
	//! Number of probability amplitudes held in stateVec by this process
	//! In the non-MPI version, this is the total number of amplitudes
	long long int numAmps;
	//! The position of the chunk of the state vector held by this process in the full state vector
	int chunkId;
	//! Number of chunks the state vector is broken up into -- the number of MPI processes used
	int numChunks;
} MultiQubit;

/** Information about the environment the program is running in.
In practice, this holds info about MPI ranks and helps to hide MPI initialization code
*/
typedef struct QUESTEnv
{
	int rank;
	int numRanks;
} QUESTEnv;


// QUEST library functions whose implementation is independent of environment (local, MPI)

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QUESTEnv env);

void destroyMultiQubit(MultiQubit multiQubit, QUESTEnv env);

void reportState(MultiQubit multiQubit);

void initStateVec(MultiQubit *multiQubit);

// QUEST library functions whose implementation depends on environment (local, MPI)

/** Initialize QUEST environment. If something needs to be done to set up the execution environment, such as 
 * initializing MPI when running in distributed mode, it is handled here
 * @param[in,out] env object representing the execution environment. A single instance is used for each program
 */
void initQUESTEnv(QUESTEnv *env);

/** Close QUEST environment. If something needs to be done to clean up the execution environment, such as 
 * finalizing MPI when running in distributed mode, it is handled here
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void closeQUESTEnv(QUESTEnv env);

/** Guarantees that all code up to the given point has been executed on all nodes. 
 */
void syncQUESTEnv(QUESTEnv env);

/** Calculate the probability of being in any state by taking the norm of the entire state vector. 
 * Should be equal to 1.
 * @param[in] multiQubit object representing a set of qubits
 * @return total probability
 */
double calcTotalProbability(MultiQubit multiQubit);

/** Rotate a single qubit in the state vector of probability amplitudes, given the angle rotation arguments.
alphaRe = cos(angle1) * cos(angle2) \n
alphaIm = cos(angle1) * sin(angle2) \n            
betaRe  = sin(angle1) * cos(angle3) \n            
betaIm  = sin(angle1) * sin(angle3) \n           

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits to be initialised
@param[in] rotQubit qubit to rotate
@param[in] alpha rotation angle
@param[in] beta rotation angle
 */
void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta);

/** Measure the probability
of a specified qubit being in the zero state.     

@param[in] multiQubit object representing the set of qubits to be initialised
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/
double findProbabilityOfZero(MultiQubit multiQubit, const int measureQubit);

# endif
