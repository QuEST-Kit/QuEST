# ifndef QUBITS
# define QUBITS

# include "precision.h"
/** @file
 * The QuEST library API and objects. 
*/

/** Represents an array of complex numbers grouped into an array of real components and an array of coressponding complex components.
*/
typedef struct ComplexArray
{
	REAL *real; 
	REAL *imag;
} ComplexArray;

/** Represents one complex number.
*/
typedef struct Complex
{
	REAL real;
	REAL imag;
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
typedef struct QuESTEnv
{
	int rank;
	int numRanks;
} QuESTEnv;

// Codes for sigmaZ phase gate variations
enum phaseGateType {SIGMA_Z=0, S_GATE=1, T_GATE=2};

// QuEST library functions whose implementation is independent of environment (local, MPI)

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env);

void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env);

void reportState(MultiQubit multiQubit);

void reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank);

void reportMultiQubitParams(MultiQubit multiQubit);

void initStateZero(MultiQubit *multiQubit);

void initStatePlus(MultiQubit *multiQubit);

void initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome);

void initStateDebug(MultiQubit *multiQubit);

void initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env);

int compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision);

void quadCPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2, 
	const int idQubit3, const int idQubit4);

void controlPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2);

void sigmaZ(MultiQubit multiQubit, const int rotQubit);

void sGate(MultiQubit multiQubit, const int rotQubit);

void tGate(MultiQubit multiQubit, const int rotQubit);


// QuEST library functions whose implementation depends on environment (local, MPI)

/** Initialize QuEST environment. If something needs to be done to set up the execution environment, such as 
 * initializing MPI when running in distributed mode, it is handled here
 * @param[in,out] env object representing the execution environment. A single instance is used for each program
 */
void initQuESTEnv(QuESTEnv *env);

/** Close QuEST environment. If something needs to be done to clean up the execution environment, such as 
 * finalizing MPI when running in distributed mode, it is handled here
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void closeQuESTEnv(QuESTEnv env);

/** Guarantees that all code up to the given point has been executed on all nodes. 
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void syncQuESTEnv(QuESTEnv env);

/** Performs a logical AND on all successCodes held by all processes. If any one process has a zero successCode
 * all processes will return a zero success code.
 * @param[in] env object representing the execution environment. A single instance is used for each program
 * @param[in] successCode 1 if process task succeeded, 0 if process task failed
 * @returns 1 if all processes succeeded, 0 if any one process failed
 */ 
int syncQuESTSuccess(QuESTEnv env, int successCode);

/** Report information about the QuEST environment
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void reportQuESTEnv(QuESTEnv env);

/** Report a list of CPU hostnames and the rank that is running on each if running with MPI enabled and an 
error message otherwise. For debugging purposes. 
@param[in] env object representing the execution environment. A single instance is used for each program
*/
void reportNodeList(QuESTEnv env);

void getEnvironmentString(QuESTEnv env, MultiQubit multiQubit, char str[200]);


/** Get the real component of the probability amplitude at an index in the state vector.
For debugging purposes.
@param[in] multiQubit object representing a set of qubits
@param[in] index index in state vector of probability amplitudes
@return real component at that index
*/
REAL getRealAmpEl(MultiQubit multiQubit, long long int index);

/** Get the imaginary component of the probability amplitude at an index in the state vector.
For debugging purposes.
@param[in] multiQubit object representing a set of qubits
@param[in] index index in state vector of probability amplitudes
@return imaginary component at that index
*/
REAL getImagAmpEl(MultiQubit multiQubit, long long int index);

/** Get the probability of the state at an index in the state vector.
@param[in] multiQubit object representing a set of qubits
@param[in] index index in state vector of probability amplitudes
@return realEl*realEl + imagEl*imagEl
*/
REAL getProbEl(MultiQubit multiQubit, long long int index);

/** Calculate the probability of being in any state by taking the norm of the entire state vector. 
 * Should be equal to 1.
 * @param[in] multiQubit object representing a set of qubits
 * @return total probability
 */
REAL calcTotalProbability(MultiQubit multiQubit);

/** Rotate a single qubit in the state vector of probability amplitudes, given the angle rotation arguments.
alphaRe = cos(angle1) * cos(angle2) \n
alphaIm = cos(angle1) * sin(angle2) \n            
betaRe  = sin(angle1) * cos(angle3) \n            
betaIm  = sin(angle1) * sin(angle3) \n           

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] alpha rotation angle
@param[in] beta rotation angle
 */
void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta);

/** Rotate a single qubit in the state vector of probability amplitudes, given the angle rotation arguments and a control qubit. Only perform 
the rotation for elements where the control qubit is 1. 
alphaRe = cos(angle1) * cos(angle2) \n
alphaIm = cos(angle1) * sin(angle2) \n            
betaRe  = sin(angle1) * cos(angle3) \n            
betaIm  = sin(angle1) * sin(angle3) \n           

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] controlQubit perform rotation if this qubit is 1
@param[in] alpha rotation angle
@param[in] beta rotation angle
 */
void controlRotateQubit(MultiQubit multiQubit, const int rotQubit, const int controlQubit, Complex alpha, Complex beta);

/** Rotate a single qubit by {{0,1},{1,0}} -- swap |0> and |1>.
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
*/
void sigmaX(MultiQubit multiQubit, const int rotQubit);

/** Rotate a single qubit by {{0,-i},{i,0}} -- swap |0> and |1> and apply
a phase of -i or i.
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
*/
void sigmaY(MultiQubit multiQubit, const int rotQubit);

/** Rotate a single qubit by {{1,0},{{0,-1}} -- apply a phase of -1 to |1>.
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
*/
void sigmaZ(MultiQubit multiQubit, const int rotQubit);

/** Rotate a single qubit by {{1,1},{1,-1}}/sqrt2 -- turn a |0> into a |+>
and a |1> into a |->.
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
*/
void hadamard(MultiQubit multiQubit, const int rotQubit);

/** Rotate a single qubit by {{0,-i},{i,0}} -- swap |0> and |1> and apply
a phase of -i or i, only for elements when control qubit is 1.
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] controlQubit perform sigmaX rotation if this qubit is 1
*/
void controlNot(MultiQubit multiQubit, const int targetQubit, const int controlQubit);

/** Measure the probability
of a specified qubit being in the zero or one state.     

@param[in] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@param[in] outcome to measure the probability of -- either zero or one
@return probability of qubit measureQubit being either zero or one
*/
REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome);

/** Update the state vector to be consistent with measuring measureQubit=0 or measureQubit=1 according to the value 
of outcome. 
Measure in Zero performs an irreversible change to the state vector: it updates the vector according
to the event that an outcome has been measured on the qubit indicated by measureQubit (where 
his label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
then renormalising based on the total probability of measuring measureQubit=0 if outcome=0 or
measureQubit=1 if outcome=1. It then returns the probability of making this measurement. 

@param[in,out] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@param[in] outcome to measure the probability of and set the state to -- either zero or one
@return probability of qubit measureQubit being either zero or one
*/
REAL measureInState(MultiQubit multiQubit, const int measureQubit, int outcome);

/** Updates the state according to this scenario: we ask "are these 3 qubits in 111" and the answer is "no".
The function returns the probability of this outcome (if zero, it will exit with error) 
@param[in,out] multiQubit object representing the set of qubits
@param[in] idQubit1, idQubit2, idQubit3 specified qubits                 
@return Total probability that the 3 qubits are not all in the 1 state. 
*/
REAL filterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3);

/** Evaluates the state according to this scenario: we ask "are these 3 qubits in 111" and the answer is "no".
The function returns the probability of this outcome (if zero, it will exit with error) 
@param[in,out] multiQubit object representing the set of qubits
@param[in] idQubit1, idQubit2, idQubit3 specified qubits                 
@return Total probability that the 3 qubits are not all in the 1 state. 
*/
REAL probOfFilterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3);

/**
Rotate a single qubit by {{1,0},{0,p}} where p is a phase term determined by the type argument
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] type the type of phase gate to apply -- one of {SIGMA_Z, S_GATE, T_GATE}
*/
void phaseGate(MultiQubit multiQubit, const int rotQubit, enum phaseGateType type);





# endif
