/** @file qubits.c
 * The core of the QuEST Library.
 */

# include <math.h>  //SCB new line
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include "precision.h"
# include "qubits.h"

# include <omp.h>

# define DEBUG 0

static int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber);

// Maihi: Where I have made changes I have marked SCB so please note those points - Simon

/** Create a MultiQubit object representing a set of qubits.
 * Allocate space for state vector of probability amplitudes, including space for temporary values to be copied from
 * one other chunk if running the distributed version. Define properties related to the size of the set of qubits.
 * @param[in,out] multiQubit object representing the set of qubits
 * @param[in] numQubits number of qubits in the system
 * @param[in] env object representing the execution environment (local, multinode etc)
 */
void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env)
{
	long long int numAmps = 1L << numQubits;
	long long int numAmpsPerRank = numAmps/env.numRanks;

        multiQubit->stateVec.real = malloc(numAmpsPerRank * sizeof(*(multiQubit->stateVec.real)));
        multiQubit->stateVec.imag = malloc(numAmpsPerRank * sizeof(*(multiQubit->stateVec.imag)));
	if (env.numRanks>1){
		multiQubit->pairStateVec.real = malloc(numAmpsPerRank * sizeof(*(multiQubit->pairStateVec.real)));
		multiQubit->pairStateVec.imag = malloc(numAmpsPerRank * sizeof(*(multiQubit->pairStateVec.imag)));
	}

        if ( (!(multiQubit->stateVec.real) || !(multiQubit->stateVec.imag))
		 && numAmpsPerRank ) {
                printf("Could not allocate memory!");
                exit (EXIT_FAILURE);
        }

	if ( env.numRanks>1 && (!(multiQubit->pairStateVec.real) || !(multiQubit->pairStateVec.imag))
		 && numAmpsPerRank ) {
                printf("Could not allocate memory!");
                exit (EXIT_FAILURE);
	}

	multiQubit->numQubits = numQubits;
	multiQubit->numAmps = numAmpsPerRank;
	multiQubit->chunkId = env.rank;
	multiQubit->numChunks = env.numRanks;

}
/** Deallocate a MultiQubit object representing a set of qubits
 * Free memory allocated to state vector of probability amplitudes, including temporary vector for
 * values copied from another chunk if running the distributed version.
 * @param[in,out] multiQubit object to be deallocated
 * @param[in] env object representing the execution environment (local, multinode etc)
 */
void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env){
	free(multiQubit.stateVec.real);
	free(multiQubit.stateVec.imag);
	if (env.numRanks>1){
		free(multiQubit.pairStateVec.real);
		free(multiQubit.pairStateVec.imag);
	}
}

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
 * @param[in,out] multiQubit object representing the set of qubits
 */
void reportState(MultiQubit multiQubit){
	FILE *state;
	char filename[100];
	long long int index;
	sprintf(filename, "state_rank_%d.csv", multiQubit.chunkId);
	state = fopen(filename, "w");
	if (multiQubit.chunkId==0) fprintf(state, "real, imag\n");

	for(index=0; index<multiQubit.numAmps; index++){
		fprintf(state, REAL_STRING_FORMAT "," REAL_STRING_FORMAT "\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
	}
	fclose(state);
}

/** Print the current state vector of probability amplitudes for a set of qubits to standard out. 
For debugging purposes. Each rank should print output serially. Only print output for systems <= 5 qubits
*/
void reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank){
	long long int index;
	int rank;
	if (multiQubit.numQubits<=5){
		for (rank=0; rank<multiQubit.numChunks; rank++){
			if (multiQubit.chunkId==rank){
				if (reportRank) {
					printf("Reporting state from rank %d [\n", multiQubit.chunkId);
					//printf("\trank, index, real, imag\n");
					printf("real, imag\n");
				} else if (rank==0) {
					printf("Reporting state [\n");
					printf("real, imag\n");
				}

				for(index=0; index<multiQubit.numAmps; index++){
					printf(REAL_STRING_FORMAT ", " REAL_STRING_FORMAT "\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
				}
				if (reportRank || rank==multiQubit.numChunks-1) printf("]\n");
			}
			syncQuESTEnv(env);
		}
	}
}

/** Report metainformation about a set of qubits: number of qubits, number of probability amplitudes.
 * @param[in,out] multiQubit object representing the set of qubits
 * @param[in] env object representing the execution environment (local, multinode etc)
 */
void reportMultiQubitParams(MultiQubit multiQubit){
	long long int numAmps = 1L << multiQubit.numQubits;
	long long int numAmpsPerRank = numAmps/multiQubit.numChunks;
	if (multiQubit.chunkId==0){
                printf("QUBITS:\n");
                printf("Number of qubits is %d.\n", multiQubit.numQubits);
                printf("Number of amps is %lld.\n", numAmps);
		printf("Number of amps per rank is %lld.\n", numAmpsPerRank);
        }
}

void getEnvironmentString(QuESTEnv env, MultiQubit multiQubit, char str[200]){
        int numThreads=1;
# ifdef _OPENMP
        numThreads=omp_get_max_threads(); 
# endif
        sprintf(str, "%dqubits_CPU_%dranksx%dthreads", multiQubit.numQubits, env.numRanks, numThreads);
}

/**
 * Initialise the state vector of probability amplitudes for a set of qubits to the zero state: |000...00>
 * @param[in,out] multiQubit object representing the set of qubits to be initialised
 */
void initStateZero (MultiQubit *multiQubit)
{
	long long int stateVecSize;
	long long int index;

	// dimension of the state vector
	stateVecSize = multiQubit->numAmps;

	// Can't use multiQubit->stateVec as a private OMP var
	REAL *stateVecReal = multiQubit->stateVec.real;
	REAL *stateVecImag = multiQubit->stateVec.imag;

	// initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecSize, stateVecReal, stateVecImag) \
	private  (index) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (index=0; index<stateVecSize; index++) {
			stateVecReal[index] = 0.0;
			stateVecImag[index] = 0.0;
		}
	}

	if (multiQubit->chunkId==0){
		// zero state |0000..0000> has probability 1
		stateVecReal[0] = 1.0;
		stateVecImag[0] = 0.0;
	}

	if (DEBUG) printf("COMPLETED INIT\n");
}

/**
 * Initialise the state vector of probability amplitudes for a set of qubits to an equal real superposition of all amplitudes: |+++...++>
 * @param[in,out] multiQubit object representing the set of qubits to be initialised
 */
void initStatePlus (MultiQubit *multiQubit)
{
	long long int chunkSize, stateVecSize;
	long long int index;

	// dimension of the state vector
	chunkSize = multiQubit->numAmps;
	stateVecSize = chunkSize*multiQubit->numChunks;
	REAL normFactor = 1.0/sqrt(stateVecSize);

	// Can't use multiQubit->stateVec as a private OMP var
	REAL *stateVecReal = multiQubit->stateVec.real;
	REAL *stateVecImag = multiQubit->stateVec.imag;

	// initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (chunkSize, stateVecReal, stateVecImag, normFactor) \
	private  (index) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (index=0; index<chunkSize; index++) {
			stateVecReal[index] = normFactor;
			stateVecImag[index] = 0.0;
		}
	}
	if (DEBUG) printf("COMPLETED INIT\n");
}

/**
 * Initialise the state vector of probability amplitudes such that one qubit is set to 'outcome' and all other qubits are in an equal superposition of zero and one.
 * @param[in,out] multiQubit object representing the set of qubits to be initialised
 * @param[in] qubitId id of qubit to set to state 'outcome'
 * @param[in] value of qubit 'qubitId'
 */
void initStateOfSingleQubit(MultiQubit *multiQubit, int qubitId, int outcome)
{
	long long int chunkSize, stateVecSize;
	long long int index;
	int bit;
	const long long int chunkId=multiQubit->chunkId;

	// dimension of the state vector
	chunkSize = multiQubit->numAmps;
	stateVecSize = chunkSize*multiQubit->numChunks;
	REAL normFactor = 1.0/sqrt(stateVecSize/2);

	// Can't use multiQubit->stateVec as a private OMP var
	REAL *stateVecReal = multiQubit->stateVec.real;
	REAL *stateVecImag = multiQubit->stateVec.imag;

	// initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (chunkSize, stateVecReal, stateVecImag, normFactor, qubitId, outcome) \
	private  (index, bit) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (index=0; index<chunkSize; index++) {
			bit = extractBit(qubitId, index+chunkId*chunkSize);
			if (bit==outcome) {
				stateVecReal[index] = normFactor;
				stateVecImag[index] = 0.0;
			} else {
				stateVecReal[index] = 0.0;
				stateVecImag[index] = 0.0;
			}
		}
	}
	if (DEBUG) printf("COMPLETED INIT\n");
}


/**
 * Initialise the state vector of probability amplitudes to an (unphysical) state with
 * each component of each probability amplitude a unique floating point value. For debugging processes
 * @param[in,out] multiQubit object representing the set of qubits to be initialised
 */
void initStateDebug (MultiQubit *multiQubit)
{
	long long int chunkSize;
	long long int index;

	// dimension of the state vector
	chunkSize = multiQubit->numAmps;

	// Can't use multiQubit->stateVec as a private OMP var
	REAL *stateVecReal = multiQubit->stateVec.real;
	REAL *stateVecImag = multiQubit->stateVec.imag;

	REAL chunkOffset = (2.0*chunkSize*multiQubit->chunkId)/10.0;

	// initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (chunkSize, stateVecReal, stateVecImag, chunkOffset) \
	private  (index) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (index=0; index<chunkSize; index++) {
			stateVecReal[index] = chunkOffset + (index*2.0)/10.0;
			stateVecImag[index] = chunkOffset + (index*2.0+1.0)/10.0;
		}
	}
}

void initializeStateFromSingleFile(MultiQubit *multiQubit, char filename[200], QuESTEnv env){
	long long int chunkSize, stateVecSize;
	long long int indexInChunk, totalIndex;

	chunkSize = multiQubit->numAmps;
	stateVecSize = chunkSize*multiQubit->numChunks;

	REAL *stateVecReal = multiQubit->stateVec.real;
	REAL *stateVecImag = multiQubit->stateVec.imag;
	
	FILE *fp;
	char line[200];

	for (int rank=0; rank<(multiQubit->numChunks); rank++){
		if (rank==multiQubit->chunkId){
			fp = fopen(filename, "r");
			indexInChunk = 0; totalIndex = 0;
			while (fgets(line, sizeof(char)*200, fp) != NULL && totalIndex<stateVecSize){
				if (line[0]!='#'){
					int chunkId = totalIndex/chunkSize;
					if (chunkId==multiQubit->chunkId){
						//! fix -- format needs to work for single precision values
						sscanf(line, "%lf, %lf", &(stateVecReal[indexInChunk]), 
								&(stateVecImag[indexInChunk]));
						indexInChunk += 1;
					}
					totalIndex += 1;
				}
			}	
			fclose(fp);
		}
		syncQuESTEnv(env);
	}
}

int compareStates(MultiQubit mq1, MultiQubit mq2, REAL precision){
	REAL diff;
	int chunkSize = mq1.numAmps;
	for (int i=0; i<chunkSize; i++){
		diff = mq1.stateVec.real[i] - mq2.stateVec.real[i];
		if (diff<0) diff *= -1;
		if (diff>precision) return 0;
		diff = mq1.stateVec.imag[i] - mq2.stateVec.imag[i];
		if (diff<0) diff *= -1;
		if (diff>precision) return 0;
	}
	return 1;
}

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
void rotateQubitLocal (MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta)
{
	long long int sizeBlock, sizeHalfBlock;
	long long int thisBlock, // current block
	     indexUp,indexLo;    // current index and corresponding index in lower half block

	REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps>>1;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// set dimensions
	sizeHalfBlock = 1LL << rotQubit;  
	sizeBlock     = 2LL * sizeHalfBlock; 

	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;
	REAL alphaImag=alpha.imag, alphaReal=alpha.real;
	REAL betaImag=beta.imag, betaReal=beta.real;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag) \
	private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {

			thisBlock   = thisTask / sizeHalfBlock;
			indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
			indexLo     = indexUp + sizeHalfBlock;

			// store current state vector values in temp variables
			stateRealUp = stateVecReal[indexUp];
			stateImagUp = stateVecImag[indexUp];

			stateRealLo = stateVecReal[indexLo];
			stateImagLo = stateVecImag[indexLo];

			// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
			stateVecReal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp 
				- betaReal*stateRealLo - betaImag*stateImagLo;
			stateVecImag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp 
				- betaReal*stateImagLo + betaImag*stateRealLo;

			// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
			stateVecReal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp 
				+ alphaReal*stateRealLo + alphaImag*stateImagLo;
			stateVecImag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp 
				+ alphaReal*stateImagLo - alphaImag*stateRealLo;
		} 
	}

} 

/** Rotate a single qubit in the state
vector of probability amplitudes, given the              
angle rotation arguments, and a subset of the state vector with upper and lower block values 
stored seperately.

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] rot1 rotation angle
@param[in] rot2 rotation angle
@param[in] stateVecUp probability amplitudes in upper half of a block
@param[in] stateVecLo probability amplitudes in lower half of a block
@param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
*/

void rotateQubitDistributed (MultiQubit multiQubit, const int rotQubit,
		Complex rot1, Complex rot2,
		ComplexArray stateVecUp,
		ComplexArray stateVecLo,
		ComplexArray stateVecOut)
{

	REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;  
	const long long int numTasks=multiQubit.numAmps;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	REAL rot1Real=rot1.real, rot1Imag=rot1.imag;
	REAL rot2Real=rot2.real, rot2Imag=rot2.imag;
	REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
	REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
	REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
			rot1Real,rot1Imag, rot2Real,rot2Imag) \
	private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			// store current state vector values in temp variables
			stateRealUp = stateVecRealUp[thisTask];
			stateImagUp = stateVecImagUp[thisTask];

			stateRealLo = stateVecRealLo[thisTask];
			stateImagLo = stateVecImagLo[thisTask];

			// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
			stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp + rot2Real*stateRealLo + rot2Imag*stateImagLo;
			stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp + rot2Real*stateImagLo - rot2Imag*stateRealLo;
		}
	}
}

/** Rotate a single qubit in the state vector of probability amplitudes, given the angle rotation arguments and 
a control qubit. Only perform the rotation for elements where the control qubit is one.
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
void controlRotateQubitLocal (MultiQubit multiQubit, const int rotQubit, const int controlQubit, 
		Complex alpha, Complex beta)
{
	long long int sizeBlock, sizeHalfBlock;
	long long int thisBlock, // current block
	     indexUp,indexLo;    // current index and corresponding index in lower half block

	REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps>>1;
	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;

	int controlBit;

	// As rotations are symmetric, we can apply rotations for all elements where
	// targetQubit==0 and controlQubit==1.  
	// However, this means we will skip the case where targetQubit==controlQubit. 
	// We check for that here. 
	// We could also choose to rotate on targetQubit==1, but are doing it this way 
	// to match the regular rotate implementation. 
	int rotateAll=(rotQubit==controlQubit);

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// set dimensions
	sizeHalfBlock = 1LL << rotQubit;  
	sizeBlock     = 2LL * sizeHalfBlock; 

	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;
	REAL alphaImag=alpha.imag, alphaReal=alpha.real;
	REAL betaImag=beta.imag, betaReal=beta.real;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag,\
			rotateAll) \
	private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {

			thisBlock   = thisTask / sizeHalfBlock;
			indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
			indexLo     = indexUp + sizeHalfBlock;

			controlBit = extractBit (controlQubit, indexUp+chunkId*chunkSize);
			if (rotateAll || controlBit){
				// store current state vector values in temp variables
				stateRealUp = stateVecReal[indexUp];
				stateImagUp = stateVecImag[indexUp];

				stateRealLo = stateVecReal[indexLo];
				stateImagLo = stateVecImag[indexLo];

				// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
				stateVecReal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp 
					- betaReal*stateRealLo - betaImag*stateImagLo;
				stateVecImag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp 
					- betaReal*stateImagLo + betaImag*stateRealLo;

				// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
				stateVecReal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp 
					+ alphaReal*stateRealLo + alphaImag*stateImagLo;
				stateVecImag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp 
					+ alphaReal*stateImagLo - alphaImag*stateRealLo;
			}
		} 
	}

} 

/** Rotate a single qubit in the state
vector of probability amplitudes, given the              
angle rotation arguments, and a subset of the state vector with upper and lower block values 
stored seperately. Only perform the rotation where the control qubit is one.

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] controlQubit qubit to determine whether or not to perform a rotation 
@param[in] rot1 rotation angle
@param[in] rot2 rotation angle
@param[in] stateVecUp probability amplitudes in upper half of a block
@param[in] stateVecLo probability amplitudes in lower half of a block
@param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
*/
void controlRotateQubitDistributed (MultiQubit multiQubit, const int rotQubit, const int controlQubit,
		Complex rot1, Complex rot2,
		ComplexArray stateVecUp,
		ComplexArray stateVecLo,
		ComplexArray stateVecOut)
{

	REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;  
	const long long int numTasks=multiQubit.numAmps;
	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;

	int controlBit;

	// As rotations are symmetric, we can apply rotations for all elements where
	// targetQubit==0 and controlQubit==1.  
	// However, this means we will skip the case where targetQubit==controlQubit. 
	// We check for that here. 
	// We could also choose to rotate on targetQubit==1, but are doing it this way 
	// to match the regular rotate implementation. 
	int rotateAll=(rotQubit==controlQubit);

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	REAL rot1Real=rot1.real, rot1Imag=rot1.imag;
	REAL rot2Real=rot2.real, rot2Imag=rot2.imag;
	REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
	REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
	REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
			rot1Real,rot1Imag, rot2Real,rot2Imag,rotateAll) \
	private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
			if (rotateAll || controlBit){
				// store current state vector values in temp variables
				stateRealUp = stateVecRealUp[thisTask];
				stateImagUp = stateVecImagUp[thisTask];

				stateRealLo = stateVecRealLo[thisTask];
				stateImagLo = stateVecImagLo[thisTask];

				// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
				stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp + rot2Real*stateRealLo + rot2Imag*stateImagLo;
				stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp + rot2Real*stateImagLo - rot2Imag*stateRealLo;
			}
		}
	}
}

/** Rotate a single qubit by {{0,1},{1,0}.
@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
 */
void sigmaXLocal(MultiQubit multiQubit, const int rotQubit)
{
	long long int sizeBlock, sizeHalfBlock;
	long long int thisBlock, // current block
	     indexUp,indexLo;    // current index and corresponding index in lower half block

	REAL stateRealUp,stateImagUp;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps>>1;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// set dimensions
	sizeHalfBlock = 1LL << rotQubit;  
	sizeBlock     = 2LL * sizeHalfBlock; 

	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
	private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			thisBlock   = thisTask / sizeHalfBlock;
			indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
			indexLo     = indexUp + sizeHalfBlock;

			stateRealUp = stateVecReal[indexUp];
			stateImagUp = stateVecImag[indexUp];

			stateVecReal[indexUp] = stateVecReal[indexLo];
			stateVecImag[indexUp] = stateVecImag[indexLo];

			stateVecReal[indexLo] = stateRealUp;
			stateVecImag[indexLo] = stateImagUp;
		} 
	}

}

/** Rotate a single qubit by {{0,1},{1,0}.
Operate on a subset of the state vector with upper and lower block values 
stored seperately. This rotation is just swapping upper and lower values, and
stateVecIn must already be the correct section for this chunk

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
@param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
*/

void sigmaXDistributed (MultiQubit multiQubit, const int rotQubit,
		ComplexArray stateVecIn,
		ComplexArray stateVecOut)
{

	long long int thisTask;  
	const long long int numTasks=multiQubit.numAmps;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	REAL *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
	REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut) \
	private  (thisTask)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			stateVecRealOut[thisTask] = stateVecRealIn[thisTask];
			stateVecImagOut[thisTask] = stateVecImagIn[thisTask];
		}
	}
} 


/** Rotate a single qubit by {{0,1},{1,0} for elements where controlQubit is one.
@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] targetQubit qubit to rotate
@param[in] controlQubit qubit to determine whether or not to perform a rotation 
 */
void controlNotLocal(MultiQubit multiQubit, const int targetQubit, const int controlQubit)
{
	long long int sizeBlock, sizeHalfBlock;
	long long int thisBlock, // current block
	     indexUp,indexLo;    // current index and corresponding index in lower half block

	REAL stateRealUp,stateImagUp;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps>>1;
	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;

	int controlBit;

	// if targetQubit==controlQubit, it is guaranteed that controlQubit==1 when
	// targetQubit==1. As rotations are symmetric, we can instead apply the rotation
	// on all amplitudes where targetQubit==0 as we do here.
	int rotateAll=(targetQubit==controlQubit);

	// test qubit valid
	assert (targetQubit >= 0 && targetQubit < multiQubit.numQubits);

	// set dimensions
	sizeHalfBlock = 1LL << targetQubit;  
	sizeBlock     = 2LL * sizeHalfBlock; 


	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag,rotateAll) \
	private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,controlBit) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			thisBlock   = thisTask / sizeHalfBlock;
			indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
			indexLo     = indexUp + sizeHalfBlock;

			controlBit = extractBit(controlQubit, indexUp+chunkId*chunkSize);
			if (rotateAll || controlBit){
				stateRealUp = stateVecReal[indexUp];
				stateImagUp = stateVecImag[indexUp];

				stateVecReal[indexUp] = stateVecReal[indexLo];
				stateVecImag[indexUp] = stateVecImag[indexLo];

				stateVecReal[indexLo] = stateRealUp;
				stateVecImag[indexLo] = stateImagUp;
			}
		} 
	}

}

/** Rotate a single qubit by {{0,1},{1,0}.
Operate on a subset of the state vector with upper and lower block values 
stored seperately. This rotation is just swapping upper and lower values, and
stateVecIn must already be the correct section for this chunk. Only perform the rotation
for elements where controlQubit is one.

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
@param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
*/

void controlNotDistributed (MultiQubit multiQubit, const int targetQubit, const int controlQubit,
		ComplexArray stateVecIn,
		ComplexArray stateVecOut)
{

	long long int thisTask;  
	const long long int numTasks=multiQubit.numAmps;
	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;
	
	// if targetQubit==controlQubit, it is guaranteed that controlQubit==1 when
	// targetQubit==1. As rotations are symmetric, we can instead apply the rotation
	// on all amplitudes where targetQubit==0 as we do here.
	int rotateAll=(targetQubit==controlQubit);

	int controlBit;

	// test qubit valid
	assert (targetQubit >= 0 && targetQubit < multiQubit.numQubits);

	REAL *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
	REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut,rotateAll) \
	private  (thisTask,controlBit)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
			if (rotateAll || controlBit){
				stateVecRealOut[thisTask] = stateVecRealIn[thisTask];
				stateVecImagOut[thisTask] = stateVecImagIn[thisTask];
			}
		}
	}
} 



/** Rotate a single qubit by {{0,-i},{i,0}.
@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
 */
void sigmaYLocal(MultiQubit multiQubit, const int rotQubit)
{
	long long int sizeBlock, sizeHalfBlock;
	long long int thisBlock, // current block
	     indexUp,indexLo;    // current index and corresponding index in lower half block

	REAL stateRealUp,stateImagUp;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps>>1;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// set dimensions
	sizeHalfBlock = 1LL << rotQubit;  
	sizeBlock     = 2LL * sizeHalfBlock; 

	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
	private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			thisBlock   = thisTask / sizeHalfBlock;
			indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
			indexLo     = indexUp + sizeHalfBlock;

			stateRealUp = stateVecReal[indexUp];
			stateImagUp = stateVecImag[indexUp];

			stateVecReal[indexUp] = stateVecImag[indexLo];
			stateVecImag[indexUp] = -stateVecReal[indexLo];

			stateVecReal[indexLo] = -stateImagUp;
			stateVecImag[indexLo] = stateRealUp;
		} 
	}
}

/** Rotate a single qubit by {{0,-i},{i,0}.
Operate on a subset of the state vector with upper and lower block values 
stored seperately. This rotation is just swapping upper and lower values, and
stateVecIn must already be the correct section for this chunk

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
@param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
@param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
*/

void sigmaYDistributed(MultiQubit multiQubit, const int rotQubit,
		ComplexArray stateVecIn,
		ComplexArray stateVecOut, 
		int updateUpper)
{

	long long int thisTask;  
	const long long int numTasks=multiQubit.numAmps;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	REAL *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
	REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

	int realSign=1, imagSign=1;
	if (updateUpper) imagSign=-1;
	else realSign = -1;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut,realSign,imagSign) \
	private  (thisTask)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			stateVecRealOut[thisTask] = realSign*stateVecImagIn[thisTask];
			stateVecImagOut[thisTask] = imagSign*stateVecRealIn[thisTask];
		}
	}
} 

/** Rotate a single qubit by {{1,1},{1,-1}}/sqrt2.
@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
 */
void hadamardLocal(MultiQubit multiQubit, const int rotQubit)
{
	long long int sizeBlock, sizeHalfBlock;
	long long int thisBlock, // current block
	     indexUp,indexLo;    // current index and corresponding index in lower half block

	REAL stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps>>1;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// set dimensions
	sizeHalfBlock = 1LL << rotQubit;  
	sizeBlock     = 2LL * sizeHalfBlock; 

	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

	REAL recRoot2 = 1.0/sqrt(2);

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, recRoot2) \
	private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) 
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			thisBlock   = thisTask / sizeHalfBlock;
			indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
			indexLo     = indexUp + sizeHalfBlock;

			stateRealUp = stateVecReal[indexUp];
			stateImagUp = stateVecImag[indexUp];

			stateRealLo = stateVecReal[indexLo];
			stateImagLo = stateVecImag[indexLo];

			stateVecReal[indexUp] = recRoot2*(stateRealUp + stateRealLo);
			stateVecImag[indexUp] = recRoot2*(stateImagUp + stateImagLo);

			stateVecReal[indexLo] = recRoot2*(stateRealUp - stateRealLo);
			stateVecImag[indexLo] = recRoot2*(stateImagUp - stateImagLo);
		} 
	}
}

/** Rotate a single qubit by {{1,1},{1,-1}}/sqrt2.
Operate on a subset of the state vector with upper and lower block values 
stored seperately. This rotation is just swapping upper and lower values, and
stateVecIn must already be the correct section for this chunk

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
@param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
@param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
*/

void hadamardDistributed(MultiQubit multiQubit, const int rotQubit,
		ComplexArray stateVecUp,
		ComplexArray stateVecLo,
		ComplexArray stateVecOut,
		int updateUpper)
{

	REAL   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
	long long int thisTask;  
	const long long int numTasks=multiQubit.numAmps;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	int sign;
	if (updateUpper) sign=1;
	else sign=-1;

	REAL recRoot2 = 1.0/sqrt(2);

	REAL *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
	REAL *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
	REAL *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
			recRoot2, sign) \
	private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			// store current state vector values in temp variables
			stateRealUp = stateVecRealUp[thisTask];
			stateImagUp = stateVecImagUp[thisTask];

			stateRealLo = stateVecRealLo[thisTask];
			stateImagLo = stateVecImagLo[thisTask];

			stateVecRealOut[thisTask] = recRoot2*(stateRealUp + sign*stateRealLo);
			stateVecImagOut[thisTask] = recRoot2*(stateImagUp + sign*stateImagLo);
		}
	}
}

/**
Rotate a single qubit by {{1,0},{0,p}} where p is a phase term determined by the type argument

@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] type the type of phase gate to apply -- one of {SIGMA_Z, S_GATE, T_GATE}
 */
void phaseGateLocal(MultiQubit multiQubit, const int rotQubit, enum phaseGateType type)
{
	long long int sizeBlock, sizeHalfBlock;
	long long int thisBlock, // current block
	     indexUp,indexLo;    // current index and corresponding index in lower half block

	REAL stateRealLo,stateImagLo;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps>>1;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// set dimensions
	sizeHalfBlock = 1LL << rotQubit;  
	sizeBlock     = 2LL * sizeHalfBlock; 

	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

	REAL recRoot2 = 1.0/sqrt(2);

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock,stateVecReal,stateVecImag,recRoot2,type) \
	private  (thisTask,thisBlock,indexUp,indexLo,stateRealLo,stateImagLo) 
# endif
	{
		if (type==SIGMA_Z){
# ifdef _OPENMP
			# pragma omp for schedule (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				//! fix -- can i rewrite this to not use mod?
				thisBlock   = thisTask / sizeHalfBlock;
				indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
				indexLo     = indexUp + sizeHalfBlock;

				stateVecReal[indexLo] = -stateVecReal[indexLo];
				stateVecImag[indexLo] = -stateVecImag[indexLo];
			} 
		} 
		
		else if (type==S_GATE){
# ifdef _OPENMP
			# pragma omp for schedule (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				//! fix -- can i rewrite this to not use mod?
				thisBlock   = thisTask / sizeHalfBlock;
				indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
				indexLo     = indexUp + sizeHalfBlock;
				stateRealLo = stateVecReal[indexLo];
				stateImagLo = stateVecImag[indexLo];

				stateVecReal[indexLo] = -stateImagLo;
				stateVecImag[indexLo] = stateRealLo;
			} 
		} else if (type==T_GATE){
# ifdef _OPENMP
			# pragma omp for schedule (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				//! fix -- can i rewrite this to not use mod?
				thisBlock   = thisTask / sizeHalfBlock;
				indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
				indexLo     = indexUp + sizeHalfBlock;
				stateRealLo = stateVecReal[indexLo];
				stateImagLo = stateVecImag[indexLo];

				stateVecReal[indexLo] = recRoot2 * (stateRealLo - stateImagLo);
				stateVecImag[indexLo] = recRoot2 * (stateRealLo + stateImagLo);
			} 
		} else printf("Type %d is an invalid phase gate\n", type);
	}
}

/**
Rotate a single qubit by {{1,0},{0,p}} where p is a phase term determined by the type argument

@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] type the type of phase gate to apply -- one of {SIGMA_Z, S_GATE, T_GATE}
 */
void phaseGateDistributed(MultiQubit multiQubit, const int rotQubit, enum phaseGateType type)
{
	REAL stateRealLo,stateImagLo;
	long long int thisTask;         
	const long long int numTasks=multiQubit.numAmps;

	// test qubit valid
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// Can't use multiQubit.stateVec as a private OMP var
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

	REAL recRoot2 = 1.0/sqrt(2);

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (stateVecReal,stateVecImag, recRoot2, type) \
	private  (thisTask,stateRealLo,stateImagLo) 
# endif
	{
		if (type==SIGMA_Z){
# ifdef _OPENMP
			# pragma omp for schedule (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				stateVecReal[thisTask] = -stateVecReal[thisTask];
				stateVecImag[thisTask] = -stateVecImag[thisTask];
			} 
		} else if (type==S_GATE){
# ifdef _OPENMP
			# pragma omp for schedule (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				stateRealLo = stateVecReal[thisTask];
				stateImagLo = stateVecImag[thisTask];

				stateVecReal[thisTask] = -stateImagLo;
				stateVecImag[thisTask] = stateRealLo;
			} 
		} else if (type==T_GATE){
# ifdef _OPENMP
			# pragma omp for schedule (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				stateRealLo = stateVecReal[thisTask];
				stateImagLo = stateVecImag[thisTask];

				stateVecReal[thisTask] = recRoot2 * (stateRealLo - stateImagLo);
				stateVecImag[thisTask] = recRoot2 * (stateRealLo + stateImagLo);
			} 
		} else printf("Type %d is an invalid phase gate\n", type);
	}
}
 
void sigmaZ(MultiQubit multiQubit, const int rotQubit)
{
	        phaseGate(multiQubit, rotQubit, SIGMA_Z);
}

void sGate(MultiQubit multiQubit, const int rotQubit)
{
	        phaseGate(multiQubit, rotQubit, S_GATE);
} 

void tGate(MultiQubit multiQubit, const int rotQubit)
{
	        phaseGate(multiQubit, rotQubit, T_GATE);
}

/** Measure the total probability of a specified qubit being in the zero state across all amplitudes in this chunk.
Size of regions to skip is less than    
the size of one chunk.                   

@param[in] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/

REAL findProbabilityOfZeroLocal (MultiQubit multiQubit,
		const int measureQubit)
{
	// ----- sizes
	long long int sizeBlock,                                           // size of blocks
	sizeHalfBlock;                                       // size of blocks halved
	// ----- indices
	long long int thisBlock,                                           // current block
	     index;                                               // current index for first half block
	// ----- measured probability
	REAL   totalProbability;                                    // probability (returned) value
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	long long int numTasks=multiQubit.numAmps>>1;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (measureQubit >= 0 && measureQubit < multiQubit.numQubits);


	// ---------------------------------------------------------------- //
	//            dimensions                                            //
	// ---------------------------------------------------------------- //
	sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
	// and then the number to skip
	sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

	// initialise returned value
	totalProbability = 0.0;

	// initialise correction for kahan summation
	if (DEBUG) printf("sizeHalfBlock=%Ld sizeBlock=%Ld numTasks=%Ld\n",sizeHalfBlock,sizeBlock,numTasks);

	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
	private   (thisTask,thisBlock,index) \
	reduction ( +:totalProbability )
# endif	
	{
# ifdef _OPENMP
		# pragma omp for schedule  (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			thisBlock = thisTask / sizeHalfBlock;
			index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;

			if (index<0){ printf("ABORTING as index=%Ld with thisBlock = %Ld  thisTask=%Ld \n", index,thisBlock,thisTask); exit(1);}

			// summation -- simple implementation
			totalProbability += stateVecReal[index]*stateVecReal[index]
				+ stateVecImag[index]*stateVecImag[index];

			/*
			// summation -- kahan correction
			y = stateVecReal[index]*stateVecReal[index]
			+ stateVecImag[index]*stateVecImag[index] - c;
			t = totalProbability + y;
			c = (t - totalProbability) - y;
			totalProbability = t;
			*/

		}
	}
	return totalProbability;
}

/** Measure the probability of a specified qubit being in the zero state across all amplitudes held in this chunk.
Size of regions to skip is a multiple of chunkSize.

@param[in] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/

REAL findProbabilityOfZeroDistributed (MultiQubit multiQubit,
		const int measureQubit)
{
	// ----- measured probability
	REAL   totalProbability;                                    // probability (returned) value
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	long long int numTasks=multiQubit.numAmps;
	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (measureQubit >= 0 && measureQubit < multiQubit.numQubits);

	// ---------------------------------------------------------------- //
	//            find probability                                      //
	// ---------------------------------------------------------------- //

	// initialise returned value
	totalProbability = 0.0;

	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	shared    (numTasks,stateVecReal,stateVecImag) \
	private   (thisTask) \
	reduction ( +:totalProbability )
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule  (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			// summation -- simple implementation
			totalProbability += stateVecReal[thisTask]*stateVecReal[thisTask]
				+ stateVecImag[thisTask]*stateVecImag[thisTask];

			/*
			// summation -- kahan correction
			y = stateVecReal[thisTask]*stateVecReal[thisTask]
			+ stateVecImag[thisTask]*stateVecImag[thisTask] - c;
			t = totalProbability + y;
			c = (t - totalProbability) - y;
			totalProbability = t;
			*/

		}
	}

	return totalProbability;
}

/** Get the value of the bit at a particular index in a number.
SCB edit: new definition of extractBit is much faster ***
 * @param[in] locationOfBitFromRight location of bit in theEncodedNumber
 * @param[in] theEncodedNumber number to search
 * @return the value of the bit in theEncodedNumber
 */
static int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber)
{
	return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

/** The control phase (the two qubit phase gate).
For each state, if both input qubits are equal to one, multiply the amplitude of that state by -1.
@param[in,out] multiQubit object representing the set of qubits
@param[in] idQubit1, idQubit2 specified qubits                 
*/

void controlPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //

	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < multiQubit.numQubits && idQubit2 < multiQubit.numQubits);

	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;

	// ---------------------------------------------------------------- //
	//            initialise the state to |0000..0>                     //
	// ---------------------------------------------------------------- //

	// dimension of the state vector
	stateVecSize = multiQubit.numAmps;
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel for \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag ) \
	private  (index,bit1,bit2)		       \
	schedule (static)
# endif
	for (index=0; index<stateVecSize; index++) {
		bit1 = extractBit (idQubit1, index+chunkId*chunkSize);
		bit2 = extractBit (idQubit2, index+chunkId*chunkSize);
		if (bit1 && bit2) {
			stateVecReal [index] = - stateVecReal [index];
			stateVecImag [index] = - stateVecImag [index];
		}
	}
}

/** The control not gate
*/
/*
void controlNotGate (MultiQubit multiQubit, const int control, const int target)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //

	assert (control >= 0 && target >= 0 && control < multiQubit.numQubits && target < multiQubit.numQubits);


	// ---------------------------------------------------------------- //
	//            initialise the state to |0000..0>                     //
	// ---------------------------------------------------------------- //

	// dimension of the state vector
	stateVecSize = multiQubit.numAmps;
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel for \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag ) \
	private  (index,bit1,bit2)		       \
	schedule (static)
# endif
	for (index=0; index<stateVecSize; index++) {
		bit1 = extractBit (control, index);
		bit2 = extractBit (target, index);
		if (bit1 && bit2) {
			stateVecReal [index] = - stateVecReal [index];
			stateVecImag [index] = - stateVecImag [index];
		}
	}
}
*/

/** The control phase (the four qubit phase gate).
For each state, if all four input qubits are equal to one, multiply the amplitude of that state by -1.
@param[in,out] multiQubit object representing the set of qubits
@param[in] idQubit1, idQubit2, idQubit3, idQubit4 specified qubits                 
*/
void quadCPhaseGate (MultiQubit multiQubit, const int idQubit1, const int idQubit2, 
		const int idQubit3, const int idQubit4)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2, bit3, bit4;
	
	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < multiQubit.numQubits && idQubit2 < multiQubit.numQubits);

	stateVecSize = multiQubit.numAmps;
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag ) \
	private  (index,bit1,bit2,bit3,bit4)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (index=0; index<stateVecSize; index++) {
			bit1 = extractBit (idQubit1, index+chunkId*chunkSize);
			bit2 = extractBit (idQubit2, index+chunkId*chunkSize);
			bit3 = extractBit (idQubit3, index+chunkId*chunkSize);
			bit4 = extractBit (idQubit4, index+chunkId*chunkSize);
			if (bit1 && bit2 && bit3 && bit4) {
				stateVecReal [index] = - stateVecReal [index];
				stateVecImag [index] = - stateVecImag [index];
			}
		}
	}
}

/** Update the state vector to be consistent with measuring measureQubit=0 if outcome=0 and measureQubit=1
if outcome=1.
Performs an irreversible change to the state vector: it updates the vector according
to the event that an outcome have been measured on the qubit indicated by measureQubit (where 
this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
then renormalising based on the total probability of measuring measureQubit=0 or 1 according to the 
value of outcome. 
In the local version, one or more blocks (with measureQubit=0 in the first half of the block and
measureQubit=1 in the second half of the block) fit entirely into one chunk. 

@param[in,out] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@param[in] totalProbability probability of qubit measureQubit being either zero or one
@param[in] outcome to measure the probability of and set the state to -- either zero or one
*/
void measureInStateLocal(MultiQubit multiQubit, int measureQubit, REAL totalProbability, int outcome)
{
	// ----- sizes
	long long int sizeBlock,                                           // size of blocks
	sizeHalfBlock;                                       // size of blocks halved
	// ----- indices
	long long int thisBlock,                                           // current block
	     index;                                               // current index for first half block
	// ----- measured probability
	REAL   renorm;                                    // probability (returned) value
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	// (good for shared memory parallelism)
	long long int numTasks=multiQubit.numAmps>>1;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (measureQubit >= 0 && measureQubit < multiQubit.numQubits);
	assert (totalProbability != 0);

	// ---------------------------------------------------------------- //
	//            dimensions                                            //
	// ---------------------------------------------------------------- //
	sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
	// and then the number to skip
	sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)
	
	renorm=1/sqrt(totalProbability);
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;


# ifdef _OPENMP
# pragma omp parallel \
	default (none) \
	shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag,renorm,outcome) \
	private   (thisTask,thisBlock,index)
# endif
	{
		if (outcome==0){
			// measure qubit is 0
# ifdef _OPENMP
			# pragma omp for schedule  (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				thisBlock = thisTask / sizeHalfBlock;
				index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
				stateVecReal[index]=stateVecReal[index]*renorm;
				stateVecImag[index]=stateVecImag[index]*renorm;

				stateVecReal[index+sizeHalfBlock]=0;
				stateVecImag[index+sizeHalfBlock]=0;
			}
		} else {
			// measure qubit is 1
# ifdef _OPENMP
			# pragma omp for schedule  (static)
# endif
			for (thisTask=0; thisTask<numTasks; thisTask++) {
				thisBlock = thisTask / sizeHalfBlock;
				index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
				stateVecReal[index]=0;
				stateVecImag[index]=0;

				stateVecReal[index+sizeHalfBlock]=stateVecReal[index+sizeHalfBlock]*renorm;
				stateVecImag[index+sizeHalfBlock]=stateVecImag[index+sizeHalfBlock]*renorm;
			}
		}
	}

}

/** Renormalise parts of the state vector where measureQubit=0 or 1, based on the total probability of that qubit being
in state 0 or 1.
Measure in Zero performs an irreversible change to the state vector: it updates the vector according
to the event that the value 'outcome' has been measured on the qubit indicated by measureQubit (where 
this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
then renormalising based on the total probability of measuring measureQubit=0 if outcome=0 and
measureQubit=1 if outcome=1.
In the distributed version, one block (with measureQubit=0 in the first half of the block and
measureQubit=1 in the second half of the block) is spread over multiple chunks, meaning that each chunks performs
only renormalisation or only setting amplitudes to 0. This function handles the renormalisation.

@param[in,out] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@param[in] totalProbability probability of qubit measureQubit being zero
*/

REAL measureInStateDistributedRenorm (MultiQubit multiQubit, const int measureQubit, const REAL totalProbability)
{
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	long long int numTasks=multiQubit.numAmps;
	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (measureQubit >= 0 && measureQubit < multiQubit.numQubits);
	assert (totalProbability != 0);

	REAL renorm=1/sqrt(totalProbability);
	
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	shared    (numTasks,stateVecReal,stateVecImag) \
	private   (thisTask)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule  (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			// summation -- simple implementation
			stateVecReal[thisTask] = stateVecReal[thisTask]*renorm;
			stateVecImag[thisTask] = stateVecImag[thisTask]*renorm;
		}
	}
	return totalProbability;
}

/** Set all amplitudes in one chunk to 0. 
Measure in Zero performs an irreversible change to the state vector: it updates the vector according
to the event that a zero have been measured on the qubit indicated by measureQubit (where 
this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
then renormalising based on the total probability of measuring measureQubit=0 or 1.
In the distributed version, one block (with measureQubit=0 in the first half of the block and
measureQubit=1 in the second half of the block) is spread over multiple chunks, meaning that each chunks performs
only renormalisation or only setting amplitudes to 0. This function handles setting amplitudes to 0.

@param[in,out] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
*/

void measureInStateDistributedSetZero(MultiQubit multiQubit, const int measureQubit)
{
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	long long int numTasks=multiQubit.numAmps;
	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (measureQubit >= 0 && measureQubit < multiQubit.numQubits);

	// ---------------------------------------------------------------- //
	//            find probability                                      //
	// ---------------------------------------------------------------- //

	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	shared    (numTasks,stateVecReal,stateVecImag) \
	private   (thisTask)
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule  (static)
# endif
		for (thisTask=0; thisTask<numTasks; thisTask++) {
			// summation -- simple implementation
			stateVecReal[thisTask] = 0;
			stateVecImag[thisTask] = 0;
		}
	}
}

/** Updates the state according to this scenario: we ask "are these 3 qubits in 111" and the answer is "no".
@param[in,out] multiQubit object representing the set of qubits
@param[in] idQubit1, idQubit2, idQubit3 specified qubits                 
@param[in] probOfFilter Total probability that the 3 qubits are not all in the 1 state. 
*/
void filterOut111Local(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3,
		const REAL probOfFilter)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2, bit3;
	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < multiQubit.numQubits && idQubit2 < multiQubit.numQubits);

	assert (probOfFilter != 0);
	stateVecSize = multiQubit.numAmps;

	if ( probOfFilter<1e-16 ){ printf("Extremely small or negative profOfFilter="REAL_STRING_FORMAT"; aborting! \n",probOfFilter); exit(1);}
	REAL myNorm=1/sqrt(probOfFilter);
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag, myNorm) \
	private  (index,bit1,bit2,bit3)		       
# endif 
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (index=0; index<stateVecSize; index++) {
			bit1 = extractBit (idQubit1, index+chunkId*chunkSize);
			bit2 = extractBit (idQubit2, index+chunkId*chunkSize);
			bit3 = extractBit (idQubit3, index+chunkId*chunkSize);
			if ((bit1 && bit2 && bit3)) {
				stateVecReal[index]=0;
				stateVecImag [index]=0;
			}else{
				stateVecReal[index] *= myNorm;
				stateVecImag[index] *= myNorm;
			}
		}
	}
}

/** Evaluates the state according to this scenario: we ask "are these 3 qubits in 111" and the answer is "no".
The function returns the probability of this outcome across all amplitudes in this chunk (if zero, it will exit with error) 
@param[in,out] multiQubit object representing the set of qubits
@param[in] idQubit1, idQubit2, idQubit3 specified qubits                 
@return Total probability that the 3 qubits are not all in the 1 state. 
*/
REAL probOfFilterOut111Local(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2, bit3;
	const long long int chunkSize=multiQubit.numAmps;
	const long long int chunkId=multiQubit.chunkId;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < multiQubit.numQubits && idQubit2 < multiQubit.numQubits);

	stateVecSize = multiQubit.numAmps;
	REAL probOfFilter=0;
	
	REAL *stateVecReal = multiQubit.stateVec.real;
	REAL *stateVecImag = multiQubit.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag) \
	private  (index,bit1,bit2,bit3)		       \
	reduction ( +:probOfFilter )
# endif
	{
# ifdef _OPENMP
		# pragma omp for schedule (static)
# endif
		for (index=0; index<stateVecSize; index++) {
			bit1 = extractBit (idQubit1, index+chunkId*chunkSize);
			bit2 = extractBit (idQubit2, index+chunkId*chunkSize);
			bit3 = extractBit (idQubit3, index+chunkId*chunkSize);
			if (!(bit1 && bit2 && bit3)) {
				probOfFilter+= stateVecReal[index]*stateVecReal[index] + stateVecImag[index]* stateVecImag [index];
			}
		}
	}
	return probOfFilter;
}

/** Get probability of the state at an index in the state vector.
For debugging purposes.
@param[in] multiQubit object representing a set of qubits
@param[in] index index in state vector of probability amplitudes
@return real component * real component + imag component * imag component
*/
REAL getProbEl(MultiQubit multiQubit, long long int index){
        REAL real;
        REAL imag;
        real = getRealAmpEl(multiQubit, index);
        imag = getImagAmpEl(multiQubit, index);
        return real*real + imag*imag;
}
