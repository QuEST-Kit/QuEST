/** @file qubits.c
 * The core of the QUEST Library.
 */

# include <math.h>  //SCB new line
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include "qubits.h"

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
void createMultiQubit(MultiQubit *multiQubit, int numQubits, QUESTEnv env)
{
	long long int numAmps = 1L << numQubits;
	long long int numAmpsPerRank = numAmps/env.numRanks;

        multiQubit->stateVec.real = malloc(numAmpsPerRank * sizeof(multiQubit->stateVec.real));
        multiQubit->stateVec.imag = malloc(numAmpsPerRank * sizeof(multiQubit->stateVec.imag));
	if (env.numRanks>1){
		multiQubit->pairStateVec.real = malloc(numAmpsPerRank * sizeof(multiQubit->pairStateVec.real));
		multiQubit->pairStateVec.imag = malloc(numAmpsPerRank * sizeof(multiQubit->pairStateVec.imag));
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
void destroyMultiQubit(MultiQubit multiQubit, QUESTEnv env){
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
		fprintf(state, "%.12f, %.12f\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
	}
	fclose(state);
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

/**
 * Initialise the state vector of probability amplitudes for a set of qubits to the zero state: |000...00>
 * @param[in,out] multiQubit object representing the set of qubits to be initialised
 */
void initStateVec (MultiQubit *multiQubit)
{
	long long int stateVecSize;
	long long int index;

	// dimension of the state vector
	stateVecSize = multiQubit->numAmps;

	// Can't use multiQubit->stateVec as a private OMP var
	double *stateVecReal = multiQubit->stateVec.real;
	double *stateVecImag = multiQubit->stateVec.imag;

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
void rotateQubitLocal (MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta)
{
	// ----- sizes
	long long int sizeBlock,                                           // size of blocks
	sizeHalfBlock;                                       // size of blocks halved
	// ----- indices
	long long int thisBlock,                                           // current block
	     indexUp,indexLo;                                     // current index and corresponding index in lower half block

	// ----- temp variables
	double   stateRealUp,stateRealLo,                             // storage for previous state values
		 stateImagUp,stateImagLo;                             // (used in updates)
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	const long long int numTasks=multiQubit.numAmps>>1;
	// (good for shared memory parallelism)


	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);


	// ---------------------------------------------------------------- //
	//            dimensions                                            //
	// ---------------------------------------------------------------- //
	sizeHalfBlock = 1LL << rotQubit;                               // size of blocks halved
	sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks


	// ---------------------------------------------------------------- //
	//            rotate                                                //
	// ---------------------------------------------------------------- //

	//
	// --- task-based shared-memory parallel implementation
	//
	
	// Can't use multiQubit.stateVec as a private OMP var
	double *stateVecReal = multiQubit.stateVec.real;
	double *stateVecImag = multiQubit.stateVec.imag;
	double alphaImag=alpha.imag, alphaReal=alpha.real;
	double betaImag=beta.imag, betaReal=beta.real;

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
			stateVecReal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo;
			stateVecImag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo;

			// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
			stateVecReal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo;
			stateVecImag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo;
		} // end for loop
	}

} // end of function definition

/** Rotate a single qubit in the state
vector of probability amplitudes, given the              
angle rotation arguments, and a subset of the state vector with upper and lower block values 
stored seperately.

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits to be initialised
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
	// ----- temp variables
	double   stateRealUp,stateRealLo,                             // storage for previous state values
	stateImagUp,stateImagLo;                             // (used in updates)
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	const long long int numTasks=multiQubit.numAmps;

	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);

	// ---------------------------------------------------------------- //
	//            rotate                                                //
	// ---------------------------------------------------------------- //

	//
	// --- task-based shared-memory parallel implementation
	//
	double rot1Real=rot1.real, rot1Imag=rot1.imag;
	double rot2Real=rot2.real, rot2Imag=rot2.imag;
	double *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
	double *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
	double *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

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
		} // end for loop
	}
} // end of function definition

/** Measure the probability
of a specified qubit being in the zero state.     
Size of regions to skip is less than    
the size of one chunk.                   

@param[in] multiQubit object representing the set of qubits to be initialised
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/

double findProbabilityOfZeroLocal (MultiQubit multiQubit,
		const int measureQubit)
{
	// ----- sizes
	long long int sizeBlock,                                           // size of blocks
	sizeHalfBlock;                                       // size of blocks halved
	// ----- indices
	long long int thisBlock,                                           // current block
	     index;                                               // current index for first half block
	// ----- measured probability
	double   totalProbability;                                    // probability (returned) value
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	long long int numTasks=multiQubit.numAmps>>1;
	// (good for shared memory parallelism)

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

	// ---------------------------------------------------------------- //
	//            find probability                                      //
	// ---------------------------------------------------------------- //

	// initialise returned value
	totalProbability = 0.0;

	// initialise correction for kahan summation
	if (DEBUG) printf("sizeHalfBlock=%Ld sizeBlock=%Ld numTasks=%Ld\n",sizeHalfBlock,sizeBlock,numTasks);

	//
	// --- task-based shared-memory parallel implementation
	//
	
	double *stateVecReal = multiQubit.stateVec.real;
	double *stateVecImag = multiQubit.stateVec.imag;

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

/** Measure the probability
of a specified qubit being in the zero state.     
Size of regions to skip is a multiple of chunkSize.

@param[in] multiQubit object representing the set of qubits to be initialised
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/

double findProbabilityOfZeroDistributed (MultiQubit multiQubit,
		const int measureQubit)
{
	// ----- measured probability
	double   totalProbability;                                    // probability (returned) value
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

	// initialise correction for kahan summation

	//
	// --- task-based shared-memory parallel implementation
	//
	
	double *stateVecReal = multiQubit.stateVec.real;
	double *stateVecImag = multiQubit.stateVec.imag;

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

// *** SCB edit: new definition of extractBit is much faster ***
/** Get the value of the bit at a particular index in a number
 * @param[in] locationOfBitFromRight location of bit in theEncodedNumber
 * @param[in] theEncodedNumber number to search
 * @return the value of the bit in theEncodedNumber
 */
static int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber)
{
	return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

/** Implement the control phase       
(the two qubit phase gate).
**** REWRITE TO USE MULTIQUBIT
     input:                                                           //
                    numQubits     -- number of qubits                 //
                    idQubit1,     -- specified qubits                 //
                    idQubit2                                          //
                    stateVecReal, -- real/imag parts of               //
                    stateVecImag     the state vector                 //
                                                                      //
     output:                                                          //
                    stateVecReal, -- real/imag parts of               //
                    stateVecImag     the state vector (overwritten)   //
                                                                      //
*/

void controlPhaseGate (const int numQubits, const int idQubit1, const int idQubit2,
		double *restrict stateVecReal, double *restrict stateVecImag)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //

	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < numQubits && idQubit2 < numQubits);


	// ---------------------------------------------------------------- //
	//            initialise the state to |0000..0>                     //
	// ---------------------------------------------------------------- //

	// dimension of the state vector
	stateVecSize = 1LL << numQubits;

# ifdef _OPENMP
# pragma omp parallel for \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag ) \
	private  (index,bit1,bit2)		       \
	schedule (static)
# endif
	for (index=0; index<stateVecSize; index++) {
		bit1 = extractBit (idQubit1, index);
		bit2 = extractBit (idQubit2, index);
		if (bit1 && bit2) {
			stateVecReal [index] = - stateVecReal [index];
			stateVecImag [index] = - stateVecImag [index];
		}
	}
}

// SCB the functions below, quadCPhaseGate and (more importantly) measureInZero are new

// SCB tripleCPhaseGate is just like controlPhaseGate except it applies the conditional phase depending on 4 qubits, not two

void quadCPhaseGate (const int numQubits, const int idQubit1, const int idQubit2, const int idQubit3, const int idQubit4, double *restrict stateVecReal, double *restrict stateVecImag)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2, bit3, bit4;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < numQubits && idQubit2 < numQubits);

	stateVecSize = 1LL << numQubits;

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
			bit1 = extractBit (idQubit1, index);
			bit2 = extractBit (idQubit2, index);
			bit3 = extractBit (idQubit3, index);
			bit4 = extractBit (idQubit4, index);
			if (bit1 && bit2 && bit3 && bit4) {
				stateVecReal [index] = - stateVecReal [index];
				stateVecImag [index] = - stateVecImag [index];
			}
		}
	}
}


// measure in Zero performs an irreversible change to the state vector: it updates the vector according
// to the event that a zero have been measured on the qubit indicated by measureQubit (where 
// this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
// then renormalising. It also returns the probability that this event would happen.

void measureInZeroLocal(MultiQubit multiQubit, int measureQubit, double totalProbability)
{
	// ----- sizes
	long long int sizeBlock,                                           // size of blocks
	sizeHalfBlock;                                       // size of blocks halved
	// ----- indices
	long long int thisBlock,                                           // current block
	     index;                                               // current index for first half block
	// ----- measured probability
	double   renorm;                                    // probability (returned) value
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	// (good for shared memory parallelism)
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

	// ---------------------------------------------------------------- //
	//            find probability                                      //
	// ---------------------------------------------------------------- //

	//
	// --- task-based shared-memory parallel implementation
	//
	renorm=1/sqrt(totalProbability);
	double *stateVecReal = multiQubit.stateVec.real;
	double *stateVecImag = multiQubit.stateVec.imag;


# ifdef _OPENMP
# pragma omp parallel \
	default (none) \
	shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag,renorm) \
	private   (thisTask,thisBlock,index)
# endif
	{
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
	}

	//SCB this is a debugging style check. It is probably useful to leave in, but it could be parallelised I guess
	//  double checkTotal=1.;
	//  for (index=0; index<2*numTasks; index++) {
	//  	checkTotal=checkTotal-(stateVecReal[index]*stateVecReal[index] + stateVecImag[index]*stateVecImag[index]);
	//  }
	//  if (checkTotal>0.00001){printf("Deviation of sum squared amps from unity is %.16f\n",checkTotal); exit(1);}
}

double measureInZeroDistributedRenorm (MultiQubit multiQubit, const int measureQubit, const double totalProbability)
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

	// initialise correction for kahan summation

	//
	// --- task-based shared-memory parallel implementation
	//
	
	double renorm=1/sqrt(totalProbability);
	
	double *stateVecReal = multiQubit.stateVec.real;
	double *stateVecImag = multiQubit.stateVec.imag;

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

void measureInZeroDistributedSetZero(MultiQubit multiQubit, const int measureQubit)
{
	// ----- measured probability
	double   totalProbability;                                    // probability (returned) value
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

	// initialise correction for kahan summation

	//
	// --- task-based shared-memory parallel implementation
	//
	
	double *stateVecReal = multiQubit.stateVec.real;
	double *stateVecImag = multiQubit.stateVec.imag;

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
			stateVecReal[thisTask] = 0;
			stateVecImag[thisTask] = 0;
		}
	}
}

// filterOut111 updates the state according to this scenario: we ask "are these 3 qubits in 111" and the answer is "no"
//              the function returns the probability of this outcome (if zero, it will exit with error) 

double filterOut111 (const int numQubits, const int idQubit1, const int idQubit2, const int idQubit3,
		double *restrict stateVecReal,
		double *restrict stateVecImag)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2, bit3;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < numQubits && idQubit2 < numQubits);

	stateVecSize = 1LL << numQubits;
	double probOfFilter=0;

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
			bit1 = extractBit (idQubit1, index);
			bit2 = extractBit (idQubit2, index);
			bit3 = extractBit (idQubit3, index);
			if (!(bit1 && bit2 && bit3)) {
				probOfFilter+= stateVecReal[index]*stateVecReal[index] + stateVecImag[index]* stateVecImag [index];
			}
		}
	}
	if ( probOfFilter<1e-16 ){ printf("Extremely small or negative profOfFilter=%.8e; aborting! \n",probOfFilter); exit(1);}
	double myNorm=1/sqrt(probOfFilter);

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
			bit1 = extractBit (idQubit1, index);
			bit2 = extractBit (idQubit2, index);
			bit3 = extractBit (idQubit3, index);
			if ((bit1 && bit2 && bit3)) {
				stateVecReal[index]=0;
				stateVecImag [index]=0;
			}else{
				stateVecReal[index] *= myNorm;
				stateVecImag[index] *= myNorm;
			}
		}
	}
	return probOfFilter;
}

// probFilterOut111 evaluates the state according to this scenario: we ask "are these 3 qubits in 111" and the answer is "no"
//              the function returns the probability of this outcome (if zero, it will exit with error) 

double probOfFilterOut111 (const int numQubits, const int idQubit1, const int idQubit2, const int idQubit3,
		double *restrict stateVecReal,
		double *restrict stateVecImag)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2, bit3;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (idQubit1 >= 0 && idQubit2 >= 0 && idQubit1 < numQubits && idQubit2 < numQubits);

	stateVecSize = 1LL << numQubits;
	double probOfFilter=0;

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
			bit1 = extractBit (idQubit1, index);
			bit2 = extractBit (idQubit2, index);
			bit3 = extractBit (idQubit3, index);
			if (!(bit1 && bit2 && bit3)) {
				probOfFilter+= stateVecReal[index]*stateVecReal[index] + stateVecImag[index]* stateVecImag [index];
			}
		}
	}
	return probOfFilter;
}
