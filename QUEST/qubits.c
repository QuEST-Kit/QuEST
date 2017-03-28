# include "math.h"  //SCB new line
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include "qubits.h"

// Maihi: Where I have made changes I have marked SCB so please note those points - Simon

void allocCircuit(Circuit *circuit, int numQubits, int rank, int numRanks)
{
	long long int numAmps = 1L << numQubits;
	long long int numAmpsPerRank = numAmps/numRanks;

        circuit->stateVec.real = malloc(numAmpsPerRank * sizeof(circuit->stateVec.real));
        circuit->stateVec.imag = malloc(numAmpsPerRank * sizeof(circuit->stateVec.imag));
        circuit->pairStateVec.real = malloc(numAmpsPerRank * sizeof(circuit->pairStateVec.real));
        circuit->pairStateVec.imag = malloc(numAmpsPerRank * sizeof(circuit->pairStateVec.imag));

        if ( (!(circuit->stateVec.real) || !(circuit->stateVec.imag)
		 || !(circuit->pairStateVec.real) || !(circuit->pairStateVec.imag))
		 && numAmpsPerRank ) {
                printf("Could not allocate memory!");
                exit (EXIT_FAILURE);
        }

	circuit->numQubits = numQubits;
	circuit->numAmps = numAmpsPerRank;
	circuit->chunkId = rank;
}

void freeCircuit(Circuit *circuit){
	free(circuit->stateVec.real);
	free(circuit->stateVec.imag);
	free(circuit->pairStateVec.real);
	free(circuit->pairStateVec.imag);
}

void reportState(Circuit circuit){
	FILE *state;
	char filename[100];
	long long int index;
	sprintf(filename, "state_rank_%d.csv", circuit.chunkId);
	state = fopen(filename, "w");
	if (circuit.chunkId==0) fprintf(state, "real, imag\n");

	for(index=0; index<circuit.numAmps; index++){
		fprintf(state, "%.12f, %.12f\n", circuit.stateVec.real[index], circuit.stateVec.imag[index]);
	}
	fclose(state);
}

// ==================================================================== //
//                                                                      //
//     initStateVec -- routine to initialise the state vector of        //
//                  probability amplitudes to the zero |0000..0000>     //
//                                                                      //
//     input:                                                           //
//                    numQubits     -- number of qubits                 //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector                 //
//                                                                      //
//     output:                                                          //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector (overwritten)   //
//                                                                      //
// ==================================================================== //

void initStateVec (Circuit *circuit)
{
	long long int stateVecSize;
	long long int index;

	// dimension of the state vector
	stateVecSize = circuit->numAmps;

	printf("stateVecSize=%Ld   now performing init with only one thread:\n",stateVecSize);

	// Can't use circuit->stateVec as a private OMP var
	double *stateVecReal = circuit->stateVec.real;
	double *stateVecImag = circuit->stateVec.imag;

	// initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel for \
	default  (none) \
	shared   (stateVecSize, stateVecReal, stateVecImag) \
	private  (index) \
	schedule (static)
# endif
	for (index=0; index<stateVecSize; index++) {
		stateVecReal[index] = 0.0;
		stateVecImag[index] = 0.0;
	}

	if (circuit->chunkId==0){
		// zero state |0000..0000> has probability 1
		stateVecReal[0] = 1.0;
		stateVecImag[0] = 0.0;
	}

	printf("COMPLETED INIT\n");
}

// ==================================================================== //
//                                                                      //
//     rotateQubitLocal -- routine to rotate a single qubit in the state//
//                    vector of probability akmplitudes, given the      //
//                    angle rotation arguments.                         //
//                                                                      //
//     input:                                                           //
//                    numTasks      -- number of amps in this chunk     //
//                    numQubits     -- number of qubits                 //
//                    rotQubit      -- qubit to rotate                  //
//                    alphaReal,    -- real/imag part of                //
//                    alphaImag        rotation angle alpha             //
//                    betaReal,     -- real/imag part of                //
//                    betaImag         rotation angle beta              //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector                 //
//                                                                      //
//     output:                                                          //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector (overwritten)   //
//                                                                      //
//     note:                                                            //
//                    qubits are zero-based and the                     //
//                    the first qubit is the rightmost                  //
//                                                                      //
//                    alphaRe = cos(angle1) * cos(angle2);              //
//                    alphaIm = cos(angle1) * sin(angle2);              //
//                    betaRe  = sin(angle1) * cos(angle3);              //
//                    betaIm  = sin(angle1) * sin(angle3);              //
//                                                                      //
// ==================================================================== //

void rotateQubitLocal (Circuit *circuit, const int rotQubit,
		double alphaReal, double alphaImag,
		double betaReal,  double betaImag)
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
	const long long int numTasks=circuit->numAmps>>1;
	// (good for shared memory parallelism)


	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (rotQubit >= 0 && rotQubit < circuit->numQubits);


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
	
	// Can't use circuit->stateVec as a private OMP var
	double *stateVecReal = circuit->stateVec.real;
	double *stateVecImag = circuit->stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
	default  (none) \
	shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag) \
	private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) 
# endif
	{
# ifdef _OPENMP
# pragma omp for \
		schedule (static)
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

// ==================================================================== 
// chunkIsUpper: returns whether a given chunk in position chunkId is in the upper or lower half of
// a block
//
// inputs: 
//      chunkId -- id of chunk in state vector
//      chunkSize -- number of amps in chunk
//      rotQubit -- qubit being rotated 
// ==================================================================== 

int chunkIsUpper(int chunkId, int chunkSize, int rotQubit)
{
	long long int sizeHalfBlock = 1LL << (rotQubit);
	long long int sizeBlock = sizeHalfBlock*2;
	long long int posInBlock = (chunkId*chunkSize) % sizeBlock;
	return posInBlock<sizeHalfBlock;
}

// ==================================================================== 
// getAlphaBeta: get rotation values for a given chunk
//
// inputs:
//      chunkIsUpper -- 1: chunk is in upper half of block, 0: chunk is in lower half
//
//      rot1Real, rot1Imag, rot2Real, rot2Imag -- rotation values to use, allocated for upper/lower such that
//      :: stateUpper = rot1 * stateUpper + conj(rot2)  * stateLower ::
//      or
//      :: stateLower = rot1 * stateUpper + conj(rot2)  * stateLower ::
//
//      aReal, aImag, bReal, bImag -- initial rotation values 
//      
// ==================================================================== 
void getAlphaBeta(int chunkIsUpper, double *rot1Real, double *rot1Imag, double *rot2Real, double *rot2Imag,
		double aReal, double aImag, double bReal, double bImag)
{
	if (chunkIsUpper){
		*rot1Real = aReal;
		*rot1Imag = aImag;
		*rot2Real = -bReal;
		*rot2Imag = -bImag;
	} else {
		*rot1Real = bReal;
		*rot1Imag = bImag;
		*rot2Real = aReal;
		*rot2Imag = aImag;
	}
}
// ==================================================================== 
// getChunkPairId: get position of corresponding chunk, holding values required to
// update values in chunk at chunkId
//
// inputs:
//      chunkIsUpper -- 1: chunk is in upper half of block, 0: chunk is in lower half
//      chunkId -- id of chunk in state vector
//      chunkSize -- number of amps in chunk
//      rotQubit -- qubit being rotated 
// ==================================================================== 

int getChunkPairId(int chunkIsUpper, int chunkId, int chunkSize, int rotQubit)
{
	long long int sizeHalfBlock = 1LL << (rotQubit);
	int chunksPerHalfBlock = sizeHalfBlock/chunkSize;
	if (chunkIsUpper){
		return chunkId + chunksPerHalfBlock;
	} else {
		return chunkId - chunksPerHalfBlock;
	}
}

// ==================================================================== 
// halfMatrixBlockFitsInChunk: return whether the current qubit rotation will use
// blocks that fit within a single chunk
// inputs:
//      chunkSize -- number of amps in chunk
//      rotQubit -- qubit being rotated 
// ==================================================================== 

int halfMatrixBlockFitsInChunk(int chunkSize, int rotQubit)
{
	long long int sizeHalfBlock = 1LL << (rotQubit);
	if (chunkSize > sizeHalfBlock) return 1;
	else return 0;
}

// ====================================================================         //
//                                                                              //
//     rotateQubitDistributed -- routine to rotate a single qubit in the state  //
//                    vector of probability akmplitudes, given the              //
//                    angle rotation arguments, for a distributed version where //
//                    upper and lower values are stored seperately              //
//                                                                              //
//     input:                                                                   //
//                    numTasks      -- num amps handled by one processor        //
//                    numQubits     -- number of qubits                         //
//                    rotQubit      -- qubit to rotate                          //
//                    alphaReal,    -- real/imag part of                        //
//                    alphaImag        rotation angle alpha                     //
//                    betaReal,     -- real/imag part of                        //
//                    betaImag         rotation angle beta                      //
//                    stateVecRealUp, -- real/imag parts of                     //
//                    stateVecImagUp     the state vector in the upper half     //
//                                       of a block                             //
//                    stateVecRealLo, -- real/imag parts of                     //
//                    stateVecImagLo     the state vector in the lower half     //
//                                       of a block                             //
//                    stateVecRealLo, -- real/imag parts of                     //
//                    stateVecImagLo     the output state vector                //
//                                                                      //
//     output:                                                          //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector (overwritten)   //
//                                                                      //
//     note:                                                            //
//                    qubits are zero-based and the                     //
//                    the first qubit is the rightmost                  //
//                                                                      //
//                    alphaRe = cos(angle1) * cos(angle2);              //
//                    alphaIm = cos(angle1) * sin(angle2);              //
//                    betaRe  = sin(angle1) * cos(angle3);              //
//                    betaIm  = sin(angle1) * sin(angle3);              //
//                                                                      //
// ==================================================================== //

void rotateQubitDistributed (Circuit *circuit, const int rotQubit,
		double rot1Real, double rot1Imag,
		double rot2Real,  double rot2Imag,
		double *stateVecRealUp, double *stateVecImagUp,
		double *stateVecRealLo, double *stateVecImagLo,
		double *stateVecRealOut, double *stateVecImagOut)
{
	// ----- temp variables
	double   stateRealUp,stateRealLo,                             // storage for previous state values
	stateImagUp,stateImagLo;                             // (used in updates)
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	const long long int numTasks=circuit->numAmps;

	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (rotQubit >= 0 && rotQubit < circuit->numQubits);

	// ---------------------------------------------------------------- //
	//            rotate                                                //
	// ---------------------------------------------------------------- //

	//
	// --- task-based shared-memory parallel implementation
	//
# pragma omp parallel \
	default  (none) \
	shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
			rot1Real,rot1Imag, rot2Real,rot2Imag) \
	private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
	{
# pragma omp for \
		schedule (static)
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

// ==================================================================== 
// isChunkToSkipInFindPZero: When calculating probability of a bit q being zero,
// sum up 2^q values, then skip 2^q values, etc. This function finds if an entire chunk
// is in the range of values to be skipped
//
// inputs:
//
//      chunkId -- id of chunk in state vector
//      chunkSize -- number of amps in chunk
//      measureQubit -- qubit being measured
//
// outputs:
//      int -- 1: skip, 0: don't skip
// ==================================================================== 

int isChunkToSkipInFindPZero(int chunkId, int chunkSize, int measureQubit){
	long long int sizeHalfBlock = 1LL << (measureQubit);
	int numChunksToSkip = sizeHalfBlock/chunkSize;
	// calculate probability by summing over numChunksToSkip, then skipping numChunksToSkip, etc
	int bitToCheck = chunkId & numChunksToSkip;
	return bitToCheck;
}


// ==================================================================== //
//                                                                      //
//     findProbabilityOfZeroLocal -- routine to measure the probabilityi//
//                              of a specified qubit in zero state.     // 
//                              Size of regions to skip is less than    //
//                              the size of one chunk                   //
//                                                                      //
//     input:                                                           //
//                    numTasks      -- number of amps in this chunk     //
//                    numQubits     -- number of qubits                 //
//                    measureQubit  -- qubit to measure                 //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector                 //
//                                                                      //
//     output:                                                          //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector (overwritten)   //
//                                                                      //
//     note:                                                            //
//                                                                      //
// ==================================================================== //

double findProbabilityOfZeroLocal (Circuit *circuit,
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
	long long int numTasks=circuit->numAmps>>1;
	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (measureQubit >= 0 && measureQubit < circuit->numQubits);


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
	printf("sizeHalfBlock=%Ld sizeBlock=%Ld numTasks=%Ld\n",sizeHalfBlock,sizeBlock,numTasks);

	//
	// --- task-based shared-memory parallel implementation
	//
	
	double *stateVecReal = circuit->stateVec.real;
	double *stateVecImag = circuit->stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel for \
	shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
	private   (thisTask,thisBlock,index) \
	schedule  (static) \
	reduction ( +:totalProbability )
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

	return totalProbability;
}

// ==================================================================== //
//                                                                      //
//     findProbabilityOfZeroDistributed                                 // 
//                              -- routine to measure the probability   //
//                              of a specified qubit in zero state.     // 
//                              Size of regions to skip is a multiple   //
//                              of chunkSize                            //
//                                                                      //
//     input:                                                           //
//                    numTasks      -- number of amps in this chunk     //
//                    numQubits     -- number of qubits                 //
//                    measureQubit  -- qubit to measure                 //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector                 //
//                                                                      //
//     output:                                                          //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector (overwritten)   //
//                                                                      //
//     note:                                                            //
//                                                                      //
// ==================================================================== //

double findProbabilityOfZeroDistributed (Circuit *circuit,
		const int measureQubit)
{
	// ----- measured probability
	double   totalProbability;                                    // probability (returned) value
	// ----- temp variables
	long long int thisTask;                                   // task based approach for expose loop with small granularity
	long long int numTasks=circuit->numAmps;
	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	assert (measureQubit >= 0 && measureQubit < circuit->numQubits);

	// ---------------------------------------------------------------- //
	//            find probability                                      //
	// ---------------------------------------------------------------- //

	// initialise returned value
	totalProbability = 0.0;

	// initialise correction for kahan summation

	//
	// --- task-based shared-memory parallel implementation
	//
	
	double *stateVecReal = circuit->stateVec.real;
	double *stateVecImag = circuit->stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel for \
	shared    (numTasks,stateVecReal,stateVecImag) \
	private   (thisTask) \
	schedule  (static) \
	reduction ( +:totalProbability )
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

	return totalProbability;
}

// ==================================================================== //
//                                                                      //
//     controlPhaseGate -- routine to implement the control phase       //
//                         (the two qubit phase gate)                   //
//                                                                      //
//     input:                                                           //
//                    numQubits     -- number of qubits                 //
//                    idQubit1,     -- specified qubits                 //
//                    idQubit2                                          //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector                 //
//                                                                      //
//     output:                                                          //
//                    stateVecReal, -- real/imag parts of               //
//                    stateVecImag     the state vector (overwritten)   //
//                                                                      //
//     note:                                                            //
//                                                                      //
// ==================================================================== //

// *** SCB edit: new definition of extractBit is much faster ***
int extractBit (const int locationOfBitFromRight, const long long int theEncodedNumber)
{
	return (theEncodedNumber & ( 1LL << locationOfBitFromRight )) >> locationOfBitFromRight;
}

void controlPhaseGate (const int numQubits, const int idQubit1, const int idQubit2,
		double *restrict stateVecReal, double *restrict stateVecImag)
{
	long long int index;
	long long int stateVecSize;
	int bit1, bit2;

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //

	// SCB -- I have made no change, but I have a question: 
	// doesn't assert terminate the program if its argument is false? And in that case, why
	// are we using OR operations between conditions that are ALL essential? This assert will
	// be passed unless ALL of these errors are present! 
	// Don't we need (condit 1) && (condit 2) && ... (condit n) ??
	assert (idQubit1 >= 0 || idQubit2 >= 0 || idQubit1 < numQubits || idQubit2 < numQubits);


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
	//SCB please see my comment on assert() in the function controlPhaseGate above. 
	//I think the line below, for the current function, should be eight conditions with && between
	//line idQubit1 >= 0 && idQubit1 < numQubits && idQubit2 >= 0 ... idQubit4 < numQubits
	//but I have not made any change
	assert (idQubit1 >= 0 || idQubit2 >= 0 || idQubit1 < numQubits || idQubit2 < numQubits);

	stateVecSize = 1LL << numQubits;

# ifdef _OPENMP
# pragma omp parallel for \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag ) \
	private  (index,bit1,bit2,bit3,bit4)		       \
	schedule (static)
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


// measure in Zero performs an irreversible change to the state vector: it updates the vector according
// to the event that a zero have been measured on the qubit indicated by measureQubit (where 
// this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
// then renormalising. It also returns the probability that this event would happen.

double measureInZero (const int numQubits, 
		const int measureQubit,
		double *restrict stateVecReal,
		double *restrict stateVecImag)
{
	// ----- sizes
	long long int numBlocks,                                           // number of blocks
	sizeBlock,                                           // size of blocks
	sizeHalfBlock;                                       // size of blocks halved
	// ----- indices
	long long int thisBlock,                                           // current block
	     index;                                               // current index for first half block
	// ----- measured probability
	double   totalProbability, renorm;                                    // probability (returned) value
	// ----- temp variables
	long long int thisTask,numTasks;                                   // task based approach for expose loop with small granularity
	// (good for shared memory parallelism)

	// ---------------------------------------------------------------- //
	//            tests                                                 //
	// ---------------------------------------------------------------- //
	//SCB please see my comment on assert() in the function controlPhaseGate above.
	assert (measureQubit >= 0 || measureQubit < numQubits);


	// ---------------------------------------------------------------- //
	//            dimensions                                            //
	// ---------------------------------------------------------------- //
	sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
	// and then the number to skip
	sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

	// ---------------------------------------------------------------- //
	//            find probability                                      //
	// ---------------------------------------------------------------- //
	numTasks = 1LL << (numQubits-1);

	// initialise returned value
	totalProbability = 0.0;

	//
	// --- task-based shared-memory parallel implementation
	//
# ifdef _OPENMP
# pragma omp parallel for \
	shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
	private   (thisTask,thisBlock,index) \
	schedule  (static) \
	reduction ( +:totalProbability )
# endif
	for (thisTask=0; thisTask<numTasks; thisTask++) {
		thisBlock = thisTask / sizeHalfBlock;
		index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;

		totalProbability += stateVecReal[index]*stateVecReal[index]
			+ stateVecImag[index]*stateVecImag[index];
	}
	renorm=1/sqrt(totalProbability);


# ifdef _OPENMP
# pragma omp parallel for \
	shared    (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag) \
	private   (thisTask,thisBlock,index) \
	schedule  (static) \
	reduction ( +:totalProbability )
# endif
	for (thisTask=0; thisTask<numTasks; thisTask++) {
		thisBlock = thisTask / sizeHalfBlock;
		index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
		stateVecReal[index]=stateVecReal[index]*renorm;
		stateVecImag[index]=stateVecImag[index]*renorm;

		stateVecReal[index+sizeHalfBlock]=0;
		stateVecImag[index+sizeHalfBlock]=0;
	}

	//SCB this is a debugging style check. It is probably useful to leave in, but it could be parallelised I guess
	//  double checkTotal=1.;
	//  for (index=0; index<2*numTasks; index++) {
	//  	checkTotal=checkTotal-(stateVecReal[index]*stateVecReal[index] + stateVecImag[index]*stateVecImag[index]);
	//  }
	//  if (checkTotal>0.00001){printf("Deviation of sum squared amps from unity is %.16f\n",checkTotal); exit(1);}

	return totalProbability;
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
	//SCB please see my comment on assert() in the function controlPhaseGate above. 
	//I think the line below, for the current function, should be eight conditions with && between
	//line idQubit1 >= 0 && idQubit1 < numQubits && idQubit2 >= 0 ... idQubit4 < numQubits
	//but I have not made any change
	assert (idQubit1 >= 0 || idQubit2 >= 0 || idQubit1 < numQubits || idQubit2 < numQubits);

	stateVecSize = 1LL << numQubits;
	double probOfFilter=0;

# ifdef _OPENMP
# pragma omp parallel for \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag ) \
	private  (index,bit1,bit2,bit3)		       \
	schedule (static)\
	reduction ( +:probOfFilter )
# endif
	for (index=0; index<stateVecSize; index++) {
		bit1 = extractBit (idQubit1, index);
		bit2 = extractBit (idQubit2, index);
		bit3 = extractBit (idQubit3, index);
		if (!(bit1 && bit2 && bit3)) {
			probOfFilter+= stateVecReal[index]*stateVecReal[index] + stateVecImag[index]* stateVecImag [index];
		}
	}
	if ( probOfFilter<1e-16 ){ printf("Extremely small or negative profOfFilter=%.8e; aborting! \n",probOfFilter); exit(1);}
	double myNorm=1/sqrt(probOfFilter);

# ifdef _OPENMP
# pragma omp parallel for \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag, myNorm ) \
	private  (index,bit1,bit2,bit3)		       \
	schedule (static)
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
	//SCB please see my comment on assert() in the function controlPhaseGate above. 
	//I think the line below, for the current function, should be eight conditions with && between
	//line idQubit1 >= 0 && idQubit1 < numQubits && idQubit2 >= 0 ... idQubit4 < numQubits
	//but I have not made any change
	assert (idQubit1 >= 0 || idQubit2 >= 0 || idQubit1 < numQubits || idQubit2 < numQubits);

	stateVecSize = 1LL << numQubits;
	double probOfFilter=0;

# ifdef _OPENMP
# pragma omp parallel for \
	default  (none)			     \
	shared   (stateVecSize, stateVecReal,stateVecImag ) \
	private  (index,bit1,bit2,bit3)		       \
	schedule (static)\
	reduction ( +:probOfFilter )
# endif
	for (index=0; index<stateVecSize; index++) {
		bit1 = extractBit (idQubit1, index);
		bit2 = extractBit (idQubit2, index);
		bit3 = extractBit (idQubit3, index);
		if (!(bit1 && bit2 && bit3)) {
			probOfFilter+= stateVecReal[index]*stateVecReal[index] + stateVecImag[index]* stateVecImag [index];
		}
	}
	return probOfFilter;
}
