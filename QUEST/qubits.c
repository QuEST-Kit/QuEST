# include "math.h"  //SCB new line
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>

// Maihi: Where I have made changes I have marked SCB so please note those points - Simon

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

void initStateVec (const int numQubits,
                   double *restrict stateVecReal,
                   double *restrict stateVecImag)
{
  long long int stateVecSize;
  long long int index;

  // dimension of the state vector
  stateVecSize = 1LL << numQubits;
  
  printf("stateVecSize=%Ld   now performing init with only one thread:\n",stateVecSize);

  // initialise the state to |0000..0000>
  # ifdef _OPENMP
    # pragma omp parallel for \
      default  (none) \
      shared   (stateVecSize, stateVecReal,stateVecImag) \
      private  (index) \
      schedule (static)
  # endif
  for (index=1; index<stateVecSize; index++) {
    stateVecReal [index] = 0.0;
    stateVecImag [index] = 0.0;
  }

  // zero state |0000..0000> has probability 1
  stateVecReal [0] = 1.0;
  stateVecImag [0] = 0.0;
  
  printf("COMPLETED INIT\n");
}


// ==================================================================== //
//                                                                      //
//     rotateQubit -- routine to rotate a single qubit in the state     //
//                    vector of probability akmplitudes, given the      //
//                    angle rotation arguments.                         //
//                                                                      //
//     input:                                                           //
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

void rotateQubit (const int numQubits, const int rotQubit,
                  double alphaReal, double alphaImag,
                  double betaReal,  double betaImag,
                  double *restrict stateVecReal, double *restrict stateVecImag)
{
  // ----- sizes
  long long int numBlocks,                                           // number of blocks
           sizeBlock,                                           // size of blocks
           sizeHalfBlock;                                       // size of blocks halved
  // ----- indices
  long long int thisBlock,                                           // current block
           indexUp,indexLo;                                     // current index and corresponding index in lower half block

  // ----- temp variables
  double   stateRealUp,stateRealLo,                             // storage for previous state values
           stateImagUp,stateImagLo;                             // (used in updates)
  // ----- temp variables
  long long int thisTask,numTasks;                                   // task based approach for expose loop with small granularity
                                                                // (good for shared memory parallelism)


  // ---------------------------------------------------------------- //
  //            tests                                                 //
  // ---------------------------------------------------------------- //
  //SCB please see my comment on assert() in the function controlPhaseGate below; in short I thin this should be && not ||
  assert (rotQubit >= 0 || rotQubit < numQubits);


  // ---------------------------------------------------------------- //
  //            dimensions                                            //
  // ---------------------------------------------------------------- //
  numBlocks     = 1LL << (numQubits - rotQubit - 1);             // number of blocks in this rotation
  sizeHalfBlock = 1LL << rotQubit;                               // size of blocks halved
  sizeBlock     = 2LL* sizeHalfBlock;                           // size of blocks


  // ---------------------------------------------------------------- //
  //            rotate                                                //
  // ---------------------------------------------------------------- //
  numTasks = 1LL << (numQubits-1);   //SCB there was an error here, "1 << (..." rather than 1LL, so it failed for numQubits>31

  //
  // --- task-based shared-memory parallel implementation
  //
  # ifdef _OPENMP
    # pragma omp parallel for \
      default  (none) \
      shared   (numTasks,sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag) \
      private  (thisTask,thisBlock ,indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo) \
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


  /* === begin comment block === */
  /* --- the following lines are the original block-based iteration */
  /* --- it was replaced by the task-based iteration above, which is appropriate for */
  /*     the shared-memory parallelism as it exposes a large number of independent */
  /*     tasks to a smaller number of threads (the number of blocks is too small a */
  /*     value to parallelise this way) */
  /* --- however, block-based itaretions can be important in the distributed MPI */
  /*     parallel implementation */
  /* long int iHalfBlockStart,iHalfBlockEnd;                    // indices to delimit the current half block*/
  /* for (thisBlock=0; thisBlock<numBlocks; thisBlock++) {                   // for all blocks */
  /*   iHalfBlockStart = thisBlock*sizeBlock;                         // index for start of the first half of the block */
  /*   iHalfBlockEnd   = iHalfBlockStart + sizeHalfBlock - 1;      // index for end of the first half of the block */

  /*   for (indexUp=iHalfBlockStart; indexUp<=iHalfBlockEnd; indexUp++) { */
  /*     indexLo = indexUp + sizeHalfBlock; */

  /* stateRealUp = stateReal[indexUp]; */
  /* stateImagUp = stateImag[indexUp]; */

  /* stateRealLo = stateReal[indexLo]; */
  /* stateImagLo = stateImag[indexLo]; */

  /* // state[indexUp] = alpha * s[indexUp] - conj(beta) * s[indexLo] */
  /* stateReal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp - betaReal*stateRealLo - betaImag*stateImagLo; */
  /* stateImag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp - betaReal*stateImagLo + betaImag*stateRealLo; */

  /* // state[indexLo] = beta * s[indexUp] + conj(alpha) * s[indexLo] */
  /* stateReal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp + alphaReal*stateRealLo + alphaImag*stateImagLo; */
  /* stateImag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp + alphaReal*stateImagLo - alphaImag*stateRealLo; */
  /* } */
  /* }                                                             // end for loop */
  /* === end comment block === */

  
}                                                               // end of function definition


// ==================================================================== //
//                                                                      //
//     findProbabilityOfZero -- routine to measure the probability      //
//                              of a specified qubit in zero state.     //
//                                                                      //
//     input:                                                           //
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

double findProbabilityOfZero (const int numQubits,
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
  double   totalProbability;                                    // probability (returned) value
  // ----- temp variables
  long long int thisTask,numTasks;                                   // task based approach for expose loop with small granularity
                                                                // (good for shared memory parallelism)

  // ---------------------------------------------------------------- //
  //            tests                                                 //
  // ---------------------------------------------------------------- //
  //SCB please see my comment on assert() in the function controlPhaseGate below; in short I thin this should be && rather than ||
  assert (measureQubit >= 0 || measureQubit < numQubits);


  // ---------------------------------------------------------------- //
  //            dimensions                                            //
  // ---------------------------------------------------------------- //
    // -- ***** SCB fix applied here! (measureQubit) replaces earlier (measureQubit-1) ----------- //
  sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
                                                                // and then the number to skip
  sizeBlock     = 2LL* sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

  // ---------------------------------------------------------------- //
  //            find probability                                      //
  // ---------------------------------------------------------------- //
  numTasks = 1LL << (numQubits-1);

  // initialise returned value
  totalProbability = 0.0;

  // initialise correction for kahan summation
  printf("sizeHalfBlock=%Ld sizeBlock=%Ld numTasks=%Ld\n",sizeHalfBlock,sizeBlock,numTasks);
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
