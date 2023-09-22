// Distributed under MIT licence. See https://github.com/QuEST-Kit/QuEST/blob/master/LICENCE.txt for details 

/** @file
 * The core of the CPU backend functionality. The CPU/MPI implementations of the pure state functions in
 * ../QuEST_ops_pure.h are in QuEST_cpu_local.c and QuEST_cpu_distributed.c which mostly wrap the core
 * functions defined here. Some additional hardware-agnostic functions are defined here
 *
 * @author Ania Brown
 * @author Tyson Jones
 * @author Balint Koczor
 */

# include "QuEST.h"
# include "QuEST_internal.h"
# include "QuEST_precision.h"
# include "QuEST_validation.h"
# include "mt19937ar.h"

# include "QuEST_cpu_internal.h"

# include <math.h>  
# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <assert.h>

# ifdef _OPENMP
# include <omp.h>
# endif

/* to support MSVC, we must remove the use of VLA in multiQubtUnitary.
 * We'll instead create stack arrays use _malloca
 */
#ifdef _WIN32
    #include <malloc.h>
#endif


/*
 * overloads for consistent API with GPU 
 */

void copyStateToGPU(Qureg qureg) {
}

void copyStateFromGPU(Qureg qureg) {
}

void statevec_copySubstateToGPU(Qureg qureg, long long int startInd, long long int numAmps) {
}

void statevec_copySubstateFromGPU(Qureg qureg, long long int startInd, long long int numAmps) {
}


/*
 * state vector and density matrix operations
 */

void densmatr_oneQubitDegradeOffDiagonal(Qureg qureg, int targetQubit, qreal retain){
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMask = 1LL << targetQubit;
    long long int outerMask = 1LL << (targetQubit + (qureg.numQubitsRepresented));

    long long int thisTask;
    long long int thisPattern;
    long long int totMask = innerMask|outerMask;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (innerMask,outerMask,totMask,qureg,retain,numTasks, targetQubit) \
    private  (thisTask,thisPattern)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++){
            thisPattern = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMask;
            if ((thisPattern==innerMask) || (thisPattern==outerMask)){
                // do dephase
                // the lines below will degrade the off-diagonal terms |..0..><..1..| and |..1..><..0..|
                qureg.stateVec.real[thisTask] = retain*qureg.stateVec.real[thisTask]; 
                qureg.stateVec.imag[thisTask] = retain*qureg.stateVec.imag[thisTask]; 
            } 
        }  
    }
}

void densmatr_mixDephasing(Qureg qureg, int targetQubit, qreal dephase) {
    qreal retain=1-dephase;
    densmatr_oneQubitDegradeOffDiagonal(qureg, targetQubit, retain);
}

void densmatr_mixTwoQubitDephasing(Qureg qureg, int qubit1, int qubit2, qreal dephase) {
    qreal retain=1-dephase;

    long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMaskQubit1 = 1LL << qubit1;
    long long int outerMaskQubit1 = 1LL << (qubit1 + (qureg.numQubitsRepresented));
    long long int innerMaskQubit2 = 1LL << qubit2;
    long long int outerMaskQubit2 = 1LL << (qubit2 + (qureg.numQubitsRepresented));
    long long int totMaskQubit1 = innerMaskQubit1|outerMaskQubit1;
    long long int totMaskQubit2 = innerMaskQubit2|outerMaskQubit2;

    long long int thisTask;
    long long int thisPatternQubit1, thisPatternQubit2;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (innerMaskQubit1,outerMaskQubit1,totMaskQubit1,innerMaskQubit2,outerMaskQubit2, \
                totMaskQubit2,qureg,retain,numTasks) \
    private  (thisTask,thisPatternQubit1,thisPatternQubit2)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++){
            thisPatternQubit1 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit1;
            thisPatternQubit2 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit2;
            
            // any mismatch |...0...><...1...| etc
            if ( (thisPatternQubit1==innerMaskQubit1) || (thisPatternQubit1==outerMaskQubit1) || 
                    (thisPatternQubit2==innerMaskQubit2) || (thisPatternQubit2==outerMaskQubit2) ){ 
                // do dephase
                // the lines below will degrade the off-diagonal terms |..0..><..1..| and |..1..><..0..|
                qureg.stateVec.real[thisTask] = retain*qureg.stateVec.real[thisTask]; 
                qureg.stateVec.imag[thisTask] = retain*qureg.stateVec.imag[thisTask]; 
            } 
        }  
    }
}

void densmatr_mixDepolarisingLocal(Qureg qureg, int targetQubit, qreal depolLevel) {
    qreal retain=1-depolLevel;

    long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMask = 1LL << targetQubit;
    long long int outerMask = 1LL << (targetQubit + (qureg.numQubitsRepresented));
    long long int totMask = innerMask|outerMask;

    long long int thisTask;
    long long int partner;
    long long int thisPattern;

    qreal realAv, imagAv;
        
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (innerMask,outerMask,totMask,qureg,retain,depolLevel,numTasks) \
    private  (thisTask,partner,thisPattern,realAv,imagAv)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++){
            thisPattern = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMask;
            if ((thisPattern==innerMask) || (thisPattern==outerMask)){
                // do dephase
                // the lines below will degrade the off-diagonal terms |..0..><..1..| and |..1..><..0..|
                qureg.stateVec.real[thisTask] = retain*qureg.stateVec.real[thisTask]; 
                qureg.stateVec.imag[thisTask] = retain*qureg.stateVec.imag[thisTask]; 
            } else {
                if ((thisTask&totMask)==0){ //this element relates to targetQubit in state 0
                    // do depolarise
                    partner = thisTask | totMask;
                    realAv =  (qureg.stateVec.real[thisTask] + qureg.stateVec.real[partner]) /2 ;
                    imagAv =  (qureg.stateVec.imag[thisTask] + qureg.stateVec.imag[partner]) /2 ;
                    
                    qureg.stateVec.real[thisTask] = retain*qureg.stateVec.real[thisTask] + depolLevel*realAv;
                    qureg.stateVec.imag[thisTask] = retain*qureg.stateVec.imag[thisTask] + depolLevel*imagAv;
                    
                    qureg.stateVec.real[partner] = retain*qureg.stateVec.real[partner] + depolLevel*realAv;
                    qureg.stateVec.imag[partner] = retain*qureg.stateVec.imag[partner] + depolLevel*imagAv;
                }
            }
        }  
    }
}

void densmatr_mixDampingLocal(Qureg qureg, int targetQubit, qreal damping) {
    qreal retain=1-damping;
    qreal dephase=sqrt(retain);

    long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMask = 1LL << targetQubit;
    long long int outerMask = 1LL << (targetQubit + (qureg.numQubitsRepresented));
    long long int totMask = innerMask|outerMask;

    long long int thisTask;
    long long int partner;
    long long int thisPattern;

    //qreal realAv, imagAv;
        
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (innerMask,outerMask,totMask,qureg,retain,damping,dephase,numTasks) \
    private  (thisTask,partner,thisPattern)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++){
            thisPattern = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMask;
            if ((thisPattern==innerMask) || (thisPattern==outerMask)){
                // do dephase
                // the lines below will degrade the off-diagonal terms |..0..><..1..| and |..1..><..0..|
                qureg.stateVec.real[thisTask] = dephase*qureg.stateVec.real[thisTask]; 
                qureg.stateVec.imag[thisTask] = dephase*qureg.stateVec.imag[thisTask]; 
            } else {
                if ((thisTask&totMask)==0){ //this element relates to targetQubit in state 0
                    // do depolarise
                    partner = thisTask | totMask;
                    //realAv =  (qureg.stateVec.real[thisTask] + qureg.stateVec.real[partner]) /2 ;
                    //imagAv =  (qureg.stateVec.imag[thisTask] + qureg.stateVec.imag[partner]) /2 ;
                    
                    qureg.stateVec.real[thisTask] = qureg.stateVec.real[thisTask] + damping*qureg.stateVec.real[partner];
                    qureg.stateVec.imag[thisTask] = qureg.stateVec.imag[thisTask] + damping*qureg.stateVec.imag[partner];
                    
                    qureg.stateVec.real[partner] = retain*qureg.stateVec.real[partner];
                    qureg.stateVec.imag[partner] = retain*qureg.stateVec.imag[partner];
                }
            }
        }  
    }
}

void densmatr_mixDepolarisingDistributed(Qureg qureg, int targetQubit, qreal depolLevel) {

    // first do dephase part. 
    // TODO -- this might be more efficient to do at the same time as the depolarise if we move to
    // iterating over all elements in the state vector for the purpose of vectorisation
    // TODO -- if we keep this split, move this function to densmatr_mixDepolarising()
    densmatr_mixDephasing(qureg, targetQubit, depolLevel);

    long long int sizeInnerBlock, sizeInnerHalfBlock;
    long long int sizeOuterColumn, sizeOuterHalfColumn;
    long long int thisInnerBlock, // current block
         thisOuterColumn, // current column in density matrix
         thisIndex,    // current index in (density matrix representation) state vector
         thisIndexInOuterColumn,
         thisIndexInInnerBlock; 
    int outerBit; 

    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeInnerHalfBlock = 1LL << targetQubit;  
    sizeInnerBlock = 2LL * sizeInnerHalfBlock; 
    sizeOuterColumn = 1LL << qureg.numQubitsRepresented;
    sizeOuterHalfColumn = sizeOuterColumn >> 1;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeInnerBlock,sizeInnerHalfBlock,sizeOuterColumn,sizeOuterHalfColumn, \
                qureg,depolLevel,numTasks,targetQubit) \
    private  (thisTask,thisInnerBlock,thisOuterColumn,thisIndex,thisIndexInOuterColumn, \
                thisIndexInInnerBlock,outerBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // thisTask iterates over half the elements in this process' chunk of the density matrix
        // treat this as iterating over all columns, then iterating over half the values
        // within one column.
        // If this function has been called, this process' chunk contains half an 
        // outer block or less
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // we want to process all columns in the density matrix,
            // updating the values for half of each column (one half of each inner block)
            thisOuterColumn = thisTask / sizeOuterHalfColumn;
            thisIndexInOuterColumn = thisTask&(sizeOuterHalfColumn-1); // thisTask % sizeOuterHalfColumn
            thisInnerBlock = thisIndexInOuterColumn/sizeInnerHalfBlock;
            // get index in state vector corresponding to upper inner block
            thisIndexInInnerBlock = thisTask&(sizeInnerHalfBlock-1); // thisTask % sizeInnerHalfBlock
            thisIndex = thisOuterColumn*sizeOuterColumn + thisInnerBlock*sizeInnerBlock 
                + thisIndexInInnerBlock;
            // check if we are in the upper or lower half of an outer block
            outerBit = extractBit(targetQubit, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBit*(sizeInnerHalfBlock);

            // NOTE: at this point thisIndex should be the index of the element we want to 
            // dephase in the chunk of the state vector on this process, in the 
            // density matrix representation. 
            // thisTask is the index of the pair element in pairStateVec

           
            // state[thisIndex] = (1-depolLevel)*state[thisIndex] + depolLevel*(state[thisIndex]
            //      + pair[thisTask])/2
            qureg.stateVec.real[thisIndex] = (1-depolLevel)*qureg.stateVec.real[thisIndex] +
                    depolLevel*(qureg.stateVec.real[thisIndex] + qureg.pairStateVec.real[thisTask])/2;
            
            qureg.stateVec.imag[thisIndex] = (1-depolLevel)*qureg.stateVec.imag[thisIndex] +
                    depolLevel*(qureg.stateVec.imag[thisIndex] + qureg.pairStateVec.imag[thisTask])/2;
        } 
    }    
}

void densmatr_mixDampingDistributed(Qureg qureg, int targetQubit, qreal damping) {
    qreal retain=1-damping;
    qreal dephase=sqrt(1-damping);

    // multiply the off-diagonal (|0><1| and |1><0|) terms by sqrt(1-damping)
    densmatr_oneQubitDegradeOffDiagonal(qureg, targetQubit, dephase);
    
    // below, we modify the diagonals terms which require |1><1| to |0><0| communication

    long long int sizeInnerBlock, sizeInnerHalfBlock;
    long long int sizeOuterColumn, sizeOuterHalfColumn;
    long long int thisInnerBlock, // current block
         thisOuterColumn, // current column in density matrix
         thisIndex,    // current index in (density matrix representation) state vector
         thisIndexInOuterColumn,
         thisIndexInInnerBlock; 
    int outerBit; 
    int stateBit;

    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeInnerHalfBlock = 1LL << targetQubit;  
    sizeInnerBlock = 2LL * sizeInnerHalfBlock; 
    sizeOuterColumn = 1LL << qureg.numQubitsRepresented;
    sizeOuterHalfColumn = sizeOuterColumn >> 1;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeInnerBlock,sizeInnerHalfBlock,sizeOuterColumn,sizeOuterHalfColumn, \
                qureg,damping, retain, dephase, numTasks,targetQubit) \
    private  (thisTask,thisInnerBlock,thisOuterColumn,thisIndex,thisIndexInOuterColumn, \
                thisIndexInInnerBlock,outerBit, stateBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // thisTask iterates over half the elements in this process' chunk of the density matrix
        // treat this as iterating over all columns, then iterating over half the values
        // within one column.
        // If this function has been called, this process' chunk contains half an 
        // outer block or less
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // we want to process all columns in the density matrix,
            // updating the values for half of each column (one half of each inner block)
            thisOuterColumn = thisTask / sizeOuterHalfColumn;
            thisIndexInOuterColumn = thisTask&(sizeOuterHalfColumn-1); // thisTask % sizeOuterHalfColumn
            thisInnerBlock = thisIndexInOuterColumn/sizeInnerHalfBlock;
            // get index in state vector corresponding to upper inner block
            thisIndexInInnerBlock = thisTask&(sizeInnerHalfBlock-1); // thisTask % sizeInnerHalfBlock
            thisIndex = thisOuterColumn*sizeOuterColumn + thisInnerBlock*sizeInnerBlock 
                + thisIndexInInnerBlock;
            // check if we are in the upper or lower half of an outer block
            outerBit = extractBit(targetQubit, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBit*(sizeInnerHalfBlock);

            // NOTE: at this point thisIndex should be the index of the element we want to 
            // dephase in the chunk of the state vector on this process, in the 
            // density matrix representation. 
            // thisTask is the index of the pair element in pairStateVec

            // Extract state bit, is 0 if thisIndex corresponds to a state with 0 in the target qubit
            // and is 1 if thisIndex corresponds to a state with 1 in the target qubit
            stateBit = extractBit(targetQubit, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId));
           
            // state[thisIndex] = (1-depolLevel)*state[thisIndex] + depolLevel*(state[thisIndex]
            //      + pair[thisTask])/2
            if(stateBit == 0){
                qureg.stateVec.real[thisIndex] = qureg.stateVec.real[thisIndex] +
                    damping*( qureg.pairStateVec.real[thisTask]);
                
                qureg.stateVec.imag[thisIndex] = qureg.stateVec.imag[thisIndex] +
                    damping*( qureg.pairStateVec.imag[thisTask]);
            } else{
                qureg.stateVec.real[thisIndex] = retain*qureg.stateVec.real[thisIndex];
            
                qureg.stateVec.imag[thisIndex] = retain*qureg.stateVec.imag[thisIndex];
            }
        } 
    }    
}

void densmatr_mixTwoQubitDepolarisingLocal(Qureg qureg, int qubit1, int qubit2, qreal delta, qreal gamma) {
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMaskQubit1 = 1LL << qubit1;
    long long int outerMaskQubit1= 1LL << (qubit1 + qureg.numQubitsRepresented);
    long long int totMaskQubit1 = innerMaskQubit1 | outerMaskQubit1;
    long long int innerMaskQubit2 = 1LL << qubit2;
    long long int outerMaskQubit2 = 1LL << (qubit2 + qureg.numQubitsRepresented);
    long long int totMaskQubit2 = innerMaskQubit2 | outerMaskQubit2;

    long long int thisTask;
    long long int partner;
    long long int thisPatternQubit1, thisPatternQubit2;

    qreal real00, imag00;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (totMaskQubit1,totMaskQubit2,qureg,delta,gamma,numTasks) \
    private  (thisTask,partner,thisPatternQubit1,thisPatternQubit2,real00,imag00)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        //--------------------------------------- STEP ONE ---------------------
        for (thisTask=0; thisTask<numTasks; thisTask++){ 
            thisPatternQubit1 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit1;
            thisPatternQubit2 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit2;
            if ((thisPatternQubit1==0) && ((thisPatternQubit2==0) 
                        || (thisPatternQubit2==totMaskQubit2))){ 
                //this element of form |...X...0...><...X...0...|  for X either 0 or 1.
                partner = thisTask | totMaskQubit1;
                real00 =  qureg.stateVec.real[thisTask];
                imag00 =  qureg.stateVec.imag[thisTask];
                    
                qureg.stateVec.real[thisTask] = qureg.stateVec.real[thisTask] 
                    + delta*qureg.stateVec.real[partner];
                qureg.stateVec.imag[thisTask] = qureg.stateVec.imag[thisTask] 
                    + delta*qureg.stateVec.imag[partner];
                    
                qureg.stateVec.real[partner] = qureg.stateVec.real[partner] + delta*real00;
                qureg.stateVec.imag[partner] = qureg.stateVec.imag[partner] + delta*imag00;
                                
            }
        }
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        //--------------------------------------- STEP TWO ---------------------
        for (thisTask=0; thisTask<numTasks; thisTask++){ 
            thisPatternQubit1 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit1;
            thisPatternQubit2 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit2;
            if ((thisPatternQubit2==0) && ((thisPatternQubit1==0) 
                        || (thisPatternQubit1==totMaskQubit1))){ 
                //this element of form |...0...X...><...0...X...|  for X either 0 or 1.
                partner = thisTask | totMaskQubit2;
                real00 =  qureg.stateVec.real[thisTask];
                imag00 =  qureg.stateVec.imag[thisTask];
                    
                qureg.stateVec.real[thisTask] = qureg.stateVec.real[thisTask] 
                    + delta*qureg.stateVec.real[partner];
                qureg.stateVec.imag[thisTask] = qureg.stateVec.imag[thisTask] 
                    + delta*qureg.stateVec.imag[partner];
                    
                qureg.stateVec.real[partner] = qureg.stateVec.real[partner] + delta*real00;
                qureg.stateVec.imag[partner] = qureg.stateVec.imag[partner] + delta*imag00;

            }
        }

# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        //--------------------------------------- STEP THREE ---------------------
        for (thisTask=0; thisTask<numTasks; thisTask++){ 
            thisPatternQubit1 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit1;
            thisPatternQubit2 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit2;
            if ((thisPatternQubit2==0) && ((thisPatternQubit1==0) 
                        || (thisPatternQubit1==totMaskQubit1))){ 
                //this element of form |...0...X...><...0...X...|  for X either 0 or 1.
                partner = thisTask | totMaskQubit2;
                partner = partner ^ totMaskQubit1;
                real00 =  qureg.stateVec.real[thisTask];
                imag00 =  qureg.stateVec.imag[thisTask];

                qureg.stateVec.real[thisTask] = gamma * (qureg.stateVec.real[thisTask] 
                        + delta*qureg.stateVec.real[partner]);
                qureg.stateVec.imag[thisTask] = gamma * (qureg.stateVec.imag[thisTask] 
                        + delta*qureg.stateVec.imag[partner]);
                    
                qureg.stateVec.real[partner] = gamma * (qureg.stateVec.real[partner] 
                        + delta*real00);
                qureg.stateVec.imag[partner] = gamma * (qureg.stateVec.imag[partner] 
                        + delta*imag00);

            }
        }
    }
}

void densmatr_mixTwoQubitDepolarisingLocalPart1(Qureg qureg, int qubit1, int qubit2, qreal delta) {
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int innerMaskQubit1 = 1LL << qubit1;
    long long int outerMaskQubit1= 1LL << (qubit1 + qureg.numQubitsRepresented);
    long long int totMaskQubit1 = innerMaskQubit1 | outerMaskQubit1;
    long long int innerMaskQubit2 = 1LL << qubit2;
    long long int outerMaskQubit2 = 1LL << (qubit2 + qureg.numQubitsRepresented);
    long long int totMaskQubit2 = innerMaskQubit2 | outerMaskQubit2;
    // correct for being in a particular chunk
    //totMaskQubit2 = totMaskQubit2&(qureg.numAmpsPerChunk-1); // totMaskQubit2 % numAmpsPerChunk
    

    long long int thisTask;
    long long int partner;
    long long int thisPatternQubit1, thisPatternQubit2;

    qreal real00, imag00;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (totMaskQubit1,totMaskQubit2,qureg,delta,numTasks) \
    private  (thisTask,partner,thisPatternQubit1,thisPatternQubit2,real00,imag00)
# endif
    {

# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        //--------------------------------------- STEP ONE ---------------------
        for (thisTask=0; thisTask<numTasks; thisTask ++){ 
            thisPatternQubit1 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit1;
            thisPatternQubit2 = (thisTask+qureg.numAmpsPerChunk*qureg.chunkId)&totMaskQubit2;
            if ((thisPatternQubit1==0) && ((thisPatternQubit2==0) 
                        || (thisPatternQubit2==totMaskQubit2))){ 
                //this element of form |...X...0...><...X...0...|  for X either 0 or 1.
                partner = thisTask | totMaskQubit1;
                real00 =  qureg.stateVec.real[thisTask];
                imag00 =  qureg.stateVec.imag[thisTask];
                    
                qureg.stateVec.real[thisTask] = qureg.stateVec.real[thisTask] 
                    + delta*qureg.stateVec.real[partner];
                qureg.stateVec.imag[thisTask] = qureg.stateVec.imag[thisTask] 
                    + delta*qureg.stateVec.imag[partner];
                    
                qureg.stateVec.real[partner] = qureg.stateVec.real[partner] + delta*real00;
                qureg.stateVec.imag[partner] = qureg.stateVec.imag[partner] + delta*imag00;
                                
            }
        }
    }
}

void densmatr_mixTwoQubitDepolarisingDistributed(Qureg qureg, int targetQubit, 
        int qubit2, qreal delta, qreal gamma) {

    long long int sizeInnerBlockQ1, sizeInnerHalfBlockQ1;
    long long int sizeInnerBlockQ2, sizeInnerHalfBlockQ2, sizeInnerQuarterBlockQ2;
    long long int sizeOuterColumn, sizeOuterQuarterColumn;
    long long int thisInnerBlockQ2,
         thisOuterColumn, // current column in density matrix
         thisIndex,    // current index in (density matrix representation) state vector
         thisIndexInOuterColumn,
         thisIndexInInnerBlockQ1, 
         thisIndexInInnerBlockQ2, 
         thisInnerBlockQ1InInnerBlockQ2;
    int outerBitQ1, outerBitQ2; 

    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>2;

    // set dimensions
    sizeInnerHalfBlockQ1 = 1LL << targetQubit;  
    sizeInnerHalfBlockQ2 = 1LL << qubit2;  
    sizeInnerQuarterBlockQ2 = sizeInnerHalfBlockQ2 >> 1;  
    sizeInnerBlockQ2 = sizeInnerHalfBlockQ2 << 1;  
    sizeInnerBlockQ1 = 2LL * sizeInnerHalfBlockQ1; 
    sizeOuterColumn = 1LL << qureg.numQubitsRepresented;
    sizeOuterQuarterColumn = sizeOuterColumn >> 2;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeInnerBlockQ1,sizeInnerHalfBlockQ1,sizeInnerBlockQ2,sizeInnerHalfBlockQ2,sizeInnerQuarterBlockQ2,\
                sizeOuterColumn,sizeOuterQuarterColumn,qureg,delta,gamma,numTasks,targetQubit,qubit2) \
    private  (thisTask,thisInnerBlockQ2,thisInnerBlockQ1InInnerBlockQ2, \
                thisOuterColumn,thisIndex,thisIndexInOuterColumn, \
                thisIndexInInnerBlockQ1,thisIndexInInnerBlockQ2,outerBitQ1,outerBitQ2)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // thisTask iterates over half the elements in this process' chunk of the density matrix
        // treat this as iterating over all columns, then iterating over half the values
        // within one column.
        // If this function has been called, this process' chunk contains half an 
        // outer block or less
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // we want to process all columns in the density matrix,
            // updating the values for half of each column (one half of each inner block)
            thisOuterColumn = thisTask / sizeOuterQuarterColumn;
            // thisTask % sizeOuterQuarterColumn
            thisIndexInOuterColumn = thisTask&(sizeOuterQuarterColumn-1); 
            thisInnerBlockQ2 = thisIndexInOuterColumn / sizeInnerQuarterBlockQ2;
            // thisTask % sizeInnerQuarterBlockQ2;
            thisIndexInInnerBlockQ2 = thisTask&(sizeInnerQuarterBlockQ2-1);
            thisInnerBlockQ1InInnerBlockQ2 = thisIndexInInnerBlockQ2 / sizeInnerHalfBlockQ1;
            // thisTask % sizeInnerHalfBlockQ1;
            thisIndexInInnerBlockQ1 = thisTask&(sizeInnerHalfBlockQ1-1);

            // get index in state vector corresponding to upper inner block
            thisIndex = thisOuterColumn*sizeOuterColumn + thisInnerBlockQ2*sizeInnerBlockQ2 
                + thisInnerBlockQ1InInnerBlockQ2*sizeInnerBlockQ1 + thisIndexInInnerBlockQ1;

            // check if we are in the upper or lower half of an outer block for Q1
            outerBitQ1 = extractBit(targetQubit, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBitQ1*(sizeInnerHalfBlockQ1);

            // check if we are in the upper or lower half of an outer block for Q2
            outerBitQ2 = extractBit(qubit2, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBitQ2*(sizeInnerQuarterBlockQ2<<1);

            // NOTE: at this point thisIndex should be the index of the element we want to 
            // dephase in the chunk of the state vector on this process, in the 
            // density matrix representation. 
            // thisTask is the index of the pair element in pairStateVec

           
            // state[thisIndex] = (1-depolLevel)*state[thisIndex] + depolLevel*(state[thisIndex]
            //      + pair[thisTask])/2
            // NOTE: must set gamma=1 if using this function for steps 1 or 2
            qureg.stateVec.real[thisIndex] = gamma*(qureg.stateVec.real[thisIndex] +
                    delta*qureg.pairStateVec.real[thisTask]);
            qureg.stateVec.imag[thisIndex] = gamma*(qureg.stateVec.imag[thisIndex] +
                    delta*qureg.pairStateVec.imag[thisTask]);
        } 
    }    
}

void densmatr_mixTwoQubitDepolarisingQ1LocalQ2DistributedPart3(Qureg qureg, int targetQubit, 
        int qubit2, qreal delta, qreal gamma) {

    long long int sizeInnerBlockQ1, sizeInnerHalfBlockQ1;
    long long int sizeInnerBlockQ2, sizeInnerHalfBlockQ2, sizeInnerQuarterBlockQ2;
    long long int sizeOuterColumn, sizeOuterQuarterColumn;
    long long int thisInnerBlockQ2,
         thisOuterColumn, // current column in density matrix
         thisIndex,    // current index in (density matrix representation) state vector
         thisIndexInPairVector,
         thisIndexInOuterColumn,
         thisIndexInInnerBlockQ1, 
         thisIndexInInnerBlockQ2, 
         thisInnerBlockQ1InInnerBlockQ2;
    int outerBitQ1, outerBitQ2; 

    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>2;

    // set dimensions
    sizeInnerHalfBlockQ1 = 1LL << targetQubit;  
    sizeInnerHalfBlockQ2 = 1LL << qubit2;  
    sizeInnerQuarterBlockQ2 = sizeInnerHalfBlockQ2 >> 1;  
    sizeInnerBlockQ2 = sizeInnerHalfBlockQ2 << 1;  
    sizeInnerBlockQ1 = 2LL * sizeInnerHalfBlockQ1; 
    sizeOuterColumn = 1LL << qureg.numQubitsRepresented;
    sizeOuterQuarterColumn = sizeOuterColumn >> 2;

//# if 0
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeInnerBlockQ1,sizeInnerHalfBlockQ1,sizeInnerBlockQ2,sizeInnerHalfBlockQ2,sizeInnerQuarterBlockQ2,\
                sizeOuterColumn,sizeOuterQuarterColumn,qureg,delta,gamma, numTasks,targetQubit,qubit2) \
    private  (thisTask,thisInnerBlockQ2,thisInnerBlockQ1InInnerBlockQ2, \
                thisOuterColumn,thisIndex,thisIndexInPairVector,thisIndexInOuterColumn, \
                thisIndexInInnerBlockQ1,thisIndexInInnerBlockQ2,outerBitQ1,outerBitQ2)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
//# endif
        // thisTask iterates over half the elements in this process' chunk of the density matrix
        // treat this as iterating over all columns, then iterating over half the values
        // within one column.
        // If this function has been called, this process' chunk contains half an 
        // outer block or less
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            // we want to process all columns in the density matrix,
            // updating the values for half of each column (one half of each inner block)
            thisOuterColumn = thisTask / sizeOuterQuarterColumn;
            // thisTask % sizeOuterQuarterColumn
            thisIndexInOuterColumn = thisTask&(sizeOuterQuarterColumn-1); 
            thisInnerBlockQ2 = thisIndexInOuterColumn / sizeInnerQuarterBlockQ2;
            // thisTask % sizeInnerQuarterBlockQ2;
            thisIndexInInnerBlockQ2 = thisTask&(sizeInnerQuarterBlockQ2-1);
            thisInnerBlockQ1InInnerBlockQ2 = thisIndexInInnerBlockQ2 / sizeInnerHalfBlockQ1;
            // thisTask % sizeInnerHalfBlockQ1;
            thisIndexInInnerBlockQ1 = thisTask&(sizeInnerHalfBlockQ1-1);

            // get index in state vector corresponding to upper inner block
            thisIndex = thisOuterColumn*sizeOuterColumn + thisInnerBlockQ2*sizeInnerBlockQ2 
                + thisInnerBlockQ1InInnerBlockQ2*sizeInnerBlockQ1 + thisIndexInInnerBlockQ1;

            // check if we are in the upper or lower half of an outer block for Q1
            outerBitQ1 = extractBit(targetQubit, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBitQ1*(sizeInnerHalfBlockQ1);
            
            // For part 3 we need to match elements such that (my Q1 != pair Q1) AND (my Q2 != pair Q2)
            // Find correct index in pairStateVector
            thisIndexInPairVector = thisTask + (1-outerBitQ1)*sizeInnerHalfBlockQ1*sizeOuterQuarterColumn -
                outerBitQ1*sizeInnerHalfBlockQ1*sizeOuterQuarterColumn;
            
            // check if we are in the upper or lower half of an outer block for Q2
            outerBitQ2 = extractBit(qubit2, (thisIndex+qureg.numAmpsPerChunk*qureg.chunkId)>>qureg.numQubitsRepresented);
            // if we are in the lower half of an outer block, shift to be in the lower half
            // of the inner block as well (we want to dephase |0><0| and |1><1| only)
            thisIndex += outerBitQ2*(sizeInnerQuarterBlockQ2<<1);
            
            
            // NOTE: at this point thisIndex should be the index of the element we want to 
            // dephase in the chunk of the state vector on this process, in the 
            // density matrix representation. 

           
            // state[thisIndex] = (1-depolLevel)*state[thisIndex] + depolLevel*(state[thisIndex]
            //      + pair[thisIndexInPairVector])/2
            qureg.stateVec.real[thisIndex] = gamma*(qureg.stateVec.real[thisIndex] +
                    delta*qureg.pairStateVec.real[thisIndexInPairVector]);
            
            qureg.stateVec.imag[thisIndex] = gamma*(qureg.stateVec.imag[thisIndex] +
                    delta*qureg.pairStateVec.imag[thisIndexInPairVector]);
        } 
    }    

}


/* Without nested parallelisation, only the outer most loops which call below are parallelised */
void zeroSomeAmps(Qureg qureg, long long int startInd, long long int numAmps) {
    long long int i;
# ifdef _OPENMP
# pragma omp parallel for schedule (static)
# endif
    for (i=startInd; i < startInd+numAmps; i++) {
        qureg.stateVec.real[i] = 0;
        qureg.stateVec.imag[i] = 0;
    }
}
void normaliseSomeAmps(Qureg qureg, qreal norm, long long int startInd, long long int numAmps) {
    long long int i;
# ifdef _OPENMP
# pragma omp parallel for schedule (static)
# endif
    for (i=startInd; i < startInd+numAmps; i++) {
        qureg.stateVec.real[i] /= norm;
        qureg.stateVec.imag[i] /= norm;
    }
}
void alternateNormZeroingSomeAmpBlocks(
    Qureg qureg, qreal norm, int normFirst, 
    long long int startAmpInd, long long int numAmps, long long int blockSize
) {     
    long long int numDubBlocks = numAmps / (2*blockSize);
    long long int blockStartInd;
    
    if (normFirst) {
        long long int dubBlockInd;
# ifdef _OPENMP
# pragma omp parallel for schedule (static) private (blockStartInd)
# endif 
        for (dubBlockInd=0; dubBlockInd < numDubBlocks; dubBlockInd++) {
            blockStartInd = startAmpInd + dubBlockInd*2*blockSize;
            normaliseSomeAmps(qureg, norm, blockStartInd,             blockSize); // |0><0|
            zeroSomeAmps(     qureg,       blockStartInd + blockSize, blockSize);
        }
    } else {
        long long int dubBlockInd;
# ifdef _OPENMP
# pragma omp parallel for schedule (static) private (blockStartInd)
# endif 
        for (dubBlockInd=0; dubBlockInd < numDubBlocks; dubBlockInd++) {
            blockStartInd = startAmpInd + dubBlockInd*2*blockSize;
            zeroSomeAmps(     qureg,       blockStartInd,             blockSize);
            normaliseSomeAmps(qureg, norm, blockStartInd + blockSize, blockSize); // |1><1|
        }
    }
}

/** Renorms (/prob) every | * outcome * >< * outcome * | state, setting all others to zero */
void densmatr_collapseToKnownProbOutcome(Qureg qureg, int measureQubit, int outcome, qreal totalStateProb) {

	// only (global) indices (as bit sequence): '* outcome *(n+q) outcome *q are spared
    // where n = measureQubit, q = qureg.numQubitsRepresented.
    // We can thus step in blocks of 2^q+n, killing every second, and inside the others,
    //  stepping in sub-blocks of 2^q, killing every second.
    // When outcome=1, we offset the start of these blocks by their size.
    long long int innerBlockSize = (1LL << measureQubit);
    long long int outerBlockSize = (1LL << (measureQubit + qureg.numQubitsRepresented));
    
    // Because there are 2^a number of nodes(/chunks), each node will contain 2^b number of blocks,
    // or each block will span 2^c number of nodes. Similarly for the innerblocks.
    long long int locNumAmps = qureg.numAmpsPerChunk;
    long long int globalStartInd = qureg.chunkId * locNumAmps;
    int innerBit = extractBit(measureQubit, globalStartInd);
    int outerBit = extractBit(measureQubit + qureg.numQubitsRepresented, globalStartInd);
    
    // If this chunk's amps are entirely inside an outer block
    if (locNumAmps <= outerBlockSize) {
        
        // if this is an undesired outer block, kill all elems
        if (outerBit != outcome) {
            zeroSomeAmps(qureg, 0, qureg.numAmpsPerChunk);
            return;
        }
        
        // othwerwise, if this is a desired outer block, and also entirely an inner block
        if (locNumAmps <= innerBlockSize) {
            
            // and that inner block is undesired, kill all elems
            if (innerBit != outcome) 
                zeroSomeAmps(qureg, 0, qureg.numAmpsPerChunk);
            // otherwise normalise all elems
            else
                normaliseSomeAmps(qureg, totalStateProb, 0, qureg.numAmpsPerChunk);
                
            return;
        }
                
        // otherwise this is a desired outer block which contains 2^a inner blocks; kill/renorm every second inner block
        alternateNormZeroingSomeAmpBlocks(
            qureg, totalStateProb, innerBit==outcome, 0, qureg.numAmpsPerChunk, innerBlockSize);
        return;
    }
    
    // Otherwise, this chunk's amps contain multiple outer blocks (and hence multiple inner blocks)
    long long int numOuterDoubleBlocks = locNumAmps / (2*outerBlockSize);
    long long int firstBlockInd;
    
    // alternate norming* and zeroing the outer blocks (with order based on the desired outcome)
    // These loops aren't parallelised, since they could have 1 or 2 iterations and will prevent
    // inner parallelisation
    if (outerBit == outcome) {

        for (long long int outerDubBlockInd = 0; outerDubBlockInd < numOuterDoubleBlocks; outerDubBlockInd++) {
            firstBlockInd = outerDubBlockInd*2*outerBlockSize;
            
            // *norm only the desired inner blocks in the desired outer block
            alternateNormZeroingSomeAmpBlocks(
                qureg, totalStateProb, innerBit==outcome, 
                firstBlockInd, outerBlockSize, innerBlockSize);
            
            // zero the undesired outer block
            zeroSomeAmps(qureg, firstBlockInd + outerBlockSize, outerBlockSize);
        }

    } else {

        for (long long int outerDubBlockInd = 0; outerDubBlockInd < numOuterDoubleBlocks; outerDubBlockInd++) {
            firstBlockInd = outerDubBlockInd*2*outerBlockSize;
            
            // same thing but undesired outer blocks come first
            zeroSomeAmps(qureg, firstBlockInd, outerBlockSize);
            alternateNormZeroingSomeAmpBlocks(
                qureg, totalStateProb, innerBit==outcome, 
                firstBlockInd + outerBlockSize, outerBlockSize, innerBlockSize);
        }
    }
    
}

qreal densmatr_calcPurityLocal(Qureg qureg) {
    
    /* sum of qureg^2, which is sum_i |qureg[i]|^2 */
    long long int index;
    long long int numAmps = qureg.numAmpsPerChunk;
        
    qreal trace = 0;
    qreal *vecRe = qureg.stateVec.real;
    qreal *vecIm = qureg.stateVec.imag;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (vecRe, vecIm, numAmps) \
    private   (index) \
    reduction ( +:trace )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0LL; index<numAmps; index++) {
                        
            trace += vecRe[index]*vecRe[index] + vecIm[index]*vecIm[index];
        }
    }
    
    return trace;
}

void densmatr_mixDensityMatrix(Qureg combineQureg, qreal otherProb, Qureg otherQureg) {
    
    /* corresponding amplitudes live on the same node (same dimensions) */
    
    // unpack vars for OpenMP
    qreal* combineVecRe = combineQureg.stateVec.real;
    qreal* combineVecIm = combineQureg.stateVec.imag;
    qreal* otherVecRe = otherQureg.stateVec.real;
    qreal* otherVecIm = otherQureg.stateVec.imag;
    long long int numAmps = combineQureg.numAmpsPerChunk;
    long long int index;
    
# ifdef _OPENMP
# pragma omp parallel \
    default (none) \
    shared  (combineVecRe,combineVecIm,otherVecRe,otherVecIm, otherProb, numAmps) \
    private (index)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0; index < numAmps; index++) {
            combineVecRe[index] *= 1-otherProb;
            combineVecIm[index] *= 1-otherProb;
            
            combineVecRe[index] += otherProb * otherVecRe[index];
            combineVecIm[index] += otherProb * otherVecIm[index];
        }
    }
}

/** computes Tr((a-b) conjTrans(a-b)) = sum of abs values of (a-b) */
qreal densmatr_calcHilbertSchmidtDistanceSquaredLocal(Qureg a, Qureg b) {
    
    long long int index;
    long long int numAmps = a.numAmpsPerChunk;
        
    qreal *aRe = a.stateVec.real;
    qreal *aIm = a.stateVec.imag;
    qreal *bRe = b.stateVec.real;
    qreal *bIm = b.stateVec.imag;
    
    qreal trace = 0;
    qreal difRe, difIm;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (aRe,aIm, bRe,bIm, numAmps) \
    private   (index,difRe,difIm) \
    reduction ( +:trace )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0LL; index<numAmps; index++) {
                        
            difRe = aRe[index] - bRe[index];
            difIm = aIm[index] - bIm[index];
            trace += difRe*difRe + difIm*difIm;
        }
    }
    
    return trace;
}

/** computes Tr(conjTrans(a) b) = sum of (a_ij^* b_ij) */
qreal densmatr_calcInnerProductLocal(Qureg a, Qureg b) {
    
    long long int index;
    long long int numAmps = a.numAmpsPerChunk;
        
    qreal *aRe = a.stateVec.real;
    qreal *aIm = a.stateVec.imag;
    qreal *bRe = b.stateVec.real;
    qreal *bIm = b.stateVec.imag;
    
    qreal trace = 0;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (aRe,aIm, bRe,bIm, numAmps) \
    private   (index) \
    reduction ( +:trace )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0LL; index<numAmps; index++) {
            trace += aRe[index]*bRe[index] + aIm[index]*bIm[index];
        }
    }
    
    return trace;
}


/** computes a few dens-columns-worth of (vec^*T) dens * vec */
qreal densmatr_calcFidelityLocal(Qureg qureg, Qureg pureState) {
        
    /* Here, elements of pureState are not accessed (instead grabbed from qureg.pair).
     * We only consult the attributes.
     *
     * qureg is a density matrix, and pureState is a statevector.
     * Every node contains as many columns of qureg as amps by pureState.
     * (each node contains an integer, exponent-of-2 number of whole columns of qureg)
     * Ergo, this node contains columns:
     * qureg.chunkID * pureState.numAmpsPerChunk  to
     * (qureg.chunkID + 1) * pureState.numAmpsPerChunk
     *
     * The first pureState.numAmpsTotal elements of qureg.pairStateVec are the
     * entire pureState state-vector
     */
    
    // unpack everything for OPENMP
    qreal* vecRe  = qureg.pairStateVec.real;
    qreal* vecIm  = qureg.pairStateVec.imag;
    qreal* densRe = qureg.stateVec.real;
    qreal* densIm = qureg.stateVec.imag;
    
    int row, col;
    int dim = (int) pureState.numAmpsTotal; 
    int colsPerNode = (int) pureState.numAmpsPerChunk;
    // using only int, because density matrix has squared as many amps so its 
    // iteration would be impossible if the pureStates numAmpsTotal didn't fit into int
    
    // starting GLOBAL column index of the qureg columns on this node
    int startCol = (int) (qureg.chunkId * pureState.numAmpsPerChunk);
    
    qreal densElemRe, densElemIm;
    qreal prefacRe, prefacIm;
    qreal rowSumRe, rowSumIm;
    qreal vecElemRe, vecElemIm;
    
    // quantity computed by this node
    qreal globalSumRe = 0;   // imag-component is assumed zero
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (vecRe,vecIm,densRe,densIm, dim,colsPerNode,startCol) \
    private   (row,col, prefacRe,prefacIm, rowSumRe,rowSumIm, densElemRe,densElemIm, vecElemRe,vecElemIm) \
    reduction ( +:globalSumRe )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        // indices of my GLOBAL row
        for (row=0; row < dim; row++) {
            
            // single element of conj(pureState)
            prefacRe =   vecRe[row];
            prefacIm = - vecIm[row];
                    
            rowSumRe = 0;
            rowSumIm = 0;
            
            // indices of my LOCAL column
            for (col=0; col < colsPerNode; col++) {
            
                // my local density element
                densElemRe = densRe[row + dim*col];
                densElemIm = densIm[row + dim*col];
            
                // state-vector element
                vecElemRe = vecRe[startCol + col];
                vecElemIm = vecIm[startCol + col];
            
                rowSumRe += densElemRe*vecElemRe - densElemIm*vecElemIm;
                rowSumIm += densElemRe*vecElemIm + densElemIm*vecElemRe;
            }
        
            globalSumRe += rowSumRe*prefacRe - rowSumIm*prefacIm;   
        }
    }
    
    return globalSumRe;
}

Complex statevec_calcInnerProductLocal(Qureg bra, Qureg ket) {
    
    qreal innerProdReal = 0;
    qreal innerProdImag = 0;
    
    long long int index;
    long long int numAmps = bra.numAmpsPerChunk;
    qreal *braVecReal = bra.stateVec.real;
    qreal *braVecImag = bra.stateVec.imag;
    qreal *ketVecReal = ket.stateVec.real;
    qreal *ketVecImag = ket.stateVec.imag;
    
    qreal braRe, braIm, ketRe, ketIm;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (braVecReal, braVecImag, ketVecReal, ketVecImag, numAmps) \
    private   (index, braRe, braIm, ketRe, ketIm) \
    reduction ( +:innerProdReal, innerProdImag )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0; index < numAmps; index++) {
            braRe = braVecReal[index];
            braIm = braVecImag[index];
            ketRe = ketVecReal[index];
            ketIm = ketVecImag[index];
            
            // conj(bra_i) * ket_i
            innerProdReal += braRe*ketRe + braIm*ketIm;
            innerProdImag += braRe*ketIm - braIm*ketRe;
        }
    }
    
    Complex innerProd;
    innerProd.real = innerProdReal;
    innerProd.imag = innerProdImag;
    return innerProd;
}



void densmatr_initClassicalState (Qureg qureg, long long int stateInd)
{
    // dimension of the state vector
    long long int densityNumElems = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    qreal *densityReal = qureg.stateVec.real;
    qreal *densityImag = qureg.stateVec.imag;

    // initialise the state to all zeros
    long long int index;
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (densityNumElems, densityReal, densityImag) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<densityNumElems; index++) {
            densityReal[index] = 0.0;
            densityImag[index] = 0.0;
        }
    }
    
    // index of the single density matrix elem to set non-zero
    long long int densityDim = 1LL << qureg.numQubitsRepresented;
    long long int densityInd = (densityDim + 1)*stateInd;

    // give the specified classical state prob 1
    if (qureg.chunkId == densityInd / densityNumElems){
        densityReal[densityInd % densityNumElems] = 1.0;
        densityImag[densityInd % densityNumElems] = 0.0;
    }
}


void densmatr_initPlusState (Qureg qureg)
{
    // |+><+| = sum_i 1/sqrt(2^N) |i> 1/sqrt(2^N) <j| = sum_ij 1/2^N |i><j|
    long long int dim = (1LL << qureg.numQubitsRepresented);
    qreal probFactor = 1.0/((qreal) dim);

    // Can't use qureg->stateVec as a private OMP var
    qreal *densityReal = qureg.stateVec.real;
    qreal *densityImag = qureg.stateVec.imag;

    long long int index;
    long long int chunkSize = qureg.numAmpsPerChunk;
    // initialise the state to |+++..+++> = 1/normFactor {1, 1, 1, ...}
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkSize, densityReal, densityImag, probFactor) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<chunkSize; index++) {
            densityReal[index] = probFactor;
            densityImag[index] = 0.0;
        }
    }
}

void densmatr_initPureStateLocal(Qureg targetQureg, Qureg copyQureg) {
    
    /* copyQureg amps aren't explicitly used - they're accessed through targetQureg.pair,
     * which contains the full pure statevector.
     * targetQureg has as many columns on node as copyQureg has amps
     */
    
    long long int colOffset = targetQureg.chunkId * copyQureg.numAmpsPerChunk;
    long long int colsPerNode = copyQureg.numAmpsPerChunk;
    long long int rowsPerNode = copyQureg.numAmpsTotal;
    
    // unpack vars for OpenMP
    qreal* vecRe = targetQureg.pairStateVec.real;
    qreal* vecIm = targetQureg.pairStateVec.imag;
    qreal* densRe = targetQureg.stateVec.real;
    qreal* densIm = targetQureg.stateVec.imag;
    
    long long int col, row, index;
    
    // a_i conj(a_j) |i><j|
    qreal ketRe, ketIm, braRe, braIm;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (colOffset, colsPerNode,rowsPerNode, vecRe,vecIm,densRe,densIm) \
    private  (col,row, ketRe,ketIm,braRe,braIm, index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // local column
        for (col=0; col < colsPerNode; col++) {
        
            // global row
            for (row=0; row < rowsPerNode; row++) {
            
                // get pure state amps
                ketRe = vecRe[row];
                ketIm = vecIm[row];
                braRe =   vecRe[col + colOffset];
                braIm = - vecIm[col + colOffset]; // minus for conjugation
            
                // update density matrix
                index = row + col*rowsPerNode; // local ind
                densRe[index] = ketRe*braRe - ketIm*braIm;
                densIm[index] = ketRe*braIm + ketIm*braRe;
            }
        }
    }
}

void statevec_setAmps(Qureg qureg, long long int startInd, qreal* reals, qreal* imags, long long int numAmps) {
    
    /* this is actually distributed, since the user's code runs on every node */
    
    // local start/end indices of the given amplitudes, assuming they fit in this chunk
    // these may be negative or above qureg.numAmpsPerChunk
    long long int localStartInd = startInd - qureg.chunkId*qureg.numAmpsPerChunk;
    long long int localEndInd = localStartInd + numAmps; // exclusive
    
    // add this to a local index to get corresponding elem in reals & imags
    long long int offset = qureg.chunkId*qureg.numAmpsPerChunk - startInd;
    
    // restrict these indices to fit into this chunk
    if (localStartInd < 0)
        localStartInd = 0;
    if (localEndInd > qureg.numAmpsPerChunk)
        localEndInd = qureg.numAmpsPerChunk;
    // they may now be out of order = no iterations
    
    // unpacking OpenMP vars
    long long int index;
    qreal* vecRe = qureg.stateVec.real;
    qreal* vecIm = qureg.stateVec.imag;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (localStartInd,localEndInd, vecRe,vecIm, reals,imags, offset) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // iterate these local inds - this might involve no iterations
        for (index=localStartInd; index < localEndInd; index++) {
            vecRe[index] = reals[index + offset];
            vecIm[index] = imags[index + offset];
        }
    }
}

void statevec_createQureg(Qureg *qureg, int numQubits, QuESTEnv env)
{
    long long int numAmps = 1LL << numQubits;
    long long int numAmpsPerRank = numAmps/env.numRanks;

    validateMemoryAllocationSize(numAmpsPerRank, __func__);

    size_t arrSize = (size_t) (numAmpsPerRank * sizeof(*(qureg->stateVec.real)));
    qureg->stateVec.real = malloc(arrSize);
    qureg->stateVec.imag = malloc(arrSize);
    if (env.numRanks>1){
        qureg->pairStateVec.real = malloc(arrSize);
        qureg->pairStateVec.imag = malloc(arrSize);
    }

    qureg->numQubitsInStateVec = numQubits;
    qureg->numAmpsTotal = numAmps;
    qureg->numAmpsPerChunk = numAmpsPerRank;
    qureg->chunkId = env.rank;
    qureg->numChunks = env.numRanks;
    qureg->isDensityMatrix = 0;

    validateQuregAllocation(qureg, env, __func__);
}

void statevec_destroyQureg(Qureg qureg, QuESTEnv env){

    qureg.numQubitsInStateVec = 0;
    qureg.numAmpsTotal = 0;
    qureg.numAmpsPerChunk = 0;

    free(qureg.stateVec.real);
    free(qureg.stateVec.imag);
    if (env.numRanks>1){
        free(qureg.pairStateVec.real);
        free(qureg.pairStateVec.imag);
    }
    qureg.stateVec.real = NULL;
    qureg.stateVec.imag = NULL;
    qureg.pairStateVec.real = NULL;
    qureg.pairStateVec.imag = NULL;
}

void statevec_applySubDiagonalOp(Qureg qureg, int* targets, SubDiagonalOp op, int conj) {
    
    // each node/chunk modifies only its values in an embarrassingly parallelisable way
    long long int numLocalAmps = qureg.numAmpsPerChunk;
    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;
    qreal* opRe = op.real;
    qreal* opIm = op.imag;
    
    long long int indPref = qureg.chunkId * numLocalAmps;
    int numTargets = op.numQubits;
    
    int conjFac = 1;
    if (conj)
        conjFac = -1;
    
    int t;
    long long int i, j, v;
    qreal elemRe, elemIm, ampRe, ampIm;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (numLocalAmps, indPref, numTargets, targets, stateRe,stateIm) \
    private   (j,i,v,t, elemRe,elemIm, ampRe,ampIm)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (j=0; j<numLocalAmps; j++) {
            
            i = indPref | j;
            v = 0;
            for (t=0; t<numTargets; t++)
                v |= extractBit(targets[t], i) << t;
                
            elemRe = opRe[v];
            elemIm = opIm[v] * conjFac;
            
            ampRe = stateRe[j];
            ampIm = stateIm[j];
            
            // (a + b i)(c + d i) = (a c - b d) + i (a d + b c)
            stateRe[j] = ampRe*elemRe - ampIm*elemIm;
            stateIm[j] = ampRe*elemIm + ampIm*elemRe;
        }
    }
}

DiagonalOp agnostic_createDiagonalOp(int numQubits, QuESTEnv env) {

    // the 2^numQubits values will be evenly split between the env.numRanks nodes
    DiagonalOp op;
    op.numQubits = numQubits;
    op.numElemsPerChunk = (1LL << numQubits) / env.numRanks;
    op.chunkId = env.rank;
    op.numChunks = env.numRanks;

    // allocate CPU memory (initialised to zero)
    op.real = (qreal*) calloc(op.numElemsPerChunk, sizeof(qreal));
    op.imag = (qreal*) calloc(op.numElemsPerChunk, sizeof(qreal));

    // check cpu memory allocation was successful
    validateDiagonalOpAllocation(&op, env, __func__);

    return op;
}

void agnostic_destroyDiagonalOp(DiagonalOp op) {
    free(op.real);
    free(op.imag);
}

void agnostic_syncDiagonalOp(DiagonalOp op) {
    // nothing to do on CPU
}

void agnostic_initDiagonalOpFromPauliHamil(DiagonalOp op, PauliHamil hamil) {
    
    /* each node modifies its op sub-partition, evaluating the full hamil 
     * for every element in the sub-partition 
     */
    
    // unpack op
    long long int offset = op.chunkId * op.numElemsPerChunk;
    long long int numElems = op.numElemsPerChunk;
    qreal* opRe = op.real;
    qreal* opIm = op.imag;
    
    // unpack hamil
    int numTerms = hamil.numSumTerms;
    int numQubits = hamil.numQubits;
    qreal* coeffs = hamil.termCoeffs;
    enum pauliOpType* codes = hamil.pauliCodes;
    
    // private OpenMP vars
    long long int i, globalInd;
    qreal elem;
    int t, q, isOddNumOnes, sign;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (offset,numElems, opRe,opIm, numTerms,numQubits,coeffs,codes) \
    private  (i,globalInd, elem, isOddNumOnes,t,q,sign) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (i=0; i<numElems; i++) {
            
            globalInd = i + offset;
            elem = 0;
            
            // add every Hamiltonian coefficient to this element, either + or -
            for (t=0; t<numTerms; t++) {
                
                // determine parity of ones (in globalInd basis state) of the current term's targets
                isOddNumOnes = 0;
                for (q=0; q<numQubits; q++)
                    if (codes[q + t*numQubits] == PAULI_Z)
                        if (extractBit(q, globalInd))
                            isOddNumOnes = !isOddNumOnes;
                
                // add +- term coeff (avoiding thread divergence)
                sign = 1 - 2*isOddNumOnes; // (-1 if isOddNumOnes, else +1)
                elem += coeffs[t] * sign;
            }
            
            opRe[i] = elem;
            opIm[i] = 0;
        }
    }
    
    // we don't synch to GPU, because in GPU mode, the GPU populates and synchs to RAM
    // agnostic_syncDiagonalOp(op);
}

void statevec_reportStateToScreen(Qureg qureg, QuESTEnv env, int reportRank){
    long long int index;
    int rank;
    if (qureg.numQubitsInStateVec<=5){
        for (rank=0; rank<qureg.numChunks; rank++){
            if (qureg.chunkId==rank){
                if (reportRank) {
                    printf("Reporting state from rank %d [\n", qureg.chunkId);
                    printf("real, imag\n");
                } else if (rank==0) {
                    printf("Reporting state [\n");
                    printf("real, imag\n");
                }

                for(index=0; index<qureg.numAmpsPerChunk; index++){
                    //printf(REAL_STRING_FORMAT ", " REAL_STRING_FORMAT "\n", qureg.pairStateVec.real[index], qureg.pairStateVec.imag[index]);
                    printf(REAL_STRING_FORMAT ", " REAL_STRING_FORMAT "\n", qureg.stateVec.real[index], qureg.stateVec.imag[index]);
                }
                if (reportRank || rank==qureg.numChunks-1) printf("]\n");
            }
            syncQuESTEnv(env);
        }
    } else printf("Error: reportStateToScreen will not print output for systems of more than 5 qubits.\n");
}

void statevec_initBlankState (Qureg qureg)
{
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

    // initialise the state-vector to all-zeroes
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
}

void statevec_initZeroState (Qureg qureg)
{
    statevec_initBlankState(qureg);
    if (qureg.chunkId==0){
        // zero state |0000..0000> has probability 1
        qureg.stateVec.real[0] = 1.0;
        qureg.stateVec.imag[0] = 0.0;
    }
}

void statevec_initPlusState (Qureg qureg)
{
    long long int chunkSize, stateVecSize;
    long long int index;

    // dimension of the state vector
    chunkSize = qureg.numAmpsPerChunk;
    stateVecSize = chunkSize*qureg.numChunks;
    qreal normFactor = 1.0/sqrt((qreal)stateVecSize);

    // Can't use qureg->stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

    // initialise the state to |+++..+++> = 1/normFactor {1, 1, 1, ...}
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
}

void statevec_initClassicalState (Qureg qureg, long long int stateInd)
{
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

    // initialise the state to vector to all zeros
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

    // give the specified classical state prob 1
    if (qureg.chunkId == stateInd/stateVecSize){
        stateVecReal[stateInd % stateVecSize] = 1.0;
        stateVecImag[stateInd % stateVecSize] = 0.0;
    }
}

void statevec_cloneQureg(Qureg targetQureg, Qureg copyQureg) {
    
    // registers are equal sized, so nodes hold the same state-vector partitions
    long long int stateVecSize;
    long long int index;

    // dimension of the state vector
    stateVecSize = targetQureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    qreal *targetStateVecReal = targetQureg.stateVec.real;
    qreal *targetStateVecImag = targetQureg.stateVec.imag;
    qreal *copyStateVecReal = copyQureg.stateVec.real;
    qreal *copyStateVecImag = copyQureg.stateVec.imag;

    // initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecSize, targetStateVecReal, targetStateVecImag, copyStateVecReal, copyStateVecImag) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            targetStateVecReal[index] = copyStateVecReal[index];
            targetStateVecImag[index] = copyStateVecImag[index];
        }
    }
}

/**
 * Initialise the state vector of probability amplitudes to an (unphysical) state with
 * each component of each probability amplitude a unique floating point value. For debugging processes
 * @param[in,out] qureg object representing the set of qubits to be initialised
 */
void statevec_initDebugState (Qureg qureg)
{
    long long int chunkSize;
    long long int index;
    long long int indexOffset;

    // dimension of the state vector
    chunkSize = qureg.numAmpsPerChunk;

    // Can't use qureg->stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

    indexOffset = chunkSize * qureg.chunkId;

    // initialise the state to |0000..0000>
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkSize, stateVecReal, stateVecImag, indexOffset) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<chunkSize; index++) {
            stateVecReal[index] = ((indexOffset + index)*2.0)/10.0;
            stateVecImag[index] = ((indexOffset + index)*2.0+1.0)/10.0;
        }
    }
}

void statevec_compactUnitaryLocal (Qureg qureg, int targetQubit, Complex alpha, Complex beta)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;
    qreal alphaImag=alpha.imag, alphaReal=alpha.real;
    qreal betaImag=beta.imag, betaReal=beta.real;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag, numTasks) \
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

void statevec_multiControlledTwoQubitUnitaryLocal(Qureg qureg, long long int ctrlMask, int q1, int q2, ComplexMatrix4 u) {

    // can't use qureg.stateVec as a private OMP var
    qreal *reVec = qureg.stateVec.real;
    qreal *imVec = qureg.stateVec.imag;
    
    // the global (between all nodes) index of this node's start index
    long long int globalIndStart = qureg.chunkId*qureg.numAmpsPerChunk; 
    
    long long int numTasks = qureg.numAmpsPerChunk >> 2; // each iteration updates 4 amplitudes
    long long int thisTask;
    long long int thisGlobalInd00;
    long long int ind00, ind01, ind10, ind11;
    qreal re00, re01, re10, re11;
    qreal im00, im01, im10, im11;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (reVec,imVec,globalIndStart,numTasks,ctrlMask,u,q2,q1) \
    private  (thisTask, thisGlobalInd00, ind00,ind01,ind10,ind11, re00,re01,re10,re11, im00,im01,im10,im11)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            
            // determine ind00 of |..0..0..>
            ind00 = insertTwoZeroBits(thisTask, q1, q2);
            
            // skip amplitude if controls aren't in 1 state (overloaded for speed)
            thisGlobalInd00 = ind00 + globalIndStart;
            if (ctrlMask && ((ctrlMask & thisGlobalInd00) != ctrlMask))
                continue;
            
            // inds of |..0..1..>, |..1..0..> and |..1..1..>
            ind01 = flipBit(ind00, q1);
            ind10 = flipBit(ind00, q2);
            ind11 = flipBit(ind01, q2);

            // extract statevec amplitudes 
            re00 = reVec[ind00]; im00 = imVec[ind00];
            re01 = reVec[ind01]; im01 = imVec[ind01];
            re10 = reVec[ind10]; im10 = imVec[ind10];
            re11 = reVec[ind11]; im11 = imVec[ind11];

            // apply u * {amp00, amp01, amp10, amp11}
            reVec[ind00] = 
                u.real[0][0]*re00 - u.imag[0][0]*im00 +
                u.real[0][1]*re01 - u.imag[0][1]*im01 +
                u.real[0][2]*re10 - u.imag[0][2]*im10 +
                u.real[0][3]*re11 - u.imag[0][3]*im11;
            imVec[ind00] =
                u.imag[0][0]*re00 + u.real[0][0]*im00 +
                u.imag[0][1]*re01 + u.real[0][1]*im01 +
                u.imag[0][2]*re10 + u.real[0][2]*im10 +
                u.imag[0][3]*re11 + u.real[0][3]*im11;
                
            reVec[ind01] = 
                u.real[1][0]*re00 - u.imag[1][0]*im00 +
                u.real[1][1]*re01 - u.imag[1][1]*im01 +
                u.real[1][2]*re10 - u.imag[1][2]*im10 +
                u.real[1][3]*re11 - u.imag[1][3]*im11;
            imVec[ind01] =
                u.imag[1][0]*re00 + u.real[1][0]*im00 +
                u.imag[1][1]*re01 + u.real[1][1]*im01 +
                u.imag[1][2]*re10 + u.real[1][2]*im10 +
                u.imag[1][3]*re11 + u.real[1][3]*im11;
                
            reVec[ind10] = 
                u.real[2][0]*re00 - u.imag[2][0]*im00 +
                u.real[2][1]*re01 - u.imag[2][1]*im01 +
                u.real[2][2]*re10 - u.imag[2][2]*im10 +
                u.real[2][3]*re11 - u.imag[2][3]*im11;
            imVec[ind10] =
                u.imag[2][0]*re00 + u.real[2][0]*im00 +
                u.imag[2][1]*re01 + u.real[2][1]*im01 +
                u.imag[2][2]*re10 + u.real[2][2]*im10 +
                u.imag[2][3]*re11 + u.real[2][3]*im11;    
                
            reVec[ind11] = 
                u.real[3][0]*re00 - u.imag[3][0]*im00 +
                u.real[3][1]*re01 - u.imag[3][1]*im01 +
                u.real[3][2]*re10 - u.imag[3][2]*im10 +
                u.real[3][3]*re11 - u.imag[3][3]*im11;
            imVec[ind11] =
                u.imag[3][0]*re00 + u.real[3][0]*im00 +
                u.imag[3][1]*re01 + u.real[3][1]*im01 +
                u.imag[3][2]*re10 + u.real[3][2]*im10 +
                u.imag[3][3]*re11 + u.real[3][3]*im11;    
        }
    }
}

int qsortComp(const void *a, const void *b) {
    return *(int*)a - *(int*)b; 
}

void statevec_multiControlledMultiQubitUnitaryLocal(Qureg qureg, long long int ctrlMask, int* targs, int numTargs, ComplexMatrixN u)
{
    // can't use qureg.stateVec as a private OMP var
    qreal *reVec = qureg.stateVec.real;
    qreal *imVec = qureg.stateVec.imag;
    
    long long int numTasks = qureg.numAmpsPerChunk >> numTargs;  // kernel called on every 1 in 2^numTargs amplitudes
    long long int numTargAmps = 1 << u.numQubits;  // num amps to be modified by each task
    
    // the global (between all nodes) index of this node's start index
    long long int globalIndStart = qureg.chunkId*qureg.numAmpsPerChunk; 
    
    long long int thisTask;
    long long int thisInd00; // this thread's index of |..0..0..> (target qubits = 0) 
    long long int thisGlobalInd00; // the global (between all nodes) index of this thread's |..0..0..> state
    long long int ind;   // each thread's iteration of amplitudes to modify
    int i, t, r, c;  // each thread's iteration of amps and targets 
    qreal reElem, imElem;  // each thread's iteration of u elements
    
    // each thread/task will record and modify numTargAmps amplitudes, privately
    // (of course, tasks eliminated by the ctrlMask won't edit their allocation)
    //
    // If we're NOT on windows, we can fortunately use the stack directly
    #ifndef _WIN32
        long long int ampInds[numTargAmps];
        qreal reAmps[numTargAmps];
        qreal imAmps[numTargAmps];

        int sortedTargs[numTargs];
    // on Windows, with no VLA, we can use _malloca to allocate on stack (must free)
    #else
        long long int* ampInds;
        qreal* reAmps;
        qreal* imAmps;
        int* sortedTargs = (int*) _malloca(numTargs * sizeof *sortedTargs);
    #endif

    // we need a sorted targets list to find thisInd00 for each task.
    // we can't modify targets, because the user-ordering of targets matters in u
    for (int t=0; t < numTargs; t++) 
        sortedTargs[t] = targs[t];
    qsort(sortedTargs, numTargs, sizeof(int), qsortComp);
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (reVec,imVec, numTasks,numTargAmps,globalIndStart, ctrlMask,targs,sortedTargs,u,numTargs) \
    private  (thisTask,thisInd00,thisGlobalInd00,ind,i,t,r,c,reElem,imElem,  ampInds,reAmps,imAmps)
# endif
    {
        // when manually allocating array memory (on Windows), this must be done in each thread
        // separately and is not performed automatically by declaring a var as omp-private
        # ifdef _WIN32
            ampInds = (long long int*) _malloca(numTargAmps * sizeof *ampInds);
            reAmps = (qreal*) _malloca(numTargAmps * sizeof *reAmps);
            imAmps = (qreal*) _malloca(numTargAmps * sizeof *imAmps);
        # endif
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            
            // find this task's start index (where all targs are 0)
            thisInd00 = thisTask;
            for (t=0; t < numTargs; t++)
                thisInd00 = insertZeroBit(thisInd00, sortedTargs[t]);
                
            // this task only modifies amplitudes if control qubits are 1 for this state
            thisGlobalInd00 = thisInd00 + globalIndStart;
            if (ctrlMask && ((ctrlMask & thisGlobalInd00) != ctrlMask))
                continue;
                
            // determine the indices and record values of this tasks's target amps
            for (i=0; i < numTargAmps; i++) {
                
                // get statevec index of current target qubit assignment
                ind = thisInd00;
                for (t=0; t < numTargs; t++)
                    if (extractBit(t, i))
                        ind = flipBit(ind, targs[t]);
                
                // update this tasks's private arrays
                ampInds[i] = ind;
                reAmps [i] = reVec[ind];
                imAmps [i] = imVec[ind];
            }
            
            // modify this tasks's target amplitudes
            for (r=0; r < numTargAmps; r++) {
                ind = ampInds[r];
                reVec[ind] = 0;
                imVec[ind] = 0;
                
                for (c=0; c < numTargAmps; c++) {
                    reElem = u.real[r][c];
                    imElem = u.imag[r][c];
                    reVec[ind] += reAmps[c]*reElem - imAmps[c]*imElem;
                    imVec[ind] += reAmps[c]*imElem + imAmps[c]*reElem;
                }
            }
        }
        // on Windows, we must explicitly free the stack structures
        #ifdef _WIN32
            _freea(ampInds);
            _freea(reAmps);
            _freea(imAmps);
        #endif
    }

    #ifdef _WIN32
        _freea(sortedTargs);
    #endif
}

void statevec_unitaryLocal(Qureg qureg, int targetQubit, ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, u,numTasks) \
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


            // state[indexUp] = u00 * state[indexUp] + u01 * state[indexLo]
            stateVecReal[indexUp] = u.real[0][0]*stateRealUp - u.imag[0][0]*stateImagUp 
                + u.real[0][1]*stateRealLo - u.imag[0][1]*stateImagLo;
            stateVecImag[indexUp] = u.real[0][0]*stateImagUp + u.imag[0][0]*stateRealUp 
                + u.real[0][1]*stateImagLo + u.imag[0][1]*stateRealLo;

            // state[indexLo] = u10  * state[indexUp] + u11 * state[indexLo]
            stateVecReal[indexLo] = u.real[1][0]*stateRealUp  - u.imag[1][0]*stateImagUp 
                + u.real[1][1]*stateRealLo  -  u.imag[1][1]*stateImagLo;
            stateVecImag[indexLo] = u.real[1][0]*stateImagUp + u.imag[1][0]*stateRealUp 
                + u.real[1][1]*stateImagLo + u.imag[1][1]*stateRealLo;

        } 
    }
} 

/** Rotate a single qubit in the state vector of probability amplitudes, 
 * given two complex numbers alpha and beta, 
 * and a subset of the state vector with upper and lower block values stored seperately.
 *                                                                       
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_compactUnitaryDistributed (Qureg qureg,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    qreal   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;

    qreal rot1Real=rot1.real, rot1Imag=rot1.imag;
    qreal rot2Real=rot2.real, rot2Imag=rot2.imag;
    qreal *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    qreal *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag,numTasks) \
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

/** Apply a unitary operation to a single qubit
 *  given a subset of the state vector with upper and lower block values 
 * stored seperately.
 *
 *  @remarks Qubits are zero-based and the first qubit is the rightmost                  
 *                                                                        
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] rot1 
 *  @param[in] rot2
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_unitaryDistributed (Qureg qureg,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    qreal   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;

    qreal rot1Real=rot1.real, rot1Imag=rot1.imag;
    qreal rot2Real=rot2.real, rot2Imag=rot2.imag;
    qreal *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    qreal *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;


# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real, rot1Imag, rot2Real, rot2Imag,numTasks) \
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

            stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp 
                + rot2Real*stateRealLo - rot2Imag*stateImagLo;
            stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp 
                + rot2Real*stateImagLo + rot2Imag*stateRealLo;
        }
    }
}

void statevec_controlledCompactUnitaryLocal (Qureg qureg, int controlQubit, int targetQubit, 
        Complex alpha, Complex beta)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;
    qreal alphaImag=alpha.imag, alphaReal=alpha.real;
    qreal betaImag=beta.imag, betaReal=beta.real;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, alphaReal,alphaImag, betaReal,betaImag, \
                numTasks,chunkId,chunkSize,controlQubit) \
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
            if (controlBit){
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

/* ctrlQubitsMask is a bit mask indicating which qubits are control Qubits
 * ctrlFlipMask is a bit mask indicating which control qubits should be 'flipped'
 * in the condition, i.e. they should have value 0 when the unitary is applied
 */
void statevec_multiControlledUnitaryLocal(
    Qureg qureg, int targetQubit, 
    long long int ctrlQubitsMask, long long int ctrlFlipMask,
    ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, u, ctrlQubitsMask,ctrlFlipMask, \
                numTasks,chunkId,chunkSize) \
    private  (thisTask,thisBlock, indexUp,indexLo, stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {

            thisBlock   = thisTask / sizeHalfBlock;
            indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
            indexLo     = indexUp + sizeHalfBlock;
            
            
            // take the basis index, flip the designated (XOR) 'control' bits, AND with the controls.
            // if this equals the control mask, the control qubits have the desired values in the basis index
            if (ctrlQubitsMask == (ctrlQubitsMask & ((indexUp+chunkId*chunkSize) ^ ctrlFlipMask))) {
                // store current state vector values in temp variables
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

                stateRealLo = stateVecReal[indexLo];
                stateImagLo = stateVecImag[indexLo];

                // state[indexUp] = u00 * state[indexUp] + u01 * state[indexLo]
                stateVecReal[indexUp] = u.real[0][0]*stateRealUp - u.imag[0][0]*stateImagUp 
                    + u.real[0][1]*stateRealLo - u.imag[0][1]*stateImagLo;
                stateVecImag[indexUp] = u.real[0][0]*stateImagUp + u.imag[0][0]*stateRealUp 
                    + u.real[0][1]*stateImagLo + u.imag[0][1]*stateRealLo;

                // state[indexLo] = u10  * state[indexUp] + u11 * state[indexLo]
                stateVecReal[indexLo] = u.real[1][0]*stateRealUp  - u.imag[1][0]*stateImagUp 
                    + u.real[1][1]*stateRealLo  -  u.imag[1][1]*stateImagLo;
                stateVecImag[indexLo] = u.real[1][0]*stateImagUp + u.imag[1][0]*stateRealUp 
                    + u.real[1][1]*stateImagLo + u.imag[1][1]*stateRealLo;
            }
        } 
    }

}

void statevec_controlledUnitaryLocal(Qureg qureg, int controlQubit, int targetQubit, 
        ComplexMatrix2 u)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, u,numTasks,chunkId,chunkSize,controlQubit) \
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
            if (controlBit){
                // store current state vector values in temp variables
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

                stateRealLo = stateVecReal[indexLo];
                stateImagLo = stateVecImag[indexLo];


                // state[indexUp] = u00 * state[indexUp] + u01 * state[indexLo]
                stateVecReal[indexUp] = u.real[0][0]*stateRealUp - u.imag[0][0]*stateImagUp 
                    + u.real[0][1]*stateRealLo - u.imag[0][1]*stateImagLo;
                stateVecImag[indexUp] = u.real[0][0]*stateImagUp + u.imag[0][0]*stateRealUp 
                    + u.real[0][1]*stateImagLo + u.imag[0][1]*stateRealLo;

                // state[indexLo] = u10  * state[indexUp] + u11 * state[indexLo]
                stateVecReal[indexLo] = u.real[1][0]*stateRealUp  - u.imag[1][0]*stateImagUp 
                    + u.real[1][1]*stateRealLo  -  u.imag[1][1]*stateImagLo;
                stateVecImag[indexLo] = u.real[1][0]*stateImagUp + u.imag[1][0]*stateRealUp 
                    + u.real[1][1]*stateImagLo + u.imag[1][1]*stateRealLo;
            }
        } 
    }

}

/** Rotate a single qubit in the state vector of probability amplitudes, given two complex 
 * numbers alpha and beta and a subset of the state vector with upper and lower block values 
 * stored seperately. Only perform the rotation where the control qubit is one.
 *                                               
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_controlledCompactUnitaryDistributed (Qureg qureg, int controlQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    qreal   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    qreal rot1Real=rot1.real, rot1Imag=rot1.imag;
    qreal rot2Real=rot2.real, rot2Imag=rot2.imag;
    qreal *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    qreal *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag,numTasks,chunkId,chunkSize,controlQubit) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
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

/** Rotate a single qubit in the state vector of probability amplitudes, given two complex 
 *  numbers alpha and beta and a subset of the state vector with upper and lower block values 
 *  stored seperately. Only perform the rotation where the control qubit is one.
 *                                                 
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] controlQubit qubit to determine whether or not to perform a rotation 
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_controlledUnitaryDistributed (Qureg qureg, int controlQubit,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    qreal   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    qreal rot1Real=rot1.real, rot1Imag=rot1.imag;
    qreal rot2Real=rot2.real, rot2Imag=rot2.imag;
    qreal *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    qreal *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag, numTasks,chunkId,chunkSize,controlQubit) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
                // store current state vector values in temp variables
                stateRealUp = stateVecRealUp[thisTask];
                stateImagUp = stateVecImagUp[thisTask];

                stateRealLo = stateVecRealLo[thisTask];
                stateImagLo = stateVecImagLo[thisTask];

                stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp 
                    + rot2Real*stateRealLo - rot2Imag*stateImagLo;
                stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp 
                    + rot2Real*stateImagLo + rot2Imag*stateRealLo;
            }
        }
    }
}

/** Apply a unitary operation to a single qubit in the state vector of probability amplitudes, given
 *  a subset of the state vector with upper and lower block values 
 stored seperately. Only perform the rotation where all the control qubits are 1.
 *                                                 
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] targetQubit qubit to rotate
 *  @param[in] ctrlQubitsMask a bit mask indicating whether each qubit is a control (1) or not (0)
 *  @param[in] ctrlFlipMask a bit mask indicating whether each qubit (only controls are relevant)
 *             should be flipped when checking the control condition
 *  @param[in] rot1 rotation angle
 *  @param[in] rot2 rotation angle
 *  @param[in] stateVecUp probability amplitudes in upper half of a block
 *  @param[in] stateVecLo probability amplitudes in lower half of a block
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_multiControlledUnitaryDistributed (
        Qureg qureg, 
        int targetQubit, 
        long long int ctrlQubitsMask, long long int ctrlFlipMask,
        Complex rot1, Complex rot2,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut)
{

    qreal   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    qreal rot1Real=rot1.real, rot1Imag=rot1.imag;
    qreal rot2Real=rot2.real, rot2Imag=rot2.imag;
    qreal *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    qreal *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            rot1Real,rot1Imag, rot2Real,rot2Imag, ctrlQubitsMask,ctrlFlipMask, numTasks,chunkId,chunkSize) \
    private  (thisTask,stateRealUp,stateImagUp,stateRealLo,stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            if (ctrlQubitsMask == (ctrlQubitsMask & ((thisTask+chunkId*chunkSize) ^ ctrlFlipMask))) {
                // store current state vector values in temp variables
                stateRealUp = stateVecRealUp[thisTask];
                stateImagUp = stateVecImagUp[thisTask];

                stateRealLo = stateVecRealLo[thisTask];
                stateImagLo = stateVecImagLo[thisTask];

                stateVecRealOut[thisTask] = rot1Real*stateRealUp - rot1Imag*stateImagUp 
                    + rot2Real*stateRealLo - rot2Imag*stateImagLo;
                stateVecImagOut[thisTask] = rot1Real*stateImagUp + rot1Imag*stateRealUp 
                    + rot2Real*stateImagLo + rot2Imag*stateRealLo;
            }
        }
    }
}

void statevec_pauliXLocal(Qureg qureg, int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateImagUp;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, numTasks) \
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
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk
 *  
 *  @remarks Qubits are zero-based and the                     
 *  the first qubit is the rightmost                  
 *                                                                        
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_pauliXDistributed (Qureg qureg,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut)
{

    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;

    qreal *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut,numTasks) \
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

void statevec_controlledNotLocal(Qureg qureg, int controlQubit, int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateImagUp;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag,numTasks,chunkId,chunkSize,controlQubit) \
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
            if (controlBit){
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
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk. Only perform the rotation
 *  for elements where controlQubit is one.
 *                                          
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] controlQubit the control qubit
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 */
void statevec_controlledNotDistributed (Qureg qureg, int controlQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut)
{

    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    qreal *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut, \
                numTasks,chunkId,chunkSize,controlQubit) \
    private  (thisTask,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
                stateVecRealOut[thisTask] = stateVecRealIn[thisTask];
                stateVecImagOut[thisTask] = stateVecImagIn[thisTask];
            }
        }
    }
}

void statevec_multiControlledMultiQubitNotLocal(Qureg qureg, int ctrlMask, int targMask) {
    long long int numAmps = qureg.numAmpsPerChunk;
    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;
    
    long long int globalOffset = qureg.chunkId * numAmps;
    
    // each amplitude is swapped with a 'mate' amplitude
    long long int ampInd, mateInd, globalInd;
    qreal mateRe, mateIm;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateRe,stateIm, numAmps, ctrlMask,targMask, globalOffset) \
    private  (ampInd, mateInd,mateRe,mateIm, globalInd)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (ampInd = 0; ampInd < numAmps; ampInd++) {
            
            /* it may be a premature optimisation to remove the seemingly wasteful continues below,
             * because the maximum skipped amplitudes is 1/2 that stored in the node 
             * (e.g. since this function is not called if all amps should be skipped via controls),
             * and since we're memory-bandwidth bottlenecked. 
             */
             
            // although amps are local, we may still be running in distributed mode, 
            // and hence need to consult the global index to determine the values of 
            // the control qubits
            globalInd = ampInd + globalOffset;
            
            // modify amplitude only if control qubits are 1 for this state
            if (ctrlMask && ((ctrlMask & globalInd) != ctrlMask))
                continue;

            mateInd = ampInd ^ targMask;
            
            // if the mate is behind, it was already processed
            if (mateInd < ampInd)
                continue;
            
            mateRe = stateRe[mateInd];
            mateIm = stateIm[mateInd];
            
            // swap amp with mate
            stateRe[mateInd] = stateRe[ampInd];
            stateIm[mateInd] = stateIm[ampInd];
            stateRe[ampInd] = mateRe;
            stateIm[ampInd] = mateIm;
        }
    }
}

void statevec_multiControlledMultiQubitNotDistributed(
    Qureg qureg, int ctrlMask, int targMask,
    ComplexArray stateVecIn,
    ComplexArray stateVecOut
) {
    long long int numAmps = qureg.numAmpsPerChunk;
    long long int globalOffset = qureg.chunkId * numAmps;
    
    /* stateVecOut is qureg's local state-vector partition, which we modify.
     * stateVecIn is the pair node's state-vector partition, in an order which 
     * does not necessarily correlate to stateVecOut's 
     */
    qreal* inReal = stateVecIn.real;
    qreal* inImag = stateVecIn.imag;
    qreal* outReal = stateVecOut.real;
    qreal* outImag = stateVecOut.imag;
    
    long long int outInd, outIndGlobal, inInd, inIndGlobal;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (inReal,inImag,outReal,outImag, numAmps,globalOffset, ctrlMask,targMask) \
    private  (outInd,outIndGlobal, inInd,inIndGlobal)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (outInd = 0; outInd < numAmps; outInd++) {
            
            // modify amplitude only if control qubits are 1 for this state
            outIndGlobal = outInd + globalOffset;
            if (ctrlMask && ((ctrlMask & outIndGlobal) != ctrlMask))
                continue;
            /* it is a premature optimisation to remove this seemingly wasteful abort above,
             * because the maximum skipped amplitudes is 1/2 that stored in the node 
             * (since this function is not called if all amps should be skipped)
             */
             
            /* unlike statevec_controlledNotDistributed(), we cannot assume stateVecOut 
             * maps contiguously/parallel into stateVecIn; we must map each amplitude, bit-wise.
             * However, the arithmetic doesn't necessitate knowing the rank of stateVecIn 
             */
            inIndGlobal = outIndGlobal ^ targMask;
            inInd = inIndGlobal % numAmps;              // = inIndGlobal - pairRank * numAmps
            
            outReal[outInd] = inReal[inInd];
            outImag[outInd] = inImag[inInd];
        }
    }
}

void statevec_pauliYLocal(Qureg qureg, int targetQubit, int conjFac)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateImagUp;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, numTasks,conjFac) \
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

            stateVecReal[indexUp] = conjFac * stateVecImag[indexLo];
            stateVecImag[indexUp] = conjFac * -stateVecReal[indexLo];
            stateVecReal[indexLo] = conjFac * -stateImagUp;
            stateVecImag[indexLo] = conjFac * stateRealUp;
        } 
    }
}

/** Rotate a single qubit by +-{{0,-i},{i,0}.
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk
 *  
 *  @remarks Qubits are zero-based and the                     
 *  the first qubit is the rightmost                  
 *                                                                        
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] stateVecIn probability amplitudes in lower or upper half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 *  @param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
 * @param[in] conjFac 1: apply conj(pauliY), 0: apply pauliY
 */
void statevec_pauliYDistributed(Qureg qureg,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, 
        int updateUpper, int conjFac)
{

    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;

    qreal *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

    int realSign=1, imagSign=1;
    if (updateUpper) imagSign=-1;
    else realSign = -1;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut, \
                realSign,imagSign, numTasks,conjFac) \
    private  (thisTask)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            stateVecRealOut[thisTask] = conjFac * realSign * stateVecImagIn[thisTask];
            stateVecImagOut[thisTask] = conjFac * imagSign * stateVecRealIn[thisTask];
        }
    }
} 




void statevec_controlledPauliYLocal(Qureg qureg, int controlQubit, int targetQubit, int conjFac)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateImagUp;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, numTasks,chunkId, \
                chunkSize,controlQubit,conjFac) \
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
            if (controlBit){
                stateRealUp = stateVecReal[indexUp];
                stateImagUp = stateVecImag[indexUp];

                // update under +-{{0, -i}, {i, 0}}
                stateVecReal[indexUp] = conjFac * stateVecImag[indexLo];
                stateVecImag[indexUp] = conjFac * -stateVecReal[indexLo];
                stateVecReal[indexLo] = conjFac * -stateImagUp;
                stateVecImag[indexLo] = conjFac * stateRealUp;
            }
        } 
    }
}


void statevec_controlledPauliYDistributed (Qureg qureg, int controlQubit,
        ComplexArray stateVecIn,
        ComplexArray stateVecOut, int conjFac)
{

    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    int controlBit;

    qreal *stateVecRealIn=stateVecIn.real, *stateVecImagIn=stateVecIn.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealIn,stateVecImagIn,stateVecRealOut,stateVecImagOut, \
                numTasks,chunkId,chunkSize,controlQubit,conjFac) \
    private  (thisTask,controlBit)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {
            controlBit = extractBit (controlQubit, thisTask+chunkId*chunkSize);
            if (controlBit){
                stateVecRealOut[thisTask] = conjFac * stateVecImagIn[thisTask];
                stateVecImagOut[thisTask] = conjFac * -stateVecRealIn[thisTask];
            }
        }
    }
} 







void statevec_hadamardLocal(Qureg qureg, int targetQubit)
{
    long long int sizeBlock, sizeHalfBlock;
    long long int thisBlock, // current block
         indexUp,indexLo;    // current index and corresponding index in lower half block

    qreal stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;         
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // set dimensions
    sizeHalfBlock = 1LL << targetQubit;  
    sizeBlock     = 2LL * sizeHalfBlock; 

    // Can't use qureg.stateVec as a private OMP var
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

    qreal recRoot2 = 1.0/sqrt(2);

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (sizeBlock,sizeHalfBlock, stateVecReal,stateVecImag, recRoot2, numTasks) \
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
 *  Operate on a subset of the state vector with upper and lower block values 
 *  stored seperately. This rotation is just swapping upper and lower values, and
 *  stateVecIn must already be the correct section for this chunk
 *                                          
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] stateVecUp probability amplitudes in upper half of a block depending on chunkId
 *  @param[in] stateVecLo probability amplitudes in lower half of a block depending on chunkId
 *  @param[out] stateVecOut array section to update (will correspond to either the lower or upper half of a block)
 *  @param[in] updateUpper flag, 1: updating upper values, 0: updating lower values in block
 */
void statevec_hadamardDistributed(Qureg qureg,
        ComplexArray stateVecUp,
        ComplexArray stateVecLo,
        ComplexArray stateVecOut,
        int updateUpper)
{

    qreal   stateRealUp,stateRealLo,stateImagUp,stateImagLo;
    long long int thisTask;  
    long long int numTasks=qureg.numAmpsPerChunk;

    int sign;
    if (updateUpper) sign=1;
    else sign=-1;

    qreal recRoot2 = 1.0/sqrt(2);

    qreal *stateVecRealUp=stateVecUp.real, *stateVecImagUp=stateVecUp.imag;
    qreal *stateVecRealLo=stateVecLo.real, *stateVecImagLo=stateVecLo.imag;
    qreal *stateVecRealOut=stateVecOut.real, *stateVecImagOut=stateVecOut.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (stateVecRealUp,stateVecImagUp,stateVecRealLo,stateVecImagLo,stateVecRealOut,stateVecImagOut, \
            recRoot2, sign, numTasks) \
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

void statevec_phaseShiftByTerm (Qureg qureg, int targetQubit, Complex term)
{       
    long long int index;
    long long int stateVecSize;
    int targetBit;
    
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;
    
    qreal stateRealLo, stateImagLo;
    qreal cosAngle = term.real;
    qreal sinAngle = term.imag;

# ifdef _OPENMP
# pragma omp parallel for \
    default  (none) \
    shared   (stateVecSize, stateVecReal,stateVecImag, cosAngle,sinAngle, \
                chunkId,chunkSize,targetQubit) \
    private  (index,targetBit,stateRealLo,stateImagLo) \
    schedule (static)
# endif
    for (index=0; index<stateVecSize; index++) {
        
        // update the coeff of the |1> state of the target qubit
        targetBit = extractBit (targetQubit, index+chunkId*chunkSize);
        if (targetBit) {
            
            stateRealLo = stateVecReal[index];
            stateImagLo = stateVecImag[index];
            
            stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
            stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;  
        }
    }
}

void statevec_controlledPhaseShift (Qureg qureg, int idQubit1, int idQubit2, qreal angle)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;
    
    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;
    
    qreal stateRealLo, stateImagLo;
    qreal cosAngle = cos(angle);
    qreal sinAngle = sin(angle);

# ifdef _OPENMP
# pragma omp parallel for \
    default  (none)              \
    shared   (stateVecSize, stateVecReal,stateVecImag, chunkId,chunkSize, \
                idQubit1,idQubit2,cosAngle,sinAngle ) \
    private  (index,bit1,bit2,stateRealLo,stateImagLo) \
    schedule (static)
# endif
    for (index=0; index<stateVecSize; index++) {
        bit1 = extractBit (idQubit1, index+chunkId*chunkSize);
        bit2 = extractBit (idQubit2, index+chunkId*chunkSize);
        if (bit1 && bit2) {
            
            stateRealLo = stateVecReal[index];
            stateImagLo = stateVecImag[index];
            
            stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
            stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;  
        }
    }
}

void statevec_multiControlledPhaseShift(Qureg qureg, int *controlQubits, int numControlQubits, qreal angle)
{
    long long int index;
    long long int stateVecSize;

    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    long long int mask = getQubitBitMask(controlQubits, numControlQubits);

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;
    
    qreal stateRealLo, stateImagLo;
    qreal cosAngle = cos(angle);
    qreal sinAngle = sin(angle);

# ifdef _OPENMP
# pragma omp parallel \
    default  (none)              \
    shared   (stateVecSize, stateVecReal, stateVecImag, mask, chunkId,chunkSize,cosAngle,sinAngle) \
    private  (index, stateRealLo, stateImagLo)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            if (mask == (mask & (index+chunkId*chunkSize)) ){
                
                stateRealLo = stateVecReal[index];
                stateImagLo = stateVecImag[index];
            
                stateVecReal[index] = cosAngle*stateRealLo - sinAngle*stateImagLo;
                stateVecImag[index] = sinAngle*stateRealLo + cosAngle*stateImagLo;  
            }
        }
    }
}

int getBitMaskParity(long long int mask) {
    int parity = 0;
    while (mask) {
        parity = !parity;
        mask = mask & (mask-1);
    }
    return parity;
}

void statevec_multiRotateZ(Qureg qureg, long long int mask, qreal angle)
{
    long long int index;
    long long int stateVecSize;

    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;
    
    qreal stateReal, stateImag;
    qreal cosAngle = cos(angle/2.0);
    qreal sinAngle = sin(angle/2.0); 
    
    // = +-1, to flip sinAngle based on target qubit parity, to effect
    // exp(-angle/2 i fac_j)|j>
    int fac; 

# ifdef _OPENMP
# pragma omp parallel \
    default  (none)              \
    shared   (stateVecSize, stateVecReal, stateVecImag, mask, chunkId,chunkSize,cosAngle,sinAngle) \
    private  (index, fac, stateReal, stateImag)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            stateReal = stateVecReal[index];
            stateImag = stateVecImag[index];
            
            // odd-parity target qubits get fac_j = -1
            fac = getBitMaskParity(mask & (index+chunkId*chunkSize))? -1 : 1;
            stateVecReal[index] = cosAngle*stateReal + fac * sinAngle*stateImag;
            stateVecImag[index] = - fac * sinAngle*stateReal + cosAngle*stateImag;  
        }
    }
}

void statevec_multiControlledMultiRotateZ(Qureg qureg, long long int ctrlMask, long long int targMask, qreal angle) 
{
    long long int offset = qureg.chunkId * qureg.numAmpsPerChunk;

    long long int stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;
    
    qreal stateReal, stateImag;
    qreal cosAngle = cos(angle/2.0);
    qreal sinAngle = sin(angle/2.0); 
    
    // = +-1, to flip sinAngle based on target qubit parity, to effect
    // exp(-angle/2 i fac_j)|j>
    int fac; 
    long long int index, globalIndex;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none)              \
    shared   (offset, stateVecSize, stateVecReal,stateVecImag, ctrlMask,targMask, cosAngle,sinAngle) \
    private  (index,globalIndex, fac, stateReal,stateImag)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            stateReal = stateVecReal[index];
            stateImag = stateVecImag[index];
            
            // states with not-all-one control qubits are unmodified
            globalIndex = index + offset;
            if (ctrlMask && ((ctrlMask & globalIndex) != ctrlMask))
                continue;
            
            // odd-parity target qubits get fac_j = -1 (avoid thread divergence)
            fac = 1-2*getBitMaskParity(targMask & globalIndex);
            stateVecReal[index] = cosAngle*stateReal + fac * sinAngle*stateImag;
            stateVecImag[index] = - fac * sinAngle*stateReal + cosAngle*stateImag;  
        }
    }    
}

qreal densmatr_findProbabilityOfZeroLocal(Qureg qureg, int measureQubit) {
    
    // computes first local index containing a diagonal element
    long long int localNumAmps = qureg.numAmpsPerChunk;
    long long int densityDim = (1LL << qureg.numQubitsRepresented);
    long long int diagSpacing = 1LL + densityDim;
    long long int maxNumDiagsPerChunk = 1 + localNumAmps / diagSpacing;
    long long int numPrevDiags = (qureg.chunkId>0)? 1+(qureg.chunkId*localNumAmps)/diagSpacing : 0;
    long long int globalIndNextDiag = diagSpacing * numPrevDiags;
    long long int localIndNextDiag = globalIndNextDiag % localNumAmps;
    
    // computes how many diagonals are contained in this chunk
    long long int numDiagsInThisChunk = maxNumDiagsPerChunk;
    if (localIndNextDiag + (numDiagsInThisChunk-1)*diagSpacing >= localNumAmps)
        numDiagsInThisChunk -= 1;
    
    long long int visitedDiags;     // number of visited diagonals in this chunk so far
    long long int basisStateInd;    // current diagonal index being considered
    long long int index;            // index in the local chunk
    
    qreal zeroProb = 0;
    qreal *stateVecReal = qureg.stateVec.real;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (localIndNextDiag, numPrevDiags, diagSpacing, stateVecReal, numDiagsInThisChunk) \
    private   (visitedDiags, basisStateInd, index) \
    reduction ( +:zeroProb )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        // sums the diagonal elems of the density matrix where measureQubit=0
        for (visitedDiags = 0; visitedDiags < numDiagsInThisChunk; visitedDiags++) {
            
            basisStateInd = numPrevDiags + visitedDiags;
            index = localIndNextDiag + diagSpacing * visitedDiags;
    
            if (extractBit(measureQubit, basisStateInd) == 0)
                zeroProb += stateVecReal[index]; // assume imag[diagonls] ~ 0

        }
    }
    
    return zeroProb;
}

/** Measure the total probability of a specified qubit being in the zero state across all amplitudes in this chunk.
 *  Size of regions to skip is less than the size of one chunk.                   
 *  
 *  @param[in] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @return probability of qubit measureQubit being zero
 */
qreal statevec_findProbabilityOfZeroLocal (Qureg qureg,
        int measureQubit)
{
    // ----- sizes
    long long int sizeBlock,                                  // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                  // current block
         index;                                               // current index for first half block
    // ----- measured probability
    qreal   totalProbability;                                  // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                         // size of blocks (pairs of measure and skip entries)

    // initialise returned value
    totalProbability = 0.0;

    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

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

            totalProbability += stateVecReal[index]*stateVecReal[index]
                + stateVecImag[index]*stateVecImag[index];
        }
    }
    return totalProbability;
}

/** Measure the probability of a specified qubit being in the zero state across all amplitudes held in this chunk.
 * Size of regions to skip is a multiple of chunkSize.
 * The results are communicated and aggregated by the caller
 *  
 *  @param[in] qureg object representing the set of qubits
 *  @return probability of qubit measureQubit being zero
 */
qreal statevec_findProbabilityOfZeroDistributed (Qureg qureg) {
    // ----- measured probability
    qreal   totalProbability;                                  // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    long long int numTasks=qureg.numAmpsPerChunk;

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    // initialise returned value
    totalProbability = 0.0;

    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

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
            totalProbability += stateVecReal[thisTask]*stateVecReal[thisTask]
                + stateVecImag[thisTask]*stateVecImag[thisTask];
        }
    }

    return totalProbability;
}

void statevec_calcProbOfAllOutcomesLocal(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {
    
    /* Below, we manually reduce amplitudes into outcomeProbs by using atomic update.
     * This maintains OpenMP 3.1 compatibility. An alternative is to use array reduction 
     * (requires OpenMP 4.5, limits #qubits since outcomeProbs must be a local stack array)
     * or a dynamic list of omp locks (duplicates memory cost of outcomeProbs).
     * Using locks was always slower than the method below. Using reduction was only 
     * faster for very few threads, or very few outcomeProbs.
     * Finally, we exclude the 'update' clause after 'atomic' to maintain MSVC compatibility 
     */

    long long int numOutcomeProbs = (1 << numQubits);
    long long int j;
    
    // clear outcomeProbs (in parallel, in case it's large)
# ifdef _OPENMP
# pragma omp parallel \
    default (none) \
    shared    (numOutcomeProbs,outcomeProbs) \
    private   (j)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (j=0; j<numOutcomeProbs; j++)
            outcomeProbs[j] = 0;
    }   
    
    long long int numTasks = qureg.numAmpsPerChunk;
    long long int offset = qureg.chunkId*qureg.numAmpsPerChunk;
    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;
    
    long long int i;
    long long int outcomeInd;
    int q;
    qreal prob;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (numTasks,offset, qubits,numQubits, stateRe,stateIm, outcomeProbs) \
    private   (i, q, outcomeInd, prob)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        // every amplitude contributes to a single element of retProbs
        for (i=0; i<numTasks; i++) {
            
            // determine index informed by qubits outcome
            outcomeInd = 0;
            for (q=0; q<numQubits; q++)
                outcomeInd += extractBit(qubits[q], i + offset) * (1LL << q);
            
            prob = stateRe[i]*stateRe[i] + stateIm[i]*stateIm[i];
            
            // atomicly update corresponding outcome array element
            # ifdef _OPENMP
            # pragma omp atomic
            # endif
            outcomeProbs[outcomeInd] += prob;
        }
    }
}

void densmatr_calcProbOfAllOutcomesLocal(qreal* outcomeProbs, Qureg qureg, int* qubits, int numQubits) {
    
    // clear outcomeProbs (in parallel, in case it's large)
    long long int numOutcomeProbs = (1 << numQubits);
    long long int j;
    
# ifdef _OPENMP
# pragma omp parallel \
    default (none) \
    shared    (numOutcomeProbs,outcomeProbs) \
    private   (j)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (j=0; j<numOutcomeProbs; j++)
            outcomeProbs[j] = 0;
    }  
    
    // compute first local index containing a diagonal element
    long long int localNumAmps = qureg.numAmpsPerChunk;
    long long int densityDim = (1LL << qureg.numQubitsRepresented);
    long long int diagSpacing = 1LL + densityDim;
    long long int maxNumDiagsPerChunk = 1 + localNumAmps / diagSpacing;
    long long int numPrevDiags = (qureg.chunkId>0)? 1+(qureg.chunkId*localNumAmps)/diagSpacing : 0;
    long long int globalIndNextDiag = diagSpacing * numPrevDiags;
    long long int localIndNextDiag = globalIndNextDiag % localNumAmps;
    
    // computes how many diagonals are contained in this chunk
    long long int numDiagsInThisChunk = maxNumDiagsPerChunk;
    if (localIndNextDiag + (numDiagsInThisChunk-1)*diagSpacing >= localNumAmps)
        numDiagsInThisChunk -= 1;
        
    long long int visitedDiags;     // number of visited diagonals in this chunk so far
    long long int basisStateInd;    // current diagonal index being considered
    long long int index;            // index in the local chunk
    
    int q;
    long long int outcomeInd;
    qreal *stateRe = qureg.stateVec.real;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (localIndNextDiag, numPrevDiags, diagSpacing, stateRe, numDiagsInThisChunk, qubits,numQubits, outcomeProbs) \
    private   (visitedDiags, basisStateInd, index, q,outcomeInd)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        // sums the diagonal elems of the density matrix where measureQubit=0
        for (visitedDiags = 0; visitedDiags < numDiagsInThisChunk; visitedDiags++) {
            
            basisStateInd = numPrevDiags + visitedDiags;
            index = localIndNextDiag + diagSpacing * visitedDiags;
            
            // determine outcome implied by basisStateInd
            outcomeInd = 0;
            for (q=0; q<numQubits; q++)
                outcomeInd += extractBit(qubits[q], basisStateInd) * (1LL << q);
    
            // atomicly update corresponding outcome array element
            # ifdef _OPENMP
            # pragma omp atomic
            # endif
            outcomeProbs[outcomeInd] += stateRe[index];
        }
    }
}

void statevec_controlledPhaseFlip (Qureg qureg, int idQubit1, int idQubit2)
{
    long long int index;
    long long int stateVecSize;
    int bit1, bit2;

    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;
    
    // dimension of the state vector
    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel for \
    default  (none) \
    shared   (stateVecSize, stateVecReal,stateVecImag, chunkId,chunkSize,idQubit1,idQubit2 ) \
    private  (index,bit1,bit2) \
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

void statevec_multiControlledPhaseFlip(Qureg qureg, int *controlQubits, int numControlQubits)
{
    long long int index;
    long long int stateVecSize;

    long long int chunkSize=qureg.numAmpsPerChunk;
    long long int chunkId=qureg.chunkId;

    long long int mask = getQubitBitMask(controlQubits, numControlQubits);

    stateVecSize = qureg.numAmpsPerChunk;
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none)              \
    shared   (stateVecSize, stateVecReal,stateVecImag, mask, chunkId,chunkSize ) \
    private  (index)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0; index<stateVecSize; index++) {
            if (mask == (mask & (index+chunkId*chunkSize)) ){
                stateVecReal [index] = - stateVecReal [index];
                stateVecImag [index] = - stateVecImag [index];
            }
        }
    }
}

/** Update the state vector to be consistent with measuring measureQubit=0 if outcome=0 and measureQubit=1
 *  if outcome=1.
 *  Performs an irreversible change to the state vector: it updates the vector according
 *  to the event that an outcome have been measured on the qubit indicated by measureQubit (where 
 *  this label starts from 0, of course). It achieves this by setting all inconsistent 
 *  amplitudes to 0 and 
 *  then renormalising based on the total probability of measuring measureQubit=0 or 1 according to the 
 *  value of outcome. 
 *  In the local version, one or more blocks (with measureQubit=0 in the first half of the block and
 *  measureQubit=1 in the second half of the block) fit entirely into one chunk. 
 *  
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @param[in] totalProbability probability of qubit measureQubit being either zero or one
 *  @param[in] outcome to measure the probability of and set the state to -- either zero or one
 */
void statevec_collapseToKnownProbOutcomeLocal(Qureg qureg, int measureQubit, int outcome, qreal totalProbability)
{
    // ----- sizes
    long long int sizeBlock,                                  // size of blocks
         sizeHalfBlock;                                       // size of blocks halved
    // ----- indices
    long long int thisBlock,                                  // current block
         index;                                               // current index for first half block
    // ----- measured probability
    qreal   renorm;                                            // probability (returned) value
    // ----- temp variables
    long long int thisTask;                                   // task based approach for expose loop with small granularity
    // (good for shared memory parallelism)
    long long int numTasks=qureg.numAmpsPerChunk>>1;

    // ---------------------------------------------------------------- //
    //            dimensions                                            //
    // ---------------------------------------------------------------- //
    sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
    // and then the number to skip
    sizeBlock     = 2LL * sizeHalfBlock;                         // size of blocks (pairs of measure and skip entries)

    renorm=1/sqrt(totalProbability);
    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;


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
 *  in state 0 or 1.
 *  Measure in Zero performs an irreversible change to the state vector: it updates the vector according
 *  to the event that the value 'outcome' has been measured on the qubit indicated by measureQubit (where 
 *  this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
 *  then renormalising based on the total probability of measuring measureQubit=0 if outcome=0 and
 *  measureQubit=1 if outcome=1.
 *  In the distributed version, one block (with measureQubit=0 in the first half of the block and
 *  measureQubit=1 in the second half of the block) is spread over multiple chunks, meaning that each chunks performs
 *  only renormalisation or only setting amplitudes to 0. This function handles the renormalisation.
 *  
 *  @param[in,out] qureg object representing the set of qubits
 *  @param[in] measureQubit qubit to measure
 *  @param[in] totalProbability probability of qubit measureQubit being zero
 */
void statevec_collapseToKnownProbOutcomeDistributedRenorm (Qureg qureg, int measureQubit, qreal totalProbability)
{
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=qureg.numAmpsPerChunk;

    qreal renorm=1/sqrt(totalProbability);

    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

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
            stateVecReal[thisTask] = stateVecReal[thisTask]*renorm;
            stateVecImag[thisTask] = stateVecImag[thisTask]*renorm;
        }
    }
}

/* Set all amplitudes in one chunk to 0. 
 *  Measure in Zero performs an irreversible change to the state vector: it updates the vector according
 *  to the event that a zero have been measured on the qubit indicated by measureQubit (where 
 *  this label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
 *  then renormalising based on the total probability of measuring measureQubit=0 or 1.
 *  In the distributed version, one block (with measureQubit=0 in the first half of the block and
 *  measureQubit=1 in the second half of the block) is spread over multiple chunks, meaning that each chunks performs
 *  only renormalisation or only setting amplitudes to 0. This function handles setting amplitudes to 0.
 *
 *  @param[in,out] qureg object representing the set of qubits
 */
void statevec_collapseToOutcomeDistributedSetZero(Qureg qureg)
{
    // ----- temp variables
    long long int thisTask;                                   
    long long int numTasks=qureg.numAmpsPerChunk;

    // ---------------------------------------------------------------- //
    //            find probability                                      //
    // ---------------------------------------------------------------- //

    qreal *stateVecReal = qureg.stateVec.real;
    qreal *stateVecImag = qureg.stateVec.imag;

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
            stateVecReal[thisTask] = 0;
            stateVecImag[thisTask] = 0;
        }
    }
}

/** It is ensured that all amplitudes needing to be swapped are on this node.
 * This means that amplitudes for |a 0..0..> to |a 1..1..> all exist on this node 
 * and each node has a different bit-string prefix "a". The prefix 'a' (and ergo,
 * the chunkID) don't enter the calculations for the offset of |a 0..1..> and 
 * |a 1..0..> from |a 0..0..> and ergo are not featured below.
 */
void statevec_swapQubitAmpsLocal(Qureg qureg, int qb1, int qb2) {

    // can't use qureg.stateVec as a private OMP var
    qreal *reVec = qureg.stateVec.real;
    qreal *imVec = qureg.stateVec.imag;
    
    long long int numTasks = qureg.numAmpsPerChunk >> 2; // each iteration updates 2 amps and skips 2 amps
    long long int thisTask;
    long long int ind00, ind01, ind10;
    qreal re01, re10;
    qreal im01, im10;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (reVec,imVec,numTasks,qb1,qb2) \
    private  (thisTask, ind00,ind01,ind10, re01,re10, im01,im10) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (thisTask=0; thisTask<numTasks; thisTask++) {    
            // determine ind00 of |..0..0..>, |..0..1..> and |..1..0..>
            ind00 = insertTwoZeroBits(thisTask, qb1, qb2);
            ind01 = flipBit(ind00, qb1);
            ind10 = flipBit(ind00, qb2);

            // extract statevec amplitudes 
            re01 = reVec[ind01]; im01 = imVec[ind01];
            re10 = reVec[ind10]; im10 = imVec[ind10];

            // swap 01 and 10 amps
            reVec[ind01] = re10; reVec[ind10] = re01;
            imVec[ind01] = im10; imVec[ind10] = im01;
        }
    }
}

/** qureg.pairStateVec contains the entire set of amplitudes of the paired node
 * which includes the set of all amplitudes which need to be swapped between
 * |..0..1..> and |..1..0..>
 */
void statevec_swapQubitAmpsDistributed(Qureg qureg, int pairRank, int qb1, int qb2) {
    
    // can't use qureg.stateVec as a private OMP var
    qreal *reVec = qureg.stateVec.real;
    qreal *imVec = qureg.stateVec.imag;
    qreal *rePairVec = qureg.pairStateVec.real;
    qreal *imPairVec = qureg.pairStateVec.imag;
    
    long long int numLocalAmps = qureg.numAmpsPerChunk;
    long long int globalStartInd = qureg.chunkId * numLocalAmps;
    long long int pairGlobalStartInd = pairRank * numLocalAmps;

    long long int localInd, globalInd;
    long long int pairLocalInd, pairGlobalInd;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (reVec,imVec,rePairVec,imPairVec,numLocalAmps,globalStartInd,pairGlobalStartInd,qb1,qb2) \
    private  (localInd,globalInd, pairLocalInd,pairGlobalInd) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (localInd=0; localInd < numLocalAmps; localInd++) { 
            
            globalInd = globalStartInd + localInd;
            if (isOddParity(globalInd, qb1, qb2)) {
                
                pairGlobalInd = flipBit(flipBit(globalInd, qb1), qb2);
                pairLocalInd = pairGlobalInd - pairGlobalStartInd;
                
                reVec[localInd] = rePairVec[pairLocalInd];
                imVec[localInd] = imPairVec[pairLocalInd];
            }
        }
    }
}

void statevec_setWeightedQureg(Complex fac1, Qureg qureg1, Complex fac2, Qureg qureg2, Complex facOut, Qureg out) {

    long long int numAmps = qureg1.numAmpsPerChunk;

    qreal *vecRe1 = qureg1.stateVec.real;
    qreal *vecIm1 = qureg1.stateVec.imag;
    qreal *vecRe2 = qureg2.stateVec.real;
    qreal *vecIm2 = qureg2.stateVec.imag;
    qreal *vecReOut = out.stateVec.real;
    qreal *vecImOut = out.stateVec.imag;

    qreal facRe1 = fac1.real; 
    qreal facIm1 = fac1.imag;
    qreal facRe2 = fac2.real;
    qreal facIm2 = fac2.imag;
    qreal facReOut = facOut.real;
    qreal facImOut = facOut.imag;

    qreal re1,im1, re2,im2, reOut,imOut;
    long long int index;

# ifdef _OPENMP
# pragma omp parallel \
    shared    (vecRe1,vecIm1, vecRe2,vecIm2, vecReOut,vecImOut, facRe1,facIm1,facRe2,facIm2, numAmps) \
    private   (index, re1,im1, re2,im2, reOut,imOut)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0LL; index<numAmps; index++) {
            re1 = vecRe1[index]; im1 = vecIm1[index];
            re2 = vecRe2[index]; im2 = vecIm2[index];
            reOut = vecReOut[index];
            imOut = vecImOut[index];

            vecReOut[index] = (facReOut*reOut - facImOut*imOut) + (facRe1*re1 - facIm1*im1) + (facRe2*re2 - facIm2*im2);
            vecImOut[index] = (facReOut*imOut + facImOut*reOut) + (facRe1*im1 + facIm1*re1) + (facRe2*im2 + facIm2*re2);
        }
    }
}

void statevec_applyDiagonalOp(Qureg qureg, DiagonalOp op) {

    // each node/chunk modifies only its values in an embarrassingly parallelisable way
    long long int numAmps = qureg.numAmpsPerChunk;

    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;
    qreal* opRe = op.real;
    qreal* opIm = op.imag;

    qreal a,b,c,d;
    long long int index;

# ifdef _OPENMP
# pragma omp parallel \
    shared    (stateRe,stateIm, opRe,opIm, numAmps) \
    private   (index, a,b,c,d)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0LL; index<numAmps; index++) {
            a = stateRe[index];
            b = stateIm[index];
            c = opRe[index];
            d = opIm[index];

            // (a + b i)(c + d i) = (a c - b d) + i (a d + b c)
            stateRe[index] = a*c - b*d;
            stateIm[index] = a*d + b*c;
        }
    }
}

void densmatr_applyDiagonalOpLocal(Qureg qureg, DiagonalOp op) {
    
    /* ALL values of op are pre-loaded into qureg.pairStateVector (on every node).
     * Furthermore, since it's gauranteed each node contains an integer number of 
     * columns of qureg (because op upperlimits the number of nodes; 1 per element),
     * then we know iteration below begins at the 'top' of a column, and there is 
     * no offset for op (pairStateVector)
     */

    long long int numAmps = qureg.numAmpsPerChunk;
    int opDim = (1 << op.numQubits);

    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;
    qreal* opRe = qureg.pairStateVec.real;
    qreal* opIm = qureg.pairStateVec.imag;
    
    qreal a,b,c,d;
    long long int index;

# ifdef _OPENMP
# pragma omp parallel \
    shared    (stateRe,stateIm, opRe,opIm, numAmps,opDim) \
    private   (index, a,b,c,d)
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0LL; index<numAmps; index++) {
            a = stateRe[index];
            b = stateIm[index];
            c = opRe[index % opDim];
            d = opIm[index % opDim];

            // (a + b i)(c + d i) = (a c - b d) + i (a d + b c)
            stateRe[index] = a*c - b*d;
            stateIm[index] = a*d + b*c;
        }
    }
}

Complex statevec_calcExpecDiagonalOpLocal(Qureg qureg, DiagonalOp op) {
    
    qreal expecRe = 0;
    qreal expecIm = 0;
    
    long long int index;
    long long int numAmps = qureg.numAmpsPerChunk;
    qreal *stateReal = qureg.stateVec.real;
    qreal *stateImag = qureg.stateVec.imag;
    qreal *opReal = op.real;
    qreal *opImag = op.imag;
    
    qreal vecRe,vecIm,vecAbs, opRe, opIm;
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (stateReal, stateImag, opReal, opImag, numAmps) \
    private   (index, vecRe,vecIm,vecAbs, opRe,opIm) \
    reduction ( +:expecRe, expecIm )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (index=0; index < numAmps; index++) {
            vecRe = stateReal[index];
            vecIm = stateImag[index];
            opRe = opReal[index];
            opIm = opImag[index];
            
            // abs(vec)^2 op
            vecAbs = vecRe*vecRe + vecIm*vecIm;
            expecRe += vecAbs*opRe;
            expecIm += vecAbs*opIm;
        }
    }
    
    Complex innerProd;
    innerProd.real = expecRe;
    innerProd.imag = expecIm;
    return innerProd;
}

Complex densmatr_calcExpecDiagonalOpLocal(Qureg qureg, DiagonalOp op) {
    
    /* since for every 1 element in \p op, there exists a column in \p qureg, 
     * we know that the elements in \p op live on the same node as the 
     * corresponding diagonal elements of \p qureg. This means, the problem is 
     * embarrassingly parallelisable, and the code below works for both 
     * serial and distributed modes.
     */
    
    // computes first local index containing a diagonal element
    long long int diagSpacing = 1LL + (1LL << qureg.numQubitsRepresented);
    long long int numPrevDiags = (qureg.chunkId>0)? 1+(qureg.chunkId*qureg.numAmpsPerChunk)/diagSpacing : 0;
    long long int globalIndNextDiag = diagSpacing * numPrevDiags;
    long long int localIndNextDiag = globalIndNextDiag % qureg.numAmpsPerChunk;
    long long int numAmps = qureg.numAmpsPerChunk;
    
    qreal* stateReal = qureg.stateVec.real;
    qreal* stateImag = qureg.stateVec.imag;
    qreal* opReal = op.real;
    qreal* opImag = op.imag;
    
    qreal expecRe = 0;
    qreal expecIm = 0;
    
    long long int stateInd;
    long long int opInd;
    qreal matRe, matIm, opRe, opIm;
    
    // visits every diagonal element with global index (2^n + 1)i for i in [0, 2^n-1]
    
# ifdef _OPENMP
# pragma omp parallel \
    shared    (stateReal,stateImag, opReal,opImag, localIndNextDiag,diagSpacing,numAmps) \
    private   (stateInd,opInd, matRe,matIm, opRe,opIm) \
    reduction ( +:expecRe, expecIm )
# endif 
    {
# ifdef _OPENMP
# pragma omp for schedule  (static)
# endif
        for (stateInd=localIndNextDiag; stateInd < numAmps; stateInd += diagSpacing) {
            
            matRe = stateReal[stateInd];
            matIm = stateImag[stateInd];
            opInd = (stateInd - localIndNextDiag) / diagSpacing;
            opRe = opReal[opInd];
            opIm = opImag[opInd];
            
            // (matRe + matIm i)(opRe + opIm i) = 
            //      (matRe opRe - matIm opIm) + i (matRe opIm + matIm opRe)
            expecRe += matRe * opRe - matIm * opIm;
            expecIm += matRe * opIm + matIm * opRe;
        }
    }
    
    Complex expecVal;
    expecVal.real = expecRe;
    expecVal.imag = expecIm;
    return expecVal;
}

void agnostic_setDiagonalOpElems(DiagonalOp op, long long int startInd, qreal* real, qreal* imag, long long int numElems) {
    
    // local start/end indices of the given amplitudes, assuming they fit in this chunk
    // these may be negative or above qureg.numAmpsPerChunk
    long long int localStartInd = startInd - op.chunkId*op.numElemsPerChunk;
    long long int localEndInd = localStartInd + numElems; // exclusive
    
    // add this to a local index to get corresponding elem in reals & imags
    long long int offset = op.chunkId*op.numElemsPerChunk - startInd;
    
    // restrict these indices to fit into this chunk
    if (localStartInd < 0)
        localStartInd = 0;
    if (localEndInd > op.numElemsPerChunk)
        localEndInd = op.numElemsPerChunk;
    // they may now be out of order = no iterations
    
    // unpacking OpenMP vars
    long long int index;
    qreal* vecRe = op.real;
    qreal* vecIm = op.imag;
    
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (localStartInd,localEndInd, vecRe,vecIm, real,imag, offset) \
    private  (index) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        // iterate these local inds - this might involve no iterations
        for (index=localStartInd; index < localEndInd; index++) {
            vecRe[index] = real[index + offset];
            vecIm[index] = imag[index + offset];
        }
    }
}

void statevec_applyPhaseFuncOverrides(
    Qureg qureg, int* qubits, int numQubits, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int numTerms, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj)
{
    // each node/chunk modifies only local values in an embarrassingly parallel way 

    // thread shared vars
    int chunkId = qureg.chunkId;
    long long int numAmps = qureg.numAmpsPerChunk;
    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;

    // thread private vars
    long long int index, globalAmpInd, phaseInd;
    int i, t, q;
    qreal phase, c, s, re, im;

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkId,numAmps, stateRe,stateIm, qubits,numQubits,encoding, coeffs,exponents,numTerms, overrideInds,overridePhases,numOverrides, conj) \
    private  (index, globalAmpInd, phaseInd, i,t,q, phase, c,s,re,im) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0LL; index<numAmps; index++) {

            // determine global amplitude index 
            globalAmpInd = chunkId * numAmps + index;

            // determine phase index of {qubits}
            phaseInd = 0LL;
            if (encoding == UNSIGNED) {
                for (q=0; q<numQubits; q++) // use significance order specified by {qubits}
                    phaseInd += (1LL << q) * extractBit(qubits[q], globalAmpInd);
            } 
            else if (encoding == TWOS_COMPLEMENT) {
                for (q=0; q<numQubits-1; q++) // use final qubit to indicate sign 
                    phaseInd += (1LL << q) * extractBit(qubits[q], globalAmpInd);
                if (extractBit(qubits[numQubits-1], globalAmpInd) == 1)
                    phaseInd -= (1LL << (numQubits-1));
            }

            // determine if this phase index has an overriden value (i < numOverrides)
            for (i=0; i<numOverrides; i++)
                if (phaseInd == overrideInds[i])
                    break;

            // determine phase from {coeffs}, {exponents} (unless overriden)
            phase = 0;
            if (i < numOverrides)
                phase = overridePhases[i];
            else
                for (t=0; t<numTerms; t++)
                    phase += coeffs[t] * pow(phaseInd, exponents[t]);
                    
            // negate phase to conjugate operator 
            if (conj)
                phase *= -1;

            // modify amp to amp * exp(i phase) 
            c = cos(phase);
            s = sin(phase);
            re = stateRe[index];
            im = stateIm[index];

            // = {re[amp] cos(phase) - im[amp] sin(phase)} + i {re[amp] sin(phase) + im[amp] cos(phase)}
            stateRe[index] = re*c - im*s;
            stateIm[index] = re*s + im*c;
        }
    }
}

void statevec_applyMultiVarPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    qreal* coeffs, qreal* exponents, int* numTermsPerReg, 
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj) 
{
    // each node/chunk modifies only local values in an embarrassingly parallel way 

    // note partitions of qubits, coeffs, exponents and overrideInds are stored flat

    // thread-shared vaes
    int chunkId = qureg.chunkId;
    long long int numAmps = qureg.numAmpsPerChunk;
    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;

    // thread-private vars
    long long int index, globalAmpInd;
    int r, q, i, t, found, flatInd;
    qreal phase, c, s, re, im;

    // each thread has a private static array of length >= numRegs (private var-length is illegal)
    long long int phaseInds[MAX_NUM_REGS_APPLY_ARBITRARY_PHASE];

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkId,numAmps, stateRe,stateIm, qubits,numQubitsPerReg,numRegs,encoding, coeffs,exponents,numTermsPerReg, overrideInds,overridePhases,numOverrides, conj) \
    private  (index,globalAmpInd, r,q,i,t,flatInd, found, phaseInds,phase, c,s,re,im) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0LL; index<numAmps; index++) {

            // determine global amplitude index 
            globalAmpInd = chunkId * numAmps + index;

            // determine phase indices
            flatInd = 0;
            for (r=0; r<numRegs; r++) {
                phaseInds[r] = 0LL;

                if (encoding == UNSIGNED) {
                    for (q=0; q<numQubitsPerReg[r]; q++)
                        phaseInds[r] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);   // qubits[flatInd] ~ qubits[r][q]
                }
                else if (encoding == TWOS_COMPLEMENT) {
                    for (q=0; q<numQubitsPerReg[r]-1; q++)  
                        phaseInds[r] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);
                    // use final qubit to indicate sign
                    if (extractBit(qubits[flatInd++], globalAmpInd) == 1)
                        phaseInds[r] -= (1LL << (numQubitsPerReg[r]-1));
                }
            }

            // determine if this phase index has an overriden value (i < numOverrides)
            for (i=0; i<numOverrides; i++) {
                found = 1;
                for (r=0; r<numRegs; r++) {
                    if (phaseInds[r] != overrideInds[i*numRegs+r]) {
                        found = 0;
                        break;
                    }
                }
                if (found)
                    break;
            }

            // compute the phase (unless overriden)
            phase = 0;
            if (i < numOverrides)
                phase = overridePhases[i];
            else {
                flatInd = 0;
                for (r=0; r<numRegs; r++) {
                    for (t=0; t<numTermsPerReg[r]; t++) {
                        phase += coeffs[flatInd] * pow(phaseInds[r], exponents[flatInd]);
                        flatInd++;
                    }
                }
            }
            
            // negate phase to conjugate operator 
            if (conj)
                phase *= -1;

            // modify amp to amp * exp(i phase) 
            c = cos(phase);
            s = sin(phase);
            re = stateRe[index];
            im = stateIm[index];

            // = {re[amp] cos(phase) - im[amp] sin(phase)} + i {re[amp] sin(phase) + im[amp] cos(phase)}
            stateRe[index] = re*c - im*s;
            stateIm[index] = re*s + im*c;
        }
    }
}

void statevec_applyParamNamedPhaseFuncOverrides(
    Qureg qureg, int* qubits, int* numQubitsPerReg, int numRegs, enum bitEncoding encoding,
    enum phaseFunc phaseFuncName,
    qreal* params, int numParams,
    long long int* overrideInds, qreal* overridePhases, int numOverrides,
    int conj
) {
    // each node/chunk modifies only local values in an embarrassingly parallel way 

    // note partitions of qubits, overrideInds are stored flat

    // thread-shared vaes
    int chunkId = qureg.chunkId;
    long long int numAmps = qureg.numAmpsPerChunk;
    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;

    // thread-private vars
    long long int index, globalAmpInd;
    int r, q, i, found, flatInd;
    qreal phase, norm, prod, dist, c, s, re, im;

    // each thread has a private static array of length >= numRegs (private var-length is illegal)
    long long int phaseInds[MAX_NUM_REGS_APPLY_ARBITRARY_PHASE];

# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (chunkId,numAmps, stateRe,stateIm, qubits,numQubitsPerReg,numRegs,encoding, phaseFuncName,params,numParams, overrideInds,overridePhases,numOverrides, conj) \
    private  (index,globalAmpInd, r,q,i,flatInd, found, phaseInds,phase,norm,prod,dist, c,s,re,im) 
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (index=0LL; index<numAmps; index++) {

            // determine global amplitude index 
            globalAmpInd = chunkId * numAmps + index;

            // determine phase indices
            flatInd = 0;
            for (r=0; r<numRegs; r++) {
                phaseInds[r] = 0LL;

                if (encoding == UNSIGNED) {
                    for (q=0; q<numQubitsPerReg[r]; q++)
                        phaseInds[r] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);   // qubits[flatInd] ~ qubits[r][q]
                }
                else if (encoding == TWOS_COMPLEMENT) {
                    for (q=0; q<numQubitsPerReg[r]-1; q++)  
                        phaseInds[r] += (1LL << q) * extractBit(qubits[flatInd++], globalAmpInd);
                    // use final qubit to indicate sign
                    if (extractBit(qubits[flatInd++], globalAmpInd) == 1)
                        phaseInds[r] -= (1LL << (numQubitsPerReg[r]-1));
                }
            }

            // determine if this phase index has an overriden value (i < numOverrides)
            for (i=0; i<numOverrides; i++) {
                found = 1;
                for (r=0; r<numRegs; r++) {
                    if (phaseInds[r] != overrideInds[i*numRegs+r]) {
                        found = 0;
                        break;
                    }
                }
                if (found)
                    break;
            }

            // compute the phase (unless overriden)
            phase = 0;
            if (i < numOverrides)
                phase = overridePhases[i];
            else {
                // compute norm related phases
                if (phaseFuncName == NORM || phaseFuncName == INVERSE_NORM ||
                    phaseFuncName == SCALED_NORM || phaseFuncName == SCALED_INVERSE_NORM ||
                    phaseFuncName == SCALED_INVERSE_SHIFTED_NORM) {

                    norm = 0;
                    if (phaseFuncName == SCALED_INVERSE_SHIFTED_NORM) {
                        for (r=0; r<numRegs; r++)
                            norm += (phaseInds[r] - params[2+r])*(phaseInds[r] - params[2+r]);
                    }
                    else
                        for (r=0; r<numRegs; r++)
                            norm += phaseInds[r]*phaseInds[r];
                    norm = sqrt(norm);
 
                    if (phaseFuncName == NORM)
                        phase = norm;
                    else if (phaseFuncName == INVERSE_NORM)
                        phase = (norm == 0.)? params[0] : 1/norm;  // smallest non-zero norm is 1
                    else if (phaseFuncName == SCALED_NORM)
                        phase = params[0] * norm;
                    else if (phaseFuncName == SCALED_INVERSE_NORM || phaseFuncName == SCALED_INVERSE_SHIFTED_NORM)
                        phase = (norm <= REAL_EPS)? params[1] : params[0] / norm; // unless shifted closer to zero
                }
                // compute product related phases
                else if (phaseFuncName == PRODUCT || phaseFuncName == INVERSE_PRODUCT ||
                         phaseFuncName == SCALED_PRODUCT || phaseFuncName == SCALED_INVERSE_PRODUCT) {

                    prod = 1;
                    for (r=0; r<numRegs; r++)
                        prod *= phaseInds[r];

                    if (phaseFuncName == PRODUCT)
                        phase = prod;
                    else if (phaseFuncName == INVERSE_PRODUCT)
                        phase = (prod == 0.)? params[0] : 1/prod;  // smallest non-zero product norm is +- 1
                    else if (phaseFuncName == SCALED_PRODUCT)
                        phase = params[0] * prod;
                    else if (phaseFuncName == SCALED_INVERSE_PRODUCT)
                        phase = (prod == 0.)? params[1] : params[0] / prod;
                }
                // compute Euclidean distance related phases 
                else if (phaseFuncName == DISTANCE || phaseFuncName == INVERSE_DISTANCE ||
                         phaseFuncName == SCALED_DISTANCE || phaseFuncName == SCALED_INVERSE_DISTANCE ||
                         phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE || phaseFuncName == SCALED_INVERSE_SHIFTED_WEIGHTED_DISTANCE) {

                    dist = 0;
                    if (phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE) {
                        for (r=0; r<numRegs; r+=2)
                            dist += (phaseInds[r] - phaseInds[r+1] - params[2+r/2])*(phaseInds[r] - phaseInds[r+1] - params[2+r/2]);
                    }
                    else if (phaseFuncName == SCALED_INVERSE_SHIFTED_WEIGHTED_DISTANCE) {
                        for (r=0; r<numRegs; r+=2)
                            dist += params[2+r] * (phaseInds[r] - phaseInds[r+1] - params[2+r+1])*(phaseInds[r] - phaseInds[r+1] - params[2+r+1]);
                    }
                    else
                        for (r=0; r<numRegs; r+=2)
                            dist += (phaseInds[r+1] - phaseInds[r])*(phaseInds[r+1] - phaseInds[r]);

                    // if sqrt() arg would be negative, set it to divergence param
                    if (dist < 0)
                        dist = 0;

                    dist = sqrt(dist);

                    if (phaseFuncName == DISTANCE)
                        phase = dist;
                    else if (phaseFuncName == INVERSE_DISTANCE)
                        phase = (dist == 0.)? params[0] : 1/dist; // smallest non-zero dist is 1
                    else if (phaseFuncName == SCALED_DISTANCE)
                        phase = params[0] * dist;
                    else if (phaseFuncName == SCALED_INVERSE_DISTANCE || phaseFuncName == SCALED_INVERSE_SHIFTED_DISTANCE || phaseFuncName == SCALED_INVERSE_SHIFTED_WEIGHTED_DISTANCE)
                        phase = (dist <= REAL_EPS)? params[1] : params[0] / dist; // unless shifted closer to 0
                }
            }
            
            // negate phase to conjugate operator 
            if (conj)
                phase *= -1;

            // modify amp to amp * exp(i phase) 
            c = cos(phase);
            s = sin(phase);
            re = stateRe[index];
            im = stateIm[index];

            // = {re[amp] cos(phase) - im[amp] sin(phase)} + i {re[amp] sin(phase) + im[amp] cos(phase)}
            stateRe[index] = re*c - im*s;
            stateIm[index] = re*s + im*c;
        }
    }
}

void densmatr_setQuregToPauliHamil(Qureg qureg, PauliHamil hamil) {

    // flattened {I,X,Y,Z} matrix elements, where [k] = [p][i][j]
    int pauliRealElems[] = {   1,0, 0,1,   0,1, 1,0,   0,0, 0,0,   1,0, 0,-1  };
    int pauliImagElems[] = {   0,0, 0,0,   0,0, 0,0,   0,-1,1,0,   0,0, 0,0   };
    
    // constants (unpacked to prevent OpenMP struct tantrum)
    int numTerms = hamil.numSumTerms;
    int numQubits = hamil.numQubits;
    enum pauliOpType* paulis = hamil.pauliCodes;
    qreal* termCoeffs = hamil.termCoeffs;
    qreal* stateRe = qureg.stateVec.real;
    qreal* stateIm = qureg.stateVec.imag;
    long long int numLocAmps = qureg.numAmpsPerChunk;
    long long int nOffset = numLocAmps * qureg.chunkId;
    
    // private threading vars
    long long int n, nGlob, r, c, t, pInd;
    int q, i, j, p, k;
    qreal elemRe, elemIm;
    int kronRe, kronIm;
    int pauliRe, pauliIm, tmp;
    
    // use a private thread for every overwritten amplitude
# ifdef _OPENMP
# pragma omp parallel \
    default  (none) \
    shared   (pauliRealElems,pauliImagElems, numTerms,numQubits,numLocAmps,nOffset,  paulis,termCoeffs, stateRe,stateIm) \
    private  (q,i,j,p,k, r,c,n,nGlob,t,pInd, elemRe,elemIm, kronRe,kronIm, pauliRe,pauliIm, tmp)
# endif
    {
# ifdef _OPENMP
# pragma omp for schedule (static)
# endif
        for (n=0; n<numLocAmps; n++) {
            
            // |nGlob> = |rank>|n>
            nGlob = nOffset + n;
            
            // |nGlob> = |c>|r>
            r = nGlob & ((1LL << numQubits) - 1);
            c = nGlob >> numQubits;

            // new amplitude of |n>
            elemRe = 0;
            elemIm = 0;
            
            for (t=0; t<numTerms; t++) {
                
                // pauliKronecker[r][c] = prod_q Pauli[q][q-th bit of r and c]
                kronRe = 1;
                kronIm = 0;
                pInd = numQubits * t;
                
                for (q=0; q<numQubits; q++) {
                    
                    // get element of Pauli matrix
                    i = (r >> q) & 1;
                    j = (c >> q) & 1;
                    p = (int) paulis[pInd++];
                    k = (p<<2) + (i<<1) + j;
                    pauliRe = pauliRealElems[k]; 
                    pauliIm = pauliImagElems[k];
                    
                    // kron *= pauli
                    tmp = (pauliRe*kronRe) - (pauliIm*kronIm);
                    kronIm = (pauliRe*kronIm) + (pauliIm*kronRe);
                    kronRe = tmp;
                }
                
                // elem = sum_t coeffs[t] pauliKronecker[r][c]
                elemRe += termCoeffs[t] * kronRe;
                elemIm += termCoeffs[t] * kronIm;
            }
            
            // overwrite the density matrix entry
            stateRe[n] = elemRe;
            stateIm[n] = elemIm;
        }
    }
}
