/** @file
 * Functions for communicating and exchanging amplitudes between compute
 * nodes, when running in distributed mode, using the C MPI standard.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/gpu/gpu_config.hpp"
#include "quest/src/comm/comm_config.hpp"

#include <vector>

#if ENABLE_DISTRIBUTION
    #include <mpi.h>
#endif



/*
 * TODO:
 *
 * - assert internal succes of all MPI calls for robustness
 * 
 * - we may wish to adjust communication design when supporting MPI+cuQuantum,
 *   wherein we must define our own Communicator type.
 *   https://github.com/NVIDIA/cuQuantum/blob/788c04862241c52a5435982e6c25b30b6e3e9324/samples/custatevec/samples_mpi/distributedIndexBitSwap.cpp#L177
 *
 * - apparently we can perform better than UCX's intra-socket inter-GPU
 *   communication, achieving 1.5-3x speedups. See this thesis:
 *   https://www.queensu.ca/academia/afsahi/pprl/thesis/Yiltan_Temucin_MASc_thesis.pdf
 * 
 * - when CUDA-aware MPI (via UCX) attempting to exchange VRAM-to-VRAM has to fall
 *   back upon routing through RAM (because e.g. there is no direct
 *   interconnect between the GPUs), it will(?) create its own RAM buffer,
 *   and cannot(?) be asked to use our existing Qureg RAM buffer. This can
 *   be catastrophic; we run out of RAM in our inter-DGX tests. Perhaps we should
 *   explicitly use UCX runtime checks to detect this scenario?
 *   *
 *   Should we just consult CUDA peer-to-peer check (either all-to-all acrosss all machines,
 *   OR specifically for the nodes we wish toc communicate between, i.e. rank & deviceId specific)
 *   and if it's false (GPUDirect definitely impossible, but perhaps another optimised mode I
 *   don't know about could be used), just route through RAM ourselves? This is at least safe;
 *   we should never have UCX (or otherwise) manually malloc CPU memory. 
 *   Must read up if there is a runtime scenario later where P2P reports false for an optimised
 *   VRAM-to-VRAM method which we'd still like to make use of. Test this on ARC!
 *   *
 *   Temp proposed implementation:
 *   > determine if inter-com GPUs are on the same machine
 *     (the assumption being all GPUs queried with seeLocalGpus() or whatever are on the same
 *      physical machine).
 *   > If not: never VRAM to VRAM; route via CPU
 *   > If so: query CUda P2P, and use VRAM-to_VRAM if so, else route through CPU
 *   *
 *   Would this temporary method preclude UCX from using shared-memory methods to exchange
 *   between same-machine GPUs when CUDA P2P reports no?? 
 * 
 * - look into UCX CUDA multi-rail:
 *   https://docs.nvidia.com/networking/display/hpcxv215/unified+communication+-+x+framework+library#src-119764120_UnifiedCommunicationXFrameworkLibrary-Multi-RailMulti-Rail 
 */



/*
 * MAXIMUM SINGLE MESSAGE LENGTH
 *
 * Messages larger than this will be divided up into smaller
 * messages and exchanged in parallel / asynchronously. Note
 * that 32 is a hard upper limit (before MPI length primitives
 * overflow), 29 (at double precision) seems to be the upper-limit 
 * for UCX exchange, so 28 (to support quad precision) seems safe
 */

qindex MAX_MESSAGE_LENGTH = powerOf2(28);



/*
 * MESSAGE SIZES
 */


int NULL_TAG = 0;

void getMessageConfig(qindex *messageSize, qindex *numMessages, qindex numAmps) {

    // determine the number of max-size messages
    *messageSize = MAX_MESSAGE_LENGTH;
    *numMessages = numAmps / *messageSize; // gauranteed to divide evenly

    // when numAmps < messageSize, we need send a single (smaller) message
    if (*numMessages == 0) {
        *messageSize = numAmps;
        *numMessages = 1;
    }
}



/*
 * PRIVATE SYNCHRONOUSE AMPLITUDE EXCHANGE
 */


void exchangeArrays(qcomp* send, qcomp* recv, qindex numElems, int pairRank) {
#if ENABLE_DISTRIBUTION

    // each message is asynchronously dispatched with a final wait, as per arxiv.org/abs/2308.07402

    // divide the data into multiple messages
    qindex messageSize, numMessages;
    getMessageConfig(&messageSize, &numMessages, numElems);

    // each asynch message below will create two requests for subsequent synch
    std::vector<MPI_Request> requests(2*numMessages, MPI_REQUEST_NULL);

    // asynchronously exchange the messages (effecting MPI_Isendrecv), exploiting orderedness gaurantee.
    // note the exploitation of orderedness means we cannot use UCX's adaptive-routing (AR).
    for (qindex m=0; m<numMessages; m++) {
        MPI_Isend(&send[m*messageSize], messageSize, MPI_QCOMP, pairRank, NULL_TAG, MPI_COMM_WORLD, &requests[2*m]);
        MPI_Irecv(&recv[m*messageSize], messageSize, MPI_QCOMP, pairRank, NULL_TAG, MPI_COMM_WORLD, &requests[2*m+1]);
    }

    // wait for all exchanges to complete (MPI willl automatically free the request memory)
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

#else
    error_commButEnvNotDistributed();
#endif
}



/*
 * PRIVATE ASYNCH SEND AND RECEIVE
 */


void asynchSendArray(qcomp* send, qindex numElems, int pairRank) {
#if ENABLE_DISTRIBUTION

    // we will not track nor wait for the asynch send; instead, the caller will later comm_sync()
    MPI_Request nullReq = MPI_REQUEST_NULL;

    // divide the data into multiple messages
    qindex messageSize, numMessages;
    getMessageConfig(&messageSize, &numMessages, numElems);

    // asynchronously send the messages; pairRank receives the same ordering
    for (qindex m=0; m<numMessages; m++)
        MPI_Isend(&send[m*messageSize], messageSize, MPI_QCOMP, pairRank, NULL_TAG, MPI_COMM_WORLD, &nullReq);

#else
    error_commButEnvNotDistributed();
#endif
}


void receiveArray(qcomp* dest, qindex numElems, int pairRank) {
#if ENABLE_DISTRIBUTION

    // expect the data in multiple messages
    qindex messageSize, numMessages;
    getMessageConfig(&messageSize, &numMessages, numElems);

    // create a request for each asynch receive below
    std::vector<MPI_Request> requests(numMessages, MPI_REQUEST_NULL);

    // listen to receive each message asynchronously (as per arxiv.org/abs/2308.07402)
    for (qindex m=0; m<numMessages; m++)
        MPI_Irecv(&dest[m*messageSize], messageSize, MPI_QCOMP, pairRank, NULL_TAG, MPI_COMM_WORLD, &requests[m]);

    // receivers wait for all messages to be received (while sender asynch proceeds)
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

#else
    error_commButEnvNotDistributed();
#endif
}



/*
 * PRIVATE GPU COMMUNICATION
 */


void exchangeGpuAmpsToGpuBuffers(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps, int pairRank) {

    // exchange GPU memory directly if possible
    if (gpu_isDirectGpuCommPossible()) {

        // ensure GPU is finished modifying gpuAmps and gpuCommBuffer
        gpu_sync();

        // communicate via GPUDirect or Peer-to-Peer
        exchangeArrays(&qureg.gpuAmps[sendInd], &qureg.gpuCommBuffer[recvInd], numAmps, pairRank);

    // otherwise route the memory through the CPU
    } else {

        // copy GPU memory (amps) into CPU memory (amps), beginning from 0
        gpu_copyGpuToCpu(qureg, &qureg.gpuAmps[sendInd], &qureg.cpuAmps[0], numAmps);

        // exchange CPU memory (amps) to other node's CPU memory (buffer), beginning from 0
        exchangeArrays(qureg.cpuAmps, qureg.cpuCommBuffer, numAmps, pairRank);

        // copy CPU memory (buffer) to GPU memory (buffer), beginning from recvInd
        gpu_copyCpuToGpu(qureg, &qureg.cpuCommBuffer[0], &qureg.gpuCommBuffer[recvInd], numAmps);
    }
}


void exchangeGpuBuffers(Qureg qureg, qindex numAmpsAndRecvInd, int pairRank) {

    // exchange GPU memory directly if possible
    if (gpu_isDirectGpuCommPossible()) {

        // ensure GPU is finished modifying gpuCommBuffer
        gpu_sync();

        // communicate via GPUDirect or Peer-to-Peer
        exchangeArrays(&qureg.gpuCommBuffer[0], &qureg.gpuCommBuffer[numAmpsAndRecvInd], numAmpsAndRecvInd, pairRank);
    
    // otherwise route the memory through the CPU
    } else {

        // copy GPU memory (buffer) into CPU memory (amps), beginning from 0
        gpu_copyGpuToCpu(qureg, &qureg.gpuCommBuffer[0], &qureg.cpuAmps[0], numAmpsAndRecvInd);

        // exchange CPU memory (amps) to other node's CPU memory (buffer), beginning from 0
        exchangeArrays(qureg.cpuAmps, qureg.cpuCommBuffer, numAmpsAndRecvInd, pairRank);

        // copy CPU memory (buffer) to GPU memory (buffer), beginning from numAmpsAndRecvInd
        gpu_copyCpuToGpu(qureg, &qureg.cpuCommBuffer[0], &qureg.gpuCommBuffer[numAmpsAndRecvInd], numAmpsAndRecvInd);
    }
}


void asynchSendGpuBuffer(Qureg qureg, qindex numElems, int pairRank) {

    // send GPU memory directly if possible
    if (gpu_isDirectGpuCommPossible()) {

        // ensure GPU is finished modifying gpuCommBuffer
        gpu_sync();

        // communicate via GPUDirect or Peer-to-Peer
        asynchSendArray(qureg.gpuCommBuffer, numElems, pairRank);

    // otherwise route the memory through the CPU
    } else {

        // copy GPU memory (buffer) into CPU memory (amps)
        gpu_copyGpuToCpu(qureg, qureg.gpuCommBuffer, qureg.cpuAmps, numElems);

        // send CPU memory (amps) to other node
        asynchSendArray(qureg.cpuAmps, numElems, pairRank);
    }
}


void receiveArrayToGpuBuffer(Qureg qureg, qindex numElems, int pairRank) {

    // receive to GPU memory directly if possible
    if (gpu_isDirectGpuCommPossible()) {

        // GPU synchronisation is not necessary; we're merely receiving to buffer

        // communicate via GPUDirect or Peer-to-Peer
        receiveArray(qureg.gpuCommBuffer, numElems, pairRank);

    // otherwise, route through CPU
    } else {

        // receive array to CPU memory (buffer)
        receiveArray(qureg.cpuCommBuffer, numElems, pairRank);

        // copy CPU memory (buffer) to GPU memory (buffer)
        gpu_copyCpuToGpu(qureg, qureg.cpuCommBuffer, qureg.gpuCommBuffer, numElems);
    }
}



/*
 * PUBLIC STATE EXCHANGE METHODS
 */


void comm_exchangeAmpsToBuffers(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps, int pairRank) {
    assert_validCommBounds(qureg, sendInd, recvInd, numAmps);
    assert_quregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        exchangeGpuAmpsToGpuBuffers(qureg, sendInd, recvInd, numAmps, pairRank);
    else 
        exchangeArrays(&qureg.cpuAmps[sendInd], &qureg.cpuCommBuffer[recvInd], numAmps, pairRank);
}


void comm_exchangeBuffers(Qureg qureg, qindex numAmpsAndRecvInd, int pairRank) {
    assert_validCommBounds(qureg, 0, numAmpsAndRecvInd, numAmpsAndRecvInd);
    assert_quregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        exchangeGpuBuffers(qureg, numAmpsAndRecvInd, pairRank);
    else 
        exchangeArrays(&qureg.cpuCommBuffer[0], &qureg.cpuCommBuffer[numAmpsAndRecvInd], numAmpsAndRecvInd, pairRank);
}


void comm_asynchSendBuffer(Qureg qureg, qindex numElems, int pairRank) {
    assert_validCommBounds(qureg, 0, 0, numElems);
    assert_quregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        asynchSendGpuBuffer(qureg, numElems, pairRank);
    else
        asynchSendArray(qureg.cpuCommBuffer, numElems, pairRank);
}


void comm_receiveArrayToBuffer(Qureg qureg, qindex numElems, int pairRank) {
    assert_validCommBounds(qureg, 0, 0, numElems);
    assert_quregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        receiveArrayToGpuBuffer(qureg, numElems, pairRank);
    else
        receiveArray(qureg.cpuCommBuffer, numElems, pairRank);
}



/*
 * PUBLIC CONVENIENCE OVERLOADS
 */


void comm_exchangeAmpsToBuffers(Qureg qureg, int pairRank) {
    assert_quregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    comm_exchangeAmpsToBuffers(qureg, 0, 0, qureg.numAmpsPerNode, pairRank);
}



/*
 * PUBLIC REDUCTION METHODS
 */


void comm_reduceAmp(qcomp* localAmp) {
#if ENABLE_DISTRIBUTION

    qcomp* globalAmp;
    MPI_Allreduce(localAmp, globalAmp, 1, MPI_QCOMP, MPI_SUM, MPI_COMM_WORLD);
    *localAmp = *globalAmp;

#else
    error_commButEnvNotDistributed();
#endif
}


bool comm_isTrueOnAllNodes(bool val) {
#if ENABLE_DISTRIBUTION

    // perform global AND and broadcast result back to all nodes
    int local = (int) val;
    int global;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return (bool) global;

#else
    error_commButEnvNotDistributed();
#endif
}