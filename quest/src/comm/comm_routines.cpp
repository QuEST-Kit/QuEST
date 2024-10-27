/** @file
 * Functions for communicating and exchanging amplitudes between compute
 * nodes, when running in distributed mode, using the C MPI standard.
 */

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_indices.hpp"

#include <vector>

#if COMPILE_MPI
    #include <mpi.h>
#endif



/*
 * TODO:
 *
 * - create an MPI error handler for custom comm internal error messages:
 *   https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPI_Comm_create_errhandler.3.html#mpi-comm-create-errhandler
 * 
 * - we may wish to adjust communication design when supporting MPI+cuQuantum,
 *   wherein we must define our own Communicator type.
 *   https://github.com/NVIDIA/cuQuantum/blob/788c04862241c52a5435982e6c25b30b6e3e9324/samples/custatevec/samples_mpi/distributedIndexBitSwap.cpp#L177
 *
 * - apparently we can perform better than UCX's intra-socket inter-GPU
 *   communication, achieving 1.5-3x speedups. See this thesis:
 *   https://www.queensu.ca/academia/afsahi/pprl/thesis/Yiltan_Temucin_MASc_thesis.pdf
 * 
 * - we may wish to think more explicitly about memory affinity, using hwloc:
 *   https://www.open-mpi.org/projects/hwloc/
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
 * MPI COMPLEX TYPE FLAG
 */


#if COMPILE_MPI

    #if (FLOAT_PRECISION == 1)
        #define MPI_QCOMP MPI_CXX_FLOAT_COMPLEX

    #elif (FLOAT_PRECISION == 2)
        #define MPI_QCOMP MPI_CXX_DOUBLE_COMPLEX

    // sometimes 'MPI_CXX_LONG_DOUBLE_COMPLEX' isn't defined
    #elif (FLOAT_PRECISION == 4) && defined(MPI_CXX_LONG_DOUBLE_COMPLEX)
        #define MPI_QCOMP MPI_CXX_LONG_DOUBLE_COMPLEX

    // in that case, fall back to the C type (identical memory layout)
    #else
        #define MPI_QCOMP MPI_C_LONG_DOUBLE_COMPLEX
    #endif

#endif



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
#if COMPILE_MPI

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

    // wait for all exchanges to complete (MPI will automatically free the request memory)
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

#else
    error_commButEnvNotDistributed();
#endif
}



/*
 * PRIVATE ASYNC SEND AND RECEIVE
 */


void asynchSendArray(qcomp* send, qindex numElems, int pairRank) {
#if COMPILE_MPI

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
#if COMPILE_MPI

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
 * PRIVATE GLOBAL COMBINATION
 */


void globallyCombineSubArrays(qcomp* recv, qcomp* send, qindex numAmpsPerRank) {
#if COMPILE_MPI

    int numNodes = comm_getNumNodes();
    int myRank = comm_getRank();
    qindex myOffset = myRank * numAmpsPerRank;

    // every node first copies their 'send' portion into a distinct part of their local 'recv',
    // which they will subsequently broadcast to the other nodes
    cpu_copyArray(&recv[myOffset], send, numAmpsPerRank);

    // determine how many broadcasts are needed per sending rank due to message-size restrictions
    // (it is almost definitely one, but we'll be rigorous for defensive design)
    qindex messageSize, numMessagesPerSender;
    getMessageConfig(&messageSize, &numMessagesPerSender, numAmpsPerRank);

    // all broadcasts will be asynch, and each involve one request per message
    qindex totalNumMessages = numMessagesPerSender * numNodes;
    std::vector<MPI_Request> requests(totalNumMessages, MPI_REQUEST_NULL);

    // each node broadcasts their partition (in-turn, but each is asynch)
    for (int senderRank=0; senderRank<numNodes; senderRank++) {
        qindex senderOffset = senderRank * numAmpsPerRank;

        // divided into many asynchronous messages
        for (int m=0; m<numMessagesPerSender; m++) {
            qindex messageOffset = senderOffset + (m * messageSize);
            qindex requestIndex = m + (senderRank * numMessagesPerSender);

            MPI_Ibcast(&recv[messageOffset], messageSize, MPI_QCOMP, senderRank, MPI_COMM_WORLD, &requests[requestIndex]);
        }
    }

    // wait for all boradcasts to complete (MPI will automatically free the request memory)
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


void exchangeGpuSubBuffers(Qureg qureg, qindex numAmps, int pairRank) {

    auto [sendInd, recvInd] = getSubBufferSendRecvInds(qureg);

    // exchange GPU memory directly if possible
    if (gpu_isDirectGpuCommPossible()) {

        // ensure GPU is finished modifying gpuCommBuffer
        gpu_sync();

        // communicate via GPUDirect or Peer-to-Peer
        exchangeArrays(&qureg.gpuCommBuffer[sendInd], &qureg.gpuCommBuffer[recvInd], numAmps, pairRank);
    
    // otherwise route the memory through the CPU
    } else {

        // copy GPU memory (buffer) into CPU memory (amps), preserving offset
        gpu_copyGpuToCpu(qureg, &qureg.gpuCommBuffer[sendInd], &qureg.cpuAmps[sendInd], numAmps);

        // exchange CPU memory (amps) to other node's CPU memory (buffer), receiving at index 0
        exchangeArrays(&qureg.cpuAmps[sendInd], &qureg.cpuCommBuffer[recvInd], numAmps, pairRank);

        // copy CPU memory (buffer) to GPU memory (buffer), receiving at index 0
        gpu_copyCpuToGpu(qureg, &qureg.cpuCommBuffer[recvInd], &qureg.gpuCommBuffer[recvInd], numAmps);
    }
}


void asynchSendGpuSubBuffer(Qureg qureg, qindex numElems, int pairRank) {

    qindex sendInd = getSubBufferSendInd(qureg);

    // send GPU memory directly if possible
    if (gpu_isDirectGpuCommPossible()) {

        // ensure GPU is finished modifying gpuCommBuffer
        gpu_sync();

        // communicate via GPUDirect or Peer-to-Peer
        asynchSendArray(&qureg.gpuCommBuffer[sendInd], numElems, pairRank);

    // otherwise route the memory through the CPU
    } else {

        // copy GPU memory (buffer) into CPU memory (amps), at offset
        gpu_copyGpuToCpu(qureg, &qureg.gpuCommBuffer[sendInd], &qureg.cpuAmps[sendInd], numElems);

        // send CPU memory (amps) to other node
        asynchSendArray(&qureg.cpuAmps[sendInd], numElems, pairRank);
    }
}


void receiveArrayToGpuBuffer(Qureg qureg, qindex numElems, int pairRank) {

    qindex recvInd = getBufferRecvInd();

    // receive to GPU memory directly if possible
    if (gpu_isDirectGpuCommPossible()) {

        // GPU synchronisation is not necessary; we're merely receiving to buffer

        // communicate via GPUDirect or Peer-to-Peer
        receiveArray(&qureg.gpuCommBuffer[recvInd], numElems, pairRank);

    // otherwise, route through CPU
    } else {

        // receive array to CPU memory (buffer), at offset
        receiveArray(&qureg.cpuCommBuffer[recvInd], numElems, pairRank);

        // copy CPU memory (buffer) to GPU memory (buffer), at offset
        gpu_copyCpuToGpu(qureg, &qureg.cpuCommBuffer[recvInd], &qureg.gpuCommBuffer[recvInd], numElems);
    }
}



/*
 * PUBLIC STATE EXCHANGE METHODS
 */


void comm_exchangeAmpsToBuffers(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps, int pairRank) {
    assert_validCommBounds(qureg, sendInd, recvInd, numAmps);
    assert_commQuregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        exchangeGpuAmpsToGpuBuffers(qureg, sendInd, recvInd, numAmps, pairRank);
    else 
        exchangeArrays(&qureg.cpuAmps[sendInd], &qureg.cpuCommBuffer[recvInd], numAmps, pairRank);
}


void comm_exchangeSubBuffers(Qureg qureg, qindex numAmps, int pairRank) {
    
    auto [sendInd, recvInd] = getSubBufferSendRecvInds(qureg);

    assert_validCommBounds(qureg, sendInd, recvInd, numAmps);
    assert_bufferSendRecvDoesNotOverlap(sendInd, recvInd, numAmps);
    assert_commQuregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        exchangeGpuSubBuffers(qureg, numAmps, pairRank);
    else 
        exchangeArrays(&qureg.cpuCommBuffer[sendInd], &qureg.cpuCommBuffer[recvInd], numAmps, pairRank);
}


void comm_asynchSendSubBuffer(Qureg qureg, qindex numElems, int pairRank) {

    auto [sendInd, recvInd] = getSubBufferSendRecvInds(qureg);

    assert_validCommBounds(qureg, sendInd, recvInd, numElems);
    assert_bufferSendRecvDoesNotOverlap(sendInd, recvInd, numElems);
    assert_commQuregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        asynchSendGpuSubBuffer(qureg, numElems, pairRank);
    else
        asynchSendArray(&qureg.cpuCommBuffer[sendInd], numElems, pairRank);
}


void comm_receiveArrayToBuffer(Qureg qureg, qindex numElems, int pairRank) {

    auto [sendInd, recvInd] = getSubBufferSendRecvInds(qureg);

    assert_validCommBounds(qureg, sendInd, recvInd, numElems);
    assert_bufferSendRecvDoesNotOverlap(sendInd, recvInd, numElems);
    assert_commQuregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        receiveArrayToGpuBuffer(qureg, numElems, pairRank);
    else
        receiveArray(&qureg.cpuCommBuffer[recvInd], numElems, pairRank);
}


void comm_combineAmpsIntoBuffer(Qureg receiver, Qureg sender) {
    assert_commQuregIsDistributed(receiver);
    assert_commQuregIsDistributed(sender);
    assert_receiverCanFitSendersEntireState(receiver, sender);

    // all configurations involve broadcasting the entirety of sender's per-node amps
    qindex numSendAmps = sender.numAmpsPerNode;
    qindex numRecvAmps = sender.numAmps;

    // note that CUDA-aware MPI permits direct GPU-to-GPU (device-to-device) exchange,
    // but does not generally permit CPU-to-GPU (host-to-device). So if only one
    // Qureg is GPU-accelerated, we have to fall back entirely to copying through host.
    // There is ergo only a single scenario possible when we can directly GPU-exchange:
    if (gpu_isDirectGpuCommPossible() && receiver.isGpuAccelerated && sender.isGpuAccelerated) {
        globallyCombineSubArrays(receiver.gpuCommBuffer, sender.gpuAmps, numSendAmps);
        return;
    }

    // otherwise, we must always transmit amps through CPU buffer memory, and merely
    // have to decide whether CPU-GPU pre- and post-copies are necessary
    if (sender.isGpuAccelerated)
        gpu_copyGpuToCpu(sender, sender.gpuAmps, sender.cpuCommBuffer, numSendAmps);

    globallyCombineSubArrays(receiver.cpuCommBuffer, sender.cpuCommBuffer, numSendAmps);

    if (receiver.isGpuAccelerated)
        gpu_copyCpuToGpu(receiver, receiver.cpuCommBuffer, receiver.gpuCommBuffer, numRecvAmps);
}


void comm_combineElemsIntoBuffer(Qureg receiver, FullStateDiagMatr sender) {
    assert_commQuregIsDistributed(receiver);
    assert_commFullStateDiagMatrIsDistributed(sender);
    assert_receiverCanFitSendersEntireElems(receiver, sender);

    // all configurations involve broadcasting the entirety of sender's per-node amps
    qindex numSendAmps = sender.numElemsPerNode;
    qindex numRecvAmps = sender.numElems;

    // like in comm_combineAmpsIntoBuffer(), direct-GPU comm only possible if both ptrs are GPU
    if (gpu_isDirectGpuCommPossible() && receiver.isGpuAccelerated && sender.isGpuAccelerated) {
        globallyCombineSubArrays(receiver.gpuCommBuffer, sender.gpuElems, numSendAmps);
        return;
    }

    // even if sender is GPU-accelerated, we safely assume its CPU elements are unchanged from GPU
    globallyCombineSubArrays(receiver.cpuCommBuffer, sender.cpuElems, numSendAmps);

    if (receiver.isGpuAccelerated)
        gpu_copyCpuToGpu(receiver, receiver.cpuCommBuffer, receiver.gpuCommBuffer, numRecvAmps);
}



/*
 * PUBLIC CONVENIENCE OVERLOADS
 */


void comm_exchangeAmpsToBuffers(Qureg qureg, int pairRank) {
    assert_commQuregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    qindex sendInd = 0;
    qindex recvInd = 0;
    comm_exchangeAmpsToBuffers(qureg, sendInd, recvInd, qureg.numAmpsPerNode, pairRank);
}



/*
 * PUBLIC MISC COMMUNICATION METHODS
 */


void comm_sendAmpsToRoot(int sendRank, qcomp* send, qcomp* recv, qindex numAmps) {
#if COMPILE_MPI

    // only the sender and root nodes need to continue
    int recvRank = 0;
    int rank = comm_getRank();
    if (rank != sendRank && rank != recvRank)
        return;

    // create an MPI_Request for every asynch MPI call
    qindex messageSize, numMessages;
    getMessageConfig(&messageSize, &numMessages, numAmps);
    std::vector<MPI_Request> requests(numMessages, MPI_REQUEST_NULL);

    // asynchronously copy 'send' in sendRank over to 'recv' in recvRank
    for (qindex m=0; m<numMessages; m++)
        if (rank == sendRank)
            MPI_Isend(&send[m*messageSize], messageSize, MPI_QCOMP, recvRank, NULL_TAG, MPI_COMM_WORLD, &requests[m]);
        else
            MPI_Irecv(&recv[m*messageSize], messageSize, MPI_QCOMP, sendRank, NULL_TAG, MPI_COMM_WORLD, &requests[m]);

    // wait for all exchanges to complete (MPI will automatically free the request memory)
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

#else
    error_commButEnvNotDistributed();
#endif
}


void comm_broadcastUnsignedsFromRoot(unsigned* arr, qindex length) {
#if COMPILE_MPI

    int sendRank = 0;
    MPI_Bcast(arr, length, MPI_UNSIGNED, sendRank, MPI_COMM_WORLD);

#else
    error_commButEnvNotDistributed();
#endif
}



/*
 * PUBLIC REDUCTION METHODS
 */


void comm_reduceAmp(qcomp* localAmp) {
#if COMPILE_MPI

    qcomp* globalAmp = nullptr;
    MPI_Allreduce(localAmp, globalAmp, 1, MPI_QCOMP, MPI_SUM, MPI_COMM_WORLD);
    *localAmp = *globalAmp;

#else
    error_commButEnvNotDistributed();
#endif
}


bool comm_isTrueOnAllNodes(bool val) {
#if COMPILE_MPI

    // perform global AND and broadcast result back to all nodes
    int local = (int) val;
    int global;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return (bool) global;

#else
    error_commButEnvNotDistributed();
    return false;
#endif
}
