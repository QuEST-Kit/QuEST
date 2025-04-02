/** @file
 * Functions for communicating and exchanging amplitudes between compute
 * nodes, when running in distributed mode, using the C MPI standard.
 * Calling these functions when COMPILE_MPI=0, or when the passed Quregs
 * are not distributed, will throw a runtime internal error. 
 * 
 * @author Tyson Jones
 * @author Jakub Adamski (sped-up large comm by asynch messages)
 * @author Oliver Brown (patched max-message inference, consulted on AR and MPICH support)
 * @author Ania (Anna) Brown (developed QuEST v1 logic)
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

#if COMPILE_MPI
    #include <mpi.h>
#endif

#include <vector>
#include <array>
#include <algorithm>

using std::vector;


/**
 * @todo
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
 * for UCX exchange, so 28 (to support quad precision) seems safe.
 * While 2^28 fits in 'int', we use 'qindex' so that arithmetic
 * never overflows, and we cast down to 'int' when safe
 */

qindex MAX_MESSAGE_LENGTH = powerOf2(28);



/*
 * MPI COMPLEX TYPE FLAG
 */


#if COMPILE_MPI

    // declare MPI types for qreal and qcomp. We always use the 
    // C macros, even when the deprecated CXX equivalents are 
    // available, to maintain compatibility with modern MPICH
    #if   (FLOAT_PRECISION == 1)
        #define MPI_QREAL MPI_FLOAT
        #define MPI_QCOMP MPI_C_FLOAT_COMPLEX
    #elif (FLOAT_PRECISION == 2)
        #define MPI_QREAL MPI_DOUBLE
        #define MPI_QCOMP MPI_C_DOUBLE_COMPLEX
    #elif (FLOAT_PRECISION == 4)
        #define MPI_QREAL MPI_LONG_DOUBLE
        #define MPI_QCOMP MPI_C_LONG_DOUBLE_COMPLEX
    #else
        #error "Something went horribly wrong in inferring the MPI types"
    #endif

#endif



/*
 * PRIVATE MESSAGE SIZES
 */


int getMaxNumMessages() {
#if COMPILE_MPI

    // the max supported tag value constrains the total number of messages 
    // we can send in a round of communication, since we uniquely tag
    // each message in a round such that we do not rely upon message-order 
    // gaurantees and ergo can safely support UCX adaptive routing (AR).
    // The MPI standard necessitates the max tag is always at least...
    int minTagUpperBound = 32767;

    // but we pedantically consult MPI in case we ever need to send MORE (smaller?)
    // messages. Beware the max is obtained via a void pointer and might be unset...
    void* tagUpperBoundPtr;
    int isAttribSet;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tagUpperBoundPtr, &isAttribSet);

    // if something went wrong with obtaining the tag bound, return the safe minimum
    if (!isAttribSet)
        return minTagUpperBound;

    // otherwise return whichever is bigger of the found bound and the minimum
    // (which is really just hiding an error; it should ALWAYS be that UB>=min)
    int tagUpperBound = *(int*) tagUpperBoundPtr;
    return std::max({minTagUpperBound, tagUpperBound});

#else
    error_commButEnvNotDistributed();
    return -1;
#endif
}


std::array<qindex,2> dividePow2PayloadIntoMessages(qindex numAmps) {
    assert_commPayloadIsPowerOf2(numAmps);

    // use single message if possible
    if (numAmps < MAX_MESSAGE_LENGTH)
        return {numAmps, 1};
    
    // else, payload divides evenly between max-size messages (always fits in int)
    qindex numMessages = numAmps / MAX_MESSAGE_LENGTH;

    // which we must be able to uniquely tag
    if (numMessages > getMaxNumMessages())
        error_commNumMessagesExceedTagMax();

    // outputs always fit in 'int' but we force them to be 'qindex' since 
    // caller will multiply them with ints and could easily overflow
    return {MAX_MESSAGE_LENGTH, numMessages};
}


std::array<qindex,3> dividePayloadIntoMessages(qindex numAmps) {

    // use single message if possible
    if (numAmps < MAX_MESSAGE_LENGTH)
        return {numAmps, 1, 0};

    // else, use as many max-size messages as possible, and one smaller msg
    qindex numMaxSizeMsgs   = numAmps / MAX_MESSAGE_LENGTH; // floors
    qindex remainingMsgSize = numAmps - numMaxSizeMsgs * MAX_MESSAGE_LENGTH;

    // all of which we must be able to uniquely tag
    if (numMaxSizeMsgs + 1 > getMaxNumMessages())
        error_commNumMessagesExceedTagMax();

    // outputs always fit in 'int' but we force them to be 'qindex' since 
    // caller will multiply them with ints and could easily overflow
    return {MAX_MESSAGE_LENGTH, numMaxSizeMsgs, remainingMsgSize};
}



/*
 * PRIVATE SYNCHRONOUSE AMPLITUDE EXCHANGE
 */


void exchangeArrays(qcomp* send, qcomp* recv, qindex numElems, int pairRank) {
#if COMPILE_MPI

    // each message is asynchronously dispatched with a final wait, as per arxiv.org/abs/2308.07402

    // we will send payload in multiple asynch messages (create two requests per msg for subsequent synch)
    auto [messageSize, numMessages] = dividePow2PayloadIntoMessages(numElems);
    vector<MPI_Request> requests(2*numMessages, MPI_REQUEST_NULL);

    // asynchronously exchange the messages (effecting MPI_Isendrecv), using unique tags 
    // so that messages are permitted to arrive out-of-order (supporting UCX adaptive-routing)
    for (qindex m=0; m<numMessages; m++) {
        int tag = static_cast<int>(m); // gauranteed int, but m*messageSize needs qindex
        MPI_Isend(&send[m*messageSize], messageSize, MPI_QCOMP, pairRank, tag, MPI_COMM_WORLD, &requests[2*m]);
        MPI_Irecv(&recv[m*messageSize], messageSize, MPI_QCOMP, pairRank, tag, MPI_COMM_WORLD, &requests[2*m+1]);
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
    auto [messageSize, numMessages] = dividePow2PayloadIntoMessages(numElems);

    // asynchronously send the uniquely-tagged messages
    for (qindex m=0; m<numMessages; m++) {
        int tag = static_cast<int>(m); // gauranteed int, but m*messageSize needs qindex
        MPI_Isend(&send[m*messageSize], messageSize, MPI_QCOMP, pairRank, tag, MPI_COMM_WORLD, &nullReq);
    }

#else
    error_commButEnvNotDistributed();
#endif
}


void receiveArray(qcomp* dest, qindex numElems, int pairRank) {
#if COMPILE_MPI

    // expect the data in multiple messages
    auto [messageSize, numMessages] = dividePow2PayloadIntoMessages(numElems);

    // create a request for each asynch receive below
    vector<MPI_Request> requests(numMessages, MPI_REQUEST_NULL);

    // listen to receive each uniquely-tagged message asynchronously (as per arxiv.org/abs/2308.07402)
    for (qindex m=0; m<numMessages; m++) {
        int tag = static_cast<int>(m); // gauranteed int, but m*messageSize needs qindex
        MPI_Irecv(&dest[m*messageSize], messageSize, MPI_QCOMP, pairRank, tag, MPI_COMM_WORLD, &requests[m]);
    }

    // receivers wait for all messages to be received (while sender asynch proceeds)
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

#else
    error_commButEnvNotDistributed();
#endif
}



/*
 * PRIVATE GLOBAL COMBINATION
 */


void globallyCombineNonUniformSubArrays(
    qcomp* recv, qcomp* send,
    vector<qindex> globalRecvIndPerRank, vector<qindex> localSendIndPerRank, vector<qindex> numSendPerRank,
    bool areGpuPtrs
) {
#if COMPILE_MPI

    int myRank = comm_getRank();
    int numNodes = comm_getNumNodes();

    if (globalRecvIndPerRank.size() != (size_t) numNodes)
        error_commGivenInconsistentNumSubArraysANodes();

    // every node first copies their 'send' portion into a distinct part of their local 'recv',
    // which they will subsequently broadcast to the other nodes, if it is non-zero in size;
    // if send is nullptr, then the caller already prepared the relevant portion of recv
    if (send != nullptr) {
        auto func = (areGpuPtrs)? gpu_copyArray : cpu_copyArray;
        func(
            &recv[globalRecvIndPerRank[myRank]], 
            &send[localSendIndPerRank[myRank]], 
            numSendPerRank[myRank]); // may be zero
    }

    // all node-broadcasts will be asynch, and each involves one request per sent message,
    // but unlikely in other routines, their payloads can differ significantly in size,
    // so we do not know the total number of requests needed in advance
    vector<MPI_Request> requests;

    // each node broadcasts their partition (in-turn, but each is asynch)...
    for (int sendRank=0; sendRank<numNodes; sendRank++) {

        // potentially using multiple messages, due to message-size restrictions
        // (they almost definitely send only one, but we divide for defensive design)
        auto [bigMsgSize, numBigMsgs, remMsgSize] = dividePayloadIntoMessages(numSendPerRank[sendRank]);

        // these involve big asynch messages (could be 'numSend', or the max message size)
        for (int m=0; m<numBigMsgs; m++) {
            qindex recvInd = globalRecvIndPerRank[sendRank] + (m * bigMsgSize);
            requests.push_back(MPI_REQUEST_NULL);
            MPI_Ibcast(&recv[recvInd], bigMsgSize, MPI_QCOMP, sendRank, MPI_COMM_WORLD, &requests.back());
        }

        // and potentially one remaining asynch message 
        if (remMsgSize > 0) {
            qindex recvInd = globalRecvIndPerRank[sendRank] + (numBigMsgs * bigMsgSize);
            requests.push_back(MPI_REQUEST_NULL);
            MPI_Ibcast(&recv[recvInd], remMsgSize, MPI_QCOMP, sendRank, MPI_COMM_WORLD, &requests.back());
        }
    }

    // wait for all broadcasts to complete (MPI will automatically free the request memory)
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

#else
    error_commButEnvNotDistributed();
#endif
}


void globallyCombineSubArrays(qcomp* recv, qcomp* send, qindex numAmpsPerRank, bool areGpuPtrs) {
#if COMPILE_MPI

    // simply wrap and call the non-uniform case has no performance penalty, 
    // and is only slightly messier than a bespoke power-of-2 msg implementation

    int numNodes = comm_getNumNodes();

    vector<qindex> recvInds(numNodes);
    vector<qindex> sendInds(numNodes);
    vector<qindex> numAmps(numNodes);

    for (int r=0; r<numNodes; r++) {
        recvInds[r] = r * numAmpsPerRank;
        sendInds[r] = 0;
        numAmps[r] = numAmpsPerRank;
    }

    globallyCombineNonUniformSubArrays(recv, send, recvInds, sendInds, numAmps, areGpuPtrs);

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
    assert_commBoundsAreValid(qureg, sendInd, recvInd, numAmps);
    assert_commQuregIsDistributed(qureg);
    assert_pairRankIsDistinct(qureg, pairRank);

    if (qureg.isGpuAccelerated)
        exchangeGpuAmpsToGpuBuffers(qureg, sendInd, recvInd, numAmps, pairRank);
    else 
        exchangeArrays(&qureg.cpuAmps[sendInd], &qureg.cpuCommBuffer[recvInd], numAmps, pairRank);
}


void comm_exchangeSubBuffers(Qureg qureg, qindex numAmps, int pairRank) {
    
    auto [sendInd, recvInd] = getSubBufferSendRecvInds(qureg);

    assert_commBoundsAreValid(qureg, sendInd, recvInd, numAmps);
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

    assert_commBoundsAreValid(qureg, sendInd, recvInd, numElems);
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

    assert_commBoundsAreValid(qureg, sendInd, recvInd, numElems);
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
    if (receiver.isGpuAccelerated && sender.isGpuAccelerated && gpu_isDirectGpuCommPossible()) {
        gpu_sync();
        globallyCombineSubArrays(receiver.gpuCommBuffer, sender.gpuAmps, numSendAmps, true);
        return;
    }

    // otherwise, we must always transmit amps through CPU memory (NOT buffer), and 
    // merely have to decide whether CPU-GPU pre- and post-copies are necessary
    if (sender.isGpuAccelerated)
        gpu_copyGpuToCpu(sender, sender.gpuAmps, sender.cpuAmps, numSendAmps);

    globallyCombineSubArrays(receiver.cpuCommBuffer, sender.cpuAmps, numSendAmps, false);

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
    if (receiver.isGpuAccelerated && sender.isGpuAccelerated && gpu_isDirectGpuCommPossible() ) {
        gpu_sync();
        globallyCombineSubArrays(receiver.gpuCommBuffer, sender.gpuElems, numSendAmps, true);
        return;
    }

    // there is no need to copy sender's elems from GPU to CPU; they should never diverge

    // even if sender is GPU-accelerated, we safely assume its CPU elements are unchanged from GPU
    globallyCombineSubArrays(receiver.cpuCommBuffer, sender.cpuElems, numSendAmps, false);

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


void comm_broadcastAmp(int sendRank, qcomp* sendAmp) {
#if COMPILE_MPI

    MPI_Bcast(sendAmp, 1, MPI_QCOMP, sendRank, MPI_COMM_WORLD);

#else
    error_commButEnvNotDistributed();
#endif
}


void comm_sendAmpsToRoot(int sendRank, qcomp* send, qcomp* recv, qindex numAmps) {
#if COMPILE_MPI

    // only the sender and root nodes need to continue
    int recvRank = ROOT_RANK;
    int myRank = comm_getRank();
    if (myRank != sendRank && myRank != recvRank)
        return;

    // create an MPI_Request for every asynch MPI call
    auto [messageSize, numMessages] = dividePow2PayloadIntoMessages(numAmps);
    vector<MPI_Request> requests(numMessages, MPI_REQUEST_NULL);

    // asynchronously copy 'send' in sendRank over to 'recv' in recvRank, using 
    // uniquely-tagged messages such that they may arrive out-of-order, enabling AR
    for (qindex m=0; m<numMessages; m++) {
        int tag = static_cast<int>(m);
        (myRank == sendRank)?
            MPI_Isend(&send[m*messageSize], messageSize, MPI_QCOMP, recvRank, tag, MPI_COMM_WORLD, &requests[m]): // sender
            MPI_Irecv(&recv[m*messageSize], messageSize, MPI_QCOMP, sendRank, tag, MPI_COMM_WORLD, &requests[m]); // root
    }

    // wait for all exchanges to complete (MPI will automatically free the request memory)
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

#else
    error_commButEnvNotDistributed();
#endif
}


void comm_broadcastIntsFromRoot(int* arr, qindex length) {
#if COMPILE_MPI

    int sendRank = ROOT_RANK;
    MPI_Bcast(arr, length, MPI_INT, sendRank, MPI_COMM_WORLD);

#else
    error_commButEnvNotDistributed();
#endif
}


void comm_broadcastUnsignedsFromRoot(unsigned* arr, qindex length) {
#if COMPILE_MPI

    int sendRank = ROOT_RANK;
    MPI_Bcast(arr, length, MPI_UNSIGNED, sendRank, MPI_COMM_WORLD);

#else
    error_commButEnvNotDistributed();
#endif
}


void comm_combineSubArrays(qcomp* recv, vector<qindex> recvInds, vector<qindex> sendInds, vector<qindex> numSend) {

    // recv has already been overwritten with local contributions, which enables
    // direct GPU-to-recv writing (avoiding a superfluous CPU-CPU copy)
    qcomp* send = nullptr;
    bool areGpuPtrs = false;
    globallyCombineNonUniformSubArrays(recv, send, recvInds, sendInds, numSend, areGpuPtrs);
}



/*
 * PUBLIC REDUCTION METHODS
 */


void comm_reduceAmp(qcomp* localAmp) {
#if COMPILE_MPI

    MPI_Allreduce(MPI_IN_PLACE, localAmp, 1, MPI_QCOMP, MPI_SUM, MPI_COMM_WORLD);

#else
    error_commButEnvNotDistributed();
#endif
}


void comm_reduceReal(qreal* localReal) {
#if COMPILE_MPI

    MPI_Allreduce(MPI_IN_PLACE, localReal, 1, MPI_QREAL, MPI_SUM, MPI_COMM_WORLD);

#else
    error_commButEnvNotDistributed();
#endif
}


void comm_reduceReals(qreal* localReals, qindex numLocalReals) {
#if COMPILE_MPI

    MPI_Allreduce(MPI_IN_PLACE, localReals, numLocalReals, MPI_QREAL, MPI_SUM, MPI_COMM_WORLD);

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


bool comm_isTrueOnRootNode(bool val) {
    #if COMPILE_MPI

    // this isn't really a reduction - it's a broadcast - but
    // it's semantically relevant to comm_isTrueOnAllNodes()

    unsigned out = (unsigned) val;
    comm_broadcastUnsignedsFromRoot(&out, 1);
    return (bool) out;

#else
    error_commButEnvNotDistributed();
    return false;
#endif
}



/*
 * PUBLIC GATHER METHODS
 */


vector<string> comm_gatherStringsToRoot(char* localChars, int maxNumLocalChars) {
#if COMPILE_MPI

    // no need to validate array sizes and memory alloc successes;
    // these are trivial O(#nodes)-size arrays containing <20 chars
    int numNodes = comm_getNumNodes();

    // root makes one big contiguous array to receive all chars contiguously
    vector<char> allChars(maxNumLocalChars * numNodes);

    // all nodes send root all their local chars
    int recvRank = ROOT_RANK;
    MPI_Gather(localChars, maxNumLocalChars, MPI_CHAR, allChars.data(),
        maxNumLocalChars, MPI_CHAR, recvRank, MPI_COMM_WORLD);

    // divide allChars into stings, delimited by each node's terminal char
    vector<string> out(numNodes);
    for (int r=0; r<numNodes; r++)
        out[r] = std::string(&allChars[r*maxNumLocalChars]);

    return out;

#else
    error_commButEnvNotDistributed();
    return {};
#endif
}
