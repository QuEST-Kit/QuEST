/** @file
 * Signatures for communicating and exchanging amplitudes between compute
 * nodes, when running in distributed mode, using the C MPI standard.
 * Calling these functions when QUEST_COMPILE_MPI=0, or when the passed Quregs
 * are not distributed, will throw a runtime internal error. 
 * 
 * @author Tyson Jones
 */

#ifndef COMM_ROUTINES_HPP
#define COMM_ROUTINES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"

#include <vector>
#include <string>

using std::vector;



/*
 * STATE EXCHANGE METHODS
 */

void comm_exchangeAmpsToBuffers(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps, int pairRank);

void comm_exchangeAmpsToBuffers(Qureg qureg, int pairRank);

void comm_exchangeSubBuffers(Qureg qureg, qindex numAmpsAndRecvInd, int pairRank);

// one disjoint sub-buffer chunk to exchange with a single pair rank, used by comm_exchangeSubBufferChunks
struct CommChunk {
    qindex sendInd; // buffer index where this chunk's amps to send begin
    qindex recvInd; // buffer index where this chunk's received amps are written
    qindex numAmps; // number of amps exchanged (a power of two)
    int pairRank;   // the partner rank for this chunk (distinct from this node)
};

void comm_exchangeSubBufferChunks(Qureg qureg, const vector<CommChunk>& chunks);

void comm_asynchSendSubBuffer(Qureg qureg, qindex numElems, int pairRank);

void comm_receiveArrayToBuffer(Qureg qureg, qindex numElems, int pairRank);

void comm_combineAmpsIntoBuffer(Qureg receiver, Qureg sender);

void comm_combineElemsIntoBuffer(Qureg receiver, FullStateDiagMatr sender);



/*
 * MISC COMMUNICATION METHODS
 */

void comm_broadcastAmp(int sendRank, qcomp* sendAmp);

void comm_sendAmpsToRoot(int sendRank, qcomp* send, qcomp* recv, qindex numAmps);

void comm_broadcastIntsFromRoot(int* arr, qindex length);

void comm_broadcastUnsignedsFromRoot(unsigned* arr, qindex length);

void comm_combineSubArrays(qcomp* recv, vector<qindex> globalRecvInds, vector<qindex> localSendInds, vector<qindex> numAmpsPerRank);



/*
 * REDUCTION METHODS
 */

void comm_reduceAmp(qcomp* localAmp);

void comm_reduceReal(qreal* localReal);

void comm_reduceReals(qreal* localReals, qindex numLocalReals);

bool comm_isTrueOnAllNodes(bool val);

bool comm_isTrueOnRootNode(bool val);



/*
 * GATHER METHODS
 */

vector<std::string> comm_gatherStringsToRoot(char* localChars, int maxNumLocalChars);



#endif // COMM_ROUTINES_HPP