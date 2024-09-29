/** @file
 * Functions for communicating and exchanging amplitudes between compute
 * nodes, when running in distributed mode, using the C MPI standard.
 */

#ifndef COMM_ROUTINES_HPP
#define COMM_ROUTINES_HPP

#include "quest/include/types.h"
#include "quest/include/qureg.h"



/*
 * MESSAGE INDICES
 */

qindex getSubBufferSendInd(Qureg qureg);

qindex getSubBufferRecvInd();



/*
 * STATE EXCHANGE METHODS
 */

void comm_exchangeAmpsToBuffers(Qureg qureg, qindex sendInd, qindex recvInd, qindex numAmps, int pairRank);

void comm_exchangeAmpsToBuffers(Qureg qureg, int pairRank);

void comm_exchangeSubBuffers(Qureg qureg, qindex numAmpsAndRecvInd, int pairRank);

void comm_asynchSendSubBuffer(Qureg qureg, qindex numElems, int pairRank);

void comm_receiveArrayToBuffer(Qureg qureg, qindex numElems, int pairRank);



/*
 * MISC COMMUNICATION METHODS
 */

void comm_sendAmpsToRoot(int sendRank, qcomp* send, qcomp* recv, qindex numAmps);

void comm_broadcastUnsignedsFromRoot(unsigned* arr, qindex length);



/*
 * REDUCTION METHODS
 */

void comm_reduceAmp(qcomp* localAmp);

bool comm_isTrueOnAllNodes(bool val);



#endif // COMM_ROUTINES_HPP