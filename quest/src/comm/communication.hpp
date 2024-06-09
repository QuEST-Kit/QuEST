/** @file
 * Functions for communicating and exchanging amplitudes between compute
 * nodes, when running in distributed mode, using the C MPI standard.
 */

#ifndef COMMUNICATION_HPP
#define COMMUNICATION_HPP



/*
 * MPI ENVIRONMENT MANAGEMENT
 */

bool comm_isMpiCompiled();

bool comm_isMpiGpuAware();

bool comm_isInit();

void comm_init();

void comm_end();

int comm_getRank();

int comm_getNumNodes();

void comm_synch();



#endif // COMMUNICATION_HPP