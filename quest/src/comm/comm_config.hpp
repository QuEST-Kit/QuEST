/** @file
 * Functions for querying the distributed configuration
 * using the MPI interface, attemptedly agnostically to
 * the implementation (like OpenMPI vs MPICH).
 */

#ifndef COMM_CONFIG_HPP
#define COMM_CONFIG_HPP


constexpr int ROOT_RANK = 0;


bool comm_isMpiCompiled();

bool comm_isMpiGpuAware();

bool comm_isInit();

void comm_init();

void comm_end();

int comm_getRank();

bool comm_isRootNode();
bool comm_isRootNode(int rank);

int comm_getNumNodes();

void comm_sync();



#endif // COMM_CONFIG_HPP