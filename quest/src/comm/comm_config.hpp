/** @file
 * Functions for querying the distributed configuration
 * using the MPI interface, agnostically to the specific
 * implementation (like OpenMPI vs MPICH). These functions
 * are callable even when MPI has not been compiled/linked.
 * 
 * @author Tyson Jones
 */

#ifndef COMM_CONFIG_HPP
#define COMM_CONFIG_HPP

constexpr int ROOT_RANK = 0;

bool comm_isMpiCompiled();
bool comm_isMpiSubCommunicatorCompiled();
bool comm_isMpiGpuAware();

void comm_init(bool userOwnsMpi);
void comm_end(bool userOwnsMpi);
void comm_sync();

int comm_getRank();
int comm_getNumNodes();

bool comm_isInit();
bool comm_isRootNode();
bool comm_isRootNode(int rank);

// Signatures containing MPI types which callers must extern:
// extern MPI_Comm comm_getMpiComm()
// extern void comm_setMpiComm(MPI_Comm newComm)

#endif // COMM_CONFIG_HPP
