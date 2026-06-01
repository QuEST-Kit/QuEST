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

// queries of MPI's global/general status (when visible)
bool comm_isMpiCompiled();
bool comm_isMpiSubCommCompiled();
bool comm_isMpiGpuAware();
bool comm_isMpiInit();

// control of QuEST's (possibly more limited) MPI env
bool comm_isActive();
void comm_init(bool userOwnsMpi);
void comm_end();
void comm_sync();

// queries of QuEST's (possibly more limited) MPI env
int comm_getRank();
int comm_getNumNodes();
bool comm_isRootNode();
bool comm_isRootNode(int rank);

// Signatures containing MPI types which callers must extern:
// extern MPI_Comm comm_getMpiComm()
// extern bool comm_setMpiComm(MPI_Comm newComm, bool userOwnsMpi)

#endif // COMM_CONFIG_HPP
