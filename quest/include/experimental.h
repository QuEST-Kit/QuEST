/** @file
 * Experimental functions which are liable to
 * API breaks within QuEST minor version releases.
 * Some optional functions require compiling this
 * file against MPI, despite being outside of /comm/, 
 * and so require opt-in macros (QUEST_COMPILE_SUBCOMM)
 * 
 * @author Oliver Brown
 * @author Tyson Jones (formatting)
 * 
 * @defgroup experimental Experimental
 * @ingroup api
 * @brief Experimental functions with tentative APIs
 * @{
 */

#ifndef EXPERIMENTAL_H
#define EXPERIMENTAL_H

#include "quest/include/config.h"

#if QUEST_COMPILE_SUBCOMM && ! QUEST_COMPILE_MPI
    #error "Macro QUEST_COMPILE_SUBCOMM was true, but QUEST_COMPILE_MPI was illegally false."
#endif

#if QUEST_COMPILE_SUBCOMM
    #include <mpi.h>
#endif

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 *
 *  Advanced initialiser which lets the user positively declare that they take responsibility for MPI.
 *  This means we assume they have called MPI_Init, and that they will call MPI_Finalize.
 * 
 * @author Oliver Brown
 */
void initCustomMpiQuESTEnv(int useDistrib, bool userOwnsMpi, int useGpuAccel, int useMultithread);


#if QUEST_COMPILE_SUBCOMM
/** @notyetdoced
 * 
 *  Advanced initialiser which allows the user to provide an MPI communicator for QuEST to use.
 *  Use of this initialiser implies userOwnsMpi = true, (exposed by initCustomMpiQuESTEnv) and 
 *  therefore that they have already initialised MPI, and they will call MPI_Finalize at the 
 *  appropriate time.
 *
 *  The user-provided MPI communicator undergoes the same validation procedure as any that QuEST
 *  would use, and so must contain a power-of-2 number of processes.
 * 
 * This function is only compiled and exposed when macro QUEST_COMPILE_SUBCOMM is 1, as is
 * defined when providing CMake option QUEST_ENABLE_SUBCOMM during building.
 *
 * @author Oliver Brown
 */
void initCustomMpiCommQuESTEnv(MPI_Comm questComm, int useGpuAccel, int useMultithread);
#endif // QUEST_COMPILE_SUBCOMM


/** @notyetdoced
 * 
 * @author Oliver Brown
 */
int getQuESTNumGpuThreadsPerBlock();


/** Overrides the number of CUDA threads per block (or @p blockDim) used by QuEST's GPU-accelerated backend.
 * 
 * This changes the GPU parallelisation granularity and can affect performance, and is useful
 * for performance tuning or diagnostics. Before this function is called, QuEST will use the
 * number as specified by the environment variable @p QUEST_DEFAULT_NUM_GPU_THREADS_PER_BLOCK,
 * if defined. Otherwise, it will use the value specified by the CMake/compile option of the
 * same name, which itself presently defaults to @p 128. After this function is called, QuEST
 * will adopt @p numThreadsPerBlock for the remainder of execution, or until this function is
 * called again.
 * 
 * Practical values of @p numThreadsPerBlock can vary with the simulation size, the user's GPU hardware,
 * and whether it is NVIDIA or AMD, which have respective warp sizes of @p 32 and @p 64.
 * 
 * @note
 * This function has no effect when QuEST is not deployed with GPU-acceleration enabled.
 *
 * @param[in] numThreadsPerBlock the new block size.
 * @throws @validationerror
 * - if the @p QuESTEnv is not initialised.
 * - if @p numThreadsPerBlock is negative.
 * - if @p numThreadsPerBlock is not a multiple of the GPU warp size.
 * - if @p numThreadsPerBlock exceeds the maximum @p blockDim imposed by the GPU hardware.
 * @see
 * - QUEST_DEFAULT_NUM_GPU_THREADS_PER_BLOCK
 * @author Oliver Brown
 * @author Tyson Jones
 */
void setQuESTNumGpuThreadsPerBlock(int numThreadsPerBlock);


// end de-mangler
#ifdef __cplusplus
}
#endif

#endif // EXPERIMENTAL_H

/** @} */ // (end file-wide doxygen defgroup)
