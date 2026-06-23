/** @file
 * Experimental functions which are liable to
 * API breaks within QuEST minor version releases.
 * Some optional functions require compiling this
 * file against MPI, despite being outside of /comm/, 
 * and so require opt-in macros (QUEST_COMPILE_SUBCOMM)
 * 
 * @author Oliver Brown
 * @author Tyson Jones (formatting)
 * @author Ashmit JaiSarita Gupta (checkpointing)
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

#include "quest/include/qureg.h"

// C++ gets string overloads
#ifdef __cplusplus
    #include <string>
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
 * > [!IMPORTANT]
 * > This function is only compiled and exposed when macro QUEST_COMPILE_SUBCOMM is 1, as is
 * > defined when providing CMake option QUEST_ENABLE_SUBCOMM during building.
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


/** Writes the contents of @p qureg to the file (or folder) @p fn, so that it may later be
 * restored with createQuregFromFile(), potentially in another process.
 * 
 * The output records the @p qureg dimension (number of qubits and whether it is a density matrix),
 * the amplitude precision, the Qureg's distribution, and the Qureg's full set of amplitudes. Other
 * deployment information, such as whether the Qureg is multithreaded or GPU-accelerated, is not
 * recorded. 
 * 
 * There is no particular file extension or folder name suffix required, though since saving is
 * performed with ADIOS2, a suffix of `.bp` is conventional.
 * 
 * > [!CAUTION]
 * > Specifying @p fn equal to an existing directory or file will cause erasure and overwriting of
 * > its contents. It is especially dangerous to pass @p fn equal to a system directory, such as
 * > @c / on Unix, and may cause system corruption. 
 * 
 * > [!IMPORTANT]
 * > This function is only callable when QuEST is compiled with CMake option @c QUEST_ENABLE_ADIOS2=1.
 *
 * @param[in] qureg the Qureg to write to disk.
 * @param[in] fn    the output file (or folder) path.
 * @throws @validationerror
 * - if @p qureg is uninitialised.
 * - if QuEST was not compiled with CMake option @c QUEST_ENABLE_ADIOS2=1.
 * - if opening or writing to @p fn fails.
 * @see
 * - createQuregFromFile() to restore a Qureg saved by this function.
 * @author Ashmit JaiSarita Gupta
 */
void saveQuregToFile(Qureg qureg, const char* fn);


/** Creates a new Qureg from a file (or folder) previously created by saveQuregToFile(),
 * with automatically chosen deployments (independent of those used when the
 * file was saved), and populates the Qureg with the saved amplitudes.
 * 
 * The chosen deployments are identical to those chosen by createQureg() and createDensityQureg().
 * 
 * > [!NOTE]
 * > The number of distributed nodes chosen by the autodeployer must agree with the
 * > number of nodes of the originally saved Qureg, else a @validationerror is thrown. Therefore,
 * > the number of MPI processes calling these functions cannot be changed between saveQuregToFile()
 * > and createQuregFromFile(), unless the Qureg was non-distributed in both settings.
 * 
 * > [!IMPORTANT]
 * > This function is only callable when QuEST is compiled with CMake option @c QUEST_ENABLE_ADIOS2=1.
 *
 * @param[in] fn the file (or folder) path previously created by saveQuregToFile().
 * @returns A new Qureg instance matching the saved dimension and amplitudes.
 * @throws @validationerror
 * - if QuEST was not compiled with CMake option @c QUEST_ENABLE_ADIOS2=1.
 * - if @p fn cannot be read (since, for example, it does not exist).
 * - if the precision of the saved Qureg differs from the current QuEST precision.
 * - if the number of distributed nodes of the saved Qureg differs from the autodeployer's chosen number.
 * - if the recorded Qureg dimensions would overflow the @c qindex type.
 * - if the recorded total Qureg memory would overflow the @c size_t type.
 * - if the system contains insufficient RAM (or VRAM) to store the Qureg in any deployment.
 * - if any Qureg memory allocation unexpectedly fails.
 * @see
 * - saveQuregToFile() to create a file readable by this function.
 * @author Ashmit JaiSarita Gupta
 * @author Tyson Jones (input validation)
 */
Qureg createQuregFromFile(const char* fn);


// end de-mangler
#ifdef __cplusplus
}
#endif



#if defined(__cplusplus)


    /** 
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - saveQuregToFile()
     */
    void saveQuregToFile(Qureg qureg, std::string);


    /** 
     * @notyetdoced
     * @cpponly
     * 
     * @see
     * - createQuregFromFile()
     */
    Qureg createQuregFromFile(std::string fn);


#endif // __cplusplus


#endif // EXPERIMENTAL_H

/** @} */ // (end file-wide doxygen defgroup)
