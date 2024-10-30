/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/precision.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/paulis.h"
#include "quest/include/channels.h"

#include "quest/src/core/validation.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/parser.hpp"
#include "quest/src/core/printer.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/cpu/cpu_config.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>

using std::string;
using std::vector;



/*
 * INVALID INPUT ERROR MESSAGES
 * which can contain variables with syntax ${VAR1} ${VAR2}, substituted at error-throw with
 * runtime parameters via assertThat(..., {{"${VAR1}",1}, {"${VAR2}",2}}, ...)
 */

namespace report {


    /*
     *  ENVIRONMENT CREATION
     */

    string QUEST_ENV_ALREADY_INIT =
        "The QuEST environment has already been initialised. This can only be performed once during program execution.";

    string QUEST_ENV_ALREADY_FINAL =
        "The QuEST environment has already been finalised, and can thereafter never be re-initialised since this leads to undefined MPI behaviour.";


    string INVALID_OPTION_FOR_ENV_IS_DISTRIB =
        "Argument 'useDistrib' must be 1 or 0 to respectively indicate whether or not to distribute the new environment, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    string INVALID_OPTION_FOR_ENV_IS_GPU_ACCEL =
        "Argument 'useGpuAccel' must be 1 or 0 to respectively indicate whether or not to GPU-accelerate the new environment, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    string INVALID_OPTION_FOR_ENV_IS_MULTITHREAD =
        "Argument 'useMultithread' must be 1 or 0 to respectively indicate whether or not to enable multithreading in the new environment, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";


    string CANNOT_MULTITHREAD_ENV_WITHOUT_OPENMP_COMPILATION =
        "Cannot create a multithreaded environment because the QuEST source was not compiled with OpenMP enabled.";

    string CANNOT_DISTRIB_ENV_WITHOUT_MPI_COMMPILATION =
        "Cannot create a distributed environment because the QuEST source was not compiled with MPI enabled.";

    string CANNOT_GPU_ACCEL_ENV_WITH_GPU_COMPILATION =
        "Cannot create a GPU-accelerated environment because the QuEST source was not compiled with GPU acceleration.";


    string CANNOT_GPU_ACCEL_ENV_WITH_NO_AVAILABLE_GPUS =
        "Cannot create a GPU-accelerated environment because there are no GPUs available.";

    string CANNOT_DISTRIB_ENV_BETWEEN_NON_POW_2_NODES =
        "Cannot distribute QuEST between ${NUM_NODES} nodes; must use a power-of-2 number of nodes.";

    string CUQUANTUM_DEPLOYED_ON_GPU_WITHOUT_MEM_POOLS =
        "Cannot use cuQuantum since your GPU does not support memory pools. Please recompile with cuQuantum disabled to fall-back to using Thrust and custom kernels.";

    
    /*
     * EXISTING QUESTENV
     */

    string QUEST_ENV_NOT_INIT =
        "The QuEST environment is not initialised. Please first call initQuESTEnv() or initCustomQuESTEnv().";


    /*
     * DEBUG UTILITIES
     */

    string INVALID_NEW_EPSILON =
        "Invalid new epsilon value (${NEW_EPS}). Must specifiy a positive number, or zero to disable numerical validation (as if the epsilon is infinity).";

    string INVALID_NUM_REPORTED_ITEMS =
        "Invalid parameter (${NUM_ITEMS}). Must specify a positive number of items to be reported, or 0 to indicate that all items should be reported.";


    /*
     * QUREG CREATION
     */

    string NON_POSITIVE_NUM_QUBITS_IN_CREATE_QUREG =
        "Cannot create Qureg of ${NUM_QUBITS} qubits; must contain one or more qubits.";


    string NEW_STATEVEC_QUREG_NUM_AMPS_WOULD_EXCEED_QINDEX = 
        "Cannot create Qureg of ${NUM_QUBITS} qubits: the statevector would contain more amplitudes (2^${NUM_QUBITS}) than the maximum which can be addressed by the qindex type (2^${MAX_QUBITS}). See reportQuESTEnv().";

    string NEW_DENSMATR_QUREG_NUM_AMPS_WOULD_EXCEED_QINDEX = 
        "Cannot create density Qureg of ${NUM_QUBITS} qubits: the density matrix would contain more amplitudes (4^${NUM_QUBITS}) than can be addressed by the qindex type (4^${MAX_QUBITS}). See reportQuESTEnv().";

    string NEW_STATEVEC_QUREG_LOCAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create Qureg of ${NUM_QUBITS} qubits distributed over ${NUM_NODES} nodes because the necessary local memory (in bytes) of each node would overflow size_t. In this deployment, the maximum number of qubits in a statevector Qureg is ${MAX_QUBITS}. See reportQuESTEnv().";

    string NEW_DENSMATR_QUREG_LOCAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create density Qureg of ${NUM_QUBITS} qubits distributed over ${NUM_NODES} nodes because the necessary local memory (in bytes) of each node would overflow size_t. In this deployment, the maximum number of qubits in a density-matrix Qureg is ${MAX_QUBITS}. See reportQuESTEnv().";

    string NEW_DISTRIB_STATEVEC_QUREG_HAS_TOO_FEW_AMPS = 
        "Cannot create a distributed Qureg of only ${NUM_QUBITS} qubits between ${NUM_NODES} nodes, because each node would contain fewer than one amplitude of the statevector. The minimum size is ${MIN_QUBITS} qubits; see reportQuESTEnv(). Consider disabling distribution for this Qureg.";

    string NEW_DISTRIB_DENSMATR_QUREG_HAS_TOO_FEW_AMPS = 
        "Cannot create a distributed density Qureg of only ${NUM_QUBITS} qubits between ${NUM_NODES} nodes, because each node would contain fewer than a column's worth of amplitudes of the density matrix. The minimum size is ${MIN_QUBITS} qubits; see reportQuESTEnv(). Consider disabling distribution for this Qureg.";


    string INVALID_OPTION_FOR_QUREG_IS_DENSMATR = 
        "Argument 'isDensityMatrix' must be 1 or 0 to respectively indicate whether the Qureg should be instantiated as a potentially-mixed density matrix or a strictly-pure state-vector.";

    string INVALID_OPTION_FOR_QUREG_IS_DISTRIB = 
        "Argument 'useDistrib' must be 1 or 0 to respectively indicate whether or not to distribute the new Qureg, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    string INVALID_OPTION_FOR_QUREG_IS_GPU_ACCEL = 
        "Argument 'useGpuAccel' must be 1 or 0 to respetively indicate whether or not to GPU-accelerate the new Qureg, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    string INVALID_OPTION_FOR_QUREG_IS_MULTITHREAD = 
        "Argument 'useMultithread' must be 1 or 0 to respectively indicate whether or not to use multithreading when processing the new Qureg, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";


    string NEW_DISTRIB_QUREG_IN_NON_DISTRIB_ENV =
        "Cannot distribute a Qureg when in a non-distributed QuEST environment.";

    string NEW_GPU_QUREG_IN_NON_GPU_ACCEL_ENV =
        "Cannot allocate a Qureg to a GPU when in a non-GPU-accelerated QuEST environment.";

    string NEW_MULTITHREAD_QUREG_IN_NON_MULTITHREAD_ENV =
        "Cannot enable multithreaded processing of a Qureg created in a non-multithreaded QuEST environment.";


    string NEW_GPU_QUREG_CANNOT_USE_MULTITHREADING = 
        "Cannot simultaneously GPU-accelerate and multithread a Qureg. Please disable multithreading, or set it to ${AUTO_DEPLOYMENT_FLAG} for QuEST to automatically disable it when deploying to GPU.";


    string NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CPU_MEM =
        "The non-distributed Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits would be too large (${QCOMP_BYTES} * ${EXP_BASE}^${NUM_QUBITS} bytes) to fit into a single node's RAM (${RAM_SIZE} bytes). See reportQuESTEnv(), and consider using distribution.";

    string NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CPU_MEM = 
        "The distributed Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits, together with its commnication buffer, would be too large (2 * ${QCOMP_BYTES} * ${EXP_BASE}^${NUM_QUBITS} bytes) to fit into the combined RAM of all ${NUM_NODES} nodes (${NUM_NODES} * ${RAM_SIZE} bytes). See reportQuESTEnv().";

    string NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CURRENT_GPU_MEM =
        "The non-distributed GPU-accelerated Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits would be too large (${QCOMP_BYTES} * ${EXP_BASE}^${NUM_QUBITS} bytes) to fit into a single node's available GPU memory (currently ${MIN_VRAM_AVAIL} bytes free). Consider additionally using distribution, or disabling GPU-acceleration (though this may greatly increase runtime).";

    string NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CURRENT_GPU_MEM =
        "The distributed GPU-accelerated Qureg (isDensity=${IS_DENS}) of ${NUM_QUBITS} qubits, together with its commnication buffer, is too large; one or more of the ${NUM_GPUS} GPUs has insufficient available memory (only ${MIN_VRAM_AVAIL} bytes) to store its Qureg partition (${QCOMP_BYTES} * 2^${LOG2_NUM_AMPS} bytes) bytes. Consider disabling GPU-acceleration.";


    string NEW_QUREG_CPU_AMPS_ALLOC_FAILED = 
        "Allocation of memory to store the CPU amplitudes failed.";

    string NEW_QUREG_GPU_AMPS_ALLOC_FAILED = 
        "Allocation of memory to store the GPU amplitudes failed.";

    string NEW_QUREG_CPU_COMM_BUFFER_ALLOC_FAILED = 
        "Allocation of memory for the distributed CPU communication buffer failed.";

    string NEW_QUREG_GPU_COMM_BUFFER_ALLOC_FAILED = 
        "Allocation of memory for the distributed GPU communication buffer failed.";


    /*
     * EXISTING QUREG
     */

    string INVALID_QUREG_FIELDS = 
        "Invalid Qureg; invalid or incompatible fields isDensityMatrix=${DENS_MATR}, numQubits=${NUM_QUBITS}, numAmps=${NUM_AMPS}. It is likely this Qureg was not initialised with createQureg().";

    string QUREG_NOT_DENSITY_MATRIX =
        "Expected a density matrix Qureg but received a statevector.";

    string QUREG_NOT_STATE_VECTOR =
        "Expected a statevector Qureg but received a density matrix.";


    /*
     * MUTABLE OBJECT FLAGS
     */

    string NEW_HEAP_FLAG_ALLOC_FAILED =
        "Attempted allocation of a heap flag (such as 'isUnitary', 'isHermitian', 'isCPTP', 'wasGpuSynced') miraculously failed, despite being a mere ${NUM_BYTES} bytes. This is unfathomably unlikely - go and have your fortune read at once!";

    string INVALID_HEAP_FLAG_PTR =
        "A flag (e.g. 'isUnitary', 'isHermitian', 'isCPTP', 'wasGpuSynced') bound to the given operator data structure (e.g. matrix, superoperator, Kraus map) was a NULL pointer, instead of an expected pointer to persistent heap memory. This may imply the structure has already been destroyed and its fields manually overwritten by the user.";

    string INVALID_HEAP_FLAG_VALUE = 
        "A flag (e.g. 'isUnitary', 'isHermitian', 'isCPTP', 'wasGpuSynced') bound to the given operator data structure (e.g. matrix, superoperator, Kraus map) had an invalid value of ${BAD_FLAG}. Allowed values are '0', '1', and (except for 'wasGpuSynced') '${UNKNOWN_FLAG}' which indicates the flag is currently unknown because its evaluation has been deferred. However, these flags should never be modified directly by the user.";


    /*
     * MATRIX CREATION
     */

    string NEW_MATRIX_NUM_QUBITS_NOT_POSITIVE = 
        "Cannot create a matrix which acts upon ${NUM_QUBITS} qubits; must target one or more qubits.";


    string NEW_DIAG_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX =
        "Cannot create a diagonal matrix of ${NUM_QUBITS} qubits: the matrix would contain more elements (2^${NUM_QUBITS}) than the maximum which can be addressed by the qindex type (2^${MAX_QUBITS}).";

    string NEW_COMP_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX =
        "Cannot create a dense matrix of ${NUM_QUBITS} qubits: the matrix would contain more elements (4^${NUM_QUBITS}) than the maximum which can be addressed by the qindex type (4^${MAX_QUBITS}).";


    string NEW_LOCAL_COMP_MATR_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a local, dense matrix of ${NUM_QUBITS} qubits because the necessary memory (in bytes) would overflow size_t. In this deployment, the maximum number of qubits in such a matrix is ${MAX_QUBITS}.";

    string NEW_LOCAL_DIAG_MATR_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a local, diagonal matrix of ${NUM_QUBITS} qubits because the necessary memory (in bytes) would overflow size_t. In this deployment, the maximum number of qubits in such a matrix is ${MAX_QUBITS}.";

    string NEW_DISTRIB_DIAG_MATR_LOCAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a diagonal matrix of ${NUM_QUBITS} qubits distributed over ${NUM_NODES} nodes because the necessary local memory (in bytes) of each node would overflow size_t. In this deployment, the maximum number of qubits in such a matrix is ${MAX_QUBITS}.";


    string NEW_DISTRIB_DIAG_MATR_HAS_TOO_FEW_AMPS =
        "Cannot create a diagonal matrix of ${NUM_QUBITS} qubits distributed between ${NUM_NODES} nodes because each node would contain fewer than one element. The minimum number of qubits in such a matrix is ${MIN_QUBITS}. Consider disabling distribution for this matrix.";


    string NEW_LOCAL_COMP_MATR_CANNOT_FIT_INTO_CPU_MEM =
        "Cannot create a local, dense matrix of ${NUM_QUBITS} qubits because the necessary memory (${QCOMP_BYTES} * 4^${NUM_QUBITS} bytes) exceeds the available RAM of ${RAM_SIZE} bytes.";

    string NEW_LOCAL_DIAG_MATR_CANNOT_FIT_INTO_CPU_MEM =
        "Cannot create a local, diagonal matrix of ${NUM_QUBITS} qubits because the necessary memory (${QCOMP_BYTES} * 2^${NUM_QUBITS} bytes) exceeds the available RAM of ${RAM_SIZE} bytes.";

    string NEW_DISTRIB_DIAG_MATR_CANNOT_FIT_INTO_CPU_MEM =
        "Cannot create a diagonal matrix of ${NUM_QUBITS} qubits distributed between ${NUM_NODES} because the necessary memory per node (${QCOMP_BYTES} * 2^${NUM_QB_MINUS_LOG_NODES} bytes) exceeds the local available RAM of ${RAM_SIZE} bytes.";


    string NEW_LOCAL_COMP_MATR_CANNOT_FIT_INTO_GPU_MEM =
        "Cannot create a local, GPU-accelerated, dense matrix of ${NUM_QUBITS} qubits because the necessary memory (${QCOMP_BYTES} * 4^${NUM_QUBITS} bytes) exceeds the available GPU memory of ${VRAM_SIZE} bytes.";

    string NEW_LOCAL_DIAG_MATR_CANNOT_FIT_INTO_GPU_MEM =
        "Cannot create a local, GPU-accelerated, diagonal matrix of ${NUM_QUBITS} qubits because the necessary memory (${QCOMP_BYTES} * 2^${NUM_QUBITS} bytes) exceeds the available GPU memory of ${VRAM_SIZE} bytes.";

    string NEW_DISTRIB_DIAG_MATR_CANNOT_FIT_INTO_GPU_MEM =
        "Cannot create a GPU-accelerated, diagonal matrix of ${NUM_QUBITS} qubits distributed between ${NUM_NODES} because the necessary memory per node (${QCOMP_BYTES} * 2^${NUM_QB_MINUS_LOG_NODES} bytes) exceeds the local available GPU memory of ${VRAM_SIZE} bytes.";


    string NEW_MATRIX_CPU_ELEMS_ALLOC_FAILED = 
        "Attempted allocation of memory for one or more rows of the matrix (a total of ${NUM_BYTES} bytes in RAM) failed.";

    string NEW_MATRIX_GPU_ELEMS_ALLOC_FAILED = 
        "Attempted allocation of the matrix's GPU memory (${NUM_BYTES} bytes in VRAM) failed.";


    string NEW_DISTRIB_MATRIX_IN_NON_DISTRIB_ENV = 
        "Cannot distribute a matrix in a non-distributed environment.";

    string INVALID_OPTION_FOR_MATRIX_IS_DISTRIB = 
        "Argument 'useDistrib' must be 1 or 0 to respectively indicate whether or not to distribute the new matrix, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";



    /*
     * MATRIX INITIALISATION
     */

    string COMP_MATR_NEW_ELEMS_WRONG_NUM_ROWS =
        "Incompatible number of rows (${NUM_GIVEN_ROWS}) of elements given to overwrite a ${NUM_QUBITS}-qubit CompMatr, which expects ${NUM_EXPECTED_ROWS} rows.";

    string COMP_MATR_NEW_ELEMS_WRONG_ROW_DIM =
        "One or more rows contained an incompatible number of elements (${NUM_GIVEN_ELEMS}). The ${NUM_QUBITS}-qubit CompMatr expects a square ${EXPECTED_DIM}x${EXPECTED_DIM} matrix.";

    string DIAG_MATR_WRONG_NUM_NEW_ELEMS = 
        "Incompatible number of elements (${NUM_GIVEN_ELEMS}) assigned to a ${NUM_QUBITS}-qubit DiagMatr, which expects ${NUM_EXPECTED_ELEMS} elements.";
    

    string FULL_STATE_DIAG_MATR_NEW_ELEMS_INVALID_START_INDEX =
        "Invalid start index (${START_IND}), which must be non-negative and smaller than the total number of diagonal elements in the matrix (${MATR_NUM_ELEMS}).";

    string FULL_STATE_DIAG_MATR_NEW_ELEMS_NUM_IS_NON_POSITIVE =
        "Invalid number of new elements (${NUM_ELEMS}). Must be greater than zero.";

    string FULL_STATE_DIAG_MATR_NEW_ELEMS_NUM_EXCEEDS_MAX_NUM =
        "The given number of new elements (${NEW_NUM_ELEMS}) exceeds the total number of diagonal elements in the matrix (${MATR_NUM_ELEMS}).";

    string FULL_STATE_DIAG_MATR_NEW_ELEMS_EXCEEDS_END_INDEX =
        "The specified range of elements to set (at indices ${START_IND} to ${END_IND_EXCL} exclusive) exceeds the bounds of the diagonal matrix (of ${MATR_NUM_ELEMS} total elements).";


    string MATR_NUM_QUBITS_MISMATCHES_INLINE_SETTER =
        "The declared number of qubits (${NUM_SETTER_QUBITS}) differs from the number of qubits of the matrix (${NUM_MATRIX_QUBITS}).";

    string MATR_NUM_ELEMS_MISMATCHES_VEC_LENGTH_IN_INLINE_SETTER =
        "The declared number of passed elements (${NUM_ELEMS}) differs from the length of the given list (${VEC_LENGTH}).";


    string MULTI_VAR_FUNC_INVALID_NUM_VARS = 
        "Invalid number of variables or dimensions (${NUM_VARS}). Must be a positive integer.";

    string MULTI_VAR_FUNC_INVALID_NUM_QUBITS_PER_VAR =
        "The variable/dimension at index ${VAR_IND} was constituted by an invalid number of qubits (${VAR_QUBITS}). Each must correspond to 1 or more qubits.";

    string MULTI_VAR_FUNC_MISMATCHING_NUM_QUBITS =
        "The total number of qubits constituting all variables/dimensions (${NUM_VAR_QUBITS}) does not match the number of qubits in the matrix (${NUM_MATR_QUBITS}).";

    string MULTI_VAR_FUNC_INVALID_ARE_SIGNED_FLAG =
        "Invalid value for the 'areSigned' flag (${ARE_SIGNED}), which must instead be 1 or 0 to indicate whether or not to interpret the variable sub-register basis states as signed integers (encoded with two's complement).";


    /*
     * EXISTING MATRIX
     */

    string INVALID_COMP_MATR_1_FIELDS =
        "Invalid CompMatr1. Targeted ${NUM_QUBITS} qubits (instead of 1) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 2x2). It is likely this matrix was not initialised with getCompMatr1().";

    string INVALID_COMP_MATR_2_FIELDS =
        "Invalid CompMatr2. Targeted ${NUM_QUBITS} qubits (instead of 2) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 4x4). It is likely this matrix was not initialised with getCompMatr2().";

    string INVALID_COMP_MATR_FIELDS =
        "Invalid CompMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createCompMatr().";


    string INVALID_DIAG_MATR_1_FIELDS =
        "Invalid DiagMatr1. Targeted ${NUM_QUBITS} qubits (instead of 1) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 2x2). It is likely this matrix was not initialised with getDiagMatr1().";

    string INVALID_DIAG_MATR_2_FIELDS =
        "Invalid DiagMatr2. Targeted ${NUM_QUBITS} qubits (instead of 2) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 4x4). It is likely this matrix was not initialised with getDiagMatr2().";

    string INVALID_DIAG_MATR_FIELDS =
        "Invalid DiagMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createDiagMatr().";


    string INVALID_FULL_STATE_DIAG_MATR_FIELDS = 
        "Invalid FullStateDiagMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createFullStateDiagMatr().";

    string FULL_STATE_DIAG_MATR_GPU_ACCEL_IN_NON_GPU_ENV =
        "Invalid FullStateDiagMatr. Was marked as GPU-accelerated in a non-GPU-accelerated environment, which is impossible. It is likely this matrix was manually modified and/or corrupted.";


    string COMP_MATR_NOT_SYNCED_TO_GPU = 
        "The CompMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please call syncCompMatr() after manually modifying elements, or overwrite all elements with setCompMatr() which automatically synchronises.";

    string DIAG_MATR_NOT_SYNCED_TO_GPU = 
        "The DiagMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please call syncDiagMatr() after manually modifying elements, or overwrite all elements with setDiagMatr() which automatically synchronises.";

    string FULL_STATE_DIAG_MATR_NOT_SYNCED_TO_GPU = 
        "The FullStateDiagMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please call syncFullStateDiagMatr() after manually modifying elements, or overwrite elements in batch with setFullStateDiagMatr() which automatically synchronises.";


    string MATRIX_NOT_UNITARY = 
        "The given matrix was not (approximately) unitary.";

    string MATRIX_NOT_HERMITIAN =
        "THe given matrix was not (approximately) Hermitian.";


    string INVALID_MATRIX_CPU_ELEMS_PTR =
        "The given matrix's outer CPU heap-memory pointer was NULL. This may imply the matrix was prior destroyed and its pointers were manually cleared.";

    string MATRIX_GPU_ELEMS_PTR_UNEXPECTEDLY_NULL =
        "The given matrix's GPU memory pointer was unexpectedly NULL. This may imply the matrix was prior destroyed and its pointer was manually cleared.";

    string MATRIX_GPU_ELEMS_PTR_UNEXPECTEDLY_NOT_NULL =
        "The given matrix's GPU memory pointer was non-NULL, implying an allocation, despite the QuEST environment being non-GPU-accelerated. This implies the matrix fields were manually modified and corrupted.";


    string FULL_STATE_DIAG_MATR_MISMATCHES_QUREG_DIM =
        "The given FullStateDiagMatr operates upon a different number of qubits (${NUM_MATR_QUBITS}) than exists in the Qureg (${NUM_QUREG_QUBITS}).";

    string FULL_STATE_DIAG_MATR_IS_DISTRIB_BUT_QUREG_ISNT =
        "The given FullStateDiagMatr is distributed but the Qureg is not, which is forbidden. Consider disabling distribution for this matrix via createCustomFullStateDiagMatr().";


    /*
     * SUPEROPERATOR CREATION
     */

    string NEW_SUPER_OP_NUM_QUBITS_NOT_POSITIVE =
        "Cannot create a superoperator of ${NUM_QUBITS} qubits. The operator must act upon one or more qubits.";


    string NEW_SUPER_OP_NUM_ELEMS_WOULD_EXCEED_QINDEX =
        "Cannot create a superoperator of ${NUM_QUBITS} qubits because the matrix would contain a total of 16^${NUM_QUBITS} elements which exceeds the maximum representable index of 'qindex' (permitting up to ${MAX_QUBITS} qubits).";
    
    string NEW_SUPER_OP_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a superoperator of ${NUM_QUBITS} qubits because the total required memory (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) exceeds the maximum representable by size_t. In this deployment, the maximum number of qubits in such a superoperator is ${MAX_QUBITS}";


    string NEW_SUPER_OP_CANNOT_FIT_INTO_CPU_MEM =
        "Cannot create a superoperator of ${NUM_QUBITS} qubits because the total memory required (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) exceeds that available (${RAM_SIZE} bytes).";

    string NEW_SUPER_OP_CANNOT_FIT_INTO_GPU_MEM =
        "Cannot create a GPU-accelerated superoperator of ${NUM_QUBITS} qubits because the total memory required (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) exceeds the available GPU memory (${VRAM_SIZE} bytes).";


    string NEW_SUPER_OP_CPU_ELEMS_ALLOC_FAILED =
        "Attempted allocation of memory for one or more rows of superoperator matrix (a total of ${NUM_BYTES} bytes in RAM) failed.";
        
    string NEW_SUPER_OP_GPU_ELEMS_ALLOC_FAILED =
        "Attempted allocation of GPU memory (a total of ${NUM_BYTES} in VRAM) for the superoperator matrix failed.";


    string NEW_INLINE_SUPER_OP_MATRIX_WRONG_NUM_ROWS =
        "The specified number of qubits (${NUM_DECLARED_QUBITS}) in the superoperator is inconsistent with the number of rows in the given matrix (expected ${NUM_DECLARED_ROWS}, received ${GIVEN_DIM}).";

    string NEW_INLINE_SUPER_OP_MATRIX_WRONG_NUM_COLS =
        "The specified number of qubits (${NUM_DECLARED_QUBITS}) in the superoperator is inconsistent with the number of columns in one or more rows of the given matrix (expected ${NUM_DECLARED_ROWS} columns, received ${GIVEN_DIM}).";


    /*
     * SUPEROPERATOR INITIALISATION
     */

    string SUPER_OP_NEW_MATRIX_ELEMS_WRONG_NUM_ROWS =
        "Incompatible number of rows (${GIVEN_DIM}) of elements given to overwrite a ${NUM_QUBITS}-qubit SuperOp, which expects ${EXPECTED_DIM} rows.";

    string SUPER_OP_NEW_MATRIX_ELEMS_WRONG_NUM_COLS =
        "Incompatible number of columns (${GIVEN_DIM}) in one or more rows of a matrix given to overwrite a ${NUM_QUBITS}-qubit SuperOp, which expects ${EXPECTED_DIM} rows.";

    string SUPER_OP_FIELDS_MISMATCH_PARAMS =
        "The specified number of qubits (${NUM_PASSED_QUBITS}) differs from the number in the SuperOp (${NUM_OP_QUBITS}).";
        

    /*
     * EXISTING SUPEROPERATOR
     */

    string INVALID_SUPER_OP_FIELDS =
        "The given SuperOp had invalid fields (${NUM_QUBITS} qubits and ${NUM_ROWS} rows), suggesting it was not initialised with createSuperOp().";

    string INVALID_SUPER_OP_CPU_MEM_PTR =
        "The given SuperOp's pointer to its CPU memory was NULL. This may imply the superoperator was already destroyed and had its memory pointers manually overwritten by the user.";

    string INVALID_SUPER_OP_GPU_MEM_PTR =
        "The given SuperOp's pointer to its GPU memory was NULL. This may imply the superoperator was already destroyed and had its memory pointers manually overwritten by the user.";

    string SUPER_OP_NOT_SYNCED_TO_GPU =
        "The SuperOp has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please call syncSuperOp() after manually modifying the superoperator elements, or overwrite all elements with setSuperOp() which will automatically syncrhonise.";


    /*
     * KRAUS MAP CREATION
     */

    string NEW_KRAUS_MAP_NUM_QUBITS_NOT_POSITIVE = 
        "Cannot create a Kraus map which acts upon ${NUM_QUBITS} qubits; must target one or more qubits.";


    string KRAUS_MAP_MATRICES_TOTAL_ELEMS_WOULD_EXCEED_QINDEX =
        "Cannot create a KrausMap with the given number of Kraus operators. The total number of elements between the ${NUM_MATRICES} given ${NUM_QUBITS}-qubit matrices would exceed the maximum number representing by the qindex type. At most, ${MAX_NUM_MATRICES} matrices of that size may be specified.";
        
    string KRAUS_MAP_MATRICES_TOTAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a KrausMap with the given number of Kraus operators. The total memory to store ${NUM_MATRICES} ${NUM_QUBITS}-qubit matrices exceeds that representable by the size_t type. At most, ${MAX_NUM_MATRICES} matrices of that size may be specified.";

    string NEW_KRAUS_MAPS_SUPER_OP_NUM_ELEMS_WOULD_EXCEED_QINDEX =
        "Cannot create a Kraus map of ${NUM_QUBITS} qubits because the corresponding superoperator would contain more elements (16^${NUM_QUBITS}) than the maximum which can be addressed by the qindex type (16^${MAX_QUBITS}).";

    string NEW_KRAUS_MAPS_SUPER_OP_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a Kraus map of ${NUM_QUBITS} qubits because the necessary memory for its corresponding superoperator (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) would overflow size_t. In this deployment, the maximum number of qubits in a Kraus map is ${MAX_QUBITS}.";


    string NEW_KRAUS_MAPS_SUPER_OP_CANNOT_FIT_INTO_CPU_MEM =
        "Cannot create a Kraus map of ${NUM_QUBITS} qubits because the total memory required by its corresponding superoperator (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) exceeds the available memory (${RAM_SIZE} bytes).";
        
    string NEW_KRAUS_MAPS_SUPER_OP_CANNOT_FIT_INTO_GPU_MEM =
        "Cannot create a GPU-accelerated Kraus map of ${NUM_QUBITS} qubits because the total memory required by its corresponding superoperator (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) exceeds the available GPU memory (${VRAM_SIZE} bytes).";


    string NEW_KRAUS_MAPS_SUPER_OP_CPU_ELEMS_ALLOC_FAILED =
        "Attempted allocation of memory for one or more rows of the Kraus map's correspondng superoperator matrix (a total of ${NUM_BYTES} bytes in RAM) failed.";
                
    string NEW_KRAUS_MAPS_SUPER_OP_GPU_ELEMS_ALLOC_FAILED =
        "Attempted allocation of GPU memory (a total of ${NUM_BYTES} bytes in VRAM) for the Kraus map's corresponding superoperator matrix failed.";

    string NEW_KRAUS_MAP_CPU_MATRICES_ALLOC_FAILED =
        "Failed to allocate memory (a total of ${NUM_BYTES} bytes) for the KrausMap's ${NUM_MATRICES} ${NUM_QUBITS}-qubit Kraus operator matrices.";


    string NEW_INLINE_KRAUS_MAP_INCOMPATIBLE_NUM_NEW_MATRICES =
        "Specified ${NUM_EXPECTED} Kraus operators, but passed ${NUM_GIVEN} matrices.";

    string NEW_INLINE_KRAUS_MAP_MATRIX_WRONG_NUM_ROWS =
        "Specified ${NUM_QUBITS} qubits, but one or more passed matrices have ${NUM_GIVEN_ROWS} rows instead of the expected ${NUM_EXPECTED_ROWS}.";

    string NEW_INLINE_KRAUS_MAP_MATRIX_WRONG_NUM_COLS =
        "The number of given columns (${NUM_GIVEN_COLS}) in one or more matrices was inconsistent with the specified number of qubits (${NUM_QUBITS}). Every Kraus operator must be specified as a ${NUM_EXPECTED_COLS}x${NUM_EXPECTED_COLS} matrix.";

    
    /*
     * KRAUS MAP INITIALISATION
     */

    string KRAUS_MAP_NUM_GIVEN_NEW_MATRICES_NOT_POSITIVE =
        "Cannot create a Kraus map from ${NUM_MATRICES} matrices; must be given a strictly positive number of matrices.";

    string KRAUS_MAP_INCOMPATIBLE_NUM_NEW_MATRICES =
        "Passed a number of matrices (${NUM_GIVEN}) which is inconsistent with the number of Kraus operators defined in the map (${NUM_EXPECTED}).";


    string KRAUS_MAP_NEW_MATRIX_ELEMS_WRONG_NUM_ROWS =
        "One or more of the given matrices had a matrix dimension (${NUM_GIVEN_ROWS} rows) incompatible with the given Kraus map of ${NUM_QUBITS} qubits, which expects matrices with ${NUM_EXPECTED_ROWS} rows.";

    string KRAUS_MAP_NEW_MATRIX_ELEMS_WRONG_ROW_DIM =
        "One or more of the given matrices had a matrix dimension (${NUM_GIVEN_COLS} columns) incompatible with the given Kraus map of ${NUM_QUBITS} qubits, which expects matrices with ${NUM_EXPECTED_COLS} columns.";

    string KRAUS_MAP_FIELDS_MISMATCH_PARAMS =
        "The specified number of Kraus operators (${NUM_PASSED_OPS}) and qubits (${NUM_PASSED_QUBITS}) differs from the number in the KrausMap (respectively ${NUM_MAP_OPS} and ${NUM_MAP_QUBITS}).";


    /*
     * EXISTING KRAUS MAP
     */

    string INVALID_KRAUS_MAP_FIELDS = 
        "Invalid KrausMap. One or more of its fields are invalid or inconsistent; it consists of ${NUM_MATRICES} Kraus operators, each upon ${NUM_QUBITS} qubits, with matrices of dimension ${NUM_ROWS}x${NUM_ROWS}. It is likely this Kraus map was not created with createKrausMap().";

    string INVALID_KRAUS_MAPS_SUPER_OP_FIELDS =
        "The given KrausMap's superoperator had invalid fields (${NUM_QUBITS} qubits and ${NUM_ROWS} rows), suggesting the map was not initialised with createKrausMap().";

    string INVALID_KRAUS_MAP_SUPER_OP_NUM_QUBITS =
        "Invalid KrausMap. The ${MAP_QUBITS}-qubit KrausMap's superoperator had an incompatible dimension of ${SUPER_OP_QUBITS}. This suggests the Kraus map was not created using createKrausMap(), or was not initialised with setKrausMap().";


    string INVALID_KRAUS_MAP_IS_CPTP_PTR = 
        "The 'isCPTP' field of the given KrausMap was a NULL pointer, instead of the expected pointer to persistent heap memory. This suggests the KrausMap was already destroyed and had its fields overwritten by the user.";

    string INVALID_KRAUS_MAP_IS_CPTP_FLAG = 
        "The 'isCPTP' field of the given Kraus map had invalid value ${BAD_FLAG}, suggesting it was manually modified. Valid values are 0, 1 and ${UNKNOWN_FLAG} (to indicate that unitarity is not yet known, deferring evaluation) although this need never be modified by the user.";


    string INVALID_KRAUS_MAPS_SUPER_OP_CPU_MEM_PTR =
        "The given KrausMap's superoperator had a NULL pointer where its CPU memory address was expected. This may imply the Kraus map was already destroyed and had its memory pointers manually overwritten by the user.";
        
    string INVALID_KRAUS_MAPS_SUPER_OP_GPU_MEM_PTR =
        "The given KrausMap's superoperator had a NULL pointer where its GPU memory address was expected. This may imply the Kraus map was already destroyed and had its memory pointers manually overwritten by the user.";

    string INVALID_KRAUS_MAP_MATRIX_LIST_MEM_PTR =
        "The KrausMap's list of Kraus operator matrices was unexpectedly a NULL pointer. This might imply the map was already destroyed and had its CPU memory pointers manually overwritten by the user.";


    string KRAUS_MAP_NOT_SYNCED_TO_GPU =
        "The KrausMap has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please call syncKrausMap() after manually modifying the superoperator elements, or overwrite all elements with setKrausMap() which will automatically syncrhonise.";

    string KRAUS_MAP_NOT_CPTP =
        "The KrausMap was not (approximately) completely positive and trace preserving (CPTP).";
    
    string KRAUS_MAP_SIZE_MISMATCHES_TARGETS =
        "The KrausMap has a size (${KRAUS_QUBITS} qubits) inconsistent with the specified number of target qubits (${TARG_QUBITS}).";
    

    /*
     * PAULI STRING CREATION
     */

    string NEW_PAULI_STR_NON_POSITIVE_NUM_PAULIS =
        "Invalid number of Paulis (${NUM_PAULIS}). The Pauli string must contain at least one Pauli operator, which is permittedly identity.";
    
    string NEW_PAULI_STR_NUM_PAULIS_EXCEEDS_TYPE =
        "Cannot make a Pauli string of ${NUM_PAULIS} Paulis since this exceeds the maximum of ${MAX_PAULIS} imposed by the typing of PauliStr's fields.";
    
    // TODO: replace BAD_CHAR ascii code with actual character, once tokenSubs is generalised to any-type
    string NEW_PAULI_STR_UNRECOGNISED_PAULI_CHAR = 
        "Given an unrecognised Pauli character (ASCII code ${BAD_CHAR}) at index ${CHAR_IND}. Each character must be one of I X Y Z (or lower case), or equivalently 0 1 2 3'.";

    string NEW_PAULI_STR_INVALID_PAULI_CODE =
        "Given an invalid Pauli code (${BAD_CODE}) at index ${CODE_IND}. Each code must be one of 0 1 2 3, corresponding respectively to I X Y Z.";

    string NEW_PAULI_STR_TERMINATION_CHAR_TOO_EARLY =
        "The given string contained fewer characters (${TERM_IND}) than the specified number of Pauli operators (${NUM_PAULIS}).";

    string NEW_PAULI_STR_INVALID_INDEX =
        "Invalid index (${BAD_IND}). Pauli indices must be non-negative and cannot equal nor exceed the maximum number of representable Pauli operators (${MAX_PAULIS}).";

    string NEW_PAULI_STR_DUPLICATED_INDEX =
        "The Pauli indices contained a duplicate. Indices must be unique.";

    string NEW_PAULI_STR_DIFFERENT_NUM_CHARS_AND_INDS = 
        "Given a different number of Pauli operators (${NUM_PAULIS}) and their qubit indices (${NUM_INDS}).";


    /*
     * PAULI STRING SUM CREATION
     */

    string NEW_PAULI_STR_SUM_NON_POSITIVE_NUM_STRINGS =
        "Cannot create a sum with ${NUM_TERMS} terms. The number of terms must be a positive integer.";

    string NEW_PAULI_STR_SUM_DIFFERENT_NUM_STRINGS_AND_COEFFS =
        "Given a different number of Pauli strings (${NUM_STRS}) and coefficients ${NUM_COEFFS}.";

    string NEW_PAULI_STR_SUM_STRINGS_ALLOC_FAILED = 
        "Attempted allocation of the PauliStrSum's ${NUM_TERMS}-term array of Pauli strings (${NUM_BYTES} bytes total) unexpectedly failed.";

    string NEW_PAULI_STR_SUM_COEFFS_ALLOC_FAILED =
        "Attempted allocation of the PauliStrSum's ${NUM_TERMS}-term array of coefficients (${NUM_BYTES} bytes total) unexpectedly failed.";


    /*
     * PAULI STRING SUM PARSING
     */

    string PARSED_PAULI_STR_SUM_UNINTERPRETABLE_LINE = 
        "Could not interpret line ${LINE_NUMBER} as a coefficient followed by a sequence of Pauli operators.";

    string PARSED_PAULI_STR_SUM_INCONSISTENT_NUM_PAULIS_IN_LINE =
        "Line ${LINE_NUMBER} specified ${NUM_LINE_PAULIS} Pauli operators which is inconsistent with the number of Paulis of the previous lines (${NUM_PAULIS}).";

    string PARSED_PAULI_STR_SUM_COEFF_IS_INVALID =
        "The coefficient of line ${LINE_NUMBER} could not be converted to a qcomp, possibly due to it exceeding the valid numerical range.";

    string PARSED_STRING_IS_EMPTY =
        "The given string was empty (contained only whitespace characters) and could not be parsed.";


    /*
     * EXISTING PAULI STRING SUM
     */

    string INVALID_PAULI_STR_SUM_FIELDS =
        "The given PauliStrSum had invalid fields (.numTerms=${NUM_TERMS}), which is likely a result from not being correctly initialised by createPauliStrSum().";

    string INVALID_PAULI_STR_HEAP_PTR =
        "One or more of the PauliStrumSum's heap pointers was unexpectedly NULL. This might imply the PauliStrSum was already destroyed and had its pointers manually overwritten by the user.";

    
    string PAULI_STR_SUM_NOT_HERMITIAN =
        "THe given PauliStrSum is not Hermitian.";


    /*
     * BASIS STATE INDICES
     */

    string INVALID_BASIS_STATE_INDEX = 
        "Classical state index ${STATE_IND} is invalid for the given ${NUM_QUBITS} qubit Qureg. Index must be greater than or equal to zero, and cannot equal nor exceed the number of unique classical states (2^${NUM_QUBITS} = ${NUM_STATES}).";

    string INVALID_BASIS_STATE_ROW_OR_COL =
        "The row and column indices (${ROW_IND}, ${COL_IND}) are invalid for the given ${NUM_QUBITS} qubit Qureg. Both indices must be greater than or equal to zero, and neither can equal nor exceed the number of unique classical states (2^${NUM_QUBITS} = ${NUM_STATES}).";


    /*
     * QUBIT INDICES
     */

    string NEGATIVE_OR_ZERO_NUM_TARGETS =
        "The specified number of target qubits (${NUM_QUBITS}) is invalid. Must specify one or more targets.";
    
    string NUM_TARGETS_EXCEEDS_QUREG_SIZE =
        "The specified number of target qubits (${NUM_QUBITS}) exceeds the number of qubits in the Qureg (${QUREG_QUBITS}).";

    string TARGET_LIST_WAS_NULL_PTR =
        "Received a NULL pointer where a list of target qubits was expected. It is only valid to pass NULL or nullptr for control qubit lists.";

    string INVALID_TARGET_QUBIT = 
        "Invalid target qubit (${QUBIT_IND}). Must be greater than or equal to zero, and less than the number of qubits in the Qureg (${NUM_QUBITS}).";

    string DUPLICATE_TARGET_QUBITS =
        "The list of target qubits contained duplicates. All qubits must be unique.";



    /*
     * FILE IO
     */

    // TODO: embed filename into error message when tokenSubs supports strings
    string CANNOT_READ_FILE = 
        "Could not load and read the given file. Make sure the file exists and is readable as plaintext.";
}



/*
 * INVALID INPUT RESPONSE BEHAVIOUR
 */

// default C/C++ compatible error response is to simply exit in fail state
void default_invalidQuESTInputError(const char* msg, const char* func) {

    // safe to call even before MPI has been setup
    print(string("")
        + "QuEST encountered a validation error during function " 
        + "'" + func + "':\n" + msg + "\n"
        + "Exiting...");

    // force a synch because otherwise non-main nodes may exit before print, and MPI
    // will then attempt to instantly abort all nodes, losing the error message.
    comm_sync();

    // finalise MPI before error-exit to avoid scaring user with giant MPI error message
    if (comm_isInit())
        comm_end();

    exit(EXIT_FAILURE);
}

// enable default error response to be user-overriden as a weak symbol (even in C, and on Windows)
extern "C" {

    #ifndef _WIN32
        #pragma weak invalidQuESTInputError
        void invalidQuESTInputError(const char* msg, const char* func) {
            default_invalidQuESTInputError(msg, func);
        }
    #elif defined(_WIN64)
        #pragma comment(linker, "/alternatename:invalidQuESTInputError=default_invalidQuESTInputError")   
    #else
        #pragma comment(linker, "/alternatename:_invalidQuESTInputError=_default_invalidQuESTInputError")
    #endif

} // end C++ de-mangler



/*
 * VALIDATION TOGGLE
 *
 * consulted by assertThat below, or earlier by more expensive
 * validation functions to avoid superfluous compute
 */

static bool isValidationEnabled = true;

void validateconfig_enable() {
    isValidationEnabled = true;
}
void validateconfig_disable() {
    isValidationEnabled = false;
}
bool validateconfig_isEnabled() {
    return isValidationEnabled;
}



/*
 * VALIDATION PRECISION
 *
 * which influences how strict unitarity, 
 * Hermiticity and CPTP checks are performed
 */

static qreal validationEpsilon = DEAULT_VALIDATION_EPSILON;

void validateconfig_setEpsilon(qreal eps) {
    validationEpsilon = eps;
}
void validateconfig_setEpsilonToDefault() {
    validationEpsilon = DEAULT_VALIDATION_EPSILON;
}
qreal validateconfig_getEpsilon() {
    return validationEpsilon;
}



/*
 * UTILITIES
 */

// map like "${X}" -> 5, with max-size signed int values to prevent overflows.
// in C++11, these can be initialised with {{"${X}", 5}, ...}
using tokenSubs = std::map<string, long long int>;

string getStringWithSubstitutedVars(string oldStr, tokenSubs varsAndVals) {

    string newStr = oldStr;

    // substitute every var,val pair into newStr
    for (auto varAndVal : varsAndVals) {

        // unpack var -> val
        string var = std::get<0>(varAndVal);
        string val = std::to_string(std::get<1>(varAndVal));

        // assert var is well-formed 
        if (var[0] != '$' || var[1] != '{' || var.back() != '}' )
            error_validationMessageVarWasIllFormed(newStr, var);

        // assert var appears at least once in string
        if (newStr.find(var) == string::npos)
            error_validationMessageVarNotSubstituted(newStr, var);

        // replace every occurrence of var with val
        size_t ind = newStr.find(var);
        while (ind != string::npos) {
            newStr = newStr.replace(ind, var.length(), val);
            ind = newStr.find(var, ind);
        }
    }

    // assert there is no $ left in the strings
    if (newStr.find("$") != string::npos)
        error_validationMessageContainedUnsubstitutedVars(newStr);

    return newStr;
}

void assertThat(bool valid, string msg, const char* func) {

    // skip validation if user has disabled
    if (!isValidationEnabled)
        return;

    // this function does not seek consensus among nodes in distributed 
    // settings in order to remain cheap (consensus requires expensive sync
    // and comm), so is suitable for validation which is gauranteed to be
    // uniform between nodes (assuming user's do not hack in rank-specific
    // arguments to the API!)

    if (!valid)
        invalidQuESTInputError(msg.c_str(), func);
}

void assertThat(bool valid, string msg, tokenSubs vars, const char* func) {

    // substitute the vars into the error message ONLY if the error message
    // is actually going to be thrown; otherwise we're needlessly increasing
    // the time to pass validation with superfluous string manipulation
    if (!valid)
        msg = getStringWithSubstitutedVars(msg, vars);

    // commenting out the above if-statement is convenient during development, 
    // because it ensures that even passed validation was going to otherwise
    // throw a well-formed error message - useful for quickly catching typos!

    assertThat(valid, msg, func);
}

void assertAllNodesAgreeThat(bool valid, string msg, const char* func) {

    // skip below consensus broadcast if user has disabled validation
    if (!isValidationEnabled)
        return;

    // this function seeks consensus among distributed nodes before validation,
    // to ensure all nodes fail together (and ergo all validly finalize MPI)
    // when performing validation that may be non-uniform between nodes. For
    // example, mallocs may succeed on one node but fail on another due to
    // inhomogeneous loads.

    if (comm_isInit())
        valid = comm_isTrueOnAllNodes(valid);

    assertThat(valid, msg, func);
}

void assertAllNodesAgreeThat(bool valid, string msg, tokenSubs vars, const char* func) {

    string newMsg = getStringWithSubstitutedVars(msg, vars);
    assertAllNodesAgreeThat(valid, newMsg, func);
}

bool isIndexListUnique(int* list, int len) {

    // use a max-size bitmask (64 bits)
    long long unsigned int mask = 0;
    int numBits = sizeof(mask) * 8;

    // check internal safety
    for (int i=0; i<len; i++)
        if (list[i] >= numBits)
            error_validationListUniquenessCheckExceededMaskSize();

    // write encountered elements to a bitmask
    for (int i=0; i<len; i++)
        if (1 & (mask >> list[i]))
            return false;
        else
            mask |= 1ULL << list[i];

    return true;
}



/*
 * ENVIRONMENT CREATION
 */

void validate_envNeverInit(bool isQuESTInit, bool isQuESTFinal, const char* caller) {

    assertThat(!isQuESTInit, report::QUEST_ENV_ALREADY_INIT, caller);
    assertThat(!isQuESTFinal, report::QUEST_ENV_ALREADY_FINAL, caller);
}

void validate_newEnvDeploymentMode(int isDistrib, int isGpuAccel, int isMultithread, const char* caller) {

    // deployment flags must be boolean or auto
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(isDistrib     == 0 || isDistrib     == 1 || isDistrib     == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_DISTRIB,     vars, caller);
    assertThat(isGpuAccel    == 0 || isGpuAccel    == 1 || isGpuAccel    == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_GPU_ACCEL,   vars, caller);
    assertThat(isMultithread == 0 || isMultithread == 1 || isMultithread == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_ENV_IS_MULTITHREAD, vars, caller);

    // if a deployment is forced (i.e. not auto), assert that the backend binaries are compiled correctly
    if (isDistrib == 1)
        assertThat(comm_isMpiCompiled(), report::CANNOT_DISTRIB_ENV_WITHOUT_MPI_COMMPILATION, caller);

    if (isMultithread == 1)
        assertThat(cpu_isOpenmpCompiled(), report::CANNOT_MULTITHREAD_ENV_WITHOUT_OPENMP_COMPILATION, caller);

    if (isGpuAccel == 1) {
        assertThat(gpu_isGpuCompiled(), report::CANNOT_GPU_ACCEL_ENV_WITH_GPU_COMPILATION, caller);

        // additionally require GPUs are runtime discoverable
        assertThat(gpu_isGpuAvailable(), report::CANNOT_GPU_ACCEL_ENV_WITH_NO_AVAILABLE_GPUS, caller);
    }

    // note we do not require that when (isDistrib && isGpuAccel), the MPI is CUDA-aware and that
    // the GPU is MPI-supporting. We have a runtime fallback to route through CPU memory in that situation.

    // note we do not here check whether, when distributed, a power-of-2 number of nodes are used,
    // because that requires we first initialise MPI, which we wish the caller to explicitly perform
}

void validate_newEnvDistributedBetweenPower2Nodes(const char* caller) {

    // note that we do NOT finalize MPI before erroring below, because that would necessitate
    // every node (launched by mpirun) serially print the error message, causing spam.
    // Instead, we permit the evil of every MPI process calling exit() and MPI aborting when
    // encountering the first non-zero exit code.

    int numNodes = comm_getNumNodes(); // callable even when not distributed
    tokenSubs vars = {{"${NUM_NODES}", numNodes}};

    assertThat(isPowerOf2(numNodes), report::CANNOT_DISTRIB_ENV_BETWEEN_NON_POW_2_NODES, vars, caller);
}

void validate_gpuIsCuQuantumCompatible(const char* caller) {

    // require node consensus in case nodes have heterogeneous GPU hardware
    assertAllNodesAgreeThat(
        gpu_doesGpuSupportMemPools(), 
        report::CUQUANTUM_DEPLOYED_ON_GPU_WITHOUT_MEM_POOLS, caller);

    // TODO:
    // check other requirements like compute-capability?
}



/*
 * EXISTING ENVIRONMENT
 */

void validate_envIsInit(const char* caller) {

    assertThat(isQuESTEnvInit(), report::QUEST_ENV_NOT_INIT, caller);
}



/*
 * DEBUG UTILITIES
 */

void validate_newEpsilonValue(qreal eps, const char* caller) {

    assertThat(eps >= 0, report::INVALID_NEW_EPSILON, {{"${NEW_EPS}", eps}}, caller);
}

void validate_newNumReportedItems(qindex num, const char* caller) {

    assertThat(num >= 0, report::INVALID_NUM_REPORTED_ITEMS, {{"${NUM_ITEMS}", num}}, caller);
}



/*
 * QUREG CREATION
 */

void assertQuregNonEmpty(int numQubits, const char* caller) {

    assertThat(numQubits >= 1, report::NON_POSITIVE_NUM_QUBITS_IN_CREATE_QUREG, {{"${NUM_QUBITS}",numQubits}}, caller);
}

void assertQuregDeployFlagsRecognised(int isDensMatr, int isDistrib, int isGpuAccel, int isMultithread, const char* caller) {

    // qureg type must be boolean
    assertThat(isDensMatr == 0 || isDensMatr == 1, report::INVALID_OPTION_FOR_QUREG_IS_DENSMATR, caller);

    // deployment flags must be boolean or auto
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(isDistrib     == 0 || isDistrib     == 1 || isDistrib     == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_QUREG_IS_DISTRIB,     vars, caller);
    assertThat(isGpuAccel    == 0 || isGpuAccel    == 1 || isGpuAccel    == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_QUREG_IS_GPU_ACCEL,   vars, caller);
    assertThat(isMultithread == 0 || isMultithread == 1 || isMultithread == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_QUREG_IS_MULTITHREAD, vars, caller);
}

void assertQuregDeploysEnabledByEnv(int isDistrib, int isGpuAccel, int isMultithread, QuESTEnv env, const char* caller) {

    // qureg cannot deploy to backend not already enabled by the environment
    if (!env.isDistributed)
        assertThat(isDistrib     == 0 || isDistrib     == modeflag::USE_AUTO, report::NEW_DISTRIB_QUREG_IN_NON_DISTRIB_ENV, caller);
    if (!env.isGpuAccelerated)
        assertThat(isGpuAccel    == 0 || isGpuAccel    == modeflag::USE_AUTO, report::NEW_GPU_QUREG_IN_NON_GPU_ACCEL_ENV, caller);
    if (!env.isMultithreaded)
        assertThat(isMultithread == 0 || isMultithread == modeflag::USE_AUTO, report::NEW_MULTITHREAD_QUREG_IN_NON_MULTITHREAD_ENV, caller);
}

void assertQuregTotalNumAmpsDontExceedMaxIndex(int numQubits, int isDensMatr, const char* caller) {

    int maxNumQubits = mem_getMaxNumQuregQubitsBeforeIndexOverflow(isDensMatr);

    // make message specific to statevector or density matrix
    string msg = (isDensMatr)? 
        report::NEW_DENSMATR_QUREG_NUM_AMPS_WOULD_EXCEED_QINDEX : 
        report::NEW_STATEVEC_QUREG_NUM_AMPS_WOULD_EXCEED_QINDEX ;

    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits}, 
        {"${MAX_QUBITS}", maxNumQubits}};

    assertThat(numQubits <= maxNumQubits, msg, vars, caller);
}

void assertQuregLocalMemDoesntExceedMaxSizeof(int numQubits, int isDensMatr, int isDistrib, QuESTEnv env, const char* caller) {

    // assume distributed (unless it is force disabled), because that reduces the memory required
    // per node and is ergo more permissive - and the auto-deployer would never choose non-distribution
    // in a distributed env if the memory would exceed the max sizeof!
    int numQuregNodes = (isDistrib == 0 || ! env.isDistributed)? 1 : env.numNodes;
    int maxNumQubits = (int) mem_getMaxNumQuregQubitsBeforeLocalMemSizeofOverflow(isDensMatr, numQuregNodes);

    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits},
        {"${NUM_NODES}",  numQuregNodes},
        {"${MAX_QUBITS}", maxNumQubits}};

    string msg = (isDensMatr)?
        report::NEW_DENSMATR_QUREG_LOCAL_MEM_WOULD_EXCEED_SIZEOF :
        report::NEW_STATEVEC_QUREG_LOCAL_MEM_WOULD_EXCEED_SIZEOF;

    assertThat(numQubits <= maxNumQubits, msg, vars, caller);
}

void assertQuregNotDistributedOverTooManyNodes(int numQubits, int isDensMatr, int isDistrib, QuESTEnv env, const char* caller) {

    // only validate when distribution is forced (auto-deployer will never over-distribute)
    if (isDistrib != 1)
        return;

    // make message specific to statevector or density matrix
    string msg = (isDensMatr)? report::NEW_DISTRIB_DENSMATR_QUREG_HAS_TOO_FEW_AMPS : report::NEW_DISTRIB_STATEVEC_QUREG_HAS_TOO_FEW_AMPS;
    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits},
        {"${NUM_NODES}",  env.numNodes},
        {"${MIN_QUBITS}", std::floor(std::log2(env.numNodes))}};

    int minQubits = mem_getMinNumQubitsForDistribution(env.numNodes);
    assertThat(numQubits >= minQubits, msg, vars, caller);
}

void assertQuregFitsInCpuMem(int numQubits, int isDensMatr, int isDistrib, QuESTEnv env, const char* caller) {

    // attempt to fetch RAM, and simply return if we fail; if we unknowingly
    // didn't have enough RAM, then alloc validation will trigger later
    size_t memPerNode = 0;
    try {
        memPerNode = mem_tryGetLocalRamCapacityInBytes();
    } catch(mem::COULD_NOT_QUERY_RAM e) {
        return;
    }

    // check whether qureg (considering if distributed) fits between node memory(s).
    // note this sets numQuregNodes=1 only when distribution is impossible/switched-off,
    // but not when it would later be automatically disabled. that's fine; the auto-deployer
    // will never disable distribution if local RAM can't store the qureg, so we don't need to
    // validate the auto-deployed-to-non-distributed scenario. we only need to ensure that
    // auto-deploying-to-distribution is permitted by memory capacity.
    int numQuregNodes = (isDistrib == 0 || ! env.isDistributed)? 1 : env.numNodes;
    bool quregFitsInMem = mem_canQuregFitInMemory(numQubits, isDensMatr, numQuregNodes, memPerNode);

    // make error message specific to whether qureg is local or distributed
    string msg = (numQuregNodes == 1)?
        report::NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CPU_MEM :
        report::NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CPU_MEM;

    tokenSubs vars = {
        {"${IS_DENS}",     isDensMatr},
        {"${NUM_QUBITS}",  numQubits},
        {"${QCOMP_BYTES}", sizeof(qcomp)},
        {"${EXP_BASE}",    (isDensMatr)? 4 : 2},
        {"${RAM_SIZE}",    memPerNode}};

    if (numQuregNodes > 1)
        vars["${NUM_NODES}"] = numQuregNodes;

    // require expensive node consensus in case of heterogeneous RAMs
    assertAllNodesAgreeThat(quregFitsInMem, report::NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CPU_MEM, vars, caller);
}

void assertQuregFitsInGpuMem(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, QuESTEnv env, const char* caller) {

    // validate when GPU-acceleration is possible; EVEN if it's auto! Because the auto-deployer will never
    // fall-back to painfully slow CPU if the GPU memory is filled. So if GPU-accel is possible, it must fit qureg.
    if (isGpuAccel == 0 || env.isGpuAccelerated == 0)
        return;

    // we consult the current available local GPU memory (being more strict than is possible for RAM)
    size_t localCurrGpuMem = gpu_getCurrentAvailableMemoryInBytes();

    // check whether qureg (considering if distributed) fits between node GPU memory(s).
    // note this sets numQuregNodes=1 only when distribution is impossible/switched-off,
    // but not when it would later be automatically disabled. that's fine; the auto-deployer
    // will never disable distribution if local GPU memory can't store the qureg, so we don't 
    // need to validate the auto-deployed-to-non-distributed scenario. we only need to ensure 
    // that auto-deploying-to-distribution is permitted by GPU memory capacity.
    int numQuregNodes = (isDistrib == 0 || ! env.isDistributed)? 1 : env.numNodes;
    bool quregFitsInMem = mem_canQuregFitInMemory(numQubits, isDensMatr, numQuregNodes, localCurrGpuMem);

    tokenSubs vars = {
        {"${IS_DENS}",        isDensMatr},
        {"${NUM_QUBITS}",     numQubits},
        {"${QCOMP_BYTES}",    sizeof(qcomp)},
        {"${MIN_VRAM_AVAIL}", localCurrGpuMem}};

    // make error message specific to whether qureg is local or distributed
    if (numQuregNodes == 1) {
        vars["${EXP_BASE}"] = (isDensMatr)? 4 : 2;

        // require expensive node consensus in case of heterogeneous GPU hardware or loads
        assertAllNodesAgreeThat(quregFitsInMem, report::NEW_QUREG_CANNOT_FIT_INTO_NON_DISTRIB_CURRENT_GPU_MEM, vars, caller);

    // when distributed, comm buffers are considered (hence +1 below)
    } else {
        vars["${LOG2_NUM_AMPS}"] = 1 + mem_getEffectiveNumStateVecQubitsPerNode(numQubits, isDensMatr, numQuregNodes);
        vars["${NUM_GPUS}"] = numQuregNodes;

        // require expensive node consensus in case of heterogeneous GPU hardware or loads
        assertAllNodesAgreeThat(quregFitsInMem, report::NEW_QUREG_CANNOT_FIT_INTO_POTENTIALLY_DISTRIB_CURRENT_GPU_MEM, vars, caller);
    }
}

void validate_newQuregNotBothMultithreadedAndGpuAccel(int useGpu, int useMultithread, const char* caller) {

    // note either or both of useGpu and useMultithread are permitted to be modeflag::USE_AUTO (=-1)
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(useGpu != 1 || useMultithread != 1, report::NEW_GPU_QUREG_CANNOT_USE_MULTITHREADING, vars, caller);
}

void validate_newQuregParams(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, int isMultithread, QuESTEnv env, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    assertQuregNonEmpty(numQubits, caller);
    assertQuregDeployFlagsRecognised(isDensMatr, isDistrib, isGpuAccel, isMultithread, caller);
    assertQuregDeploysEnabledByEnv(isDistrib, isGpuAccel, isMultithread, env, caller);
    assertQuregTotalNumAmpsDontExceedMaxIndex(numQubits, isDensMatr, caller);
    assertQuregLocalMemDoesntExceedMaxSizeof(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregNotDistributedOverTooManyNodes(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInCpuMem(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInGpuMem(numQubits, isDensMatr, isDistrib, isGpuAccel, env, caller);

    validate_newQuregNotBothMultithreadedAndGpuAccel(isGpuAccel, isMultithread, caller);
}

void validate_newQuregAllocs(Qureg qureg, const char* caller) {

    // we get node consensus in case mallocs fail on some nodes but not others, as may occur
    // in heterogeneous settings, or where nodes may have other processes and loads hogging RAM. 
    assertAllNodesAgreeThat(mem_isAllocated(qureg.cpuAmps), report::NEW_QUREG_CPU_AMPS_ALLOC_FAILED, caller);

    if (qureg.isGpuAccelerated)
        assertAllNodesAgreeThat(mem_isAllocated(qureg.gpuAmps), report::NEW_QUREG_GPU_AMPS_ALLOC_FAILED, caller);

    if (qureg.isDistributed)
        assertAllNodesAgreeThat(mem_isAllocated(qureg.cpuCommBuffer), report::NEW_QUREG_CPU_COMM_BUFFER_ALLOC_FAILED, caller);

    if (qureg.isDistributed && qureg.isGpuAccelerated)
        assertAllNodesAgreeThat(mem_isAllocated(qureg.gpuCommBuffer), report::NEW_QUREG_GPU_COMM_BUFFER_ALLOC_FAILED, caller);
}



/*
 * EXISTING QUREG
 */

void validate_quregFields(Qureg qureg, const char* caller) {

    // attempt to detect the Qureg was not initialised with createQureg by the 
    // struct fields being randomised, and ergo being dimensionally incompatible
    bool valid = true;
    valid &= (qureg.numQubits > 0);
    valid &= (qureg.isDensityMatrix == 0 || qureg.isDensityMatrix == 1);
    valid &= (qureg.numAmps == powerOf2(((qureg.isDensityMatrix)? 2:1) * qureg.numQubits));

    // we do not bother checking slightly more involved fields like numAmpsPerNode

    tokenSubs vars = {
        {"${DENS_MATR}", qureg.isDensityMatrix},
        {"${NUM_QUBITS}", qureg.numQubits},
        {"${NUM_AMPS}", qureg.numAmps}};
        
    assertThat(valid, report::INVALID_QUREG_FIELDS, vars, caller);

    // In theory, we could check qureg's malloc'd pointers are not-NULL.
    // However, this wouldn't catch when Qureg was un-initialised because most
    // compilers will not automatically set struct pointers to NULL (though
    // that scenario will likely have been caught by above checks). It also
    // will not catch that the Qureg has already been validly created then
    // destroyed because the struct pointers will not set to NULL (because
    // C passes a struct copy). So NULL checks could only check for the specific
    // scenario of the user explicitly overwriting valid pointers with NULL - 
    // this is not worth catching (the eventual NULL segfault might be better)
}

void validate_quregIsStateVector(Qureg qureg, const char* caller) {

    assertThat(!qureg.isDensityMatrix, report::QUREG_NOT_STATE_VECTOR, caller);
}

void validate_quregIsDensityMatrix(Qureg qureg, const char* caller) {

    assertThat(qureg.isDensityMatrix, report::QUREG_NOT_DENSITY_MATRIX, caller);
}



/*
 * MATRIX CREATION
 */

void assertMatrixDeployFlagsRecognised(int isDistrib, const char* caller) {

    // deployment flags must be boolean or auto
    assertThat(
        isDistrib == 0 || isDistrib == 1 || isDistrib == modeflag::USE_AUTO, 
        report::INVALID_OPTION_FOR_MATRIX_IS_DISTRIB,
        {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}}, caller);
}

void assertMatrixDeploysEnabledByEnv(int isDistrib, int envIsDistrib, const char* caller) {

    // cannot deploy to backend not already enabled by the environment
    if (!envIsDistrib)
        assertThat(isDistrib == 0 || isDistrib == modeflag::USE_AUTO, report::NEW_DISTRIB_MATRIX_IN_NON_DISTRIB_ENV, caller);
}

void assertMatrixNonEmpty(int numQubits, const char* caller) {

    assertThat(numQubits >= 1, report::NEW_MATRIX_NUM_QUBITS_NOT_POSITIVE, {{"${NUM_QUBITS}",numQubits}}, caller);
}

void assertMatrixTotalNumElemsDontExceedMaxIndex(int numQubits, bool isDense, const char* caller) {

    int maxNumQubits = mem_getMaxNumMatrixQubitsBeforeIndexOverflow(isDense);

    string msg = (isDense)?
        report::NEW_DIAG_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX :
        report::NEW_COMP_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX ;

    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits}, 
        {"${MAX_QUBITS}", maxNumQubits}};

    assertThat(numQubits <= maxNumQubits, msg, vars, caller);
}

void assertMatrixLocalMemDoesntExceedMaxSizeof(int numQubits, bool isDense, int isDistrib, int numEnvNodes, const char* caller) {

    // 'isDistrib' can be 0 (user-disabled, or matrix is a local type), 1 (user-forced)
    // or modeflag::USE_AUTO=-1 (automatic; the type can be distributed but the user has 
    // not forced it). Currently, only distributed diagonal matrices are supported, so 
    // isDistrib must be specifically 0 for dense matrices.
    if (isDense && isDistrib != 0)
        error_validationEncounteredUnsupportedDistributedDenseMatrix();

    // assume distributed (unless it is force disabled), because that reduces the memory required
    // per node and is ergo more permissive - and the auto-deployer would never choose non-distribution
    // in a distributed env if the memory would exceed the max sizeof!
    int numMatrNodes = (isDistrib == 0)? 1 : numEnvNodes;
    int maxNumQubits = mem_getMaxNumMatrixQubitsBeforeLocalMemSizeofOverflow(isDense, numMatrNodes);

    // make error message specific to whether the matrix is distributed or non-distributed type;
    // non-distributed matrices (e.g. CompMatr) should only ever cause the local error message
    string msg = (numMatrNodes > 1)?
        report::NEW_DISTRIB_DIAG_MATR_LOCAL_MEM_WOULD_EXCEED_SIZEOF :
        ((isDense)?
            report::NEW_LOCAL_COMP_MATR_MEM_WOULD_EXCEED_SIZEOF :
            report::NEW_LOCAL_DIAG_MATR_MEM_WOULD_EXCEED_SIZEOF);
    
    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits},
        {"${MAX_QUBITS}", maxNumQubits}};
    if (numMatrNodes > 1)
        vars["${NUM_NODES}"] = numMatrNodes;

    assertThat(numQubits <= maxNumQubits, msg, vars, caller);
}

void assertMatrixNotDistributedOverTooManyNodes(int numQubits, bool isDense, int isDistrib, int numEnvNodes, const char* caller) {

    // 'isDistrib' can be 0 (user-disabled, or matrix is a local type), 1 (user-forced)
    // or modeflag::USE_AUTO=-1 (automatic; the type can be distributed but the user has 
    // not forced it). Currently, only distributed diagonal matrices are supported, so 
    // isDistrib must be specifically 0 for dense matrices.
    if (isDense && isDistrib != 0)
        error_validationEncounteredUnsupportedDistributedDenseMatrix();

    // only need to validate when distribution is forced (auto-deployer will never over-distribute,
    // and non-distributed types (like CompMatr) will pass isDistrib=0
    if (isDistrib != 1)
        return;

    // distributed diagonal matrices require at least 1 element per node, and while we
    // don't yet support distributed complex matrices, let's for now assume they would
    // require having at leat 1 column per node, like density matrix Quregs do.
    int minQubits = mem_getMinNumQubitsForDistribution(numEnvNodes);

    string msg = report::NEW_DISTRIB_DIAG_MATR_HAS_TOO_FEW_AMPS;
    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits},
        {"${MIN_QUBITS}", minQubits},
        {"${NUM_NODES}",  numEnvNodes}};
        
    assertThat(numQubits >= minQubits, msg, vars, caller);
}

void assertMatrixFitsInCpuMem(int numQubits, bool isDense, int isDistrib, int numEnvNodes, const char* caller) {

    // 'isDistrib' can be 0 (user-disabled, or matrix is a local type), 1 (user-forced)
    // or modeflag::USE_AUTO=-1 (automatic; the type can be distributed but the user has 
    // not forced it). Currently, only distributed diagonal matrices are supported, so 
    // isDistrib must be specifically 0 for dense matrices.
    if (isDense && isDistrib != 0)
        error_validationEncounteredUnsupportedDistributedDenseMatrix();

    // attempt to fetch RAM, and simply return if we fail; if we unknowingly
    // didn't have enough RAM, then alloc validation will trigger later
    size_t memPerNode = 0;
    try {
        memPerNode = mem_tryGetLocalRamCapacityInBytes();
    } catch(mem::COULD_NOT_QUERY_RAM e) {
        return;
    }

    // check whether matrix (considering if distributed) fits between node memory(s).
    // note this sets numMatrNodes=1 only when distribution is impossible/switched-off,
    // but not when it would later be automatically disabled. That's fine; the auto-deployer
    // will never disable distribution if the RAM can't store the entire matrix, so we 
    // don't need to validate the auto-deployed-to-non-distributed scenario. We only need to 
    // ensure that auto-deploying-to-distribution is permitted by memory capacity. Note too
    // that the distinction between (env.isDistributed) and (env.numNodes>1) is unimportant
    // for matrix structs because they never store communication buffers.
    int numMatrNodes = (isDistrib == 0)? 1 : numEnvNodes;
    bool matrFitsInMem = mem_canMatrixFitInMemory(numQubits, isDense, numMatrNodes, memPerNode);

    // specialise error message to whether matrix is distributed and dense or diag
    string msg = (isDense)?
        report::NEW_LOCAL_COMP_MATR_CANNOT_FIT_INTO_CPU_MEM :
        ((numMatrNodes == 1)?
            report::NEW_LOCAL_DIAG_MATR_CANNOT_FIT_INTO_CPU_MEM :
            report::NEW_DISTRIB_DIAG_MATR_CANNOT_FIT_INTO_CPU_MEM );

    tokenSubs vars = {
        {"${NUM_QUBITS}",  numQubits},
        {"${QCOMP_BYTES}", sizeof(qcomp)},
        {"${RAM_SIZE}",    memPerNode}};

    if (numMatrNodes > 1) {
        vars["${NUM_NODES}"] = numMatrNodes;
        vars["${NUM_QB_MINUS_LOG_NODES}"] = numQubits - logBase2(numMatrNodes);
    }
    
    // TODO:
    // seek expensive node consensus in case of heterogeneous RAM - alas this may induce
    // unnecessary slowdown (due to sync and broadcast) in applications allocating many
    // small matrices in the heap. If this turns out to be the case, we could opt to
    // enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    // chance of it fitting into some node RAM but not others isn't negligible.
    assertAllNodesAgreeThat(matrFitsInMem, msg, vars, caller);
}

void assertMatrixFitsInGpuMem(int numQubits, bool isDense, int isDistrib, int isGpu, int numEnvNodes, const char* caller) {

    // 'isDistrib' can be 0 (user-disabled, or matrix is a local type), 1 (user-forced)
    // or modeflag::USE_AUTO=-1 (automatic; the type can be distributed but the user has 
    // not forced it). Currently, only distributed diagonal matrices are supported, so 
    // isDistrib must be specifically 0 for dense matrices.
    if (isDense && isDistrib != 0)
        error_validationEncounteredUnsupportedDistributedDenseMatrix();

    // 'isGpu' can be 0 (user-disabled, or matrix is a local type), 1 (user-forced), or
    // modeflag::USE_AUTO=-1 (automatic; the type can be GPU-accel but the user has not forced it).
    // Whenever GPU isn't forced, we do not need to check there's sufficient GPU memory, because
    // if there's insufficient, the auto-deployer will never choose it. So proceed only if GPU is forced.
    if (isGpu != 1)
        return;

    // we consult the current available local GPU memory (being more strict than is possible for RAM)
    size_t localCurrGpuMem = gpu_getCurrentAvailableMemoryInBytes();

    // check whether matrix (considering if distributed) fits between node memory(s).
    // note this sets numMatrNodes=1 only when distribution is impossible/switched-off,
    // but not when it would later be automatically disabled. That's fine; the auto-deployer
    // will never disable distribution if the local GPU can't store the entire matrix, so we 
    // don't need to validate the auto-deployed-to-non-distributed scenario. We only need to 
    // ensure that auto-deploying-to-distribution is permitted by memory capacity. Note too
    // that the distinction between (env.isDistributed) and (env.numNodes>1) is unimportant
    // for matrix structs because they never store communication buffers.
    int numMatrNodes = (isDistrib == 0)? 1 : numEnvNodes;
    bool matrFitsInMem = mem_canMatrixFitInMemory(numQubits, isDense, numMatrNodes, localCurrGpuMem);

    // specialise error message to whether matrix is distributed and dense or diag
    string msg = (isDense)?
        report::NEW_LOCAL_COMP_MATR_CANNOT_FIT_INTO_GPU_MEM :
        ((numMatrNodes == 1)?
            report::NEW_LOCAL_DIAG_MATR_CANNOT_FIT_INTO_GPU_MEM :
            report::NEW_DISTRIB_DIAG_MATR_CANNOT_FIT_INTO_GPU_MEM );

    tokenSubs vars = {
        {"${NUM_QUBITS}",  numQubits},
        {"${QCOMP_BYTES}", sizeof(qcomp)},
        {"${VRAM_SIZE}",   localCurrGpuMem}};
        
    if (numMatrNodes > 1) {
        vars["${NUM_NODES}"] = numMatrNodes;
        vars["${NUM_QB_MINUS_LOG_NODES}"] = numQubits - logBase2(numMatrNodes);
    }
    
    // TODO:
    // seek expensive node consensus in case of heterogeneous GPU hardware - alas this may 
    // induce unnecessary slowdown (due to sync and broadcast) in applications allocating many
    // small matrices in the GPU. If this turns out to be the case, we could opt to
    // enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    // chance of it fitting into some GPU's memory but not others isn't negligible.
    assertAllNodesAgreeThat(matrFitsInMem, msg, vars, caller);
}

void assertNewMatrixParamsAreValid(int numQubits, int useDistrib, int useGpu, bool isDenseType, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    QuESTEnv env = getQuESTEnv();
    assertMatrixNonEmpty(numQubits, caller);
    assertMatrixTotalNumElemsDontExceedMaxIndex(numQubits, isDenseType, caller);
    assertMatrixLocalMemDoesntExceedMaxSizeof(numQubits,  isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixNotDistributedOverTooManyNodes(numQubits, isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixFitsInCpuMem(numQubits, isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixFitsInGpuMem(numQubits, isDenseType, useDistrib, useGpu, env.numNodes, caller);
}

void validate_newCompMatrParams(int numQubits, const char* caller) {
    validate_envIsInit(caller);

    // CompMatr can never be distributed
    int useDistrib = 0;

    // CompMatr is always GPU accelerated whenever enabled by the environment
    bool useGpu = getQuESTEnv().isGpuAccelerated;

    // CompMatr stores 2^(2*numQubits) elems
    bool isDenseType = true;

    assertNewMatrixParamsAreValid(numQubits, useDistrib, useGpu, isDenseType, caller);
}
void validate_newDiagMatrParams(int numQubits, const char* caller) {
    validate_envIsInit(caller);

    // DiagMatr can never be distributed
    int useDistrib = 0;

    // DiagMatr is always GPU accelerated whenever enabled by the environment
    bool useGpu = getQuESTEnv().isGpuAccelerated;

    // DiagMatr stores only the diagonals; 2^numQubits elems
    bool isDenseType = false;

    assertNewMatrixParamsAreValid(numQubits, useDistrib, useGpu, isDenseType, caller);
}
void validate_newFullStateDiagMatrParams(int numQubits, int useDistrib, int useGpu, const char* caller) {
    validate_envIsInit(caller);

    // FullStateDiagMatr stores only the diagonals
    bool isDenseType = false;

    assertNewMatrixParamsAreValid(numQubits, useDistrib, useGpu, isDenseType, caller);
}

// T can be CompMatr, DiagMatr, FullStateDiagMatr (i.e. heap-based matrices)
template <typename T>
void assertNewMatrixAllocsSucceeded(T matr, qindex numBytes, const char* caller) {

    // we expensively get node consensus about malloc failure, in case of heterogeneous hardware/loads,
    // but we avoid this potetially expensive synchronisation if validation is anyway disabled
    // (which also avoids us enumerating the matrix rows)
    if (!isValidationEnabled)
        return;

    tokenSubs vars = {{"${NUM_BYTES}", numBytes}};

    // assert CPU array (which may be nested arrays) all allocated successfully
    bool isAlloc;
    if constexpr (util_isDenseMatrixType<T>())
        isAlloc = mem_isAllocated(matr.cpuElems, matr.numRows);
    else
        isAlloc = mem_isAllocated(matr.cpuElems);
    assertAllNodesAgreeThat(isAlloc, report::NEW_MATRIX_CPU_ELEMS_ALLOC_FAILED, vars, caller);

    // optionally assert GPU memory was malloc'd successfully
    bool gpuShouldBeAlloc = getQuESTEnv().isGpuAccelerated;
    if constexpr (util_isFullStateDiagMatr<T>())
        gpuShouldBeAlloc &= matr.isGpuAccelerated;
        
    if (gpuShouldBeAlloc)
        assertAllNodesAgreeThat(mem_isAllocated(util_getGpuMemPtr(matr)), report::NEW_MATRIX_GPU_ELEMS_ALLOC_FAILED, vars, caller);

    // assert the teeny-tiny heap flags are alloc'd
    vars["${NUM_BYTES}"] = sizeof(*(matr.isUnitary));
    assertAllNodesAgreeThat(mem_isAllocated(matr.isUnitary),    report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
    assertAllNodesAgreeThat(mem_isAllocated(matr.isHermitian),  report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
    assertAllNodesAgreeThat(mem_isAllocated(matr.wasGpuSynced), report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
}

void validate_newMatrixAllocs(CompMatr matr, const char* caller) {

    bool isDenseMatrix = true;
    int numNodes = 1;
    qindex numBytes = mem_getLocalMatrixMemoryRequired(matr.numQubits, isDenseMatrix, numNodes);
    assertNewMatrixAllocsSucceeded(matr, numBytes, caller);
}
void validate_newMatrixAllocs(DiagMatr matr, const char* caller) {

    bool isDenseMatrix = false;
    int numNodes = 1;
    qindex numBytes = mem_getLocalMatrixMemoryRequired(matr.numQubits, isDenseMatrix, numNodes);
    assertNewMatrixAllocsSucceeded(matr, numBytes, caller);
}
void validate_newMatrixAllocs(FullStateDiagMatr matr, const char* caller) {

    bool isDenseMatrix = false;
    int numNodes = (matr.isDistributed)? getQuESTEnv().numNodes : 1;
    qindex numBytes = mem_getLocalMatrixMemoryRequired(matr.numQubits, isDenseMatrix, numNodes);
    assertNewMatrixAllocsSucceeded(matr, numBytes, caller);
}



/*
 * MATRIX INITIALISATION
 */

void validate_matrixNumNewElems(int numQubits, vector<vector<qcomp>> elems, const char* caller) {

    // CompMatr accept 2D elems   
    qindex dim = powerOf2(numQubits);
    tokenSubs vars = {
        {"${NUM_QUBITS}",        numQubits},
        {"${NUM_EXPECTED_ROWS}", dim},
        {"${NUM_GIVEN_ROWS}",    elems.size()}};

    assertThat(dim == (qindex) elems.size(), report::COMP_MATR_NEW_ELEMS_WRONG_NUM_ROWS, vars, caller);

    for(auto & row : elems) {

        vars = {
            {"${NUM_QUBITS}",      numQubits},
            {"${EXPECTED_DIM}",    dim},
            {"${NUM_GIVEN_ELEMS}", row.size()}};

        assertThat(dim == (qindex) row.size(), report::COMP_MATR_NEW_ELEMS_WRONG_ROW_DIM, vars, caller);
    }
}
void validate_matrixNumNewElems(int numQubits, vector<qcomp> elems, const char* caller) {

    // DiagMatr accept 1D elems
    qindex dim = powerOf2(numQubits);
    tokenSubs vars = {
        {"${NUM_QUBITS}",        numQubits},
        {"${NUM_EXPECTED_ELEMS}", dim},
        {"${NUM_GIVEN_ELEMS}",    elems.size()}};

    assertThat(dim == (qindex) elems.size(), report::DIAG_MATR_WRONG_NUM_NEW_ELEMS, vars, caller);
}

void validate_fullStateDiagMatrNewElems(FullStateDiagMatr matr, qindex startInd, qindex numElems, const char* caller) {

    assertThat(
        startInd >= 0 && startInd < matr.numElems, 
        report::FULL_STATE_DIAG_MATR_NEW_ELEMS_INVALID_START_INDEX, 
        {{"${START_IND}", startInd}, {"${MATR_NUM_ELEMS}", matr.numElems}},
        caller);

    assertThat(
        numElems > 0,
        report::FULL_STATE_DIAG_MATR_NEW_ELEMS_NUM_IS_NON_POSITIVE, 
        {{"${NUM_ELEMS}", numElems}}, caller);

    assertThat(
        numElems <= matr.numElems,
        report::FULL_STATE_DIAG_MATR_NEW_ELEMS_NUM_EXCEEDS_MAX_NUM, 
        {{"${NEW_NUM_ELEMS}", numElems}, {"${MATR_NUM_ELEMS}", matr.numElems}},
        caller);

    qindex endIndExcl = startInd + numElems;
    tokenSubs vars = {
        {"${START_IND}",      startInd},
        {"${MATR_NUM_ELEMS}", matr.numElems},
        {"${END_IND_EXCL}",   endIndExcl}};
       
    assertThat(
        endIndExcl <= matr.numElems, 
        report::FULL_STATE_DIAG_MATR_NEW_ELEMS_EXCEEDS_END_INDEX,  
        vars, caller);
}

void validate_matrixNumQubitsMatchesParam(int numMatrQubits, int numSetterQubits, const char* caller) {

    tokenSubs vars = {
        {"${NUM_SETTER_QUBITS}", numSetterQubits},
        {"${NUM_MATRIX_QUBITS}", numMatrQubits}};

    assertThat(numMatrQubits == numSetterQubits, report::MATR_NUM_QUBITS_MISMATCHES_INLINE_SETTER, vars, caller);
}

void validate_declaredNumElemsMatchesVectorLength(qindex numElems, qindex vecLength, const char* caller) {

    tokenSubs vars = {
        {"${NUM_ELEMS}", numElems},
        {"${VEC_LENGTH}", vecLength}};

    assertThat(numElems == vecLength, report::MATR_NUM_ELEMS_MISMATCHES_VEC_LENGTH_IN_INLINE_SETTER, vars, caller);
}

void validate_multiVarFuncQubits(int numMatrQubits, int* numQubitsPerVar, int numVars, const char* caller) {

    assertThat(numVars > 0, report::MULTI_VAR_FUNC_INVALID_NUM_VARS, {{"${NUM_VARS}", numVars}}, caller);

    for (int v=0; v<numVars; v++)
        assertThat(numQubitsPerVar[v] > 0, report::MULTI_VAR_FUNC_INVALID_NUM_QUBITS_PER_VAR, {{"${VAR_IND}",v},{"${VAR_QUBITS}",numQubitsPerVar[v]}}, caller);

    int numVarQubits = 0;
    for (int v=0; v<numVars; v++)
        numVarQubits += numQubitsPerVar[v];

    assertThat(numMatrQubits == numVarQubits, report::MULTI_VAR_FUNC_MISMATCHING_NUM_QUBITS, 
        {{"${NUM_MATR_QUBITS}", numMatrQubits}, {"${NUM_VAR_QUBITS}", numVarQubits}}, caller);
}

void validate_funcVarSignedFlag(int areSigned, const char* caller) {

    assertThat(areSigned == 0 || areSigned == 1, report::MULTI_VAR_FUNC_INVALID_ARE_SIGNED_FLAG, {{"${ARE_SIGNED}", areSigned}}, caller);
}



/*
 * EXISTING MATRICES
 */

// T can be CompMatr, DiagMatr, FullStateDiagMatr
template <class T>
void assertAdditionalHeapMatrixFieldsAreValid(T matr, const char* caller) {

    // assert heap pointers are not NULL
    assertThat(mem_isAllocated(matr.isUnitary),    report::INVALID_HEAP_FLAG_PTR, caller);
    assertThat(mem_isAllocated(matr.isHermitian),  report::INVALID_HEAP_FLAG_PTR, caller);
    assertThat(mem_isAllocated(matr.wasGpuSynced), report::INVALID_HEAP_FLAG_PTR, caller);

    tokenSubs vars = {{"${BAD_FLAG}", 0}, {"${UNKNOWN_FLAG}", validate_STRUCT_PROPERTY_UNKNOWN_FLAG}};

    // assert isUnitary has valid value
    int flag = *matr.isUnitary;
    vars["${BAD_FLAG}"] = flag;
    assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);

    // assert isHermitian has valid value
    flag = *matr.isHermitian;
    vars["${BAD_FLAG}"] = flag;
    assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);

    // assert wasGpuSynced has a valid value
    flag = *matr.wasGpuSynced;
    vars["${BAD_FLAG}"] = flag;
    assertThat(flag == 0 || flag == 1, report::INVALID_HEAP_FLAG_VALUE, vars, caller);

    // checks whether users have, after destroying their struct, manually set the outer
    // heap-memory pointers to NULL. We do not check inner pointers of 2D structures (which may
    // be too expensive to enumerate), and we cannot determine whether the struct was not validly
    // created because un-initialised structs will have random non-NULL pointers.
    assertThat(mem_isOuterAllocated(matr.cpuElems), report::INVALID_MATRIX_CPU_ELEMS_PTR, caller);

    // check whether GPU status/memory pointers are consistent with env
    validate_envIsInit(caller);
    bool envIsGpuAccel = getQuESTEnv().isGpuAccelerated;
    bool matrHasGpuAlloc = mem_isOuterAllocated(util_getGpuMemPtr(matr));

    // FullStateDiagMatr can disable GPU-accel even in GPU-accelerated environments
    if constexpr (util_isFullStateDiagMatr<T>()) {
        if (matr.isGpuAccelerated) {
            assertThat(envIsGpuAccel, report::FULL_STATE_DIAG_MATR_GPU_ACCEL_IN_NON_GPU_ENV, caller);
            assertThat(matrHasGpuAlloc, report::MATRIX_GPU_ELEMS_PTR_UNEXPECTEDLY_NULL, caller);
        } else
            assertThat( ! matrHasGpuAlloc, report::MATRIX_GPU_ELEMS_PTR_UNEXPECTEDLY_NOT_NULL, caller);

    // all other heap-matrices always GPU-accel in GPU-accelerated envs
    } else {
        if (envIsGpuAccel)
            assertThat(matrHasGpuAlloc, report::MATRIX_GPU_ELEMS_PTR_UNEXPECTEDLY_NULL, caller);
        else
            assertThat( ! matrHasGpuAlloc, report::MATRIX_GPU_ELEMS_PTR_UNEXPECTEDLY_NOT_NULL, caller);
    }
}

// T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template <class T>
void assertMatrixFieldsAreValid(T matr, int expectedNumQb, string badFieldMsg, const char* caller) {

    qindex dim = util_getMatrixDim(matr);
    tokenSubs vars = {
        {"${NUM_QUBITS}", matr.numQubits},
        {"${NUM_ROWS}",   dim}};

    // assert correct fixed-size numQubits (caller gaurantees this passes for dynamic-size),
    // where the error message string already contains the expected numQb
    assertThat(matr.numQubits == expectedNumQb, badFieldMsg, vars, caller);

    // validate .numQubits and .numRows or .numElems
    qindex expectedDim = powerOf2(matr.numQubits);
    assertThat(matr.numQubits >= 1, badFieldMsg, vars, caller);
    assertThat(dim == expectedDim,  badFieldMsg, vars, caller);

    if constexpr (util_isHeapMatrixType<T>())
        assertAdditionalHeapMatrixFieldsAreValid(matr, caller);

    // we do not bother checking slightly more involved fields like numAmpsPerNode - there's
    // no risk that they're wrong (because they're const so users cannot modify them) unless 
    // the struct was unitialised, which we have already validated against
}

void validate_matrixFields(CompMatr1 matr, const char* caller) { assertMatrixFieldsAreValid(matr, 1,              report::INVALID_COMP_MATR_1_FIELDS, caller); }
void validate_matrixFields(CompMatr2 matr, const char* caller) { assertMatrixFieldsAreValid(matr, 2,              report::INVALID_COMP_MATR_2_FIELDS, caller); }
void validate_matrixFields(CompMatr  matr, const char* caller) { assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_COMP_MATR_FIELDS,   caller); }
void validate_matrixFields(DiagMatr1 matr, const char* caller) { assertMatrixFieldsAreValid(matr, 1,              report::INVALID_DIAG_MATR_1_FIELDS, caller); }
void validate_matrixFields(DiagMatr2 matr, const char* caller) { assertMatrixFieldsAreValid(matr, 2,              report::INVALID_DIAG_MATR_2_FIELDS, caller); }
void validate_matrixFields(DiagMatr  matr, const char* caller) { assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_DIAG_MATR_FIELDS,   caller); }
void validate_matrixFields(FullStateDiagMatr matr, const char* caller) { assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_FULL_STATE_DIAG_MATR_FIELDS, caller); }

// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void assertMatrixIsSynced(T matr, string errMsg, const char* caller) {

    // we don't need to perform any sync check in CPU-only mode
    if (!mem_isAllocated(util_getGpuMemPtr(matr)))
        return;

    // check if GPU amps have EVER been overwritten; we sadly cannot check the LATEST changes were pushed though
    assertThat(*(matr.wasGpuSynced) == 1, errMsg, caller);
}

void validate_matrixIsSynced(CompMatr matr, const char* caller) {
    assertMatrixIsSynced(matr, report::COMP_MATR_NOT_SYNCED_TO_GPU, caller);
}
void validate_matrixIsSynced(DiagMatr matr, const char* caller) {
    assertMatrixIsSynced(matr, report::DIAG_MATR_NOT_SYNCED_TO_GPU, caller);
}
void validate_matrixIsSynced(FullStateDiagMatr matr, const char* caller) {
    assertMatrixIsSynced(matr, report::FULL_STATE_DIAG_MATR_NOT_SYNCED_TO_GPU, caller);
}

// type T can be CompMatr, DiagMatr, FullStateDiagMatr
template <class T> 
void ensureMatrixUnitarityIsKnown(T matr) {

    // do nothing if we already know unitarity
    if (*(matr.isUnitary) != validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        return;

    // determine local unitarity, modifying matr.isUnitary. This will
    // involve MPI communication if matr is a distributed type
    *(matr.isUnitary) = util_isUnitary(matr, validationEpsilon);
}

// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template <class T> 
void assertMatrixIsUnitary(T matr, const char* caller) {

    // avoid expensive unitarity check if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    // unitarity is determined differently depending on matrix type
    bool isUnitary = false;

    // fixed-size matrices have their unitarity calculated afresh (since cheap)
    if constexpr (util_isFixedSizeMatrixType<T>())
        isUnitary = util_isUnitary(matr, validationEpsilon);

    // dynamic matrices have their field consulted, which may invoke lazy eval and global synchronisation
    else {
        ensureMatrixUnitarityIsKnown(matr);
        isUnitary = *(matr.isUnitary);
    }

    assertThat(isUnitary, report::MATRIX_NOT_UNITARY, caller);
}

void validate_matrixIsUnitary(CompMatr1 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(CompMatr2 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(CompMatr matr, const char* caller) {
    validate_matrixFields(matr, caller);
    validate_matrixIsSynced(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(DiagMatr1 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(DiagMatr2 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(DiagMatr matr, const char* caller) {
    validate_matrixFields(matr, caller);
    validate_matrixIsSynced(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}
void validate_matrixIsUnitary(FullStateDiagMatr matr, const char* caller) {
    validate_matrixIsSynced(matr, caller);
    validate_matrixFields(matr, caller);
    assertMatrixIsUnitary(matr, caller);
}

// type T can be CompMatr, DiagMatr, FullStateDiagMatr
template <class T> 
void ensureMatrHermiticityIsKnown(T matr) {

    // do nothing if we already know hermiticity
    if (*(matr.isHermitian) != validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        return;

    // determine local unitarity, modifying matr.isHermitian
    *(matr.isHermitian) = util_isHermitian(matr, validationEpsilon);
}

// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template <class T> 
void assertMatrixIsHermitian(T matr, const char* caller) {

    // avoid expensive unitarity check if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    // unitarity is determined differently depending on matrix type
    bool isHermitian = false;

    // fixed-size matrices have their hermiticity calculated afresh (since cheap)
    if constexpr (util_isFixedSizeMatrixType<T>())
        isHermitian = util_isHermitian(matr, validationEpsilon);

    // dynamic matrices have their field consulted, which may invoke lazy eval and global synchronisation
    else {
        ensureMatrHermiticityIsKnown(matr);
        isHermitian = *(matr.isHermitian);
    }

    assertThat(isHermitian, report::MATRIX_NOT_HERMITIAN, caller);
}

void validate_matrixIsHermitian(CompMatr1 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsHermitian(matr, caller);
}
void validate_matrixIsHermitian(CompMatr2 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsHermitian(matr, caller);
}
void validate_matrixIsHermitian(CompMatr matr, const char* caller) {
    validate_matrixFields(matr, caller);
    validate_matrixIsSynced(matr, caller);
    assertMatrixIsHermitian(matr, caller);
}
void validate_matrixIsHermitian(DiagMatr1 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsHermitian(matr, caller);
}
void validate_matrixIsHermitian(DiagMatr2 matr, const char* caller) {
    validate_matrixFields(matr, caller);
    assertMatrixIsHermitian(matr, caller);
}
void validate_matrixIsHermitian(DiagMatr matr, const char* caller) {
    validate_matrixFields(matr, caller);
    validate_matrixIsSynced(matr, caller);
    assertMatrixIsHermitian(matr, caller);
}
void validate_matrixIsHermitian(FullStateDiagMatr matr, const char* caller) {
    validate_matrixIsSynced(matr, caller);
    validate_matrixFields(matr, caller);
    assertMatrixIsHermitian(matr, caller);
}

void validate_matrixAndQuregAreCompatible(FullStateDiagMatr matr, Qureg qureg, const char* caller) {

    // we do not need to define this function for the other matrix types,
    // since their validation will happen through validation of the
    // user-given list of target qubits. But we do need to define it for
    // FullStatedDiagMatr to check both distribution compatibility, and
    // that dimensions match

    tokenSubs vars = {
        {"${NUM_MATR_QUBITS}",  matr.numQubits},
        {"${NUM_QUREG_QUBITS}", qureg.numQubits}};

    // dimensions must match
    assertThat(matr.numQubits == qureg.numQubits, report::FULL_STATE_DIAG_MATR_MISMATCHES_QUREG_DIM, vars, caller);

    // when matrix is duplicated on every node, its application is trivial
    if (!matr.isDistributed)
        return;

    // but when it's distributed, so too must be the qureg; the precise reason why is 
    // specific to whether qureg is a statevector or density matrix, but boils down
    // to there being no communication buffers available to broadcast matr
    assertThat(qureg.isDistributed, report::FULL_STATE_DIAG_MATR_IS_DISTRIB_BUT_QUREG_ISNT, caller); // did not pass vars
}



/*
 * SUPEROPERATOR CREATION
 */

void assertSuperOpNonEmpty(int numQubits, const char* caller) {

    assertThat(numQubits >= 1, report::NEW_SUPER_OP_NUM_QUBITS_NOT_POSITIVE, {{"${NUM_QUBITS}",numQubits}}, caller);
}

void assertSuperOpTotalNumElemsDontExceedMaxIndex(int numQubits, bool isInKrausMap, const char* caller) {

    int maxQubits = mem_getMaxNumSuperOpQubitsBeforeIndexOverflow();

    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits}, 
        {"${MAX_QUBITS}", maxQubits}};

    auto msg = (isInKrausMap)? 
        report::NEW_KRAUS_MAPS_SUPER_OP_NUM_ELEMS_WOULD_EXCEED_QINDEX: 
        report::NEW_SUPER_OP_NUM_ELEMS_WOULD_EXCEED_QINDEX;
    assertThat(numQubits <= maxQubits, msg, vars, caller);
}

void assertSuperOpLocalMemDoesntExceedMaxSizeof(int numQubits, bool isInKrausMap, const char* caller) {

    int maxNumQubits = mem_getMaxNumSuperOpQubitsBeforeLocalMemSizeofOverflow();

    tokenSubs vars = {
        {"${NUM_QUBITS}", numQubits},
        {"${MAX_QUBITS}", maxNumQubits},
        {"${QCOMP_BYTES}", sizeof(qcomp)}};

    auto msg = (isInKrausMap)?
        report::NEW_KRAUS_MAPS_SUPER_OP_MEM_WOULD_EXCEED_SIZEOF:
        report::NEW_SUPER_OP_MEM_WOULD_EXCEED_SIZEOF;
    assertThat(numQubits <= maxNumQubits, msg, vars, caller);
}

void assertSuperOpFitsInCpuMem(int numQubits, bool isInKrausMap, const char* caller) {

    // attempt to fetch RAM, and simply return if we fail; if we unknowingly
    // didn't have enough RAM, then alloc validation will trigger later
    size_t memPerNode = 0;
    try {
        memPerNode = mem_tryGetLocalRamCapacityInBytes();
    } catch(mem::COULD_NOT_QUERY_RAM e) {
        return;
    }

    bool matrFitsInMem = mem_canSuperOpFitInMemory(numQubits, memPerNode);

    tokenSubs vars = {
        {"${NUM_QUBITS}",  numQubits},
        {"${QCOMP_BYTES}", sizeof(qcomp)},
        {"${RAM_SIZE}",    memPerNode}};
    
    // TODO:
    // seek expensive node consensus in case of heterogeneous RAM - alas this may induce
    // unnecessary slowdown (due to sync and broadcast) in applications allocating many
    // small matrices in the heap. If this turns out to be the case, we could opt to
    // enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    // chance of it fitting into some node RAM but not others isn't negligible.
    auto msg = (isInKrausMap)?
        report::NEW_KRAUS_MAPS_SUPER_OP_CANNOT_FIT_INTO_CPU_MEM:
        report::NEW_SUPER_OP_CANNOT_FIT_INTO_CPU_MEM;
    assertAllNodesAgreeThat(matrFitsInMem, msg, vars, caller);
}

void assertSuperOpFitsInGpuMem(int numQubits, int isEnvGpuAccel, bool isInKrausMap, const char* caller) {

    // kraus map GPU memory will always be allocated when env is GPU-accelerated
    if (!isEnvGpuAccel)
        return;

    // we consult the current available local GPU memory (being more strict than is possible for RAM)
    size_t localCurrGpuMem = gpu_getCurrentAvailableMemoryInBytes();
    bool matrFitsInMem = mem_canSuperOpFitInMemory(numQubits, localCurrGpuMem);

    tokenSubs vars = {
        {"${NUM_QUBITS}",  numQubits},
        {"${QCOMP_BYTES}", sizeof(qcomp)},
        {"${VRAM_SIZE}",   localCurrGpuMem}};

    auto msg = (isInKrausMap)?
        report::NEW_KRAUS_MAPS_SUPER_OP_CANNOT_FIT_INTO_GPU_MEM:
        report::NEW_SUPER_OP_CANNOT_FIT_INTO_GPU_MEM;

    // TODO:
    // seek expensive node consensus in case of heterogeneous GPU hardware - alas this may 
    // induce unnecessary slowdown (due to sync and broadcast) in applications allocating many
    // small matrices in the GPU. If this turns out to be the case, we could opt to
    // enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    // chance of it fitting into some GPU's memory but not others isn't negligible.
    assertAllNodesAgreeThat(matrFitsInMem, msg, vars, caller);
}

void validate_newSuperOpParams(int numQubits, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    validate_envIsInit(caller);
    QuESTEnv env = getQuESTEnv();

    bool isInKrausMap = false;
    assertSuperOpNonEmpty(numQubits, caller);
    assertSuperOpTotalNumElemsDontExceedMaxIndex(numQubits, isInKrausMap, caller);
    assertSuperOpLocalMemDoesntExceedMaxSizeof(numQubits, isInKrausMap, caller);
    assertSuperOpFitsInCpuMem(numQubits, isInKrausMap, caller);
    assertSuperOpFitsInGpuMem(numQubits, env.isGpuAccelerated, isInKrausMap, caller);
}

void assertNewSuperOpAllocs(SuperOp op, bool isInKrausMap, const char* caller) {
    tokenSubs vars = {{"${NUM_BYTES}", mem_getLocalSuperOpMemoryRequired(op.numQubits)}};

    // we expensively get node consensus about malloc failure, in case of heterogeneous hardware/loads,
    // but we avoid this if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    // assert CPU array of rows was alloc'd successfully
    auto msg = (isInKrausMap)?
        report::NEW_KRAUS_MAPS_SUPER_OP_CPU_ELEMS_ALLOC_FAILED:
        report::NEW_SUPER_OP_CPU_ELEMS_ALLOC_FAILED;
    assertAllNodesAgreeThat(mem_isAllocated(op.cpuElems, op.numRows), msg, vars, caller);

    // optionally assert GPU memory was malloc'd successfully
    msg = (isInKrausMap)?
        report::NEW_KRAUS_MAPS_SUPER_OP_GPU_ELEMS_ALLOC_FAILED:
        report::NEW_SUPER_OP_GPU_ELEMS_ALLOC_FAILED;
    if (getQuESTEnv().isGpuAccelerated)
        assertAllNodesAgreeThat(mem_isAllocated(util_getGpuMemPtr(op)), msg, vars, caller);

    // assert the teeny-tiny heap flag was alloc'd
    vars["${NUM_BYTES}"] = sizeof(*(op.wasGpuSynced));
    assertAllNodesAgreeThat(mem_isAllocated(op.wasGpuSynced), report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
}

void validate_newSuperOpAllocs(SuperOp op, const char* caller) {

    bool isInKrausMap = false;
    assertNewSuperOpAllocs(op, isInKrausMap, caller);
}

void validate_newInlineSuperOpDimMatchesVectors(int numDeclaredQubits, vector<vector<qcomp>> matrix, const char* caller) {

    // avoid potentially expensive matrix enumeration if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    qindex numDeclaredRows = powerOf2(2*numDeclaredQubits);
    tokenSubs vars = {
        {"${NUM_DECLARED_QUBITS}", numDeclaredQubits},
        {"${NUM_DECLARED_ROWS}",   numDeclaredRows},
        {"${GIVEN_DIM}",           matrix.size()}};

    // assert the given matrix has the correct number of rows
    assertThat(numDeclaredRows == (qindex) matrix.size(), report::NEW_INLINE_SUPER_OP_MATRIX_WRONG_NUM_ROWS, vars, caller);

    // and that each row has the correct length
    for (qindex r=0; r<numDeclaredRows; r++) {
        qindex numGivenCols = matrix[r].size();
        vars["${GIVEN_DIM}"] = numGivenCols;
        assertThat(numDeclaredRows == numGivenCols, report::NEW_INLINE_SUPER_OP_MATRIX_WRONG_NUM_COLS, vars, caller);
    }
}



/*
 * SUPEROPERATOR INITIALISATION
 */

void validate_superOpNewMatrixDims(SuperOp op, vector<vector<qcomp>> matrix, const char* caller) {

    // avoid potentially expensive matrix enumeration if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    tokenSubs vars = {
        {"${NUM_QUBITS}",   op.numQubits},
        {"${EXPECTED_DIM}", op.numRows},
        {"${GIVEN_DIM}",    matrix.size()}};

    // assert the matrix has the correct number of rows
    assertThat(op.numRows == (qindex) matrix.size(), report::SUPER_OP_NEW_MATRIX_ELEMS_WRONG_NUM_ROWS, vars, caller);

    // and that each row has the correct length
    for (qindex r=0; r<op.numRows; r++) {
        qindex numCols = matrix[r].size();
        vars["${GIVEN_DIM}"] = numCols;
        assertThat(op.numRows == numCols, report::SUPER_OP_NEW_MATRIX_ELEMS_WRONG_NUM_COLS, vars, caller);
    }
}

void validate_superOpFieldsMatchPassedParams(SuperOp op, int numQb, const char* caller) {

    tokenSubs vars = {
        {"${NUM_PASSED_QUBITS}", numQb},
        {"${NUM_OP_QUBITS}",     op.numQubits}};

    assertThat(op.numQubits == numQb, report::SUPER_OP_FIELDS_MISMATCH_PARAMS, vars, caller);
}



/*
 * EXISTING SUPEROPERATORS
 */

void assertSuperOpFieldsAreValid(SuperOp op, bool isInKrausMap, const char* caller) {

    tokenSubs vars = {
        {"${NUM_QUBITS}", op.numQubits},
        {"${NUM_ROWS}",   op.numRows}
    };

    // assert valid fields
    auto msg = (isInKrausMap)? report::INVALID_KRAUS_MAPS_SUPER_OP_FIELDS : report::INVALID_SUPER_OP_FIELDS;
    assertThat(op.numQubits >= 1, msg, vars, caller);
    assertThat(op.numRows == powerOf2(2 * op.numQubits), msg, vars, caller);
    
    // only check outer point is allocated, to avoid inefficiently enumerating matrix rows
    msg = (isInKrausMap)? report::INVALID_KRAUS_MAPS_SUPER_OP_CPU_MEM_PTR : report::INVALID_SUPER_OP_CPU_MEM_PTR;
    assertThat(mem_isOuterAllocated(op.cpuElems), msg, caller);

    validate_envIsInit(caller);
    msg = (isInKrausMap)? report::INVALID_KRAUS_MAPS_SUPER_OP_GPU_MEM_PTR : report::INVALID_SUPER_OP_GPU_MEM_PTR;
    if (getQuESTEnv().isGpuAccelerated)
        assertThat(mem_isOuterAllocated(util_getGpuMemPtr(op)), msg, caller);

    // check the teeny tiny heap pointer is not NULL
    assertThat(mem_isAllocated(op.wasGpuSynced), report::INVALID_HEAP_FLAG_PTR, caller);

    // and that its value is a boolean
    int flag = *op.wasGpuSynced;
    tokenSubs moreVars = {{"${BAD_FLAG}", flag}, {"${UNKNOWN_FLAG}", validate_STRUCT_PROPERTY_UNKNOWN_FLAG}};
    assertThat(flag == 0 || flag == 1, report::INVALID_HEAP_FLAG_VALUE, moreVars, caller);

}

void validate_superOpFields(SuperOp op, const char* caller) {

    bool isInKrausMap = false;
    assertSuperOpFieldsAreValid(op, isInKrausMap, caller);
}

void validate_superOpIsSynced(SuperOp op, const char* caller) {

    // we don't need to perform any sync check in CPU-only mode
    if (!mem_isAllocated(util_getGpuMemPtr(op)))
        return;

    // check if GPU amps have EVER been overwritten; we sadly cannot check the LATEST changes were pushed though
    assertThat(*(op.wasGpuSynced), report::SUPER_OP_NOT_SYNCED_TO_GPU, caller);
}



/*
 * KRAUS MAP CREATION
 */

void assertKrausMapNonEmpty(int numQubits, const char* caller) {

    assertThat(numQubits >= 1, report::NEW_KRAUS_MAP_NUM_QUBITS_NOT_POSITIVE, {{"${NUM_QUBITS}",numQubits}}, caller);
}

void assertKrausMapValidNumMatrices(int numQubits, int numMatrices, const char* caller) {

    assertThat(
        numMatrices >= 1, report::KRAUS_MAP_NUM_GIVEN_NEW_MATRICES_NOT_POSITIVE, {{"${NUM_MATRICES}", numMatrices}}, caller);

    // this won't overflow given huge 'numQubits', because we prior validate the implied superoperator
    qindex maxNumBeforeIndOverflow = mem_getMaxNumKrausMapMatricesBeforeIndexOverflow(numQubits);

    assertThat(
        numMatrices <= maxNumBeforeIndOverflow, report::KRAUS_MAP_MATRICES_TOTAL_ELEMS_WOULD_EXCEED_QINDEX, 
        {{"${NUM_MATRICES}", numMatrices}, {"${NUM_QUBITS}", numQubits}, {"${MAX_NUM_MATRICES}", maxNumBeforeIndOverflow}}, caller);

    qindex maxNumBeforeMemOverflow = mem_getMaxNumKrausMapMatricesBeforeLocalMemSizeofOverflow(numQubits);
    assertThat(
        numMatrices <= maxNumBeforeMemOverflow, report::KRAUS_MAP_MATRICES_TOTAL_MEM_WOULD_EXCEED_SIZEOF, 
        {{"${NUM_MATRICES}", numMatrices}, {"${NUM_QUBITS}", numQubits}, {"${MAX_NUM_MATRICES}", maxNumBeforeMemOverflow}}, caller);
}

void validate_newKrausMapParams(int numQubits, int numMatrices, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    validate_envIsInit(caller);
    QuESTEnv env = getQuESTEnv();

    // first ensure that numQubits is >0 so below validation algebra is correct/safe
    assertKrausMapNonEmpty(numQubits, caller);
    
    // then ensure the Kraus map's superoperator has safe non-overflowing dimensions
    bool isInKrausMap = true;
    assertSuperOpTotalNumElemsDontExceedMaxIndex(numQubits, isInKrausMap, caller);
    assertSuperOpLocalMemDoesntExceedMaxSizeof(numQubits, isInKrausMap, caller);

    // ensure the superoperator can fit into current memory - we do NOT bother checking
    // whether the Kraus operator matrix list fits in memory since that's a very rare
    // scenario (the memory is linear with the user's numMatrices parameter, so they will 
    // not be asonisthed) and will be posteriori caught after memory allocation failure
    assertSuperOpFitsInCpuMem(numQubits, isInKrausMap, caller);
    assertSuperOpFitsInGpuMem(numQubits, env.isGpuAccelerated, isInKrausMap, caller);

    // ensure the number of Kraus operators isn't invalid, nor will cause overflows/seg-faults
    assertKrausMapValidNumMatrices(numQubits, numMatrices, caller);
}

void validate_newKrausMapAllocs(KrausMap map, const char* caller) {

    // we expensively get node consensus about malloc failure, in case of heterogeneous hardware/loads,
    // but we avoid this if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    // prior validation gaurantees this will not overflow
    qindex matrListMem = map.numMatrices * mem_getLocalMatrixMemoryRequired(map.numQubits, true, 1);
    tokenSubs vars = {
        {"${NUM_BYTES}", matrListMem},
        {"${NUM_MATRICES}", map.numMatrices},
        {"${NUM_QUBITS}", map.numQubits}};

    // assert the list of Kraus operator matrices, and all matrices nad rows therein, were allocated
    bool krausAreAlloc = mem_isAllocated(map.matrices, map.numMatrices, map.numQubits);
    assertAllNodesAgreeThat(krausAreAlloc, report::NEW_KRAUS_MAP_CPU_MATRICES_ALLOC_FAILED, vars, caller);

    // assert the teeny-tiny heap flag was alloc'd
    assertAllNodesAgreeThat(mem_isAllocated(map.isCPTP), report::NEW_HEAP_FLAG_ALLOC_FAILED, {{"${NUM_BYTES}", sizeof(*(map.isCPTP))}}, caller);

    // assert that the superoperator itself was allocated (along with its own heap fields)
    bool isInKrausMap = true;
    assertNewSuperOpAllocs(map.superop, isInKrausMap, caller);
}

void validate_newInlineKrausMapDimMatchesVectors(int numQubits, int numOperators, vector<vector<vector<qcomp>>> matrices, const char* caller) {

    // avoid potentially expensive matrix enumeration if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    qindex numRows = powerOf2(numQubits);
    
    assertThat(numOperators == (int) matrices.size(), report::NEW_INLINE_KRAUS_MAP_INCOMPATIBLE_NUM_NEW_MATRICES,
        {{"${NUM_GIVEN}", matrices.size()}, {"${NUM_EXPECTED}", numOperators}}, caller);

    // check each given matrix...
    for (int i=0; i<numOperators; i++) {
        
        // has a correct number of rows
        assertThat(numRows == (qindex) matrices[i].size(), report::NEW_INLINE_KRAUS_MAP_MATRIX_WRONG_NUM_ROWS, 
            {{"${NUM_QUBITS}", numQubits}, {"${NUM_EXPECTED_ROWS}", numRows}, {"${NUM_GIVEN_ROWS}", matrices[i].size()}}, caller);

        // and that each row has a correct number of elements/columns
        for (qindex r=0; r<numRows; r++)
            assertThat(numRows == (qindex) matrices[i][r].size(), report::NEW_INLINE_KRAUS_MAP_MATRIX_WRONG_NUM_COLS,
                {{"${NUM_QUBITS}", numQubits}, {"${NUM_EXPECTED_COLS}", numRows}, {"${NUM_GIVEN_COLS}", matrices[i][r].size()}}, caller);
    }
}



/*
 * KRAUS MAP INITIALISATION
 */

void validate_krausMapNewMatrixDims(KrausMap map, vector<vector<vector<qcomp>>> matrices, const char* caller) {

    // avoid potentially expensive matrix enumeration if validation is anyway disabled
    if (!isValidationEnabled)
        return;
    
    assertThat(map.numMatrices == (int) matrices.size(), report::KRAUS_MAP_INCOMPATIBLE_NUM_NEW_MATRICES,
        {{"${NUM_GIVEN}", matrices.size()}, {"${NUM_EXPECTED}", map.numMatrices}}, caller);

    // check each given matrix...
    for (int i=0; i<map.numMatrices; i++) {
        
        // has a correct number of rows
        assertThat(map.numRows == (qindex) matrices[i].size(), report::KRAUS_MAP_NEW_MATRIX_ELEMS_WRONG_NUM_ROWS, 
            {{"${NUM_QUBITS}", map.numQubits}, {"${NUM_EXPECTED_ROWS}", map.numRows}, {"${NUM_GIVEN_ROWS}", matrices[i].size()}}, caller);

        // and that each row has a correct number of elements/columns
        for (qindex r=0; r<map.numRows; r++)
            assertThat(map.numRows == (qindex) matrices[i][r].size(), report::KRAUS_MAP_NEW_MATRIX_ELEMS_WRONG_ROW_DIM,
                {{"${NUM_QUBITS}", map.numQubits}, {"${NUM_EXPECTED_COLS}", map.numRows}, {"${NUM_GIVEN_COLS}", matrices[i][r].size()}}, caller);
    }
}

void validate_krausMapFieldsMatchPassedParams(KrausMap map, int numQb, int numOps, const char* caller) {

    tokenSubs vars = {
        {"${NUM_MAP_QUBITS}",    map.numQubits},
        {"${NUM_MAP_OPS}",       map.numMatrices},
        {"${NUM_PASSED_QUBITS}", numQb},
        {"${NUM_PASSED_OPS}",    numOps}};

    bool valid = (map.numQubits == numQb) && (map.numMatrices == numOps);
    assertThat(valid, report::KRAUS_MAP_FIELDS_MISMATCH_PARAMS, vars, caller);
}



/*
 * EXISTING KRAUS MAPS
 */

void validate_krausMapFields(KrausMap map, const char* caller) {

    tokenSubs vars = {
        {"${NUM_QUBITS}",   map.numQubits},
        {"${NUM_MATRICES}", map.numMatrices},
        {"${NUM_ROWS}",     map.numRows}};

    // assert valid fields
    auto msg = report::INVALID_KRAUS_MAP_FIELDS;
    assertThat(map.numQubits >= 1, msg, vars, caller);
    assertThat(map.numMatrices >= 1, msg, vars, caller);
    assertThat(map.numRows == powerOf2(map.numQubits), msg, vars, caller);

    // check only outer CPU matrix list is allocated, to avoid expensive enumerating of matrices/rows
    assertThat(mem_isOuterAllocated(map.matrices), report::INVALID_KRAUS_MAP_MATRIX_LIST_MEM_PTR, caller);

    // assert isCPTP heap flag allocated, and that is has a valid value
    assertThat(mem_isAllocated(map.isCPTP), report::INVALID_HEAP_FLAG_PTR, caller);

    // and that its value is a boolean
    int flag = *map.isCPTP;
    bool valid = flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG;
    tokenSubs moreVars = {{"${BAD_FLAG}", flag}, {"${UNKNOWN_FLAG}", validate_STRUCT_PROPERTY_UNKNOWN_FLAG}};
    assertThat(valid, report::INVALID_HEAP_FLAG_VALUE, moreVars, caller);

    // assert the superoperator dimension matches the map's
    assertThat(map.numQubits == map.superop.numQubits, report::INVALID_KRAUS_MAP_SUPER_OP_NUM_QUBITS, 
        {{"${MAP_QUBITS}", map.numQubits}, {"${SUPER_OP_QUBITS}", map.superop.numQubits}}, caller);

    // assert the superoperator fields are valid (e.g. dims match, ptrs are valid)
    bool isInKrausMap = true;
    assertSuperOpFieldsAreValid(map.superop, isInKrausMap, caller);
}

void validate_krausMapIsSynced(KrausMap map, const char* caller) {

    // we don't need to perform any sync check in CPU-only mode
    if (!mem_isAllocated(util_getGpuMemPtr(map.superop)))
        return;

    // assert the map's superoperator has been synced
    assertThat(*(map.superop.wasGpuSynced), report::KRAUS_MAP_NOT_SYNCED_TO_GPU, caller);
}

void validate_krausMapIsCPTP(KrausMap map, const char* caller) {
    validate_krausMapFields(map, caller);
    validate_krausMapIsSynced(map, caller);

    // avoid expensive CPTP check if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    // evaluate CPTPness if it isn't already known 
    if (*(map.isCPTP) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(map.isCPTP) = util_isCPTP(map, validationEpsilon);

    assertThat(*(map.isCPTP), report::KRAUS_MAP_NOT_CPTP, caller);
}

void validate_krausMapMatchesTargets(KrausMap map, int numTargets, const char* caller) {

    tokenSubs vars = {{"${KRAUS_QUBITS}", map.numQubits}, {"${TARG_QUBITS}", numTargets}};
    assertThat(map.numQubits == numTargets, report::KRAUS_MAP_SIZE_MISMATCHES_TARGETS, vars, caller);
}



/*
 * PAULI STRING CREATION
 */

void assertCorrectNumPauliCharsBeforeTerminationChar(const char* paulis, int numPaulis, const char* caller) {

    char termChar = '\0';
    int numCharsBeforeTerm = 0;

    // the termination char can actually be AFTER numPaulis; that's fine
    for (int i=0; i<numPaulis && paulis[i] != termChar; i++)
        numCharsBeforeTerm++;

    tokenSubs vars = {{"${TERM_IND}", numCharsBeforeTerm}, {"${NUM_PAULIS}", numPaulis}};
    assertThat(numCharsBeforeTerm >= numPaulis, report::NEW_PAULI_STR_TERMINATION_CHAR_TOO_EARLY, vars, caller);
}

void assertRecognisedNewPaulis(const char* paulis, int numPaulis, const char* caller) {

    // paulis might also contain '\0' char (string termination),  
    // but not before numPaulis (as prior validated)

    for (int i=0; i<numPaulis; i++) {

        // TODO: we can only display the ascii code of unrecognised characters,
        // because tokenSubs only accepts integers (not chars/substrings). Fix this!
        char ch = paulis[i];
        int ascii = (int) ch;

        assertThat(
            parser_RECOGNISED_PAULI_CHARS.find(ch) != string::npos,
            report::NEW_PAULI_STR_UNRECOGNISED_PAULI_CHAR,
            {{"${BAD_CHAR}", ascii}, {"${CHAR_IND}", i}}, caller);
    }
}

void assertValidNewPauliCodes(int* paulis, int numPaulis, const char* caller) {

    for (int i=0; i<numPaulis; i++) {
        int code = paulis[i];

        assertThat(
            code>=0 && code<=3, 
            report::NEW_PAULI_STR_INVALID_PAULI_CODE, 
            {{"${BAD_CODE}", code}, {"${CODE_IND}", i}}, caller);
    }
}

void assertValidNewPauliIndices(int* indices, int numInds, int maxIndExcl, const char* caller) {

    // check each index is valid
    for (int i=0; i<numInds; i++) {
        int ind = indices[i];
        assertThat(
            ind >= 0 && ind < maxIndExcl, report::NEW_PAULI_STR_INVALID_INDEX,
            {{"${BAD_IND}", ind}, {"${MAX_PAULIS}", maxIndExcl}}, caller);
    }

    // check no index is duplicated
    assertThat(isIndexListUnique(indices, numInds), report::NEW_PAULI_STR_DUPLICATED_INDEX, caller);
}

void validate_newPauliStrNumPaulis(int numPaulis, int maxNumPaulis, const char* caller) {

    tokenSubs vars = {{"${NUM_PAULIS}", numPaulis}};
    assertThat(numPaulis > 0, report::NEW_PAULI_STR_NON_POSITIVE_NUM_PAULIS, vars, caller);

    vars["${MAX_PAULIS}"] = maxNumPaulis;
    assertThat(numPaulis <= maxNumPaulis, report::NEW_PAULI_STR_NUM_PAULIS_EXCEEDS_TYPE, vars, caller);
}

void validate_newPauliStrParams(const char* paulis, int* indices, int numPaulis, int maxNumPaulis, const char* caller) {

    validate_newPauliStrNumPaulis(numPaulis, maxNumPaulis, caller);
    assertCorrectNumPauliCharsBeforeTerminationChar(paulis, numPaulis, caller);
    assertRecognisedNewPaulis(paulis, numPaulis, caller);
    assertValidNewPauliIndices(indices, numPaulis, maxNumPaulis, caller);
}
void validate_newPauliStrParams(int* paulis, int* indices, int numPaulis, int maxNumPaulis, const char* caller) {

    validate_newPauliStrNumPaulis(numPaulis, maxNumPaulis, caller);
    assertValidNewPauliCodes(paulis, numPaulis, caller);
    assertValidNewPauliIndices(indices, numPaulis, maxNumPaulis, caller);
}

void validate_newPauliStrNumChars(int numPaulis, int numIndices, const char* caller) {

    // this is a C++-only validation, because only std::string gaurantees we can know
    // the passed string length (C char arrays might not contain termination char)
    tokenSubs vars = {{"${NUM_PAULIS}", numPaulis}, {"${NUM_INDS}", numIndices}};
    assertThat(numPaulis == numIndices, report::NEW_PAULI_STR_DIFFERENT_NUM_CHARS_AND_INDS, vars, caller);
}



/*
 * PAULI STRING SUM CREATION
 */

void validate_newPauliStrSumParams(qindex numTerms, const char* caller) {

    // note we do not bother checking whether RAM has enough memory to contain
    // the new Pauli sum, because the caller to this function has already
    // been passed data of the same size (and it's unlikely the user is about
    // to max RAM), and the memory requirements scale only linearly with the
    // parameters (e.g. numTerms), unlike the exponential scaling of the memory
    // of Qureg and CompMatr, for example

    assertThat(numTerms > 0, report::NEW_PAULI_STR_SUM_NON_POSITIVE_NUM_STRINGS, {{"${NUM_TERMS}", numTerms}}, caller);
}

void validate_newPauliStrSumMatchingListLens(qindex numStrs, qindex numCoeffs, const char* caller) {

    tokenSubs vars = {{"${NUM_STRS}", numStrs}, {"${NUM_COEFFS}", numCoeffs}};
    assertThat(numStrs == numCoeffs, report::NEW_PAULI_STR_SUM_DIFFERENT_NUM_STRINGS_AND_COEFFS, vars, caller);
}

void validate_newPauliStrSumAllocs(PauliStrSum sum, qindex numBytesStrings, qindex numBytesCoeffs, const char* caller) {

    assertThat(
        mem_isAllocated(sum.strings), report::NEW_PAULI_STR_SUM_STRINGS_ALLOC_FAILED, 
        {{"${NUM_TERMS}", sum.numTerms}, {"${NUM_BYTES}", numBytesStrings}}, caller);

    assertThat(
        mem_isAllocated(sum.coeffs), report::NEW_PAULI_STR_SUM_COEFFS_ALLOC_FAILED, 
        {{"${NUM_TERMS}", sum.numTerms}, {"${NUM_BYTES}", numBytesCoeffs}}, caller);

    assertThat(
        mem_isAllocated(sum.isHermitian), report::NEW_HEAP_FLAG_ALLOC_FAILED, 
        {{"${NUM_BYTES}", sizeof(*(sum.isHermitian))}}, caller);
}



/*
 * PAULI STRING SUM PARSING
 */

void validate_parsedPauliStrSumLineIsInterpretable(bool isInterpretable, string line, qindex lineIndex, const char* caller) {

    // TODO: we cannot yet report 'line' because tokenSubs so far only accepts integers :(

    tokenSubs vars = {{"${LINE_NUMBER}", lineIndex + 1}}; // line numbers begin at 1
    assertThat(isInterpretable, report::PARSED_PAULI_STR_SUM_UNINTERPRETABLE_LINE, vars, caller);
}

void validate_parsedPauliStrSumLineHasConsistentNumPaulis(int numPaulis, int numLinePaulis, string line, qindex lineIndex, const char* caller) {

    // TODO: we cannot yet report 'line' because tokenSubs so far only accepts integers :(

    tokenSubs vars = {
        {"${NUM_PAULIS}",      numPaulis},
        {"${NUM_LINE_PAULIS}", numLinePaulis},
        {"${LINE_NUMBER}",      lineIndex + 1}}; // line numbers begin at 1
    assertThat(numPaulis == numLinePaulis, report::PARSED_PAULI_STR_SUM_INCONSISTENT_NUM_PAULIS_IN_LINE, vars, caller);
}

void validate_parsedPauliStrSumCoeffIsValid(bool isCoeffValid, string line, qindex lineIndex, const char* caller) {

    // TODO: we cannot yet report 'line' because tokenSubs so far only accepts integers :(

    tokenSubs vars = {{"${LINE_NUMBER}", lineIndex + 1}}; // lines begin at 1
    assertThat(isCoeffValid, report::PARSED_PAULI_STR_SUM_COEFF_IS_INVALID, vars, caller);
}

void validate_parsedStringIsNotEmpty(bool stringIsNotEmpty, const char* caller) {

    assertThat(stringIsNotEmpty, report::PARSED_STRING_IS_EMPTY, caller);
}



/*
 * EXISTING PAULI STRING SUMS
 */

void validate_pauliStrSumFields(PauliStrSum sum, const char* caller) {

    assertThat(sum.numTerms > 0, report::INVALID_PAULI_STR_SUM_FIELDS, {{"${NUM_TERMS}", sum.numTerms}}, caller);

    assertThat(mem_isAllocated(sum.coeffs),  report::INVALID_PAULI_STR_HEAP_PTR, caller);
    assertThat(mem_isAllocated(sum.strings), report::INVALID_PAULI_STR_HEAP_PTR, caller);

    assertThat(mem_isAllocated(sum.isHermitian), report::INVALID_HEAP_FLAG_PTR, caller);

    // assert isHermitian has valid value
    int flag = *sum.isHermitian;
    tokenSubs vars = {
        {"${BAD_FLAG}", flag}, 
        {"${UNKNOWN_FLAG}", validate_STRUCT_PROPERTY_UNKNOWN_FLAG}};
    assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);
}

void valdidate_pauliStrSumIsHermitian(PauliStrSum sum, const char* caller) {

    // ensure hermiticity is known (if not; compute it)
    if (*(sum.isHermitian) == validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        *(sum.isHermitian) = util_isHermitian(sum, validationEpsilon);

    assertThat(*(sum.isHermitian), report::PAULI_STR_SUM_NOT_HERMITIAN, caller);
}



/*
 * BASIS STATE INDICES
 */

void validate_basisStateIndex(Qureg qureg, qindex ind, const char* caller) {

    qindex maxIndExcl = powerOf2(qureg.numQubits);

    tokenSubs vars = {
        {"${STATE_IND}",  ind},
        {"${NUM_QUBITS}", qureg.numQubits},
        {"${NUM_STATES}", maxIndExcl}};

    assertThat(ind >= 0 && ind < maxIndExcl, report::INVALID_BASIS_STATE_INDEX, vars, caller);
}

void validate_basisStateRowCol(Qureg qureg, qindex row, qindex col, const char* caller) {

    qindex maxIndExcl = powerOf2(qureg.numQubits);

    tokenSubs vars = {
        {"${ROW_IND}",  row},
        {"${COL_IND}",  col},
        {"${NUM_QUBITS}", qureg.numQubits},
        {"${NUM_STATES}", maxIndExcl}};

    bool valid = row >= 0 && row < maxIndExcl && col >= 0 && col < maxIndExcl;
    assertThat(valid, report::INVALID_BASIS_STATE_ROW_OR_COL, vars, caller);
}



/*
 * QUBIT INDICES
 */

bool areQubitsUnique(int* qubits, int numQubits) {

    // assumes all elemtns of qubits are < 64
    qindex mask = 0;

    for (int n=0; n<numQubits; n++)
        if (getBit(mask, qubits[n]))
            return false;
        else
            mask = setBit(mask, qubits[n], 1);

    return true;
}

bool areQubitsDisjoint(int* qubitsA, int numQubitsA, int* qubitsB, int numQubitsB) {

    // assumes all elemtns of qubits are < 64
    qindex maskA = getBitMask(qubitsA, numQubitsA);

    for (int n=0; n<numQubitsB; n++)
        if (getBit(maskA, qubitsB[n]))
            return false;
    
    return true;
}

void assertValidQubit(Qureg qureg, int qubitInd, string msg, const char* caller) {

    tokenSubs vars = {
        {"${QUBIT_IND}",  qubitInd},
        {"${NUM_QUBITS}", qureg.numQubits}};

    assertThat(qubitInd >= 0 && qubitInd < qureg.numQubits, msg, vars, caller);
}

void assertValidQubits(
    Qureg qureg, int* qubits, int numQubits, bool canNumIncludeZero, const char* caller,
    string msgNegNum, string msgNumExceedsQureg, string msgNullPtr, string msgBadInd, string msgDuplicates
) {
    assertThat(numQubits > (canNumIncludeZero? -1 : 0), msgNegNum, {{"${NUM_QUBITS}", numQubits}}, caller);
    assertThat(numQubits <= qureg.numQubits, msgNumExceedsQureg, {{"${NUM_QUBITS}", numQubits}, {"${QUREG_QUBITS}", qureg.numQubits}}, caller);

    if (numQubits > 0)
        assertThat(qubits != nullptr, msgNullPtr, caller);

    for (int n=0; n<numQubits; n++)
        assertValidQubit(qureg, qubits[n], msgBadInd, caller);
    
    assertThat(areQubitsUnique(qubits, numQubits), msgDuplicates, caller);
}

void validate_target(Qureg qureg, int target, const char* caller) {

    assertValidQubit(qureg, target, report::INVALID_TARGET_QUBIT, caller);
}

void validate_targets(Qureg qureg, int* targets, int numTargets, const char* caller) {

    // must always have at least 1 target
    bool numCanBeZero = false;

    assertValidQubits(
        qureg, targets, numTargets, numCanBeZero, caller, 
        report::NEGATIVE_OR_ZERO_NUM_TARGETS, report::NUM_TARGETS_EXCEEDS_QUREG_SIZE, report::TARGET_LIST_WAS_NULL_PTR,
        report::INVALID_TARGET_QUBIT,         report::DUPLICATE_TARGET_QUBITS);
}
/*
 * MEASUREMENT PARAMETERS
 */

void validate_measurementOutcomeIsValid(int outcome, const char* caller) {

    assertThat(outcome == 0 || outcome == 1, report::ONE_QUBIT_MEASUREMENT_OUTCOME_INVALID, {{"${OUTCOME}", outcome}}, caller);
}

void validate_measurementOutcomesAreValid(int* outcomes, int numOutcomes, const char* caller) {

    // no need to validate numOutcomes; it is already validated by caller (e.g. through numTargets)

    for (int i=0; i<numOutcomes; i++)
        assertThat(
            outcomes[i] == 0 || outcomes[i]== 1, report::MANY_QUBIT_MEASUREMENTS_OUTCOME_INVALID, 
            {{"${INDEX}", i}, {"${OUTCOME}", outcomes[i]}}, caller);
}

void validate_measurementOutcomesFitInGpuMem(Qureg qureg, int numQubits, const char* caller) {

    // only GPU backend needs temp memory
    if (!qureg.isGpuAccelerated)
        return;

    qindex numOutcomes = powerOf2(numQubits);
    size_t memAvail = gpu_getCurrentAvailableMemoryInBytes();
    size_t memNeeded = sizeof(qreal) * numOutcomes;

    tokenSubs vars = {
        {"{NUM_QUBITS}",    numQubits}, 
        {"${NUM_OUTCOMES}", numOutcomes}, 
        {"${MEM_NEEDED}",   memNeeded}, 
        {"${MEM_AVAIL}",    memAvail}
    };

    assertThat(memAvail > memNeeded, report::GPU_CANNOT_FIT_TEMP_MEASUREMENT_OUTCOME_PROBS, vars, caller);
}





/*
 * FILE IO
 */

void validate_canReadFile(string fn, const char* caller) {

    // TODO: embed filename into error message when tokenSubs is updated to permit strings
    assertThat(parser_canReadFile(fn), report::CANNOT_READ_FILE, caller);
}
