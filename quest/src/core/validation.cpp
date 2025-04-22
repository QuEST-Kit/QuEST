/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 * 
 * @author Tyson Jones
 * @author Richard Meister (patched v3 mem-leak and MSVC build)
 * @author Kshitij Chhabra (patched v3 overflow bug)
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

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>

using std::string;
using std::vector;



/**
 * @todo
 *   Invalid input error messages are currently limited to containing only integer
 *   variables, whereas the ability to embed floats and strings would be useful.
 */



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

    string MULTIPLE_NODES_BOUND_TO_SAME_GPU =
        "Multiple MPI processes (nodes) were bound to the same GPU which is detrimental to performance and almost never intended. Please re-deploy QuEST with no more MPI processes than there are total GPUs. Alternatively, recompile QuEST with macro PERMIT_NODES_TO_SHARE_GPU=1.";

    string CUQUANTUM_DEPLOYED_ON_BELOW_CC_GPU =
        "Cannot use cuQuantum on a GPU with compute-capability ${OUR_CC}; a compute-capability of ${MIN_CC} or above is required. Recompile with cuQuantum disabled to fall-back to using Thrust and custom kernels.";

    string CUQUANTUM_DEPLOYED_ON_GPU_WITHOUT_MEM_POOLS =
        "Cannot use cuQuantum since your GPU does not support memory pools. Recompile with cuQuantum disabled to fall-back to using Thrust and custom kernels.";

    
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

    string INVALID_NUM_REPORTED_SCALARS =
        "Invalid parameter (${NUM_ITEMS}). Must specify a positive number of scalars to be reported, or 0 to indicate that all scalars should be reported.";

    string INVALID_NUM_REPORTED_SIG_FIGS =
        "Invalid number of significant figures (${NUM_SIG_FIGS}). Cannot be less than one.";

    string INVALID_NUM_RANDOM_SEEDS =
        "Invalid number of random seeds (${NUM_SEEDS}). Must specify one or more. In distributed settings, only the root node needs to pass a valid number of seeds (other node arguments are ignored).";
    
    string INVALID_NUM_REPORTED_NEWLINES =
        "Invalid number of trailing newlines (${NUM_NEWLINES}). Cannot generally be less than zero, and must not be zero when calling multi-line reporting functions like reportQureg().";

    string INSUFFICIENT_NUM_REPORTED_NEWLINES =
        "The number of trailing newlines (set by setNumReportedNewlines()) is zero which is not permitted when calling multi-line reporters.";

    string INVALID_NUM_NEW_PAULI_CHARS =
        "Given an invalid number of Pauli characters. Must specify precisely four to respectively replace IXYZ.";

    string INVALID_REPORTED_PAULI_STR_STYLE_FLAG =
        "Given an unrecognised style flag (${FLAG}). Legal flags are 0 and 1.";


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
        "Cannot create Qureg of ${NUM_QUBITS} qubits distributed over ${NUM_NODES} nodes because the size (in bytes) of the necessary global memory would overflow size_t. In this deployment, the maximum number of qubits in a statevector Qureg is ${MAX_QUBITS}. See reportQuESTEnv().";

    string NEW_DENSMATR_QUREG_LOCAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create density Qureg of ${NUM_QUBITS} qubits distributed over ${NUM_NODES} nodes because the size (in bytes) of the necessary global memory would overflow size_t. In this deployment, the maximum number of qubits in a density-matrix Qureg is ${MAX_QUBITS}. See reportQuESTEnv().";

    string NEW_DISTRIB_STATEVEC_QUREG_HAS_TOO_FEW_AMPS = 
        "Cannot create a distributed Qureg of only ${NUM_QUBITS} qubits between ${NUM_NODES} nodes, because each node would contain fewer than one amplitude of the statevector. The minimum size in this deployment is ${MIN_QUBITS} qubits; see reportQuESTEnv(). Consider disabling distribution for this Qureg.";

    string NEW_DISTRIB_DENSMATR_QUREG_HAS_TOO_FEW_AMPS = 
        "Cannot create a distributed density Qureg of only ${NUM_QUBITS} qubits between ${NUM_NODES} nodes, because each node would contain fewer than a column's worth of amplitudes of the density matrix. The minimum size in this deployment is ${MIN_QUBITS} qubits; see reportQuESTEnv(). Consider disabling distribution for this Qureg.";


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
        "Received an invalid Qureg, which had invalid or incompatible fields isDensityMatrix=${DENS_MATR}, numQubits=${NUM_QUBITS}, numAmps=${NUM_AMPS}. It is likely this Qureg was not initialised with createQureg().";

    string QUREG_NOT_DENSITY_MATRIX =
        "Expected a density matrix Qureg but received a statevector.";

    string QUREG_NOT_STATE_VECTOR =
        "Expected a statevector Qureg but received a density matrix.";


    /*
     * MUTABLE OBJECT FLAGS
     */

    string NEW_HEAP_FLAG_ALLOC_FAILED =
        "Attempted allocation of a heap flag (such as 'isApproxUnitary') miraculously failed, despite being a mere ${NUM_BYTES} bytes. This is unfathomably unlikely - go and have your fortune read at once!";

    string INVALID_HEAP_FLAG_PTR =
        "A flag (such as 'isApproxUnitary') bound to the given operator data structure (e.g. matrix, superoperator, Kraus map) was a NULL pointer, instead of an expected pointer to persistent heap memory. This may imply the structure was not created by its proper function (e.g. createCompMatr()).";

    string INVALID_HEAP_FLAG_VALUE = 
        "A flag (such as 'isApproxUnitary') bound to the given operator data structure (e.g. matrix, superoperator, Kraus map) had an invalid value of ${BAD_FLAG}. Allowed values are '0', '1', and (except for 'wasGpuSynced') '${UNKNOWN_FLAG}', though these flags should not be modified directly by the user. It is likely the structure was not created by its proper function (e.g. createCompMatr()).";


    /*
     * MATRIX CREATION
     */

    string INVALID_OPTION_FOR_MATRIX_IS_DISTRIB = 
        "Argument 'useDistrib' must be 1 or 0 to respectively indicate whether or not to distribute the new FullStateDiagMatr, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    string INVALID_OPTION_FOR_MATRIX_IS_GPU_ACCEL = 
        "Argument 'useGpuAccel' must be 1 or 0 to respetively indicate whether or not to GPU-accelerate the new FullStateDiagMatr, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";

    string INVALID_OPTION_FOR_MATRIX_IS_MULTITHREAD = 
        "Argument 'useMultithread' must be 1 or 0 to respectively indicate whether or not to use multithreading when initialising the new FullStateDiagMatr, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";


    string NEW_DISTRIB_MATRIX_IN_NON_DISTRIB_ENV =
        "Cannot distribute a FullStateDiagMatr when in a non-distributed QuEST environment.";

    string NEW_GPU_MATRIX_IN_NON_GPU_ACCEL_ENV =
        "Cannot allocate a FullStateDiagMatr to a GPU when in a non-GPU-accelerated QuEST environment.";

    string NEW_MULTITHREAD_MATRIX_IN_NON_MULTITHREAD_ENV =
        "Cannot enable multithreaded processing of a FullStateDiagMatr created in a non-multithreaded QuEST environment.";

    
    string NEW_MATRIX_NUM_QUBITS_NOT_POSITIVE = 
        "Cannot create a matrix which acts upon ${NUM_QUBITS} qubits; must target one or more qubits.";


    string NEW_DIAG_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX =
        "Cannot create a diagonal matrix of ${NUM_QUBITS} qubits: the matrix would contain more elements (2^${NUM_QUBITS}) than the maximum which can be addressed by the qindex type (2^${MAX_QUBITS}).";

    string NEW_COMP_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX =
        "Cannot create a dense matrix of ${NUM_QUBITS} qubits: the matrix would contain more elements (4^${NUM_QUBITS}) than the maximum which can be addressed by the qindex type (4^${MAX_QUBITS}).";


    string NEW_LOCAL_COMP_MATR_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a local, dense matrix of ${NUM_QUBITS} qubits because the size (in bytes) of the necessary memory would overflow size_t. In this deployment, the maximum number of qubits in such a matrix is ${MAX_QUBITS}.";

    string NEW_LOCAL_DIAG_MATR_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a local, diagonal matrix of ${NUM_QUBITS} qubits because the size (in bytes) of the necessary memory would overflow size_t. In this deployment, the maximum number of qubits in such a matrix is ${MAX_QUBITS}.";

    string NEW_DISTRIB_DIAG_MATR_LOCAL_MEM_WOULD_EXCEED_SIZEOF =
        "Cannot create a diagonal matrix of ${NUM_QUBITS} qubits distributed over ${NUM_NODES} nodes because the size (in bytes) of the necessary global memory would overflow size_t. In this deployment, the maximum number of qubits in such a matrix is ${MAX_QUBITS}.";


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
        "Attempted allocation of the matrix memory (a total of ${NUM_BYTES} bytes in RAM) failed.";

    string NEW_MATRIX_GPU_ELEMS_ALLOC_FAILED = 
        "Attempted allocation of the matrix's GPU memory (${NUM_BYTES} bytes in VRAM) failed.";


    /*
     * MATRIX INITIALISATION
     */

    string COMP_MATR_NEW_ELEMS_WRONG_NUM_ROWS =
        "Incompatible number of rows (${NUM_GIVEN_ROWS}) of elements given to overwrite a ${NUM_QUBITS}-qubit CompMatr, which expects ${NUM_EXPECTED_ROWS} rows.";

    string COMP_MATR_NEW_ELEMS_WRONG_ROW_DIM =
        "One or more rows contained an incompatible number of elements (${NUM_GIVEN_ELEMS}). The ${NUM_QUBITS}-qubit CompMatr expects a square ${EXPECTED_DIM}x${EXPECTED_DIM} matrix.";

    string DIAG_MATR_WRONG_NUM_NEW_ELEMS = 
        "Incompatible number of elements (${NUM_GIVEN_ELEMS}) assigned to a ${NUM_QUBITS}-qubit DiagMatr, which expects ${NUM_EXPECTED_ELEMS} elements.";

    string DIAG_MATR_NEW_ELEMS_NULL_PTR =
        "The given list of new elements was a null pointer.";

    string DENSE_MATR_NEW_ELEMS_OUTER_NULL_PTR =
        "The given matrix of new elements was a null pointer.";

    string DENSE_MATR_NEW_ELEMS_INNER_NULL_PTR =
        "The given matrix of new elements contained a null pointer, suggesting an issue with initialising one or more rows.";
    

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


    string NESTED_VECTOR_MATRIX_HAS_INCONSISTENT_NUM_COLUMNS =
        "The given matrix (a nested list of vectors) had an inconsistent number of columns. The first row had ${FIRST_LEN} columns, while the row at index ${OUTLIER_IND} had ${OUTLIER_LEN} columns.";


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


    string MATRIX_SIZE_MISMATCHES_NUM_TARGETS =
        "The given matrix has an inconsistent size (${MATRIX_NUM_QUBITS}) with the specified number of target qubits (${NUM_TARGS}).";


    string MATRIX_NOT_UNITARY = 
        "The given matrix was not (approximately) unitary.";
    
    string UNITARY_DIAG_MATR_EXPONENT_NOT_APPROX_REAL =
        "The given exponent was not approximately real (i.e. the imaginary component exceeded epsilon) such that the given diagonal matrix raised to the exponent was no longer approximately unitary. Consider changing the validation epsilon.";

    string MATRIX_NOT_HERMITIAN =
        "The given matrix was not (approximately) Hermitian.";

    string DIAG_MATR_APPROX_ZERO_WHILE_EXPONENT_REAL_AND_NEGATIVE =
        "The given exponent was (real and) negative while one or more elements of the diagonal matrix had magnitudes near (within epsilon) to zero, which would cause divergences or divison-by-zero errors.";

    string HERMITIAN_DIAG_MATR_NEGATIVE_WHILE_EXPONENT_NOT_INTEGER =
        "The given exponent was a non-integer while one or more real components of the given Hermitian diagonal matrix were negative. The exponentiated matrix would ergo contain complex elements which violate Hermiticity. This validation is strict and not affected by the validation epsilon since it is necessary to avoid NaN output from a strictly-real pow() function, used for its improved numerical accuracy over complex pow().";


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
        "Attempted allocation of memory for the superoperator matrix (a total of ${NUM_BYTES} bytes in RAM) failed.";
        
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


    string SUPER_OP_SIZE_MISMATCHES_NUM_TARGETS =
        "The SuperOp dimension (${OP_QUBITS} qubits) is inconsistent with the specified number of target qubits (${NUM_TARGS}).";


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
        "Cannot create a Kraus map of ${NUM_QUBITS} qubits because the size (in bytes) of the necessary memory for its corresponding superoperator (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) would overflow size_t. In this deployment, the maximum number of qubits in a Kraus map is ${MAX_QUBITS}.";


    string NEW_KRAUS_MAPS_SUPER_OP_CANNOT_FIT_INTO_CPU_MEM =
        "Cannot create a Kraus map of ${NUM_QUBITS} qubits because the total memory required by its corresponding superoperator (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) exceeds the available memory (${RAM_SIZE} bytes).";
        
    string NEW_KRAUS_MAPS_SUPER_OP_CANNOT_FIT_INTO_GPU_MEM =
        "Cannot create a GPU-accelerated Kraus map of ${NUM_QUBITS} qubits because the total memory required by its corresponding superoperator (${QCOMP_BYTES} * 16^${NUM_QUBITS} bytes) exceeds the available GPU memory (${VRAM_SIZE} bytes).";


    string NEW_KRAUS_MAPS_SUPER_OP_CPU_ELEMS_ALLOC_FAILED =
        "Attempted allocation of memory for the Kraus map's correspondng superoperator matrix (a total of ${NUM_BYTES} bytes in RAM) failed.";

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
        "The 'isApproxCPTP' field of the given KrausMap was a NULL pointer, instead of the expected pointer to persistent heap memory. This suggests the KrausMap was already destroyed and had its fields overwritten by the user.";

    string INVALID_KRAUS_MAP_IS_CPTP_FLAG = 
        "The 'isApproxCPTP' field of the given Kraus map had invalid value ${BAD_FLAG}, suggesting it was manually modified. Valid values are 0, 1 and ${UNKNOWN_FLAG} (to indicate that unitarity is not yet known, deferring evaluation) although this need never be modified by the user.";


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
    
    /// @todo replace BAD_CHAR ascii code with actual character, once tokenSubs is generalised to any-type
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

    string NEW_PAULI_STR_DIFFERENT_NUM_CHARS_AND_INDICES = 
        "Given a different number of Pauli operators (${NUM_PAULIS}) and their qubit indices (${NUM_INDS}).";


    /*
     * EXISTING PAULI STRING
     */

    string PAULI_STR_TARGETS_EXCEED_QUREG =
        "The highest-index non-identity Pauli operator (index ${MAX_TARG}) in the given PauliStr exceeds the maximum target of the ${QUREG_QUBITS}-qubit Qureg.";

    string PAULI_STR_OVERLAPS_CONTROLS =
        "A control qubit overlaps a non-identity Pauli operator in the given PauliStr.";


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
        "The given PauliStrSum is not Hermitian.";

    string PAULI_STR_SUM_NOT_ALL_I_Z =
        "The given PauliStrSum contained X and/or Y operators and is ergo not diagonal, so cannot be used to initialise the given FullStateDiagMatr.";


    string PAULI_STR_SUM_EXCEEDS_QUREG_NUM_QUBITS =
        "The given PauliStrSum includes non-identity upon a qubit of index ${MAX_IND} and so is only compatible with Quregs of at least ${NUM_PSS_QUBITS} qubits. It cannot act upon, nor initialise, the given ${NUM_QUREG_QUBITS}-qubit Qureg.";

    string PAULI_STR_SUM_EXCEEDS_MATR_NUM_QUBITS =
        "The given PauliStrSum includes non-identity upon a qubit of index ${MAX_IND} and so is only compatible with FullStateDiagMatr containing at least ${NUM_PSS_QUBITS}. It cannot initialise the given ${NUM_MATR_QUBITS}-qubit FullStateDiagMatr.";


    /*
     * BASIS STATE INDICES
     */

    string INVALID_BASIS_STATE_INDEX = 
        "Basis state index ${STATE_IND} is invalid for the given ${NUM_QB} qubit Qureg. Index must be greater than or equal to zero, and cannot equal nor exceed the number of unique classical states (2^${NUM_QB} = ${NUM_STATES}).";

    string INVALID_BASIS_STATE_ROW_OR_COL =
        "The row and column indices (${ROW_IND}, ${COL_IND}) are invalid for the given ${NUM_QB} qubit Qureg. Both indices must be greater than or equal to zero, and neither can equal nor exceed the number of unique classical states (2^${NUM_QB} = ${NUM_STATES}).";


    string INVALID_STARTING_BASIS_STATE_INDEX =
        "The starting basis state index ${START_IND} is invalid for the given ${NUM_QB} qubit Qureg. Index must be greater than or equal to zero, and cannot equal nor exceed the number of unique classical states (2^${NUM_QB} = ${MAX_IND_EXCL}).";

    string INVALID_NUM_BASIS_STATE_INDICES =
        "The specified number of amplitudes or basis states (${NUM_INDS}) is invalid for the given ${NUM_QB} qubit Qureg, which contains ${MAX_NUM_INDS_INCL} amplitudes.";

    string INVALID_ENDING_BASIS_STATE_INDEX =
        "The combination of the starting basis state index (${START_IND}) and the number of amplitudes (${NUM_INDS}) implies an end index of ${END_IND_EXCL} (exclusive) which is invalid for the ${NUM_QB}-qubit ${MAX_IND_EXCL}-amplitude Qureg.";


    string INVALID_STARTING_BASIS_ROW_OR_COL =
        "Either or both of the starting row and column (${START_ROW} and ${START_COL} respectively) are invalid for the given ${NUM_QB} qubit Qureg. Both indices must be greater than or equal to zero, and cannot exceed equal nor exceed the density matrix dimension of 2^${NUM_QB} = ${MAX_IND_EXCL}.";

    string INVALID_NUM_BASIS_ROWS_OR_COLS =
        "Either or both of the number of rows and columns (${NUM_ROWS} and ${NUM_COLS} respectively) are invalid for the given ${NUM_QB} qubit Qureg, which has a dimension of ${MAX_NUM_INCL}.";

    string INVALID_ENDING_BASIS_ROW_OR_COL =
        "The combination of starting (row, column) = (${START_ROW}, ${START_COL}) and number of rows and columns (${NUM_ROWS} and ${NUM_COLS} respectively) implies final (row, column) = (${END_ROW_EXCL}, ${END_COL_EXCL}), which is invalid for the given ${NUM_QB}-qubit ${MAX_END_INCL}-dimension density matrix.";


    string INVALID_STARTING_LOCAL_AMP_INDEX =
        "The starting local amplitude index ${START_IND} is invalid for the given ${NUM_QB}-qubit ${NUM_AMPS_TOTAL}-amplitude Qureg, distributed as ${MAX_IND_EXCL} amplitudes in each of ${NUM_NODES} nodes. Valid local indices are between 0 (inclusive) and ${MAX_IND_EXCL} (exclusive).";

    string INVALID_NUM_LOCAL_AMP_INDICES =
        "The specified number of local amplitudes (${NUM_INDS}) to synchronise is invalid for the given ${NUM_QB}-qubit ${NUM_AMPS_TOTAL}-amplitude Qureg, distributed as ${MAX_IND_EXCL} amplitudes in each of ${NUM_NODES} nodes. A valid number of amplitudes is between 0 (inclusive) and ${MAX_IND_EXCL} (inclusive).";

    string INVALID_ENDING_LOCAL_AMP_INDEX =
        "The combination of the starting local amplitude index (${START_IND}) and the number of amplitudes (${NUM_INDS}) implies an end local index of ${END_IND_EXCL} (exclusive) which is invalid for the ${NUM_QB}-qubit ${NUM_AMPS_TOTAL}-amplitude Qureg distributed over ${NUM_NODES} nodes. Valid end indices must be less than the nmber of local amplitudes, i.e. ${MAX_IND_EXCL}.";


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
        "The given target qubits contained duplicates. All qubits must be unique.";


    string NEGATIVE_NUM_CONTROLS =
        "The specified number of control qubits (${NUM_QUBITS}) is invalid. Must specify zero or more controls.";

    string NUM_CONTROLS_EXCEEDS_QUREG_SIZE =
        "The specified number of control qubits (${NUM_QUBITS}) exceeds the number of qubits in the Qureg (${QUREG_QUBITS}).";

    string NON_EMPTY_CONTROL_LIST_WAS_NULL_PTR =
        "Received a NULL pointer for a list of control qubits which was declared to have a non-zero length.";

    string INVALID_CONTROL_QUBIT = 
        "Invalid control qubit (${QUBIT_IND}). Must be greater than or equal to zero, and less than the number of qubits in the Qureg (${NUM_QUBITS}).";

    string DUPLICATE_CONTROL_QUBITS =
        "The control qubits contained duplicates. All qubits must be unique.";

    
    string CONTROLS_OVERLAP_TARGETS =
        "A qubit appeared among both the control and target qubits, which cannot overlap.";

    string INVALID_CONTROL_STATE =
        "The control qubit at index ${INDEX} has an invalid control-state of ${STATE}. Valid states are 0 and 1.";

    string DIFFERENT_NUM_CTRLS_AND_STATES =
        "A differing number of control qubits (${NUM_CTRLS}) and control states (${NUM_STATES}) was given.";


    /*
    * MEASUREMENT PARAMETERS
    */

    string ONE_QUBIT_MEASUREMENT_OUTCOME_INVALID =
        "The given qubit measurement outcome (${OUTCOME}) is invalid. Valid outcomes are 0 and 1.";

    string ONE_QUBIT_MEASUREMENT_OUTCOME_IMPOSSIBLY_UNLIKELY =
        "The specified measurement outcome (${OUTCOME}) is impossibly unlikely (i.e. has probability less than epsilon), so the post-measurement state cannot be reliably renormalised.";

    string MANY_QUBIT_MEASUREMENTS_OUTCOME_INVALID =
        "The given qubit measurement outcome (${OUTCOME}) at index ${INDEX} is invalid. Valid outcomes are 0 and 1.";

    string MANY_QUBIT_MEASUREMENT_OUTCOME_IMPOSSIBLY_UNLIKELY =
        "The specified multi-qubit measurement outcome (with binary value ${OUTCOME_VALUE}) is impossibly unlikely (i.e. has probability less than epsilon), so the post-measurement state cannot be reliably renormalised.";

    string OUTCOME_PROBS_DO_NOT_SUM_TO_ONE =
        "The state was unnormalised such that the probabilities of possible outcomes did not sum to one, and ergo cannot be meaningfully sampled.";

    string GPU_CANNOT_FIT_TEMP_MEASUREMENT_OUTCOME_PROBS =
        "The GPU has less available memory (${MEM_AVAIL} bytes) than that needed (${MEM_NEEDED} bytes) to temporarily store the ${NUM_OUTCOMES} outcome probabilities of the specified ${NUM_QUBITS} qubits.";

    string MEASUREMENT_OUTCOMES_MISMATCH_NUM_TARGETS =
        "The given number of measurement outcomes (${NUM_OUTCOMES}) is inconsistent with the given number of qubits (${NUM_QUBITS}).";


    /*
    * MISC GATE PARAMETERS
    */

    string ROTATION_AXIS_VECTOR_IS_ZERO =
        "The rotation axis vector was all zero, or within epsilion magnitude to the zero vector.";


    string CANNOT_FIT_MIXED_STATEVEC_AMPS_INTO_SINGLE_NODE =
        "Cannot perform this ${NUM_TARGS}-target operation upon a ${NUM_QUREG_QUBITS}-qubit statevector distributed between ${NUM_NODES} nodes, since each node's communication buffer (with capacity for ${NUM_QUREG_AMPS_PER_NODE} amps) cannot simultaneously store the ${NUM_TARG_AMPS} mixed remote amplitudes.";

    string CANNOT_FIT_MIXED_DENSMATR_AMPS_INTO_SINGLE_NODE =
        "Cannot perform this ${NUM_TARGS}-target operation upon a ${NUM_QUREG_QUBITS}-qubit density-matrix distributed between ${NUM_NODES} nodes, since each node's communication buffer (with capacity for ${NUM_QUREG_AMPS_PER_NODE} amps) cannot simultaneously store the ${NUM_TARG_AMPS} mixed remote amplitudes.";


    string INVALID_TROTTER_ORDER =
        "Invalid Trotter order (${ORDER}). The order parameter must be positive and even, or unity.";

    string INVALID_TROTTER_REPS =
        "Invalid number of Trotter repetitions (${REPS}). The number of repetitions must be positive.";


    /*
     * CHANNEL PARAMETERS 
     */

    string INVALID_PROBABILITY =
        "The given probability is invalid, and must instead be between 0 and 1 (both inclusive).";


    string ONE_QUBIT_DEPHASING_PROB_EXCEEDS_MAXIMAL_MIXING =
        "The given one-qubit dephasing probability exceeds that which induces maximal mixing, i.e. 1/2.";

    string TWO_QUBIT_DEPHASING_PROB_EXCEEDS_MAXIMAL_MIXING =
        "The given two-qubit dephasing probability exceeds that which induces maximal mixing, i.e. 3/4.";

    string ONE_QUBIT_DEPOLARISING_PROB_EXCEEDS_MAXIMAL_MIXING =
        "The given one-qubit depolarising probability exceeds that which induces maximal mixing, i.e. 3/4.";

    string TWO_QUBIT_DEPOLARISING_PROB_EXCEEDS_MAXIMAL_MIXING =
        "The given two-qubit depolarising probability exceeds that which induces maximal mixing, i.e. 15/16.";

    string ONE_QUBIT_PAULI_CHANNEL_TOTAL_PROBS_EXCEED_ONE =
        "The given Pauli error probabilities combine to exceed one.";

    string ONE_QUBIT_PAULI_CHANNEL_PROBS_EXCEED_MAXIMAL_MIXING =
        "The given Pauli error probabilities exceed that which induce maximal mixing. That is, the probability of any particular X, Y or Z error exceeds that of I (no error).";


    /*
     * QUREG COMBINATION
     */

    string MIXED_QUREG_NOT_DENSITY_MATRIX =
        "The first Qureg, which will undergo mixing, must be a density matrix.";

    string MIXED_QUREGS_HAVE_DIFFERENT_NUM_QUBITS =
        "The given Quregs contain an inconsistent number of qubits (${NUM_A} and ${NUM_B}) and cannot be mixed.";

    string MIXED_DENSITY_MATRICES_ARE_DIFFERENTLY_DISTRIBUED =
        "The given density matrix Quregs are differently distributed and cannot be mixed.";

    string MIXED_DENSITY_MATRIX_LOCAL_BUT_STATEVEC_DISTRIBUTED =
        "The given density matrix was local, but the statevector was distributed; this configuration is unsupported (and is ridiculous!).";


    string SUPERPOSED_QUREGS_ARE_NOT_ALL_STATEVECTORS =
        "Cannot superpose a density matrix. All quregs must be statevectors.";

    string SUPERPOSED_QUREGS_HAVE_INCONSISTENT_NUM_QUBITS =
        "Cannot superpose Quregs with differing numbers of qubits.";

    string SUPERPOSED_QUREGS_HAVE_INCONSISTENT_GPU_DEPLOYMENT =
        "Cannot superpose Quregs with inconsistent GPU deployments. All or no Quregs must be GPU-accelerated.";

    string SUPERPOSED_QUREGS_HAVE_INCONSISTENT_DISTRIBUTION =
        "Cannot superpose Quregs which are inconsistently distributed. All or no Quregs must be distributed.";


    string INIT_PURE_STATE_IS_DENSMATR =
        "The pure-state Qureg (the second argument) must be a statevector, not a density matrix.";

    string INIT_QUREG_HAS_DIFFERENT_NUM_QUBITS_TO_PURE_STATE =
        "The pure-state Qureg (the second argument) contained a differing number of qubits (${NUM_PURE_QUBITS}) as the modified Qureg (the first argument) which contained ${NUM_QUREG_QUBITS} qubits.";

    string INIT_DENSMATR_LOCAL_BUT_PURE_STATE_DISTRIBUTED = 
        "The pure-state Qureg cannot be distributed if the larger density-matrix Qureg is not.";

    string INIT_STATEVEC_DIFFERING_GPU_DEPLOYMENT_TO_PURE_STATE =
        "When the modified Qureg (the first argument) is a statevector, it must have the same GPU deployment as the pure-state Qureg (the second argument).";

    string INIT_STATEVEC_DIFFERING_DISTRIBUTION_TO_PURE_STATE =
        "When the modified Qureg (the first argument) is a statevector, it must be identically distributed to the pure-state Qureg (the second argument).";


    string CLONED_QUREGS_DIFFER_IN_NUM_QUBITS =
        "The cloned Qureg has a different number of qubits (${NUM_COPY_QUBITS}) than the target Qureg (${NUM_TARGET_QUBITS}).";

    string CLONED_QUREGS_INCONSISTENT_TYPES = 
        "The cloned and target Quregs must both be statevectors, or both be density matrices.";

    string CLONED_QUREGS_HAVE_DIFFERENT_DISTRIBUTIONS =
        "The cloned and target Quregs cannot be differently distributed.";

    string CLONED_STATEVECS_HAVE_DIFFERENT_GPU_DEPLOYMENTS =
        "The cloned and target Quregs (when both are statevectors) must have the same GPU deployment.";


    string PRODUCTED_QUREGS_HAVE_DIFFERENT_NUM_QUBITS =
        "The given Quregs contained differing numbers of qubits (${NUM_A} and ${NUM_B} qubits respectively) and are incompatible.";

    string PRODUCTED_SAME_TYPE_QUREGS_HAVE_DIFFERING_GPU_ACCELS =
        "The given Quregs were both statevectors or both density matrices but had differing GPU-accelerations. This is an illegal (and nonsensical) configuration which would invoke unexpectedly poor performance. Consider enabling GPU-aceleration for both registers.";

    string PRODUCTED_STATEVEC_DISTRIB_BUT_DENSMATR_LOCAL =
        "The given statevector Qureg was distributed while the larger density matrix Qureg was not. This is an illegal (and nonsensical) configuration. Consider distributing the density matrix.";

    
    string QUREG_IS_INCOMPATIBLE_WITH_WORKSPACE =
        "The primary Qureg is incompatible with the given workspace Qureg. The Quregs must have the same dimensions and be identically distributed and GPU-accelerated.";


    string CALC_FIDELITY_OF_DENSITY_MATRICES_NOT_YET_SUPPORTED =
        "Quregs cannot both be density matrices. Calculation of the fidelity between two mixed states is not currently supported.";

    string CALC_FIDELITY_NOT_APPROX_REAL =
        "The calculated fidelity between the given statevector and density matrix was not approximately real (the imaginary component was beyond epsilon), suggesting the density matrix was unnormalised and/or not Hermitian. Use calcInnerProduct() to find the unnormalised fidelity.";

    string CALC_BURES_DISTANCE_MAG_EXCEEDED_ONE =
        "The calculation of the Bures distance between statevectors failed, because the magnitude of their inner-product exceeded one (suggesting a complex distance), indicating one or both Quregs were unnormalised.";

    string CALC_PURIFIED_DISTANCE_NOT_APPROX_REAL =
        "The fidelity between the given statevector and density matrix was not approximately real (suggesting a complex purified distance), indicating the density matrix was unnormalised and/or not Hermitian.";

    string CALC_PURIFIED_DISTANCE_REAL_EXCEEDED_ONE =
        "The fidelity between the given statevector and density matrix exceed one (suggesting a complex purified distance), indicating either or both Quregs were unnormalised.";


    /*
     * QUREG MODIFICATION
     */

    string QUREG_RENORM_PROB_IS_ZERO =
        "Could not renormalise the Qureg because the current total probability is zero, or within epsilon to zero.";

    string INVALID_NUM_INIT_PURE_STATES =
        "The specified number of random pure states to mix (${NUM_STATES}) is invalid. Must specify 1 or more states.";


    /*
     * EXPECTATION VALUES
     */

    string CALC_STATEVEC_EXPECTED_PAULI_STR_VALUE_WAS_NOT_APPROX_REAL =
        "The calculated statevector expectation value was not approximately real (i.e. within epsilon). This should not result even from an unnormalised state, and is the result of unexpected arithmetic errors; please notify the QuEST developers!";

    string CALC_DENSMATR_EXPECTED_PAULI_STR_VALUE_WAS_NOT_APPROX_REAL =
        "The calculated density-matrix expectation value was not approximately real (i.e. within epsilon). This suggests the density matrix was unnormalised and/or not Hermitian.";


    string CALC_STATEVEC_EXPECTED_PAULI_STR_SUM_VALUE_WAS_NOT_APPROX_REAL =
        "The calculated statevector expectation value was not approximately real (i.e. was not within epsilon). This suggests that the PauliStrSum, despite being validated as (approximately) Hermitian, contained coefficients with sub-epsilon but non-negligible imaginary components which accumulated in the output value.";

    string CALC_DENSMATR_EXPECTED_PAULI_STR_SUM_VALUE_WAS_NOT_APPROX_REAL =
        "The calculated density-matrix expectation value was not approximately real (i.e. within epsilon). This suggests the density matrix was unnormalised and/or not (sufficiently close to) Hermitian.";


    string CALC_DENSMATR_EXPECTED_DIAG_MATR_VALUE_WAS_NOT_APPROX_REAL =
        "The calculated density-matrix expectation value was not approximately real (i.e. within epsilon). This suggests either the Qureg was incorrectly normalised (i.e. contained diagonal elements with non-negligible imaginary components), or that the negligible imaginary components of the FullStateDiagMatr (despite being validated as approximately Hermitian) accumulated non-negligibly in the output value.";

    string  CALC_DENSMATR_EXPECTED_DIAG_MATR_POWER_VALUE_WAS_NOT_APPROX_REAL =
        "The calculated density-matrix expectation value was not approximately real (i.e. within epsilon). This suggests either the Qureg was incorrectly normalised (i.e. contained diagonal elements with non-negligible imaginary components), or that the FullStateDiagMatr raised to the given power was no longer approximately Hermitian, or that the negligible imaginary components of the FullStateDiagMatr accumulated non-negligibly in the output value.";


    /*
     * PARTIAL TRACE
     */

    string NUM_TRACE_QUBITS_EQUALS_QUREG_SIZE =
        "Cannot trace out every qubit in the register. The reduced Qureg must contain at least one qubit.";

    string NUM_TRACE_QUBITS_EXCEEDS_DISTRIBUTED_MAX =
        "Cannot trace out ${NUM_TRACE_QUBITS} qubits of a ${NUM_QUREG_QUBITS}-qubit Qureg distributed between ${NUM_NODES}, because the reduced Qureg (which is necessarily also distributed) would have fewer than one column's worth of amplitudes per node (an illegal configuration). The maximum number of qubits that can be traced out of this Qureg is ${MAX_TRACE_QUBITS}.";

    string NUM_TRACE_QUBITS_INCONSISTENT_WITH_REDUCED_QUREG =
        "The given reduced Qureg contains an incorrect number of qubits (${OUT_QUBITS}) to be the result of tracing out ${TRACE_QUBITS} from the ${IN_QUBITS}-qubit input Qureg. Please instead pass a Qureg with ${RETAIN_QUBITS} qubits.";

    string REDUCED_QUREG_DIFFERING_DISTRIBUTION_TO_IN_QUREG =
        "The input and output Quregs must have the same GPU deployment (i.e. both or neither GPU-accelerated), despite that the output Qureg is smaller.";

    string REDUCED_QUREG_DIFFERING_GPU_DEPLOYMENT_TO_IN_QUREG =
        "The input and output Quregs must have the same distribution (i.e. both or neither distributed), despite that the output Qureg is smaller.";


    /*
     * FILE IO
     */

    string CANNOT_READ_FILE = 
        "Could not load and read the given file. Make sure the file exists and is readable as plaintext.";


    /*
     * TEMPORARY ALLOCATIONS
     */

    string TEMP_ALLOC_FAILED =
        "A temporary allocation of ${NUM_ELEMS} elements (each of ${NUM_BYTES_PER_ELEM} bytes) failed, possibly because of insufficient memory.";

}



/*
 * INVALID INPUT RESPONSE BEHAVIOUR
 */

void default_inputErrorHandler(const char* func, const char* msg) {

    // safe to call even before MPI has been setup, and ignores user-set trailing newlines.
    // It begins with \n to interrupt half-printed lines (when trailing newlines are set to
    // 0 via setNumReportedNewlines(0)), for visual clarity. Note that user's overriding
    // functions might not think to print an initial newline but oh well!
    print(string("\n")
        + "QuEST encountered a validation error during function " 
        + "'" + func + "':\n" + msg + "\n"
        + "Exiting...\n");

    // force a synch because otherwise non-main nodes may exit before print, and MPI
    // will then attempt to instantly abort all nodes, losing the error message.
    comm_sync();

    // finalise MPI before error-exit to avoid scaring user with giant MPI error message
    if (comm_isInit())
        comm_end();

    // simply exit, interrupting any other process (potentially leaking)
    exit(EXIT_FAILURE);
}

void (*global_inputErrorHandler)(const char*, const char*) = default_inputErrorHandler;

void validateconfig_setErrorHandler(void (*callback)(const char*, const char*)) {

    global_inputErrorHandler = callback;
}



/*
 * VALIDATION TOGGLE
 *
 * which disables all forms of validation (e.g. numerical
 * Hermiticity of operators, and whether qubit indices are
 * valid, etc), such that global_validationEpsilon (defined
 * below) is ignored. Note that many validation functions
 * will still actually perform the validation checks/compute,
 * and thereafter decide whether to throw an error within
 * assertThat(). Ergo expensive validators should check that
 * global_isValidationEnabled=true before even performing
 * their checks, to avoid superfluous computation.
 */

static bool global_isValidationEnabled = true;

void validateconfig_enable() {
    global_isValidationEnabled = true;
}
void validateconfig_disable() {
    global_isValidationEnabled = false;
}
bool validateconfig_isEnabled() {
    return global_isValidationEnabled;
}

// this file validates that the outputs of reductions (like 
// calcFidelity) are within numerical bounds (e.g. approx
// real). Because these reductions combine exponentially
// many decimals, they are much less precise than API inputs,
// so should be validated with a more tolerant epsilon. We
// here hackily specify the factor by which to expand eps.
qreal REDUCTION_EPSILON_FACTOR = 100;



/*
 * VALIDATION PRECISION
 *
 * which influences how strict numerical checks are performed,
 * e.g. checking unitarity/hermiticity/CPTP-ness. This is only
 * consulted when global_isValidationEnabled=true, and can be
 * separately disabled by setting epsilon=0, in which case 
 * permanent properties of structs (like .isApproxCPTP) will not be 
 * overwritten (so will stay validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
 */

static qreal global_validationEpsilon = DEFAULT_VALIDATION_EPSILON;

void validateconfig_setEpsilon(qreal eps) {
    global_validationEpsilon = eps;
}
void validateconfig_setEpsilonToDefault() {
    global_validationEpsilon = DEFAULT_VALIDATION_EPSILON;
}
qreal validateconfig_getEpsilon() {
    return global_validationEpsilon;
}

bool isNumericalValidationDisabled() {
    return (
        global_isValidationEnabled == 0 || 
        global_validationEpsilon   == 0);
}



/*
 * UTILITIES
 */


/// @todo
/// this method of stringification is terrible; it precludes us
/// from displaying non-integral types, like floating-point or
/// strings. But further, it requires casting integral types
/// (like size_t) to long long int. But alas size_t can be bigger
/// than long long int, causing overflow - this happens when
/// reporting failed max-smemory allocations. We cannot simply
/// switch to 'size_t' which is unsigned, losing the ability to
/// show small negative numbers. Grr!

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
    if (!global_isValidationEnabled)
        return;

    // this function does not seek consensus among nodes in distributed 
    // settings in order to remain cheap (consensus requires expensive sync
    // and comm), so is suitable for validation which is gauranteed to be
    // uniform between nodes (assuming user's do not hack in rank-specific
    // arguments to the API!)

    // invoke the potentially user-overriden error function
    if (!valid)
        global_inputErrorHandler(func, msg.c_str());
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

void assertAllNodesAgreeThat(bool valid, string msg, tokenSubs vars, const char* func) {

    // avoid communication if validation not enabled anyway
    if (!global_isValidationEnabled)
        return;

    // this function seeks consensus among distributed nodes before validation,
    // to ensure all nodes fail together (and ergo all validly finalize MPI)
    // when performing validation that may be non-uniform between nodes. For
    // example, mallocs may succeed on one node but fail on another due to
    // inhomogeneous loads.
    if (comm_isInit())
        valid = comm_isTrueOnAllNodes(valid);

    // prepare error message only if validation will fail
    if (!valid)
        if (!vars.empty())
            msg = getStringWithSubstitutedVars(msg, vars);

    // commenting out the top if-statement above (preserving the bottom one)
    // is convenient during development, since it ensures that even validation
    // which passed has a correctly formatted error message, else errors.

    assertThat(valid, msg, func);
}
void assertAllNodesAgreeThat(bool valid, string msg, const char* func) {

    assertAllNodesAgreeThat(valid, msg, {}, func);
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

void validate_newEnvNodesEachHaveUniqueGpu(const char* caller) {

    // this validation can be disabled for debugging/dev purposes
    // (caller should explicitly check this preprocessor too for clarity)
    if (PERMIT_NODES_TO_SHARE_GPU)
        return;

    bool uniqueGpus = ! gpu_areAnyNodesBoundToSameGpu();
    assertAllNodesAgreeThat(uniqueGpus, report::MULTIPLE_NODES_BOUND_TO_SAME_GPU, caller);
}

void validate_gpuIsCuQuantumCompatible(const char* caller) {

    int minCC = 70;
    int ourCC = gpu_getComputeCapability();
    tokenSubs vars = {
        {"${MIN_CC}", minCC},
        {"${OUR_CC}", ourCC}
    };
    assertAllNodesAgreeThat(ourCC >= minCC, report::CUQUANTUM_DEPLOYED_ON_BELOW_CC_GPU, vars, caller);

    bool hasMemPools = gpu_doesGpuSupportMemPools();
    assertAllNodesAgreeThat(hasMemPools, report::CUQUANTUM_DEPLOYED_ON_GPU_WITHOUT_MEM_POOLS, caller);
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

void validate_randomSeeds(unsigned* seeds, int numSeeds, const char* caller) {

    // only the root node's seeds are consulted, so we permit all non-root
    // nodes to have invalid parameters. All nodes however must know/agree
    // when the root node's seeds are invalid, to synchronise validation

    int numRootSeeds = numSeeds;
    if (getQuESTEnv().isDistributed)
        comm_broadcastIntsFromRoot(&numRootSeeds, 1);

    assertThat(numRootSeeds > 0, report::INVALID_NUM_RANDOM_SEEDS, {{"${NUM_SEEDS}", numSeeds}}, caller);
}

void validate_newEpsilonValue(qreal eps, const char* caller) {

    assertThat(eps >= 0, report::INVALID_NEW_EPSILON, {{"${NEW_EPS}", eps}}, caller);
}

void validate_newMaxNumReportedScalars(qindex numRows, qindex numCols, const char* caller) {

    assertThat(numRows >= 0, report::INVALID_NUM_REPORTED_SCALARS, {{"${NUM_ITEMS}", numRows}}, caller);
    assertThat(numCols >= 0, report::INVALID_NUM_REPORTED_SCALARS, {{"${NUM_ITEMS}", numCols}}, caller);
}

void validate_newMaxNumReportedSigFigs(int numSigFigs, const char* caller) {

    assertThat(numSigFigs >= 1, report::INVALID_NUM_REPORTED_SIG_FIGS, {{"${NUM_SIG_FIGS}", numSigFigs}}, caller);
}

void validate_newNumReportedNewlines(int numNewlines, const char* caller) {

    assertThat(numNewlines >= 0, report::INVALID_NUM_REPORTED_NEWLINES, {{"${NUM_NEWLINES}", numNewlines}}, caller);
}

void validate_numReportedNewlinesAboveZero(const char* caller) {

    assertThat(printer_getNumTrailingNewlines() > 0, report::INSUFFICIENT_NUM_REPORTED_NEWLINES, caller);
}

void validate_numPauliChars(const char* paulis, const char* caller) {

    // check position of terminal char, else default to numChars=5 (illegal)
    int numChars = 0;
    for (int i=0; i<5 && paulis[i] != '\0'; i++)
        numChars++;
    
    assertThat(numChars==4, report::INVALID_NUM_NEW_PAULI_CHARS, caller);
}

void validate_reportedPauliStrStyleFlag(int flag, const char* caller) {

    assertThat(flag==0 || flag==1, report::INVALID_REPORTED_PAULI_STR_STYLE_FLAG, {{"${FLAG}",flag}}, caller);
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
    int maxNumQubits = (int) mem_getMaxNumQuregQubitsBeforeGlobalMemSizeofOverflow(isDensMatr, numQuregNodes);

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

void validate_newQuregParams(int numQubits, int isDensMatr, int isDistrib, int isGpuAccel, int isMultithread, QuESTEnv env, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!global_isValidationEnabled)
        return;

    assertQuregNonEmpty(numQubits, caller);
    assertQuregDeployFlagsRecognised(isDensMatr, isDistrib, isGpuAccel, isMultithread, caller);
    assertQuregDeploysEnabledByEnv(isDistrib, isGpuAccel, isMultithread, env, caller);
    assertQuregTotalNumAmpsDontExceedMaxIndex(numQubits, isDensMatr, caller);
    assertQuregLocalMemDoesntExceedMaxSizeof(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregNotDistributedOverTooManyNodes(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInCpuMem(numQubits, isDensMatr, isDistrib, env, caller);
    assertQuregFitsInGpuMem(numQubits, isDensMatr, isDistrib, isGpuAccel, env, caller);
}

void validate_newQuregAllocs(Qureg qureg, const char* caller) {

    // this validation is called AFTER the caller has checked for failed
    // allocs and (in that scenario) freed every pointer, but does not 
    // overwrite any pointers to nullptr, so the failed alloc is known.
    // This is only safe to do so (rather than making the caller set ptrs
    // to nullptr) because the struct contains only 1D pointers (no nesting)

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

void assertMatrixDeployFlagsRecognised(int isDistrib, int isGpuAccel, int isMultithread, const char* caller) {

    // deployment flags must be boolean or auto
    tokenSubs vars = {{"${AUTO_DEPLOYMENT_FLAG}", modeflag::USE_AUTO}};
    assertThat(isDistrib     == 0 || isDistrib     == 1 || isDistrib     == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_MATRIX_IS_DISTRIB,     vars, caller);
    assertThat(isGpuAccel    == 0 || isGpuAccel    == 1 || isGpuAccel    == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_MATRIX_IS_GPU_ACCEL,   vars, caller);
    assertThat(isMultithread == 0 || isMultithread == 1 || isMultithread == modeflag::USE_AUTO, report::INVALID_OPTION_FOR_MATRIX_IS_MULTITHREAD, vars, caller);
}

void assertMatrixDeploysEnabledByEnv(int isDistrib, int isGpuAccel, int isMultithread, QuESTEnv env, const char* caller) {

    // cannot deploy to backend not already enabled by the environment
    if (!env.isDistributed)
        assertThat(isDistrib     == 0 || isDistrib     == modeflag::USE_AUTO, report::NEW_DISTRIB_MATRIX_IN_NON_DISTRIB_ENV, caller);
    if (!env.isGpuAccelerated)
        assertThat(isGpuAccel    == 0 || isGpuAccel    == modeflag::USE_AUTO, report::NEW_GPU_MATRIX_IN_NON_GPU_ACCEL_ENV, caller);
    if (!env.isMultithreaded)
        assertThat(isMultithread == 0 || isMultithread == modeflag::USE_AUTO, report::NEW_MULTITHREAD_MATRIX_IN_NON_MULTITHREAD_ENV, caller);
}

void assertMatrixNonEmpty(int numQubits, const char* caller) {

    assertThat(numQubits >= 1, report::NEW_MATRIX_NUM_QUBITS_NOT_POSITIVE, {{"${NUM_QUBITS}",numQubits}}, caller);
}

void assertMatrixTotalNumElemsDontExceedMaxIndex(int numQubits, bool isDense, const char* caller) {

    int maxNumQubits = mem_getMaxNumMatrixQubitsBeforeIndexOverflow(isDense);

    string msg = (isDense)?
        report::NEW_COMP_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX:
        report::NEW_DIAG_MATR_NUM_ELEMS_WOULD_EXCEED_QINDEX;

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
    int maxNumQubits = mem_getMaxNumMatrixQubitsBeforeGlobalMemSizeofOverflow(isDense, numMatrNodes);

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
    
    /// @todo
    /// seek expensive node consensus in case of heterogeneous RAM - alas this may induce
    /// unnecessary slowdown (due to sync and broadcast) in applications allocating many
    /// small matrices in the heap. If this turns out to be the case, we could opt to
    /// enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    /// chance of it fitting into some node RAM but not others isn't negligible.
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
    
    /// @todo
    /// seek expensive node consensus in case of heterogeneous GPU hardware - alas this may 
    /// induce unnecessary slowdown (due to sync and broadcast) in applications allocating many
    /// small matrices in the GPU. If this turns out to be the case, we could opt to
    /// enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    /// chance of it fitting into some GPU's memory but not others isn't negligible.
    assertAllNodesAgreeThat(matrFitsInMem, msg, vars, caller);
}

void assertNewMatrixParamsAreValid(int numQubits, int useDistrib, int useGpu, int useMultithread, bool isDenseType, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!global_isValidationEnabled)
        return;

    QuESTEnv env = getQuESTEnv();
    assertMatrixDeployFlagsRecognised(useDistrib, useGpu, useMultithread, caller);
    assertMatrixDeploysEnabledByEnv(useDistrib, useGpu, useMultithread, env, caller);
    assertMatrixNonEmpty(numQubits, caller);
    assertMatrixTotalNumElemsDontExceedMaxIndex(numQubits, isDenseType, caller);
    assertMatrixLocalMemDoesntExceedMaxSizeof(numQubits,  isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixNotDistributedOverTooManyNodes(numQubits, isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixFitsInCpuMem(numQubits, isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixFitsInGpuMem(numQubits, isDenseType, useDistrib, useGpu, env.numNodes, caller);
}

void validate_newCompMatrParams(int numQubits, const char* caller) {
    validate_envIsInit(caller);

    // CompMatr can never be distributed nor multithreaded
    int useDistrib = 0;
    int useMultithread = 0;

    // CompMatr is always GPU accelerated whenever enabled by the environment
    bool useGpu = getQuESTEnv().isGpuAccelerated;

    // CompMatr stores 2^(2*numQubits) elems
    bool isDenseType = true;

    assertNewMatrixParamsAreValid(numQubits, useDistrib, useGpu, useMultithread, isDenseType, caller);
}
void validate_newDiagMatrParams(int numQubits, const char* caller) {
    validate_envIsInit(caller);

    // DiagMatr can never be distributed nor multithreaded
    int useDistrib = 0;
    int useMultithread = 0;

    // DiagMatr is always GPU accelerated whenever enabled by the environment
    bool useGpu = getQuESTEnv().isGpuAccelerated;

    // DiagMatr stores only the diagonals; 2^numQubits elems
    bool isDenseType = false;

    assertNewMatrixParamsAreValid(numQubits, useDistrib, useGpu, useMultithread, isDenseType, caller);
}
void validate_newFullStateDiagMatrParams(int numQubits, int useDistrib, int useGpu, int useMultithread, const char* caller) {
    validate_envIsInit(caller);

    // FullStateDiagMatr stores only the diagonals
    bool isDenseType = false;

    assertNewMatrixParamsAreValid(numQubits, useDistrib, useGpu, useMultithread, isDenseType, caller);
}

// T can be CompMatr, DiagMatr, FullStateDiagMatr (i.e. heap-based matrices)
template <typename T>
void assertNewMatrixAllocsSucceeded(T matr, size_t numBytes, const char* caller) {

    // this validation is called AFTER the caller has checked for failed
    // allocs and (in that scenario) freed every pointer, but does not 
    // overwrite any pointers to nullptr, so the failed alloc is known.
    // This is only safe to do so (rather than making the caller set ptrs
    // to nullptr) because the structs contains only 1D pointer; even 
    // CompMatr which "fakes" a 2D ptr via offsets of a contiguous array.

    // we expensively get node consensus about malloc failure, in case of heterogeneous hardware/loads,
    // but we avoid this potetially expensive synchronisation if validation is anyway disabled
    // (which also avoids us enumerating the matrix rows)
    if (!global_isValidationEnabled)
        return;

    /// @todo
    /// fix this abhorrent hackiness! Presently, tokenSubs accepts only qindex (in lieu of 
    /// unsigned size_t) in order to be able to report negative numbers. But alas, the max
    /// size_t is bigger than max qindex, due to the sign-bit. So trying to report numBytes
    /// can cause a size_t -> qindex overflow. This is a realistic scenario, occurring when
    /// when the user tries to allocate the max-size memory for which malloc incidentally
    /// fails. This is when numBytes is +1 too big to be a qindex; we simply reduce by 2!
    /// We under-report the memory by 2 bytes, instead of 1, just to avoid an odd number 
    /// which an astute user would immediately notice is not a power-of-2 and be confused by.
    /// This is a hacky evil, but it is better than reporting a negative memory size!
    tokenSubs vars;
    vars["${NUM_BYTES}"] = ((qindex) numBytes <= 0)?
        numBytes - 2 : numBytes;

    // assert CPU array (which may be nested arrays) all allocated successfully
    bool isAlloc;
    if constexpr (util_isDenseMatrixType<T>()) {
        // size of .cpuElems isn't included in numBytes report which is fine; it's
        // quadratically smaller than .cpuElemsFlat so quickly negligible
        isAlloc = mem_isAllocated(matr.cpuElemsFlat) && mem_isOuterAllocated(matr.cpuElems);
    } else
        isAlloc = mem_isAllocated(matr.cpuElems);
    assertAllNodesAgreeThat(isAlloc, report::NEW_MATRIX_CPU_ELEMS_ALLOC_FAILED, vars, caller);

    // optionally assert GPU memory was malloc'd successfully
    bool gpuShouldBeAlloc = getQuESTEnv().isGpuAccelerated;
    if constexpr (util_isFullStateDiagMatr<T>())
        gpuShouldBeAlloc &= matr.isGpuAccelerated;
        
    if (gpuShouldBeAlloc)
        assertAllNodesAgreeThat(mem_isAllocated(util_getGpuMemPtr(matr)), report::NEW_MATRIX_GPU_ELEMS_ALLOC_FAILED, vars, caller);

    // assert the teeny-tiny heap flags are alloc'd
    vars["${NUM_BYTES}"] = sizeof(*(matr.isApproxUnitary)); // all fields are same-size
    assertAllNodesAgreeThat(mem_isAllocated(matr.isApproxUnitary),     report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
    assertAllNodesAgreeThat(mem_isAllocated(matr.isApproxHermitian),   report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
    assertAllNodesAgreeThat(mem_isAllocated(matr.wasGpuSynced),  report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);

    // only diagonal matrices (which can be exponentiated) have these additional flags
    if constexpr (!util_isDenseMatrixType<T>()) {
        assertAllNodesAgreeThat(mem_isAllocated(matr.isApproxNonZero),       report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
        assertAllNodesAgreeThat(mem_isAllocated(matr.isStrictlyNonNegative), report::NEW_HEAP_FLAG_ALLOC_FAILED, vars, caller);
    }
}

void validate_newMatrixAllocs(CompMatr matr, const char* caller) {

    bool isDenseMatrix = true;
    int numNodes = 1;
    size_t numBytes = mem_getLocalMatrixMemoryRequired(matr.numQubits, isDenseMatrix, numNodes);
    assertNewMatrixAllocsSucceeded(matr, numBytes, caller);
}
void validate_newMatrixAllocs(DiagMatr matr, const char* caller) {

    bool isDenseMatrix = false;
    int numNodes = 1;
    size_t numBytes = mem_getLocalMatrixMemoryRequired(matr.numQubits, isDenseMatrix, numNodes);
    assertNewMatrixAllocsSucceeded(matr, numBytes, caller);
}
void validate_newMatrixAllocs(FullStateDiagMatr matr, const char* caller) {

    bool isDenseMatrix = false;
    int numNodes = (matr.isDistributed)? comm_getNumNodes() : 1;
    size_t numBytes = mem_getLocalMatrixMemoryRequired(matr.numQubits, isDenseMatrix, numNodes);
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

void validate_matrixNewElemsPtrNotNull(qcomp* elems, const char* caller) {

    assertThat(mem_isAllocated(elems), report::DIAG_MATR_NEW_ELEMS_NULL_PTR, caller);
}

void validate_matrixNewElemsPtrNotNull(qcomp** elems, qindex numRows, const char* caller) {

    // messages are suitable for all dense matrices, including SuperOp

    assertThat(mem_isOuterAllocated(elems), report::DENSE_MATR_NEW_ELEMS_OUTER_NULL_PTR, caller);

    for (qindex i=0; i<numRows; i++)
        assertThat(mem_isAllocated(elems[i]), report::DENSE_MATR_NEW_ELEMS_INNER_NULL_PTR, caller);
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

void validate_matrixRowsAllSameSize(vector<vector<qcomp>> matrix, const char* caller) {

    if (matrix.empty())
        return;

    size_t dim = matrix[0].size();

    size_t row=0;
    bool allSame = true;
    for (row=0; row<matrix.size() && allSame; row++)
        allSame = (matrix[row].size() == dim);

    // extremely lazily avoiding seg-fault from premature matrix[row==end])
    if (allSame)
        return;

    tokenSubs vars = {
        {"{FIRST_LEN}",    dim},
        {"${OUTLIER_IND}", row},
        {"${OUTLIER_LEN}", matrix[row].size()}};

    assertThat(allSame, report::NESTED_VECTOR_MATRIX_HAS_INCONSISTENT_NUM_COLUMNS, vars, caller);
}



/*
 * EXISTING MATRICES
 */

// T can be CompMatr, DiagMatr, FullStateDiagMatr
template <class T>
void assertAdditionalHeapMatrixFieldsAreValid(T matr, const char* caller) {

    // assert heap pointers are not NULL
    assertThat(mem_isAllocated(matr.isApproxUnitary),     report::INVALID_HEAP_FLAG_PTR, caller);
    assertThat(mem_isAllocated(matr.isApproxHermitian),   report::INVALID_HEAP_FLAG_PTR, caller);
    assertThat(mem_isAllocated(matr.wasGpuSynced),  report::INVALID_HEAP_FLAG_PTR, caller);

    // only diagonal matrices (which can be exponentiated) have these additional flags
    if constexpr (!util_isDenseMatrixType<T>()) {
        assertThat(mem_isAllocated(matr.isApproxNonZero),       report::INVALID_HEAP_FLAG_PTR, caller);
        assertThat(mem_isAllocated(matr.isStrictlyNonNegative), report::INVALID_HEAP_FLAG_PTR, caller);
    }

    tokenSubs vars = {{"${BAD_FLAG}", 0}, {"${UNKNOWN_FLAG}", validate_STRUCT_PROPERTY_UNKNOWN_FLAG}};

    // assert isApproxUnitary has valid value
    int flag = *matr.isApproxUnitary;
    vars["${BAD_FLAG}"] = flag;
    assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);

    // assert isApproxHermitian has valid value
    flag = *matr.isApproxHermitian;
    vars["${BAD_FLAG}"] = flag;
    assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);

    // assert wasGpuSynced has a valid value
    flag = *matr.wasGpuSynced;
    vars["${BAD_FLAG}"] = flag;
    assertThat(flag == 0 || flag == 1, report::INVALID_HEAP_FLAG_VALUE, vars, caller);

    // assert isApproxNonZero and isStrictlyNonNegative have valid values (only bound to diagonal matrices)
    if constexpr (!util_isDenseMatrixType<T>()) {
        flag = *matr.isApproxNonZero;
        vars["${BAD_FLAG}"] = flag;
        assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);

        flag = *matr.isStrictlyNonNegative;
        vars["${BAD_FLAG}"] = flag;
        assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);
    }

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
void validate_matrixFields(CompMatr1 m, const char* caller) { assertMatrixFieldsAreValid(m, 1,           report::INVALID_COMP_MATR_1_FIELDS, caller); }
void validate_matrixFields(CompMatr2 m, const char* caller) { assertMatrixFieldsAreValid(m, 2,           report::INVALID_COMP_MATR_2_FIELDS, caller); }
void validate_matrixFields(CompMatr  m, const char* caller) { assertMatrixFieldsAreValid(m, m.numQubits, report::INVALID_COMP_MATR_FIELDS,   caller); }
void validate_matrixFields(DiagMatr1 m, const char* caller) { assertMatrixFieldsAreValid(m, 1,           report::INVALID_DIAG_MATR_1_FIELDS, caller); }
void validate_matrixFields(DiagMatr2 m, const char* caller) { assertMatrixFieldsAreValid(m, 2,           report::INVALID_DIAG_MATR_2_FIELDS, caller); }
void validate_matrixFields(DiagMatr  m, const char* caller) { assertMatrixFieldsAreValid(m, m.numQubits, report::INVALID_DIAG_MATR_FIELDS,   caller); }
void validate_matrixFields(FullStateDiagMatr m, const char* caller) { assertMatrixFieldsAreValid(m, m.numQubits, report::INVALID_FULL_STATE_DIAG_MATR_FIELDS, caller); }

// type T can be CompMatr, DiagMatr or FullStateDiagMatr
template <class T>
void assertMatrixIsSynced(T matr, string errMsg, const char* caller) {

    // we don't need to perform any sync check in CPU-only mode
    if (!mem_isAllocated(util_getGpuMemPtr(matr)))
        return;

    // check if GPU amps have EVER been overwritten; we sadly cannot check the LATEST changes were pushed though.
    // note we check this whenever the matrix has GPU memory, even if it is being applied upon a Qureg which is
    // NOT GPU-accelerated and ergo the GPU memory is not consulted. It's best to build the habit in the user!
    assertThat(*(matr.wasGpuSynced) == 1, errMsg, caller);
}
void validate_matrixIsSynced(CompMatr matr, const char* caller) { assertMatrixIsSynced(matr, report::COMP_MATR_NOT_SYNCED_TO_GPU, caller);}
void validate_matrixIsSynced(DiagMatr matr, const char* caller) { assertMatrixIsSynced(matr, report::DIAG_MATR_NOT_SYNCED_TO_GPU, caller); }
void validate_matrixIsSynced(FullStateDiagMatr matr, const char* caller) { assertMatrixIsSynced(matr, report::FULL_STATE_DIAG_MATR_NOT_SYNCED_TO_GPU, caller); }

// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template <class T> 
void assertMatrixIsUnitary(T matr, const char* caller) {

    // validate both stack and heap matrices have been correctly initialised
    validate_matrixFields(matr, caller);

    // validate heap matrices have ever written to their GPU memories (if exists)
    if constexpr (util_isHeapMatrixType<T>())
        validate_matrixIsSynced(matr, caller);

    // avoid superfluous expensive unitarity check below (do not overwrite .isApproxUnitary)
    if (isNumericalValidationDisabled())
        return;

    // may overwrite matr.isApproxUnitary of heap matrices, otherwise ignores epsilon
    assertThat(util_isUnitary(matr, global_validationEpsilon), report::MATRIX_NOT_UNITARY, caller);
}
void validate_matrixIsUnitary(CompMatr1 m, const char* caller) { assertMatrixIsUnitary(m, caller); }
void validate_matrixIsUnitary(CompMatr2 m, const char* caller) { assertMatrixIsUnitary(m, caller); }
void validate_matrixIsUnitary(CompMatr  m, const char* caller) { assertMatrixIsUnitary(m, caller); }
void validate_matrixIsUnitary(DiagMatr1 m, const char* caller) { assertMatrixIsUnitary(m, caller); }
void validate_matrixIsUnitary(DiagMatr2 m, const char* caller) { assertMatrixIsUnitary(m, caller); }
void validate_matrixIsUnitary(DiagMatr  m, const char* caller) { assertMatrixIsUnitary(m, caller); }
void validate_matrixIsUnitary(FullStateDiagMatr m, const char* caller) { assertMatrixIsUnitary(m, caller); }

void validate_unitaryExponentIsReal(qcomp exponent, const char* caller) {

    // this validation always checks that 'exponent' is only
    // approximately real, permitting non-zero imaginary
    // component. Functions which require a strictly real
    // exponent never call this; they accept 'qreal' exponents.

    if (isNumericalValidationDisabled())
        return;

    // assesses diag^expo unitarity for any unitary diag
    assertThat(util_isApproxReal(exponent, global_validationEpsilon), report::UNITARY_DIAG_MATR_EXPONENT_NOT_APPROX_REAL, caller);
}

// type T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template <class T> 
void assertMatrixIsHermitian(T matr, const char* caller) {

    // validate both stack and heap matrices have been correctly initialised
    validate_matrixFields(matr, caller);

    // validate heap matrices have ever written to their GPU memories (if exists)
    if constexpr (util_isHeapMatrixType<T>())
        validate_matrixIsSynced(matr, caller);

    // avoid superfluous expensive hermiticity check below (do not overwrite  matr.isApproxHermitian)
    if (isNumericalValidationDisabled())
        return;

    // may overwrite matr.isApproxHermitian of heap matrices, otherwise ignores epsilon
    assertThat(util_isHermitian(matr, global_validationEpsilon), report::MATRIX_NOT_HERMITIAN, caller);
}
void validate_matrixIsHermitian(CompMatr1 m, const char* caller) { assertMatrixIsHermitian(m, caller); }
void validate_matrixIsHermitian(CompMatr2 m, const char* caller) { assertMatrixIsHermitian(m, caller); }
void validate_matrixIsHermitian(CompMatr  m, const char* caller) { assertMatrixIsHermitian(m, caller); }
void validate_matrixIsHermitian(DiagMatr1 m, const char* caller) { assertMatrixIsHermitian(m, caller); }
void validate_matrixIsHermitian(DiagMatr2 m, const char* caller) { assertMatrixIsHermitian(m, caller); }
void validate_matrixIsHermitian(DiagMatr  m, const char* caller) { assertMatrixIsHermitian(m, caller); }
void validate_matrixIsHermitian(FullStateDiagMatr m, const char* caller) { assertMatrixIsHermitian(m, caller); }

// type T can be DiagMatr, FullStateDiagMatr
template <class T> 
void assertMatrExpIsNonDiverging(T matr, qcomp exponent, const char* caller) {

    validate_matrixFields(matr, caller);
    validate_matrixIsSynced(matr, caller);

    // avoid exepensive and epsilon-dependent validation below (do not overwrite matr.isApproxNonZero)
    if (isNumericalValidationDisabled())
        return;

    // divergences are only validated when the imaginary component is strictly 
    // zero, otherwise alternate complex exponentiation is sometimes performed
    // with a more complicated numerical stability
    if (std::imag(exponent) != 0)
        return;

    // when the real exponent is STRICTLY less than zero, it is required that every
    // matrix element's magnitude is APPROX non-zero, to avoid 1/0 divergences. 
    // We do this independent of the size of exponent, even despite that exponents
    // really close to 0 (from below) can "counteract" the blowing up, because
    // precision is too poor near epsilon for this to be implicitly relied upon.
    if (std::real(exponent) < 0)
        assertThat(util_isApproxNonZero(matr, global_validationEpsilon), report::DIAG_MATR_APPROX_ZERO_WHILE_EXPONENT_REAL_AND_NEGATIVE, caller);
}
void validate_matrixExpIsNonDiverging(DiagMatr          m, qcomp p, const char* caller) { assertMatrExpIsNonDiverging(m, p, caller); }
void validate_matrixExpIsNonDiverging(FullStateDiagMatr m, qcomp p, const char* caller) { assertMatrExpIsNonDiverging(m, p, caller); }

// type T can be DiagMatr, FullStateDiagMatr
template <class T> 
void assertMatrExpIsHermitian(T matr, qreal exponent, const char* caller) {

    // below validation can invoke communication and expensive data processing,
    // even when numerical-epsilon is zero, which we avoid when all validation is off.
    // we always proceed however if merely numerical-validation is disabled (via zero
    // epsilon) because some checks below are epsilon independent
    if (!global_isValidationEnabled)
        return;

    // the calling function will use the std::pow(qreal,qreal) overload, rather than
    // std::pow(qcomp,qcomp), passing in the real components of matr and the given
    // exponent. As such, the result must never be complex which instead becomes NaN.
    // The validation below ergo ensures that pow(a,b) is always real and stable. 
    // All validations upon 'matr' consult existing properties (like .isHermitian),
    // computing them fresh and recording them if not already known

    // the matrix itself must be approximately real, since we consult only its reals.
    // this also validates the matrix fields, and whether GPU-matrices are synced
    validate_matrixIsHermitian(matr, caller);

    // when the exponent is not STRICTLY an integer, it is required that every matrix
    // elem's real component is STRICTLY positive, so that real-pow doesn't create NaNs.
    // this can overwrite matr.isStrictlyNonNegative even when validation epsilon=0
    if (!util_isStrictlyInteger(exponent))
        assertThat(util_isStrictlyNonNegative(matr), report::HERMITIAN_DIAG_MATR_NEGATIVE_WHILE_EXPONENT_NOT_INTEGER, caller);

    // divergences don't break Hermiticity per se, but do sabotage numerical accuracy.
    validate_matrixExpIsNonDiverging(matr, qcomp(exponent,0), caller);

    // the final plausible scenario we have not checked is when both the matrix elem
    // and exponent are strictly positive but very close to zero. In that case, the
    // result tends to 1 so does not vanish or blow up unexpectedly. All fine!
}

void validate_matrixExpIsHermitian(DiagMatr          m, qreal p, const char* caller) { assertMatrExpIsHermitian(m, p, caller); }
void validate_matrixExpIsHermitian(FullStateDiagMatr m, qreal p, const char* caller) { assertMatrExpIsHermitian(m, p, caller); }

template <class T>
void assertMatrixDimMatchesTargs(T matr, int numTargs, const char* caller) {

    validate_matrixFields(matr, caller);

    if constexpr (util_isHeapMatrixType<T>())
        validate_matrixIsSynced(matr, caller);

    int numMatrQubits = 1; // CompMatr1 or DiagMatr1
    if constexpr (util_isCompMatr2<T>() || util_isDiagMatr2<T>())
        numMatrQubits = 2;
    if constexpr (util_isHeapMatrixType<T>())
        numMatrQubits = matr.numQubits;

    // note (numMatrQubits <= qureg.numQubits) was prior validated by numTargs <= qureg.numQubits
    tokenSubs vars = {{"${MATRIX_NUM_QUBITS}", numMatrQubits}, {"${NUM_TARGS}", numTargs}};
    assertThat(numMatrQubits == numTargs, report::MATRIX_SIZE_MISMATCHES_NUM_TARGETS, vars, caller);
}

void validate_matrixDimMatchesTargets(CompMatr1 matr, int numTargs, const char* caller) { assertMatrixDimMatchesTargs(matr, numTargs, caller); }
void validate_matrixDimMatchesTargets(CompMatr2 matr, int numTargs, const char* caller) { assertMatrixDimMatchesTargs(matr, numTargs, caller); }
void validate_matrixDimMatchesTargets(CompMatr  matr, int numTargs, const char* caller) { assertMatrixDimMatchesTargs(matr, numTargs, caller); }
void validate_matrixDimMatchesTargets(DiagMatr1 matr, int numTargs, const char* caller) { assertMatrixDimMatchesTargs(matr, numTargs, caller); }
void validate_matrixDimMatchesTargets(DiagMatr2 matr, int numTargs, const char* caller) { assertMatrixDimMatchesTargs(matr, numTargs, caller); }
void validate_matrixDimMatchesTargets(DiagMatr  matr, int numTargs, const char* caller) { assertMatrixDimMatchesTargs(matr, numTargs, caller); }

void validate_matrixAndQuregAreCompatible(FullStateDiagMatr matr, Qureg qureg, bool expecOnly, const char* caller) {

    // we do not need to define this function for the other matrix types,
    // since their validation will happen through validation of the
    // user-given list of target qubits. But we do need to define it for
    // FullStatedDiagMatr to check both distribution compatibility, and
    // that dimensions match. When expecOnly=true, we relax the necessity
    // that the distributions match; one or both can be distributed

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
    if (!expecOnly)
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

    int maxNumQubits = mem_getMaxNumSuperOpQubitsBeforeGlobalMemSizeofOverflow();

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
    
    /// @todo
    /// seek expensive node consensus in case of heterogeneous RAM - alas this may induce
    /// unnecessary slowdown (due to sync and broadcast) in applications allocating many
    /// small matrices in the heap. If this turns out to be the case, we could opt to
    /// enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    /// chance of it fitting into some node RAM but not others isn't negligible.
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

    /// @todo
    /// seek expensive node consensus in case of heterogeneous GPU hardware - alas this may 
    /// induce unnecessary slowdown (due to sync and broadcast) in applications allocating many
    /// small matrices in the GPU. If this turns out to be the case, we could opt to
    /// enforce consensus only when the needed memory is large (e.g. >1GB) and ergo the 
    /// chance of it fitting into some GPU's memory but not others isn't negligible.
    assertAllNodesAgreeThat(matrFitsInMem, msg, vars, caller);
}

void validate_newSuperOpParams(int numQubits, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!global_isValidationEnabled)
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

    // this validation is called AFTER the caller has checked for failed
    // allocs and (in that scenario) freed every pointer, but does not 
    // overwrite any pointers to nullptr, so the failed alloc is known.
    // This is only safe to do so (rather than making the caller set ptrs
    // to nullptr) because the struct contains only 1D pointers (no nesting)

    tokenSubs vars = {{"${NUM_BYTES}", mem_getLocalSuperOpMemoryRequired(op.numQubits)}};

    // we expensively get node consensus about malloc failure, in case of heterogeneous hardware/loads,
    // but we avoid this if validation is anyway disabled
    if (!global_isValidationEnabled)
        return;

    // assert CPU array of rows was alloc'd successfully
    auto msg = (isInKrausMap)?
        report::NEW_KRAUS_MAPS_SUPER_OP_CPU_ELEMS_ALLOC_FAILED:
        report::NEW_SUPER_OP_CPU_ELEMS_ALLOC_FAILED;

    // note .cpuElems size is not included in error msg; fine since quadratically smaller than .cpuElemsFlat
    bool isAlloc = mem_isAllocated(op.cpuElemsFlat) && mem_isOuterAllocated(op.cpuElems);
    assertAllNodesAgreeThat(isAlloc, msg, vars, caller);

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
    if (!global_isValidationEnabled)
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
    if (!global_isValidationEnabled)
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

void validate_superOpDimMatchesTargs(SuperOp op, int numTargets, const char* caller) {

    tokenSubs vars = {{"${OP_QUBITS}", op.numQubits}, {"{NUM_TARGS}", numTargets}};
    assertThat(op.numQubits == numTargets, report::SUPER_OP_SIZE_MISMATCHES_NUM_TARGETS, vars, caller);
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
    if (!global_isValidationEnabled)
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

    // unlike other post-creation allocation validation, this function
    // expects that when allocation failed and the heap fields have already
    // been cleared, that any nested field (like map.matrices) has had the
    // outer pointer set to null. Otherwise, we would illegally attempt to
    // enumerate the outer pointer to check non-null-ness of inner pointers,
    // which would segmentation fault after the outer pointer was freed!
    // Ergo, we know map.matrices=nullptr whenever anything else failed
    // (and is nullptr), so we must check it last so as not to false report 
    // it as the cause of the failure!

    // we expensively get node consensus about malloc failure, in case of heterogeneous hardware/loads,
    // but we avoid this if validation is anyway disabled
    if (!global_isValidationEnabled)
        return;

    // prior validation gaurantees this will not overflow
    qindex matrListMem = map.numMatrices * mem_getLocalMatrixMemoryRequired(map.numQubits, true, 1);
    tokenSubs vars = {
        {"${NUM_BYTES}",    matrListMem},
        {"${NUM_MATRICES}", map.numMatrices},
        {"${NUM_QUBITS}",   map.numQubits}};

    // assert the teeny-tiny heap flag was alloc'd
    assertAllNodesAgreeThat(mem_isAllocated(map.isApproxCPTP), report::NEW_HEAP_FLAG_ALLOC_FAILED, {{"${NUM_BYTES}", sizeof(*(map.isApproxCPTP))}}, caller);

    // assert that the superoperator itself was allocated (along with its own heap fields)
    bool isInKrausMap = true;
    assertNewSuperOpAllocs(map.superop, isInKrausMap, caller);

    // assert the list of Kraus operator matrices, and all matrices nad rows therein, were allocated
    // (this must be done last, since caller sets .matrices=nullptr) whenever an inner alloc failed
    bool krausAreAlloc = mem_isOuterAllocated(map.matrices);
    assertAllNodesAgreeThat(krausAreAlloc, report::NEW_KRAUS_MAP_CPU_MATRICES_ALLOC_FAILED, vars, caller);
}

void validate_newInlineKrausMapDimMatchesVectors(int numQubits, int numOperators, vector<vector<vector<qcomp>>> matrices, const char* caller) {

    // avoid potentially expensive matrix enumeration if validation is anyway disabled
    if (!global_isValidationEnabled)
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
    if (!global_isValidationEnabled)
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

    // assert isApproxCPTP heap flag allocated, and that is has a valid value
    assertThat(mem_isAllocated(map.isApproxCPTP), report::INVALID_HEAP_FLAG_PTR, caller);

    // and that its value is a boolean
    int flag = *map.isApproxCPTP;
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

    // avoid expensive CPTP check (and do not overwrite .isApproxCPTP) if validation is anyway disabled
    if (isNumericalValidationDisabled())
        return;

    // use existing CPTPness or calculate afresh
    assertThat(util_isCPTP(map, global_validationEpsilon), report::KRAUS_MAP_NOT_CPTP, caller);
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

        /// @todo we can only display the ascii code of unrecognised characters,
        /// because tokenSubs only accepts integers (not chars/substrings). Fix this!
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
    assertThat(numPaulis == numIndices, report::NEW_PAULI_STR_DIFFERENT_NUM_CHARS_AND_INDICES, vars, caller);
}



/*
 * EXISTING PAULI STRING
 */

extern int paulis_getPauliAt(PauliStr str, int ind);
extern int paulis_getIndOfLefmostNonIdentityPauli(PauliStr str);

void validate_pauliStrTargets(Qureg qureg, PauliStr str, const char* caller) {

    // avoid producing a list of targets which requires enumerating all bits
    int maxTarg = paulis_getIndOfLefmostNonIdentityPauli(str);

    tokenSubs vars = {{"${MAX_TARG}", maxTarg}, {"${QUREG_QUBITS}", qureg.numQubits}};
    assertThat(maxTarg < qureg.numQubits, report::PAULI_STR_TARGETS_EXCEED_QUREG, vars, caller);
}

void validate_controlsAndPauliStrTargets(Qureg qureg, int* ctrls, int numCtrls, PauliStr str, const char* caller) {

    // validate targets and controls in isolation
    validate_pauliStrTargets(qureg, str, caller);
    validate_controls(qureg, ctrls, numCtrls, caller);

    // validate that they do not overlap (i.e. str has only I at ctrls, never X Y Z)
    const int I = 0;
    for (int n=0; n<numCtrls; n++)
        assertThat(paulis_getPauliAt(str, ctrls[n]) == I, report::PAULI_STR_OVERLAPS_CONTROLS, caller);
}

void validate_controlAndPauliStrTargets(Qureg qureg, int ctrl, PauliStr str, const char* caller) {

    validate_controlsAndPauliStrTargets(qureg, &ctrl, 1, str, caller);
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

    // this validation is called AFTER the caller has checked for failed
    // allocs and (in that scenario) freed every pointer, but does not 
    // overwrite any pointers to nullptr, so the failed alloc is known.
    // This is only safe to do so (rather than making the caller set ptrs
    // to nullptr) because the struct contains only 1D pointers (no nesting)

    assertThat(
        mem_isAllocated(sum.strings), report::NEW_PAULI_STR_SUM_STRINGS_ALLOC_FAILED, 
        {{"${NUM_TERMS}", sum.numTerms}, {"${NUM_BYTES}", numBytesStrings}}, caller);

    assertThat(
        mem_isAllocated(sum.coeffs), report::NEW_PAULI_STR_SUM_COEFFS_ALLOC_FAILED, 
        {{"${NUM_TERMS}", sum.numTerms}, {"${NUM_BYTES}", numBytesCoeffs}}, caller);

    assertThat(
        mem_isAllocated(sum.isApproxHermitian), report::NEW_HEAP_FLAG_ALLOC_FAILED, 
        {{"${NUM_BYTES}", sizeof(*(sum.isApproxHermitian))}}, caller);
}



/*
 * PAULI STRING SUM PARSING
 */

void validate_parsedPauliStrSumLineIsInterpretable(bool isInterpretable, string line, qindex lineIndex, const char* caller) {

    /// @todo we cannot yet report 'line' because tokenSubs so far only accepts integers :(

    tokenSubs vars = {{"${LINE_NUMBER}", lineIndex + 1}}; // line numbers begin at 1
    assertThat(isInterpretable, report::PARSED_PAULI_STR_SUM_UNINTERPRETABLE_LINE, vars, caller);
}

void validate_parsedPauliStrSumLineHasConsistentNumPaulis(int numPaulis, int numLinePaulis, string line, qindex lineIndex, const char* caller) {

    /// @todo we cannot yet report 'line' because tokenSubs so far only accepts integers :(

    tokenSubs vars = {
        {"${NUM_PAULIS}",      numPaulis},
        {"${NUM_LINE_PAULIS}", numLinePaulis},
        {"${LINE_NUMBER}",      lineIndex + 1}}; // line numbers begin at 1
    assertThat(numPaulis == numLinePaulis, report::PARSED_PAULI_STR_SUM_INCONSISTENT_NUM_PAULIS_IN_LINE, vars, caller);
}

void validate_parsedPauliStrSumCoeffIsValid(bool isCoeffValid, string line, qindex lineIndex, const char* caller) {

    /// @todo we cannot yet report 'line' because tokenSubs so far only accepts integers :(

    tokenSubs vars = {{"${LINE_NUMBER}", lineIndex + 1}}; // lines begin at 1
    assertThat(isCoeffValid, report::PARSED_PAULI_STR_SUM_COEFF_IS_INVALID, vars, caller);
}

void validate_parsedStringIsNotEmpty(bool stringIsNotEmpty, const char* caller) {

    assertThat(stringIsNotEmpty, report::PARSED_STRING_IS_EMPTY, caller);
}



/*
 * EXISTING PAULI STRING SUMS
 */

extern bool paulis_containsXOrY(PauliStrSum sum);
extern int paulis_getIndOfLefmostNonIdentityPauli(PauliStrSum sum);

void validate_pauliStrSumFields(PauliStrSum sum, const char* caller) {

    assertThat(sum.numTerms > 0, report::INVALID_PAULI_STR_SUM_FIELDS, {{"${NUM_TERMS}", sum.numTerms}}, caller);

    assertThat(mem_isAllocated(sum.coeffs),  report::INVALID_PAULI_STR_HEAP_PTR, caller);
    assertThat(mem_isAllocated(sum.strings), report::INVALID_PAULI_STR_HEAP_PTR, caller);

    assertThat(mem_isAllocated(sum.isApproxHermitian), report::INVALID_HEAP_FLAG_PTR, caller);

    // assert isApproxHermitian has valid value
    int flag = *sum.isApproxHermitian;
    tokenSubs vars = {
        {"${BAD_FLAG}", flag}, 
        {"${UNKNOWN_FLAG}", validate_STRUCT_PROPERTY_UNKNOWN_FLAG}};
    assertThat(flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_HEAP_FLAG_VALUE, vars, caller);
}

void validate_pauliStrSumIsHermitian(PauliStrSum sum, const char* caller) {

    // avoid expensive hermiticity check (and do not overwrite .isApproxHermitian) if validation is anyway disabled
    if (isNumericalValidationDisabled())
        return;

    // consult existing Hermiticity or compute it afresh
    assertThat(util_isHermitian(sum, global_validationEpsilon), report::PAULI_STR_SUM_NOT_HERMITIAN, caller);
}

void validate_pauliStrSumTargets(PauliStrSum sum, Qureg qureg, const char* caller) {

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(sum);
    int minNumQb = maxInd + 1;

    tokenSubs vars = {
        {"${NUM_QUREG_QUBITS}", qureg.numQubits},
        {"${MAX_IND}", maxInd}, 
        {"${NUM_PSS_QUBITS}", minNumQb}};

    assertThat(qureg.numQubits >= minNumQb, report::PAULI_STR_SUM_EXCEEDS_QUREG_NUM_QUBITS, vars, caller);
}

void validate_pauliStrSumCanInitMatrix(FullStateDiagMatr matr, PauliStrSum sum, const char* caller) {

    assertThat(!paulis_containsXOrY(sum), report::PAULI_STR_SUM_NOT_ALL_I_Z, caller);

    int maxInd = paulis_getIndOfLefmostNonIdentityPauli(sum);
    int minNumQb = maxInd + 1;

    tokenSubs vars = {
        {"${NUM_MATR_QUBITS}", matr.numQubits},
        {"${MAX_IND}", maxInd}, 
        {"${NUM_PSS_QUBITS}", minNumQb}};

    assertThat(matr.numQubits >= minNumQb, report::PAULI_STR_SUM_EXCEEDS_MATR_NUM_QUBITS, vars, caller);
}



/*
 * BASIS STATE INDICES
 */

void validate_basisStateIndex(Qureg qureg, qindex ind, const char* caller) {

    qindex maxIndExcl = powerOf2(qureg.numQubits);

    tokenSubs vars = {
        {"${STATE_IND}",  ind},
        {"${NUM_QB}", qureg.numQubits},
        {"${NUM_STATES}", maxIndExcl}};

    assertThat(ind >= 0 && ind < maxIndExcl, report::INVALID_BASIS_STATE_INDEX, vars, caller);
}

void validate_basisStateRowCol(Qureg qureg, qindex row, qindex col, const char* caller) {

    qindex maxIndExcl = powerOf2(qureg.numQubits);

    tokenSubs vars = {
        {"${ROW_IND}",  row},
        {"${COL_IND}",  col},
        {"${NUM_QB}", qureg.numQubits},
        {"${NUM_STATES}", maxIndExcl}};

    bool valid = row >= 0 && row < maxIndExcl && col >= 0 && col < maxIndExcl;
    assertThat(valid, report::INVALID_BASIS_STATE_ROW_OR_COL, vars, caller);
}

void validate_basisStateIndices(Qureg qureg, qindex startInd, qindex numInds, const char* caller) {

    assertThat(
        startInd >= 0 && startInd < qureg.numAmps, 
        report::INVALID_STARTING_BASIS_STATE_INDEX, 
        {{"${START_IND}", startInd}, {"${MAX_IND_EXCL}", qureg.numAmps}, {"${NUM_QB}", qureg.numQubits}},
        caller);

    // permit numInds=0
    assertThat(
        numInds >= 0 && numInds <= qureg.numAmps,
        report::INVALID_NUM_BASIS_STATE_INDICES, 
        {{"${NUM_INDS}", numInds}, {"${MAX_NUM_INDS_INCL}", qureg.numAmps}, {"${NUM_QB}", qureg.numQubits}}, caller);

    qindex endIndExcl = startInd + numInds;
    tokenSubs vars = {
        {"${NUM_QB}",       qureg.numQubits},
        {"${START_IND}",    startInd},
        {"${NUM_INDS}",     numInds},
        {"${MAX_IND_EXCL}", qureg.numAmps},
        {"${END_IND_EXCL}", endIndExcl}};
       
    assertThat(
        endIndExcl <= qureg.numAmps, 
        report::INVALID_ENDING_BASIS_STATE_INDEX,  
        vars, caller);
}

void validate_basisStateRowCols(Qureg qureg, qindex startRow, qindex startCol, qindex numRows, qindex numCols, const char* caller) {

    qindex maxRowOrColExcl = powerOf2(qureg.numQubits);

    assertThat(
        (startRow >= 0 && startRow < maxRowOrColExcl) && 
        (startCol >= 0 && startCol < maxRowOrColExcl), 
        report::INVALID_STARTING_BASIS_ROW_OR_COL, 
        {{"${START_ROW}", startRow}, {"${START_COL}", startCol}, {"${MAX_IND_EXCL}", maxRowOrColExcl}, {"${NUM_QB}", qureg.numQubits}}, 
        caller);

    assertThat(
        (numRows > 0 && numRows <= maxRowOrColExcl) &&
        (numCols > 0 && numCols <= maxRowOrColExcl),
        report::INVALID_NUM_BASIS_ROWS_OR_COLS, 
        {{"${NUM_ROWS}", numRows}, {"${NUM_COLS}", numCols}, {"${MAX_NUM_INCL}", maxRowOrColExcl}, {"${NUM_QB}", qureg.numQubits}}, 
        caller);

    qindex endRowExcl = startRow + numRows;
    qindex endColExcl = startCol + numCols;

    tokenSubs vars = {
        {"${NUM_QB}", qureg.numQubits},
        {"${START_ROW}", startRow}, 
        {"${START_COL}", startCol},
        {"${NUM_ROWS}", numRows},
        {"${NUM_COLS}", numCols},
        {"${END_ROW_EXCL}", endRowExcl},
        {"${END_COL_EXCL}", endColExcl},
        {"${MAX_END_INCL}", maxRowOrColExcl}};
       
    assertThat(
        endRowExcl <= maxRowOrColExcl && endColExcl <= maxRowOrColExcl,
        report::INVALID_ENDING_BASIS_ROW_OR_COL,  
        vars, caller);
}

void validate_localAmpIndices(Qureg qureg, qindex localStartInd, qindex numInds, const char* caller) {

    // note that localStartInd and numInds can validly DIFFER between nodes,
    // so we use assertAllNodesAgreeThat() in lieu of assertThat()

    tokenSubs baseVars = {
        {"${NUM_QB}", qureg.numQubits},
        {"${NUM_AMPS_TOTAL}", qureg.numAmps},
        {"${MAX_IND_EXCL}", qureg.numAmpsPerNode},
        {"${NUM_NODES}", qureg.numNodes}
    };

    // when numInds=0, we permit startInd to be anything (even something invalid)
    if (numInds == 0)
        return;

    // unlike validate_basisStateIndices(), we limit to #amps per node
    tokenSubs firstVars = baseVars;
    firstVars["${START_IND}"] = localStartInd;
    assertAllNodesAgreeThat(localStartInd >= 0 && localStartInd < qureg.numAmpsPerNode, report::INVALID_STARTING_LOCAL_AMP_INDEX, firstVars, caller);

    // unlike validate_basisStateIndices(), we permit numInds == 0
    tokenSubs secondVars = baseVars;
    secondVars["${NUM_INDS}"] = numInds;
    assertAllNodesAgreeThat(numInds >= 0 && numInds <= qureg.numAmps, report::INVALID_NUM_LOCAL_AMP_INDICES, secondVars, caller);

    qindex endIndExcl = localStartInd + numInds;
    baseVars["${START_IND}"] = localStartInd;
    baseVars["${NUM_INDS}"] = numInds;
    baseVars["${END_IND_EXCL}"] = endIndExcl;
    assertAllNodesAgreeThat(endIndExcl <= qureg.numAmpsPerNode, report::INVALID_ENDING_LOCAL_AMP_INDEX, baseVars, caller);
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
    assertThat(numQubits >= (canNumIncludeZero? 0 : 1), msgNegNum, {{"${NUM_QUBITS}", numQubits}}, caller);
    assertThat(numQubits <= qureg.numQubits, msgNumExceedsQureg, {{"${NUM_QUBITS}", numQubits}, {"${QUREG_QUBITS}", qureg.numQubits}}, caller);

    if (numQubits > 0)
        assertThat(mem_isAllocated(qubits), msgNullPtr, caller);

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
void validate_twoTargets(Qureg qureg, int target1, int target2, const char* caller) {

    int targs[] = {target1, target2};
    validate_targets(qureg, targs, 2, caller);
}

void validate_controls(Qureg qureg, int* ctrls, int numCtrls, const char* caller) {

    // it is fine to have zero controls
    bool numCanBeZero = true;

    assertValidQubits(
        qureg, ctrls, numCtrls, numCanBeZero, caller, 
        report::NEGATIVE_NUM_CONTROLS, report::NUM_CONTROLS_EXCEEDS_QUREG_SIZE, report::NON_EMPTY_CONTROL_LIST_WAS_NULL_PTR,
        report::INVALID_CONTROL_QUBIT, report::DUPLICATE_CONTROL_QUBITS);
}

void validate_controlsAndTargets(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, const char* caller) {

    // validate controls and targets in isolation
    validate_targets(qureg, targs, numTargs, caller);
    validate_controls(qureg, ctrls, numCtrls, caller);

    // validate that they do not intersect
    assertThat(areQubitsDisjoint(ctrls, numCtrls, targs, numTargs), report::CONTROLS_OVERLAP_TARGETS, caller);
}
void validate_controlAndTarget(Qureg qureg, int ctrl, int targ, const char* caller) {

    validate_controlsAndTargets(qureg, &ctrl, 1, &targ, 1, caller);
}
void validate_controlAndTargets(Qureg qureg, int ctrl, int* targs, int numTargs, const char* caller) {

    validate_controlsAndTargets(qureg, &ctrl, 1, targs, numTargs, caller);
}
void validate_controlsAndTarget(Qureg qureg, int* ctrls, int numCtrls, int targ, const char* caller) {

    validate_controlsAndTargets(qureg, ctrls, numCtrls, &targ, 1, caller);
}
void validate_controlAndTwoTargets(Qureg qureg, int ctrl, int targ1, int targ2, const char* caller) {

    int targs[] = {targ1, targ2};
    validate_controlsAndTargets(qureg, &ctrl, 1, targs, 2, caller);
}
void validate_controlsAndTwoTargets(Qureg qureg, int* ctrls, int numCtrls, int targ1, int targ2, const char* caller) {

    int targs[] = {targ1, targ2};
    validate_controlsAndTargets(qureg, ctrls, numCtrls, targs, 2, caller);
}

void validate_controlStates(int* states, int numCtrls, const char* caller) {

    // states is permittedly unallocated (nullptr) even when numCtrls != 0
    if (!mem_isAllocated(states))
        return;

    for (int n=0; n<numCtrls; n++)
        assertThat(states[n] == 0 || states[n] == 1, report::INVALID_CONTROL_STATE, {{"${INDEX}", n}, {"${STATE}", states[n]}}, caller);
}

void validate_controlsMatchStates(int numCtrls, int numStates, const char* caller) {

    // only invocable by the C++ interface
    tokenSubs vars = {
        {"${NUM_CTRLS}",  numCtrls},
        {"${NUM_STATES}", numStates}};

    assertThat(numCtrls == numStates, report::DIFFERENT_NUM_CTRLS_AND_STATES, vars, caller);
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

void validate_measurementOutcomeProbNotZero(int outcome, qreal prob, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo report 'prob' once validation reporting can handle floats

    assertThat(prob >= global_validationEpsilon, report::ONE_QUBIT_MEASUREMENT_OUTCOME_IMPOSSIBLY_UNLIKELY, {{"${OUTCOME}", outcome}}, caller);
}

void validate_measurementOutcomesProbNotZero(int* outcomes, int numQubits, qreal prob, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo report 'prob' and 'outcomes' (as binary sequence) once validation reporting can handle floats
    qindex outcomeValue = getIntegerFromBits(outcomes, numQubits);

    assertThat(prob >= global_validationEpsilon, report::MANY_QUBIT_MEASUREMENT_OUTCOME_IMPOSSIBLY_UNLIKELY, {{"${OUTCOME_VALUE}", outcomeValue}}, caller);
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

void validate_measurementProbsAreNormalised(vector<qreal> probs, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo include 'total' in validation msg once supported

    qreal total = util_getSum(probs);
    qreal dist = std::abs(total - 1);
    assertThat(dist <= global_validationEpsilon, report::OUTCOME_PROBS_DO_NOT_SUM_TO_ONE, caller);
}

void validate_measurementOutcomesMatchTargets(int numQubits, int numOutcomes, const char* caller) {

    // invoked only by the C++ user interface
    tokenSubs vars = {
        {"${NUM_QUBITS}",    numQubits},
        {"${NUM_OUTCOMES}",  numOutcomes}};

    assertThat(numQubits == numOutcomes, report::MEASUREMENT_OUTCOMES_MISMATCH_NUM_TARGETS, caller);
}



/*
 * MISC GATE PARAMETERS
 */

void validate_rotationAxisNotZeroVector(qreal x, qreal y, qreal z, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    qreal norm = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));

    assertThat(norm > global_validationEpsilon, report::ROTATION_AXIS_VECTOR_IS_ZERO, caller);
}

void validate_mixedAmpsFitInNode(Qureg qureg, int numTargets, const char* caller) {

    // only relevant to distributed quregs
    if (!qureg.isDistributed)
        return;

    // note the number of mixed amplitudes is independent of whether
    // qureg is a density matrix or not (consider unitaries)
    qindex numTargAmps = powerOf2(numTargets);

    tokenSubs vars = {
        {"${NUM_TARGS}",        numTargets},
        {"${NUM_TARG_AMPS}",    numTargAmps},
        {"${NUM_NODES}",               qureg.numNodes},
        {"${NUM_QUREG_QUBITS}",        qureg.numQubits},
        {"${NUM_QUREG_AMPS_PER_NODE}", qureg.numAmpsPerNode}
    };

    string msg = (qureg.isDensityMatrix)?
        report::CANNOT_FIT_MIXED_DENSMATR_AMPS_INTO_SINGLE_NODE:
        report::CANNOT_FIT_MIXED_STATEVEC_AMPS_INTO_SINGLE_NODE;

    assertThat(qureg.numAmpsPerNode >= numTargAmps, msg, vars, caller);
}

void validate_trotterParams(Qureg qureg, int order, int reps, const char* caller) {

    bool isEven = (order % 2) == 0;
    assertThat(order > 0 && (isEven || order==1), report::INVALID_TROTTER_ORDER, {{"${ORDER}", order}}, caller);
    assertThat(reps > 0, report::INVALID_TROTTER_REPS, {{"${REPS}", reps}}, caller);
}



/*
 * CHANNEL PARAMETERS 
 */

void validate_probability(qreal prob, const char* caller) {

    /// @todo report 'prob' once validation reporting can handle floats

    assertThat(prob >= 0 && prob <= 1, report::INVALID_PROBABILITY, caller);
}

void validate_oneQubitDepashingProb(qreal prob, const char* caller) {

    /// @todo report 'prob' once validation reporting can handle floats

    validate_probability(prob, caller);
    assertThat(
        prob <= util_getMaxProbOfOneQubitDephasing(), 
        report::ONE_QUBIT_DEPHASING_PROB_EXCEEDS_MAXIMAL_MIXING, caller);
}

void validate_twoQubitDepashingProb(qreal prob, const char* caller) {

    /// @todo report 'prob' once validation reporting can handle floats

    validate_probability(prob, caller);
    assertThat(
        prob <= util_getMaxProbOfTwoQubitDephasing(), 
        report::TWO_QUBIT_DEPHASING_PROB_EXCEEDS_MAXIMAL_MIXING, caller);
}

void validate_oneQubitDepolarisingProb(qreal prob, const char* caller) {

    /// @todo report 'prob' once validation reporting can handle floats

    validate_probability(prob, caller);
    assertThat(
        prob <= util_getMaxProbOfOneQubitDepolarising(), 
        report::ONE_QUBIT_DEPOLARISING_PROB_EXCEEDS_MAXIMAL_MIXING, caller);
}

void validate_twoQubitDepolarisingProb(qreal prob, const char* caller) {

    /// @todo report 'prob' once validation reporting can handle floats

    validate_probability(prob, caller);
    assertThat(
        prob <= util_getMaxProbOfTwoQubitDepolarising(), 
        report::TWO_QUBIT_DEPOLARISING_PROB_EXCEEDS_MAXIMAL_MIXING, caller);
}

void validate_oneQubitDampingProb(qreal prob, const char* caller) {

    /// @todo report 'prob' once validation reporting can handle floats

    // permit one-qubit amplitude damping of any valid probability, 
    // so that it can surpass maximal mixing and induce purity
    validate_probability(prob, caller);
}

void validate_oneQubitPauliChannelProbs(qreal pX, qreal pY, qreal pZ, const char* caller) {

    validate_probability(pX, caller);
    validate_probability(pY, caller);
    validate_probability(pZ, caller);

    /// @todo report 'prob' once validation reporting can handle floats

    qreal pXYZ = pX + pY + pZ;
    assertThat(pXYZ <= 1, report::ONE_QUBIT_PAULI_CHANNEL_TOTAL_PROBS_EXCEED_ONE, caller);

    qreal probI = 1 - pX - pY - pZ;
    qreal probM = 0;
    if (pX > probM) probM = pX;
    if (pY > probM) probM = pY;
    if (pZ > probM) probM = pZ;
    assertThat(probM <= probI, report::ONE_QUBIT_PAULI_CHANNEL_PROBS_EXCEED_MAXIMAL_MIXING, caller);
}



/*
 * QUREG COMBINATION
 */

void validate_quregCanBeWorkspace(Qureg qureg, Qureg workspace, const char* caller) {

    assertThat(
        (qureg.numQubits        == workspace.numQubits       ) &&
        (qureg.isDensityMatrix  == workspace.isDensityMatrix ) &&
        (qureg.isDistributed    == workspace.isDistributed   ) &&
        (qureg.isGpuAccelerated == workspace.isGpuAccelerated),
        report::QUREG_IS_INCOMPATIBLE_WITH_WORKSPACE, caller);
}

void validate_quregsCanBeMixed(Qureg quregOut, Qureg quregIn, const char* caller) {

    // mixing must be mathematically possible; dims are compatible, but quregIn can be a statevector
    assertThat(quregOut.isDensityMatrix, report::MIXED_QUREG_NOT_DENSITY_MATRIX, caller);
    assertThat(
        quregOut.numQubits == quregIn.numQubits, report::MIXED_QUREGS_HAVE_DIFFERENT_NUM_QUBITS, 
        {{"${NUM_A}", quregOut.numQubits}, {"${NUM_B}", quregIn.numQubits}}, caller);

    // density matrices must be equally distributed (because there is insufficient buffer to broadcast),
    // but they may differ in GPU deployment (because accelerator will copy memory as is necessary)
    if (quregIn.isDensityMatrix)
        assertThat(quregOut.isDistributed == quregIn.isDistributed, report::MIXED_DENSITY_MATRICES_ARE_DIFFERENTLY_DISTRIBUED, caller);

    // the statevector can only be distributed if the density matrix is (because there is otherwise no buffer to receive broadcast)
    if (!quregOut.isDistributed)
        assertThat(!quregIn.isDistributed, report::MIXED_DENSITY_MATRIX_LOCAL_BUT_STATEVEC_DISTRIBUTED, caller);
}

void validate_quregsCanBeSuperposed(Qureg qureg1, Qureg qureg2, Qureg qureg3, const char* caller) {

    // all quregs must be statevectors
    assertThat(
        !qureg1.isDensityMatrix && !qureg2.isDensityMatrix && !qureg3.isDensityMatrix,
        report::SUPERPOSED_QUREGS_ARE_NOT_ALL_STATEVECTORS, caller);

    // and the same dimension
    int nQb = qureg1.numQubits;
    assertThat(
        qureg2.numQubits == nQb && qureg3.numQubits == nQb, 
        report::SUPERPOSED_QUREGS_HAVE_INCONSISTENT_NUM_QUBITS, caller);

    // and all the same deployment (GPU & distribution; multithreading doesn't matter)
    int isGpu = qureg1.isGpuAccelerated;
    assertThat(
        qureg2.isGpuAccelerated == isGpu && qureg3.isGpuAccelerated == isGpu, 
        report::SUPERPOSED_QUREGS_HAVE_INCONSISTENT_GPU_DEPLOYMENT, caller);

    int isDis = qureg1.isDistributed;
    assertThat(
        qureg2.isDistributed == isDis && qureg3.isDistributed == isDis, 
        report::SUPERPOSED_QUREGS_HAVE_INCONSISTENT_DISTRIBUTION, caller);
}

void validateDensMatrCanBeInitialisedToPureState(Qureg qureg, Qureg pure, const char* caller) {

    // initPureState calls mixQureg which only additionally
    // constrains that pure.isDistributed only if qureg.isDistributed

    if (pure.isDistributed)
        assertThat(qureg.isDistributed, report::INIT_DENSMATR_LOCAL_BUT_PURE_STATE_DISTRIBUTED, caller);
}

void validateStateVecCanBeInitialisedToPureState(Qureg qureg, Qureg pure, const char* caller) {

    // statevectors must be identically distributed and GPU-accelerated

    assertThat(
        qureg.isGpuAccelerated == pure.isGpuAccelerated,
        report::INIT_STATEVEC_DIFFERING_GPU_DEPLOYMENT_TO_PURE_STATE, caller);

    assertThat(
        qureg.isDistributed == pure.isDistributed,
        report::INIT_STATEVEC_DIFFERING_DISTRIBUTION_TO_PURE_STATE, caller);
}

void validate_quregCanBeInitialisedToPureState(Qureg qureg, Qureg pure, const char* caller) {

    assertThat(!pure.isDensityMatrix, report::INIT_PURE_STATE_IS_DENSMATR, caller);

    // quregs must have the same number of qubits, regardless of dimension
    assertThat(
        qureg.numQubits == pure.numQubits,
        report::INIT_QUREG_HAS_DIFFERENT_NUM_QUBITS_TO_PURE_STATE, 
        {{"${NUM_QUREG_QUBITS}", qureg.numQubits}, {"${NUM_PURE_QUBITS}", pure.numQubits}},
        caller);

    (qureg.isDensityMatrix)?
        validateDensMatrCanBeInitialisedToPureState(qureg, pure, caller):
        validateStateVecCanBeInitialisedToPureState(qureg, pure, caller);
}

void validate_quregsCanBeCloned(Qureg quregA, Qureg quregB, const char* caller) {

    // quregs must have identical sizes... 
    assertThat(
        quregA.numQubits == quregB.numQubits, report::CLONED_QUREGS_DIFFER_IN_NUM_QUBITS, 
        {{"${NUM_TARGET_QUBITS}", quregA.numQubits}, {"${NUM_COPY_QUBITS}", quregB.numQubits}}, caller);
    
    // and types...
    assertThat(quregA.isDensityMatrix == quregB.isDensityMatrix, report::CLONED_QUREGS_INCONSISTENT_TYPES, caller);

    // and distributions...
    assertThat(quregA.isDistributed == quregB.isDistributed, report::CLONED_QUREGS_HAVE_DIFFERENT_DISTRIBUTIONS, caller);

    // but density matrices may differ in GPU-deployments, because cloning 
    // invokes mixQureg which supports heterogeneous deployment
    if (quregA.isDensityMatrix)
        return;

    // statevectors however MUST be identically GPU-accelerated
    assertThat(quregA.isGpuAccelerated == quregB.isGpuAccelerated, report::CLONED_STATEVECS_HAVE_DIFFERENT_GPU_DEPLOYMENTS, caller);
}

void validate_quregsCanBeProducted(Qureg quregA, Qureg quregB, const char* caller) {

   // number of qubits must always match
    assertThat(
        quregA.numQubits == quregB.numQubits, 
        report::PRODUCTED_QUREGS_HAVE_DIFFERENT_NUM_QUBITS,
        {{"${NUM_A}", quregA.numQubits}, {"${NUM_B}", quregB.numQubits}}, caller);

    // if both are statevecs or both are densmatr...
    if (quregA.isDensityMatrix == quregB.isDensityMatrix) {

        // then they must have the same GPU accel (because falling back to CPU is slow)
        assertThat(
            quregA.isGpuAccelerated == quregB.isGpuAccelerated,
            report::PRODUCTED_SAME_TYPE_QUREGS_HAVE_DIFFERING_GPU_ACCELS, caller);

        // but their distributions may differ, so finish
        return;
    }

    // when one is a statevector and the other is a density matrix...
    Qureg dm = (quregA.isDensityMatrix)? quregA : quregB;
    Qureg sv = (quregA.isDensityMatrix)? quregB : quregA;

    // then statevec can only be distributed if density matrix is
    if (sv.isDistributed)
        assertThat(dm.isDistributed, report::PRODUCTED_STATEVEC_DISTRIB_BUT_DENSMATR_LOCAL, caller); 

    // their GPU accelerations may differ however
}

void validate_throwErrorBecauseCalcFidOfDensMatrNotYetImplemented(const char* caller) {

    assertThat(false, report::CALC_FIDELITY_OF_DENSITY_MATRICES_NOT_YET_SUPPORTED, caller);
}

void validate_fidelityIsReal(qcomp fid, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo include imag(fid) in error message when non-integers are supported

    qreal eps = REDUCTION_EPSILON_FACTOR * global_validationEpsilon;
    assertThat(util_isApproxReal(fid, eps), report::CALC_FIDELITY_NOT_APPROX_REAL, caller);
}

void validate_buresDistanceInnerProdIsNormalised(qreal mag, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo include mag in error message when non-integers are supported

    qreal eps = REDUCTION_EPSILON_FACTOR * global_validationEpsilon;
    assertThat(mag <= 1 + eps, report::CALC_BURES_DISTANCE_MAG_EXCEEDED_ONE, caller);
}

void validate_purifiedDistanceIsNormalised(qcomp fid, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo include scalars in error message when non-integers are supported
    
    qreal eps = REDUCTION_EPSILON_FACTOR * global_validationEpsilon;
    assertThat(util_isApproxReal(fid, eps), report::CALC_PURIFIED_DISTANCE_NOT_APPROX_REAL, caller);
    assertThat(std::real(fid) <= 1 + eps, report::CALC_PURIFIED_DISTANCE_REAL_EXCEEDED_ONE, caller);
}



/*
 * QUREG MODIFICATION
 */

void validate_quregRenormProbIsNotZero(qreal prob, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo include 'prob' in error message when non-integers are supported

    // note use abs(prob) of in lieu of prob; we permit the probability
    // to be negative as can happen during setQuregToRenormalized() when
    // given an invalid density-matrix. We only require the magnitude is
    // non-zero so that division doesn't numerically diverge
    assertThat(std::abs(prob) > global_validationEpsilon, report::QUREG_RENORM_PROB_IS_ZERO, caller);
}

void validate_numInitRandomPureStates(qindex numPureStates,  const char* caller) {

    assertThat(numPureStates >= 1, report::INVALID_NUM_INIT_PURE_STATES, {{"${NUM_STATES}", numPureStates}}, caller);
}



/*
 * EXPECTATION VALUES
 */

void validate_expecPauliStrValueIsReal(qcomp value, bool isDensMatr, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    /// @todo include imag(value) in error message when non-integers are supported

    string msg = (isDensMatr)?
        report::CALC_DENSMATR_EXPECTED_PAULI_STR_VALUE_WAS_NOT_APPROX_REAL:
        report::CALC_STATEVEC_EXPECTED_PAULI_STR_VALUE_WAS_NOT_APPROX_REAL;

    qreal eps = REDUCTION_EPSILON_FACTOR * global_validationEpsilon;
    assertThat(util_isApproxReal(value, eps), msg, caller);
}

void validate_expecPauliStrSumValueIsReal(qcomp value, bool isDensMatr, const char* caller) {

    if (isNumericalValidationDisabled())
        return;

    string msg = (isDensMatr)?
        report::CALC_DENSMATR_EXPECTED_PAULI_STR_SUM_VALUE_WAS_NOT_APPROX_REAL:
        report::CALC_STATEVEC_EXPECTED_PAULI_STR_SUM_VALUE_WAS_NOT_APPROX_REAL;

    qreal eps = REDUCTION_EPSILON_FACTOR * global_validationEpsilon;
    assertThat(util_isApproxReal(value, eps), msg, caller);
}

void validate_densMatrExpecDiagMatrValueIsReal(qcomp value, qcomp exponent, const char* caller) {

    // this function is only ever called to validate the output of
    // expected value calculations of hermitian diagonal matrices
    // upon density matrices (NOT statevectors) because errors in
    // the latter (due to state unnormalisation, or exponent doma-
    // in, or unintended imaginary components of the diagonal)
    // never damage the real component, which we always safely return

    if (isNumericalValidationDisabled())
        return;

    // precise comparison since non-power overload passes epxonent=1
    string msg = (exponent == qcomp(1,0))?
        report::CALC_DENSMATR_EXPECTED_DIAG_MATR_VALUE_WAS_NOT_APPROX_REAL:
        report::CALC_DENSMATR_EXPECTED_DIAG_MATR_POWER_VALUE_WAS_NOT_APPROX_REAL;

    qreal eps = REDUCTION_EPSILON_FACTOR * global_validationEpsilon;
    assertThat(util_isApproxReal(value, eps), msg, caller);
}



/*
 * PARTIAL TRACE
 */

void validate_quregCanBeReduced(Qureg qureg, int numTraceQubits, const char* caller) {

    // 0 < numTraceQubits <= numQubits is assured by validate_targets(), but
    // numTraceQubits == numQubtis is permitted there though forbidden here
    assertThat(numTraceQubits < qureg.numQubits, report::NUM_TRACE_QUBITS_EQUALS_QUREG_SIZE, caller);

    // when not distributed, there are no further restrictions
    if (!qureg.isDistributed)
        return;

    // if the reduced qureg were hypothetically permitted to be non-distributed 
    // despite qureg being distributed, then maxTr <= numQb - ceil(logNumNodes/2).
    // we presently do not support this however, and insist the reduced qureg
    // matches the input qureg's distribution, such that the maximum traced qubits
    // is simply that which still retains at least one column per node
    int maxNumTraceQubits = qureg.numQubits - qureg.logNumNodes;

    tokenSubs vars = {
        {"${NUM_TRACE_QUBITS}", numTraceQubits},
        {"${MAX_TRACE_QUBITS}", maxNumTraceQubits},
        {"${NUM_QUREG_QUBITS}", qureg.numQubits},
        {"${NUM_NODES}",        qureg.numNodes}};

    assertThat(numTraceQubits <= maxNumTraceQubits, report::NUM_TRACE_QUBITS_EXCEEDS_DISTRIBUTED_MAX, vars, caller);
}

void validate_quregCanBeSetToReducedDensMatr(Qureg out, Qureg in, int numTraceQubits, const char* caller) {

    int numRemainingQubits = in.numQubits - numTraceQubits;

    tokenSubs vars = {
        {"${IN_QUBITS}",  in.numQubits},
        {"${OUT_QUBITS}", out.numQubits},
        {"${TRACE_QUBITS}",  numTraceQubits},
        {"${RETAIN_QUBITS}", numRemainingQubits}};

    assertThat(out.numQubits == numRemainingQubits, report::NUM_TRACE_QUBITS_INCONSISTENT_WITH_REDUCED_QUREG, vars, caller);

    assertThat(out.isDistributed    == in.isDistributed,    report::REDUCED_QUREG_DIFFERING_DISTRIBUTION_TO_IN_QUREG,   caller);
    assertThat(out.isGpuAccelerated == in.isGpuAccelerated, report::REDUCED_QUREG_DIFFERING_GPU_DEPLOYMENT_TO_IN_QUREG, caller);
}



/*
 * FILE IO
 */

void validate_canReadFile(string fn, const char* caller) {

    /// @todo embed filename into error message when tokenSubs is updated to permit strings
    assertThat(parser_canReadFile(fn), report::CANNOT_READ_FILE, caller);
}



/*
 * TEMPORARY ALLOCATIONS
 */

void validate_tempAllocSucceeded(bool succeeded, qindex numElems, qindex numBytesPerElem, const char* caller) {

    // avoid showing total bytes in case it overflows
    tokenSubs vars = {
        {"${NUM_ELEMS}", numElems},
        {"${NUM_BYTES_PER_ELEM}", numBytesPerElem}};

    assertThat(succeeded, report::TEMP_ALLOC_FAILED, vars, caller);
}
