/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/environment.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/paulis.h"

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
     * QUREG INITIALISATIONS
     */

    string INVALID_STATE_INDEX = 
        "Classical state index ${STATE_IND} is invalid for the given ${NUM_QUBITS} qubit Qureg. Index must be greater than or equal to zero, and cannot equal nor exceed the number of unique classical states (2^${NUM_QUBITS} = ${NUM_STATES}).";


    /*
     * EXISTING QUREG
     */

    string INVALID_QUREG_FIELDS = 
        "Invalid Qureg; invalid or incompatible fields isDensityMatrix=${DENS_MATR}, numQubits=${NUM_QUBITS}, numAmps=${NUM_AMPS}. It is likely this Qureg was not initialised with createQureg().";


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
        "Cannot create a diagonal matrix of ${NUM_QUBITS} distributed between ${NUM_NODES} nodes because each node would contain fewer than one element. The minimum number of qubits in such a matrix is ${MIN_QUBITS}. Consider disabling distribution for this matrix.";


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
        "Attempted allocation of memory (${NUM_BYTES} bytes in RAM) failed.";

    string NEW_MATRIX_GPU_ELEMS_ALLOC_FAILED = 
        "Attempted allocation of GPU memory (${NUM_BYTES} bytes in VRAM) failed.";

    string NEW_MATRIX_IS_UNITARY_FLAG_ALLOC_FAILED = 
        "Attempted allocation of the 'isUnitary' heap field (a mere ${NUM_BYTES} bytes) failed! This is unprecedentedly unlikely - go and have your fortune read at once.";


    string NEW_DISTRIB_MATRIX_IN_NON_DISTRIB_ENV = 
        "Cannot distribute a matrix in a non-distributed environment.";

    string INVALID_OPTION_FOR_MATRIX_IS_DISTRIB = 
        "Argument 'useDistrib' must be 1 or 0 to respectively indicate whether or not to distribute the new matrix, or ${AUTO_DEPLOYMENT_FLAG} to let QuEST choose automatically.";



    /*
     * MATRIX INITIALISATION
     */

    string MATRIX_NEW_ELEMS_CONTAINED_GPU_SYNC_FLAG = 
        "The new matrix elements contained a reserved, forbidden value as the first element, used internally to detect that whether GPU memory has not synchronised. The value was intended to be extremely unlikely to be used by users - go buy a lottery ticket! If you insist on using this value in the first element, add a numerically insignificant perturbation.";


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


    /*
     * EXISTING MATRIX
     */

    string INVALID_COMP_MATR_1_FIELDS =
        "Invalid CompMatr1. Targeted ${NUM_QUBITS} qubits (instead of 1) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 2x2). It is likely this matrix was not initialised with getCompMatr1().";

    string INVALID_COMP_MATR_2_FIELDS =
        "Invalid CompMatr2. Targeted ${NUM_QUBITS} qubits (instead of 2) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 4x4). It is likely this matrix was not initialised with getCompMatr2().";

    string INVALID_COMP_MATR_FIELDS =
        "Invalid CompMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createCompMatr().";

    string INVALID_COMP_MATR_CPU_MEM_ALLOC =
        "Invalid CompMatr. One or more rows of the 2D CPU memory (RAM) was seemingly unallocated. It is likely this matrix was not initialised with createCompMatr().";

    string INVALID_COMP_MATR_GPU_MEM_ALLOC =
        "Invalid CompMatr. The GPU memory (VRAM) was seemingly unallocated. It is likely this matrix was not initialised with createCompMatr().";


    string INVALID_DIAG_MATR_1_FIELDS =
        "Invalid DiagMatr1. Targeted ${NUM_QUBITS} qubits (instead of 1) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 2x2). It is likely this matrix was not initialised with getDiagMatr1().";

    string INVALID_DIAG_MATR_2_FIELDS =
        "Invalid DiagMatr2. Targeted ${NUM_QUBITS} qubits (instead of 2) and had a dimension of ${NUM_ROWS}x${NUM_ROWS} (instead of 4x4). It is likely this matrix was not initialised with getDiagMatr2().";

    string INVALID_DIAG_MATR_FIELDS =
        "Invalid DiagMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createDiagMatr().";

    string INVALID_DIAG_MATR_CPU_MEM_ALLOC =
        "Invalid DiagMatr. The CPU memory (RAM) was seemingly unallocated. It is likely this matrix was not initialised with createDiagMatr().";

    string INVALID_DIAG_MATR_GPU_MEM_ALLOC =
        "Invalid DiagMatr. The GPU memory (VRAM) was seemingly unallocated. It is likely this matrix was not initialised with createDiagMatr().";


    string INVALID_FULL_STATE_DIAG_MATR_FIELDS = 
        "Invalid FullStateDiagMatr. Targeted ${NUM_QUBITS} qubits and had a dimension of ${NUM_ROWS}x${NUM_ROWS}. It is likely this matrix was not created with createFullStateDiagMatr().";

    string INVALID_FULL_STATE_DIAG_MATR_CPU_MEM_ALLOC =
        "Invalid FullStateDiagMatr. The CPU memory (RAM) was seemingly unallocated. It is likely this matrix was not initialised with createFullStateDiagMatr().";

    string INVALID_FULL_STATE_DIAG_MATR_GPU_MEM_ALLOC =
        "Invalid FullStateDiagMatr. The GPU memory (VRAM) was seemingly unallocated. It is likely this matrix was not initialised with createFullStateDiagMatr().";


    string COMP_MATR_NOT_SYNCED_TO_GPU = 
        "The CompMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please call syncCompMatr() after manually modifying elements, or overwrite all elements with setCompMatr() which automatically synchronises.";

    string DIAG_MATR_NOT_SYNCED_TO_GPU = 
        "The DiagMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please first call syncDiagMatr() after manually modifying elements, or overwrite all elements with setDiagMatr().";

    string FULL_STATE_DIAG_MATR_NOT_SYNCED_TO_GPU = 
        "The FullStateDiagMatr has yet not been synchronised with its persistent GPU memory, so potential changes to its elements are being ignored. Please first call syncFullStateDiagMatr() after manually modifying elements, or overwrite elements in batch with setFullStateDiagMatr().";


    string INVALID_MATRIX_IS_UNITARY_FLAG = 
        "The 'isUnitary' field of the given matrix had invalid value ${BAD_FLAG}, suggesting it was manually modified. Valid values are 0, 1 and ${UNKNOWN_FLAG} (to indicate that unitarity is not yet known, deferring evaluation) although this need never be modified by the user.";

    string MATRIX_NOT_UNITARY = 
        "The given matrix was not (approximately) unitary.";


    string FULL_STATE_DIAG_MATR_MISMATCHES_QUREG_DIM =
        "The given FullStateDiagMatr operates upon a different number of qubits (${NUM_MATR_QUBITS}) than exists in the Qureg (${NUM_QUREG_QUBITS}).";

    string FULL_STATE_DIAG_MATR_IS_DISTRIB_BUT_QUREG_ISNT =
        "The given FullStateDiagMatr is distributed but the Qureg is not, which is forbidden. Consider disabling distribution for this matrix via createCustomFullStateDiagMatr().";


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

    
    /*
     * OPERATOR PARAMETERS
     */
    
    string INVALID_TARGET_QUBIT = 
        "Invalid target qubit (${TARGET}). Must be greater than or equal to zero, and less than the number of qubits in the Qureg (${NUM_QUBITS}).";


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

void validate_enable() {
    isValidationEnabled = true;
}
void validate_disable() {
    isValidationEnabled = false;
}
bool validate_isEnabled() {
    return isValidationEnabled;
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

    string newMsg = getStringWithSubstitutedVars(msg, vars);
    assertThat(valid, newMsg, func);
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

void validate_numReportedItems(qindex num, const char* caller) {

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
    assertAllNodesAgreeThat(qureg.cpuAmps != nullptr, report::NEW_QUREG_CPU_AMPS_ALLOC_FAILED, caller);

    if (qureg.isGpuAccelerated)
        assertAllNodesAgreeThat(qureg.gpuAmps != nullptr, report::NEW_QUREG_GPU_AMPS_ALLOC_FAILED, caller);

    if (qureg.isDistributed)
        assertAllNodesAgreeThat(qureg.cpuCommBuffer != nullptr, report::NEW_QUREG_CPU_COMM_BUFFER_ALLOC_FAILED, caller);

    if (qureg.isDistributed && qureg.isGpuAccelerated)
        assertAllNodesAgreeThat(qureg.gpuCommBuffer != nullptr, report::NEW_QUREG_GPU_COMM_BUFFER_ALLOC_FAILED, caller);
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

void assertMatrixFitsInGpuMem(int numQubits, bool isDense, int isDistrib, int isEnvGpuAccel, int numEnvNodes, const char* caller) {

    // 'isDistrib' can be 0 (user-disabled, or matrix is a local type), 1 (user-forced)
    // or modeflag::USE_AUTO=-1 (automatic; the type can be distributed but the user has 
    // not forced it). Currently, only distributed diagonal matrices are supported, so 
    // isDistrib must be specifically 0 for dense matrices.
    if (isDense && isDistrib != 0)
        error_validationEncounteredUnsupportedDistributedDenseMatrix();

    // matrix GPU memory will always be allocated when env is GPU-accelerated
    if (!isEnvGpuAccel)
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

void assertNewMatrixParamsAreValid(int numQubits, int useDistrib, bool isDenseType, const char* caller) {

    // some of the below validation involves getting distributed node consensus, which
    // can be an expensive synchronisation, which we avoid if validation is anyway disabled
    if (!isValidationEnabled)
        return;

    validate_envIsInit(caller);
    QuESTEnv env = getQuESTEnv();

    assertMatrixNonEmpty(numQubits, caller);
    assertMatrixTotalNumElemsDontExceedMaxIndex(numQubits, isDenseType, caller);
    assertMatrixLocalMemDoesntExceedMaxSizeof(numQubits,  isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixNotDistributedOverTooManyNodes(numQubits, isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixFitsInCpuMem(numQubits, isDenseType, useDistrib, env.numNodes, caller);
    assertMatrixFitsInGpuMem(numQubits, isDenseType, useDistrib, env.isGpuAccelerated, env.numNodes, caller);
}

void validate_newCompMatrParams(int numQubits, const char* caller) {

    int useDistrib = 0;
    bool isDenseType = true;
    assertNewMatrixParamsAreValid(numQubits, useDistrib, isDenseType, caller);
}
void validate_newDiagMatrParams(int numQubits, const char* caller) {

    int useDistrib = 0;
    bool isDenseType = false;
    assertNewMatrixParamsAreValid(numQubits, useDistrib, isDenseType, caller);
}
void validate_newFullStateDiagMatrParams(int numQubits, int useDistrib, const char* caller) {

    bool isDenseType = false;
    assertNewMatrixParamsAreValid(numQubits, useDistrib, isDenseType, caller);
}

// T can be CompMatr, DiagMatr, FullStateDiagMatr (i.e. heap-based matrices)
template <typename T>
void assertNewMatrixAllocsSucceeded(T matr, qindex numBytes, const char* caller) {

    // we expensively get node consensus about malloc failure, in case of heterogeneous hardware/loads,
    // but we avoid this potetially expensive synchronisation if validation is anyway disabled
    // (which also avoids us enumerating the matrix rows)
    if (!isValidationEnabled)
        return;

    // assert CPU array of rows was malloc'd successfully
    tokenSubs vars = {{"${NUM_BYTES}", numBytes}};
    assertAllNodesAgreeThat(matr.cpuElems != nullptr, report::NEW_MATRIX_CPU_ELEMS_ALLOC_FAILED, vars, caller);

    // assert each CPU row was calloc'd successfully (if matrix is dense)
    if constexpr (util_isDenseMatrixType<T>())
        for (qindex r=0; r<matr.numRows; r++)
            assertAllNodesAgreeThat(matr.cpuElems[r] != nullptr, report::NEW_MATRIX_CPU_ELEMS_ALLOC_FAILED, vars, caller);
    
    // optionally assert GPU memory was malloc'd successfully
    if (getQuESTEnv().isGpuAccelerated)
        assertAllNodesAgreeThat(util_getGpuMemPtr(matr) != nullptr, report::NEW_MATRIX_GPU_ELEMS_ALLOC_FAILED, vars, caller);

    // assert the teeny-tiny isUnitary flag was alloc'd
    vars["${NUM_BYTES}"] = sizeof(*(matr.isUnitary));
    assertAllNodesAgreeThat(matr.isUnitary != nullptr, report::NEW_MATRIX_IS_UNITARY_FLAG_ALLOC_FAILED, vars, caller);
}

void validate_newMatrixAllocs(CompMatr matr, qindex numBytes, const char* caller) {
    assertNewMatrixAllocsSucceeded(matr, numBytes, caller);
}
void validate_newMatrixAllocs(DiagMatr matr, qindex numBytes, const char* caller) {
    assertNewMatrixAllocsSucceeded(matr, numBytes, caller);
}
void validate_newMatrixAllocs(FullStateDiagMatr matr, qindex numBytes, const char* caller) {
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

void validate_matrixNewElemsDontContainUnsyncFlag(qcomp firstElem, const char* caller) {

    // avoid potentially expensive node-consensus if validation anyway disabled
    if (!isValidationEnabled)
        return;

    // we permit the matrix to contain the GPU-mem-unsync flag in CPU-only mode,
    // to avoid astonishing a CPU-only user with a GPU-related error message
    if (!getQuESTEnv().isGpuAccelerated)
        return;

    // to work with distributed FullStateDiagMatr, whereby we wish to check that
    // every node's GPU elems have been overwritten (not just e.g. the root node's),
    // we need to gather expensive consensus on the validity. 
    assertAllNodesAgreeThat(!gpu_doCpuAmpsHaveUnsyncMemFlag(firstElem), report::MATRIX_NEW_ELEMS_CONTAINED_GPU_SYNC_FLAG, caller);
}



/*
 * EXISTING MATRICES
 */

// T can be CompMatr1, CompMatr2, CompMatr, DiagMatr1, DiagMatr2, DiagMatr, FullStateDiagMatr
template <class T>
void assertMatrixFieldsAreValid(T matr, int expectedNumQb, string errMsg, const char* caller) {

    qindex dim = util_getMatrixDim(matr);
    tokenSubs vars = {
        {"${NUM_QUBITS}", matr.numQubits},
        {"${NUM_ROWS}",   dim}};

    // assert correct fixed-size numQubits (caller gaurantees this passes for dynamic-size),
    // where the error message string already contains the expected numQb
    assertThat(matr.numQubits == expectedNumQb, errMsg, vars, caller);

    // validate .numQubits and .numRows or .numElems
    qindex expectedDim = powerOf2(matr.numQubits);
    assertThat(matr.numQubits >= 1, errMsg, vars, caller);
    assertThat(dim == expectedDim,  errMsg, vars, caller);

    // assert isUnitary flag is valid; it is user-modifiable (if they are naughty)!
    if constexpr (util_doesMatrixHaveIsUnitaryFlag<T>()) {
        int flag = (matr.isUnitary == nullptr)? 0 : *matr.isUnitary; // lazy segfault if failed allocation (caught below)
        assertThat(
            flag == 0 || flag == 1 || flag == validate_STRUCT_PROPERTY_UNKNOWN_FLAG, report::INVALID_MATRIX_IS_UNITARY_FLAG,
            {{"${BAD_FLAG}", flag}, {"${UNKNOWN_FLAG}", validate_STRUCT_PROPERTY_UNKNOWN_FLAG}}, caller);
    }

    // we do not bother checking slightly more involved fields like numAmpsPerNode - there's
    // no risk that they're wrong (because they're constant) unless the struct was unitialised,
    // which we have already validated against
}

// T can be CompMatr, DiagMatr, FullStateDiagMatr (i.e. matrix structs with heap pointers)
template <class T>
void assertMatrixAllocsAreValid(T matr, string cpuErrMsg, string gpuErrMsg, const char* caller) {

    // TODO:
    // I'm fairly sure this validation is completely pointless - the memory being unallocated
    // being the struct was manually initialised cannot be determined reliably by a NULL check,
    // because the pointer fields will be un-initialised i.e. random. Only malloc will sensibly
    // set the fields to NULL, which immediately post-allocation validation would capture. Eh!

    // assert CPU memory is allocated
    assertThat(matr.cpuElems != nullptr, cpuErrMsg, caller);

    // we do not check that each CPU memory row (of CompMatr; irrelevant to DiagMatr)
    // is not-NULL because it's really unlikely that inner memory wasn't allocaed but
    // outer was and this wasn't somehow caught by post-creation validation (i.e. the 
    // user manually malloc'd memory after somehow getting around const fields). Further,
    // since this function is called many times (i.e. each time the user passes a matrix
    // to a simulation function like multiQubitUnitary()), it may be inefficient to 
    // serially process each row pointer.

    // optionally assert GPU memory is allocated
    validate_envIsInit(caller);
    if (getQuESTEnv().isGpuAccelerated)
        assertThat(util_getGpuMemPtr(matr) != nullptr, gpuErrMsg, caller);
}

void validate_matrixFields(CompMatr1 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 1, report::INVALID_COMP_MATR_1_FIELDS, caller);
}
void validate_matrixFields(CompMatr2 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 2, report::INVALID_COMP_MATR_2_FIELDS, caller);
}
void validate_matrixFields(CompMatr matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_COMP_MATR_FIELDS, caller);
    assertMatrixAllocsAreValid(matr, report::INVALID_COMP_MATR_CPU_MEM_ALLOC, report::INVALID_COMP_MATR_GPU_MEM_ALLOC, caller);
}
void validate_matrixFields(DiagMatr1 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 1, report::INVALID_DIAG_MATR_1_FIELDS, caller);
}
void validate_matrixFields(DiagMatr2 matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, 2, report::INVALID_DIAG_MATR_2_FIELDS, caller);
}
void validate_matrixFields(DiagMatr matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_DIAG_MATR_FIELDS, caller);
    assertMatrixAllocsAreValid(matr, report::INVALID_DIAG_MATR_CPU_MEM_ALLOC, report::INVALID_DIAG_MATR_GPU_MEM_ALLOC, caller);
}
void validate_matrixFields(FullStateDiagMatr matr, const char* caller) {
    assertMatrixFieldsAreValid(matr, matr.numQubits, report::INVALID_FULL_STATE_DIAG_MATR_FIELDS, caller);
    assertMatrixAllocsAreValid(matr, report::INVALID_FULL_STATE_DIAG_MATR_CPU_MEM_ALLOC, report::INVALID_FULL_STATE_DIAG_MATR_GPU_MEM_ALLOC, caller);
}

// type T can be CompMatr or DiagMatr
template <class T>
void assertMatrixIsSynced(T matr, string errMsg, const char* caller) {

    // we don't need to perform any sync check in CPU-only mode
    if (util_getGpuMemPtr(matr) == nullptr)
        return;

    // check if GPU amps have EVER been overwritten; we sadly cannot check the LATEST changes were pushed though
    assertThat(gpu_haveGpuAmpsBeenSynced(util_getGpuMemPtr(matr)), errMsg, caller);
}

// type T can be CompMatr, DiagMatr, FullStateDiagMatr
template <class T> 
void ensureMatrixUnitarityIsKnown(T matr) {

    // do nothing if we already know unitarity
    if (*(matr.isUnitary) != validate_STRUCT_PROPERTY_UNKNOWN_FLAG)
        return;

    // determine local unitarity, modifying matr.isUnitary
    *(matr.isUnitary) = util_isUnitary(matr);

    // get node consensus on global unitarity, if necessary
    if constexpr (util_isDistributableMatrixType<T>())
        if (matr.isDistributed)
            *(matr.isUnitary) = comm_isTrueOnAllNodes(*(matr.isUnitary));
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
        isUnitary = util_isUnitary(matr);

    // dynamic matrices have their field consulted, which may invoke lazy eval and global synchronisation
    else {
        ensureMatrixUnitarityIsKnown(matr);
        isUnitary = *(matr.isUnitary);
    }

    assertThat(isUnitary, report::MATRIX_NOT_UNITARY, caller);
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

void validate_matrixIsCompatibleWithQureg(FullStateDiagMatr matr, Qureg qureg, const char* caller) {

    // we do not need to define this function for the other matrix types,
    // since their validation will happen through validaiton of the
    // user-given list of target qubits. But we do need to define it for
    // FullStatedDiagMatr to check both distribution compatibility, and
    // that dimensions match

    tokenSubs vars = {
        {"${MATR_NUM_QUBITS}",  matr.numQubits},
        {"${QUREG_NUM_QUBITS}", qureg.numQubits}};

    // dimensions must match
    assertThat(matr.numQubits == qureg.numQubits, report::FULL_STATE_DIAG_MATR_MISMATCHES_QUREG_DIM, vars, caller);

    // when matrix is duplicated on every node, its application is trivial
    if (!matr.isDistributed)
        return;

    // but when it's distributed, so too must be the qureg so that comm isn't necessary
    assertThat(qureg.isDistributed, report::FULL_STATE_DIAG_MATR_IS_DISTRIB_BUT_QUREG_ISNT, caller); // did not pass vars
}


/*
 * QUREG INITIALISATIONS
 */

void validate_initClassicalStateIndex(Qureg qureg, qindex ind, const char* caller) {

    qindex maxIndExcl = powerOf2(qureg.numQubits);

    tokenSubs vars = {
        {"${STATE_IND}",  ind},
        {"${NUM_QUBITS}", qureg.numQubits},
        {"${NUM_STATES}", maxIndExcl}};

    assertThat(ind >= 0 && ind < maxIndExcl, report::INVALID_STATE_INDEX, vars, caller);
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

    // note we do not bother checking whether RAM has enough memory to contain
    // the new Pauli string, because the caller to this function has already
    // been passed data of the same size (and it's unlikely the user is about
    // to max RAM), and the memory requirements scale only linearly with the
    // parameters (e.g. numTerms), unlike the exponential scaling of the memory
    // of Qureg and CompMatr, for example

    validate_newPauliStrNumPaulis(numPaulis, maxNumPaulis, caller);
    assertCorrectNumPauliCharsBeforeTerminationChar(paulis, numPaulis, caller);
    assertRecognisedNewPaulis(paulis, numPaulis, caller);
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

    assertThat(numTerms > 0, report::NEW_PAULI_STR_SUM_NON_POSITIVE_NUM_STRINGS, {{"${NUM_TERMS}", numTerms}}, caller);
}

void validate_newPauliStrSumMatchingListLens(qindex numStrs, qindex numCoeffs, const char* caller) {

    tokenSubs vars = {{"${NUM_STRS}", numStrs}, {"${NUM_COEFFS}", numCoeffs}};
    assertThat(numStrs == numCoeffs, report::NEW_PAULI_STR_SUM_DIFFERENT_NUM_STRINGS_AND_COEFFS, vars, caller);
}

void validate_newPauliStrSumAllocs(PauliStrSum sum, qindex numBytesStrings, qindex numBytesCoeffs, const char* caller) {

    assertThat(
        sum.strings != nullptr, report::NEW_PAULI_STR_SUM_STRINGS_ALLOC_FAILED, 
        {{"${NUM_TERMS}", sum.numTerms}, {"${NUM_BYTES}", numBytesStrings}}, caller);

    assertThat(
        sum.coeffs != nullptr, report::NEW_PAULI_STR_SUM_COEFFS_ALLOC_FAILED, 
        {{"${NUM_TERMS}", sum.numTerms}, {"${NUM_BYTES}", numBytesCoeffs}}, caller);
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

    // no point checking sum's pointers are not nullptr; see the explanation in validate_quregFields()
}



/*
 * OPERATOR PARAMETERS
 */

void validate_target(Qureg qureg, int target, const char* caller) {

    tokenSubs vars = {
        {"${TARGET}",  target},
        {"${NUM_QUBITS}", qureg.numQubits}};

    assertThat(target >= 0 && target < qureg.numQubits, report::INVALID_TARGET_QUBIT, vars, caller);
}



/*
 * FILE IO
 */

void validate_canReadFile(string fn, const char* caller) {

    // TODO: embed filename into error message when tokenSubs is updated to permit strings
    assertThat(parser_canReadFile(fn), report::CANNOT_READ_FILE, caller);
}