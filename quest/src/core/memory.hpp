/** @file
 * Internal functions which query available CPU memory (in an
 * attemptedly OS-agnostic way) and determine the maximum
 * number of qubits which can be simualted. Note GPU memory
 * querying is performed by the dedicated GPU backend, 
 * though this file is always compiled (even in GPU mode) 
 * because GPU-acceleration still requires accompanying
 * CPU memory arrays.
 */