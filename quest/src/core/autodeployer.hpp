/** @file
 * Functions which automatically choose QuESTEnv, Qureg and
 * FullStateDiagMatr deployment parameters, replacing flag 
 * modeflag::USE_AUTO with 0 or 1, depending on the compiled 
 * facilities, available hardware, and object dimensions.
 * 
 * @author Tyson Jones
 */

#ifndef AUTODEPLOYER_HPP
#define AUTODEPLOYER_HPP

#include "quest/include/environment.h"



#define MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_MULTITHREADING 8

#define MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_GPU_ACCELERATION 12

#define MIN_NUM_LOCAL_QUBITS_FOR_AUTO_QUREG_DISTRIBUTION 26



void autodep_chooseQuESTEnvDeployment(int &useDistrib, int &useGpuAccel, int &useMultithread);

void autodep_chooseQuregDeployment(int numQubits, int isDensMatr, int &useDistrib, int &useGpuAccel, int &useMultithread, QuESTEnv env);

void autodep_chooseFullStateDiagMatrDeployment(int numQubits, int &useDistrib, int &useGpuAccel, int &useMultithread, QuESTEnv env);



#endif // AUTODEPLOYER_HPP