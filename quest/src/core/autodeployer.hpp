/** @file
 * Functions which automatically choose QuESTEnv and Qureg deployment
 * parameters, replacing flag modeflag::USE_AUTO with 0 or 1,
 * depending on the compiled facilities, available hardware, and
 * Qureg dimensions.
 */

#ifndef AUTODEPLOYER_HPP
#define AUTODEPLOYER_HPP

#include "quest/include/environment.h"



void autodep_chooseQuESTEnvDeployment(int &useDistrib, int &useGpuAccel, int &useMultithread);

void autodep_chooseQuregDeployment(int numQubits, int isDensMatr, int &useDistrib, int &useGpuAccel, int &useMultithread, QuESTEnv env);



#endif // AUTODEPLOYER_HPP