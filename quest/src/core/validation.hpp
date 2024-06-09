/** @file
 * Validation of user inputs which check all preconditions of the API
 * functions are satisfied, and otherwise throws a user-readable error.
 */

#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include "quest/include/environment.h"
#include "quest/include/qureg.h"



/*
 * ENVIRONMENT CREATION
 */

void validate_existingEnv(QuESTEnv env, const char* caller);

void validate_envNotYetInit(const char* caller);

void validate_envDeploymentMode(int isDistrib, int isGpuAccel, int isMultithread, const char* caller);

void validate_envDistributedBetweenPower2Nodes(int numNodes, const char* caller);



/*
 * QUREG CREATION
 */

void validate_quregAllocs(Qureg qureg, bool isNewQureg, const char* caller);



#endif // VALIDATION_HPP