/** @file
 * Functions for loading environment variables, useful for
 * configuring QuEST ahead of calling initQuESTEnv(), after
 * compilation.
 * 
 * @author Tyson Jones
 */

#ifndef ENVVARS_HPP
#define ENVVARS_HPP

#include <string>


namespace envvar_names { 
    extern std::string QUEST_PERMIT_NODES_TO_SHARE_GPU;
    extern std::string QUEST_DEFAULT_VALIDATION_EPSILON;
    extern std::string QUEST_DEFAULT_NUM_GPU_THREADS_PER_BLOCK;
}


/*
 * LOAD VARS
 */

void envvars_validateAndLoadEnvVars(const char* caller);


/*
 * GET VAR
 */

bool envvars_getWhetherGpuSharingIsPermitted();

qreal envvars_getDefaultValidationEpsilon();

int envvars_getDefaultNumGpuThreadsPerBlock();


#endif // ENVVARS_HPP
