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
    extern std::string PERMIT_NODES_TO_SHARE_GPU;
    extern std::string DEFAULT_VALIDATION_EPSILON;
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


#endif // ENVVARS_HPP
