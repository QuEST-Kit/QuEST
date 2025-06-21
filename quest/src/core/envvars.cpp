/** @file
 * Functions for loading environment variables, useful for
 * configuring QuEST ahead of calling initQuESTEnv(), after
 * compilation.
 * 
 * @author Tyson Jones
 */

#include "quest/include/precision.h"
#include "quest/include/types.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/parser.hpp"
#include "quest/src/core/validation.hpp"

#include <string>
#include <cstdlib>

using std::string;



/*
 * FIXED ENV-VAR NAMES
 */


namespace envvar_names {
    string PERMIT_NODES_TO_SHARE_GPU = "PERMIT_NODES_TO_SHARE_GPU";
    string DEFAULT_VALIDATION_EPSILON = "DEFAULT_VALIDATION_EPSILON";
}



/*
 * USER-OVERRIDABLE DEFAULT ENV-VAR VALUES
 */


namespace envvar_values {

    // by default, do not permit GPU sharing since it sabotages performance
    // and should only ever be carefully, deliberately enabled
    bool PERMIT_NODES_TO_SHARE_GPU = false;

    // by default, the initial validation epsilon (before being overriden
    // by users at runtime) should depend on qreal (i.e. FLOAT_PRECISION)
    qreal DEFAULT_VALIDATION_EPSILON = UNSPECIFIED_DEFAULT_VALIDATION_EPSILON;
}


// indicates whether envvars_validateAndLoadEnvVars() has been called
bool global_areEnvVarsLoaded = false;



/*
 * PRIVATE UTILITIES
 */


bool isEnvVarSpecified(string name) {

    // note var="" is considered unspecified, but var=" " is specified
    const char* ptr = std::getenv(name.c_str());
    return (ptr != nullptr) && (ptr[0] != '\0');
}


string getSpecifiedEnvVarValue(string name) {

    // assumes isEnvVarSpecified returned true
    // (calling getenv() a second time is fine)
    return std::string(std::getenv(name.c_str()));
}


void assertEnvVarsAreLoaded() {

    if (!global_areEnvVarsLoaded)
        error_envVarsNotYetLoaded();
}



/*
 * PRIVATE BESPOKE ENV-VAR LOADERS
 *
 * which we have opted to not-yet make generic 
 * (e.g. for each type) since YAGNI
 */


void validateAndSetWhetherGpuSharingIsPermitted(const char* caller) {

    // permit unspecified, falling back to default value
    string name = envvar_names::PERMIT_NODES_TO_SHARE_GPU;
    if (!isEnvVarSpecified(name))
        return;

    // otherwise ensure value == '0' or '1' precisely (no whitespace)
    string value = getSpecifiedEnvVarValue(name);
    validate_envVarPermitNodesToShareGpu(value, caller);

    // overwrite default env-var value
    envvar_values::PERMIT_NODES_TO_SHARE_GPU = (value[0] == '1');
}


void validateAndSetDefaultValidationEpsilon(const char* caller) {

    // permit unspecified, falling back to the hardcoded precision-specific default
    string name = envvar_names::DEFAULT_VALIDATION_EPSILON;
    if (!isEnvVarSpecified(name))
        return;
    
    // otherwise, validate user passed a positive real integer (or zero)
    string value = getSpecifiedEnvVarValue(name);
    validate_envVarDefaultValidationEpsilon(value, caller);

    // overwrite default env-var value
    envvar_values::DEFAULT_VALIDATION_EPSILON = parser_parseReal(value);    
}



/*
 * PUBLIC
 */


void envvars_validateAndLoadEnvVars(const char* caller) {

    // error if loaded twice since this indicates spaghetti
    if (global_areEnvVarsLoaded)
        error_envVarsAlreadyLoaded();

    // load all env-vars
    validateAndSetWhetherGpuSharingIsPermitted(caller);
    validateAndSetDefaultValidationEpsilon(caller);

    // ensure no re-loading
    global_areEnvVarsLoaded = true;
}


bool envvars_getWhetherGpuSharingIsPermitted() {
    assertEnvVarsAreLoaded();

    return envvar_values::PERMIT_NODES_TO_SHARE_GPU;
}


qreal envvars_getDefaultValidationEpsilon() {
    assertEnvVarsAreLoaded();

    return envvar_values::DEFAULT_VALIDATION_EPSILON;
}
