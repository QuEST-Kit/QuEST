/** @file
 * Testing utilities for loading environment variables
 * which configure the unit tests, independent of QuEST's 
 * internal environment variable facilities
 *
 * @author Tyson Jones
 */

#include <string>
#include <cstdlib>
#include <stdexcept>

using std::string;


/*
 * PRIVATE
 */

string getEnvVarValue(string name) {

    // unspecified var returns empty string
    const char* ptr = std::getenv(name.c_str());
    return (ptr == nullptr)? "" : std::string(ptr);
}

int getIntEnvVarValueOrDefault(string name, int defaultValue) {

    string strValue = getEnvVarValue(name);
    int intValue = defaultValue;

    // overwrite default only when passed variable is interpretable
    try {
        intValue = std::stoi(strValue);
    } 
    catch (const std::out_of_range&) { } 
    catch (const std::invalid_argument&) { }
    return intValue;
}


/*
 * PUBLIC TEST ENV VARS
 */

int getNumQubitsInUnitTestedQuregs() {

    static int value = getIntEnvVarValueOrDefault("QUEST_TEST_NUM_QUBITS_IN_QUREG", 6);
    return value;
}

int getMaxNumTestedQubitPermutations() {

    static int value = getIntEnvVarValueOrDefault("QUEST_TEST_MAX_NUM_QUBIT_PERMUTATIONS", 0);
    return value;
}

int getMaxNumTestedSuperoperatorTargets() {

    static int value = getIntEnvVarValueOrDefault("QUEST_TEST_MAX_NUM_SUPEROP_TARGETS", 4);
    return value;
}

int getNumTestedMixedDeploymentRepetitions() {

    static int value = getIntEnvVarValueOrDefault("QUEST_TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS", 10);
    return value;
}

bool getWhetherToTestAllDeployments() {

    static bool value = getIntEnvVarValueOrDefault("QUEST_TEST_TRY_ALL_DEPLOYMENTS", 1);
    return value;
}



/*
 * PUBLIC QUEST ENV VARS
 */

int getDefaultNumGpuThreadsPerBlock() {

    // when the env-var is not present, we MUST return the default assumed by the QuEST src code,
    // which at the time of writing, is a fixed 128 (rather than hardware-specific value)
    const int compileTimeDefaultTPB = 128;

    // when the env-var is present, we consult that, just like QuEST
    static int value = getIntEnvVarValueOrDefault("QUEST_NUM_GPU_THREADS_PER_BLOCK", compileTimeDefaultTPB);
    return value;
}
