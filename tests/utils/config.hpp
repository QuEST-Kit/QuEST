/** @file
 * Testing utilities for loading environment variables
 * which configure the unit tests, independent of QuEST's 
 * internal environment variable facilities
 *
 * @author Tyson Jones
 */


/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilsconfig Config
 * @ingroup testutils
 * @brief
 * Testing utilities for loading environment variables
 * which configure the unit tests, independent of QuEST's 
 * internal environment variable facilities
 * @{
 */

#ifndef CONFIG_HPP
#define CONFIG_HPP


/*
 * SPECIFYING ENV-VARS 
 */

// spoofing as macros to doc; beware that the values below
// merely duplicate but do not change the default values
// which are hardcoded in config.cpp
#if 0

    /// @envvardoc
    const int TEST_NUM_QUBITS_IN_QUREG = 6;

    /** @envvardoc
     * 
     * Specifies the maximum number of control and target qubit permutations for which to unit test each relevant 
     * API function.
     * 
     * Many QuEST functions accept a varying number of target qubits (like applyCompMatr()) and/or control qubits
     * (like applyMultiControlledCompMatr()). The unit tests will run these functions, passing every possible number
     * of target qubits (alongside every possible number of control qubits, if possible), from one (zero) up to the
     * number contained within the tested `Qureg` (minus the number of target qubits).
     * 
     * For each of these tested number-of-targets and number-of-controls combinations, there are factorially-many 
     * possible choices of the arbitrarily-ordered qubit indices, i.e. sub-permutations of all Qureg qubits. 
     * By default, the unit tests deterministically check every permutation in-turn. This can become prohibitively 
     * slow when the tested `Qureg` are large. For example, there are `604,800` unique, non-overlapping choices of 
     * `4`  targets and `3` controls in a Qureg containing `10` qubits.
     * 
     * When this environment variable is set to a non-zero value, the unit tests will forego testing every permutation
     * and instead perform only the number specified, randomising the involved qubits. This can significantly speed up the
     * tests though risks missing esoteric edge-cases. The runtime of the tests are approximately linearly proportional
     * to the specified number of permutations. When the specified non-zero value exceeds the number of unique 
     * permutations, the tests will revert to deterministically evaluating each once.
     * 
     * @envvarvalues
     * 
     * - set to `0` (default) to systematically test all permutations.
     * - set to a positive integer (e.g. `50`) to test (at most) that many random permutations and accelerate the tests.
     * 
     * @author Tyson Jones
     */
    const int TEST_MAX_NUM_QUBIT_PERMUTATIONS = 0;

    /// @envvardoc
    const int TEST_MAX_NUM_SUPEROP_TARGETS = 4;

    /// @envvardoc
    const int TEST_ALL_DEPLOYMENTS = 1;

    /// @envvardoc
    const int TEST_NUM_MIXED_DEPLOYMENT_REPETITIONS = 10;

#endif


/*
 * ACCESSING ENV-VARS 
 */

int getNumQubitsInUnitTestedQuregs();
int getMaxNumTestedQubitPermutations();
int getMaxNumTestedSuperoperatorTargets();
int getNumTestedMixedDeploymentRepetitions();
bool getWhetherToTestAllDeployments();


#endif // CONFIG_PP

/** @} (end defgroup) */
