/** @file
 * API signatures for debugging QuEST behaviour, 
 * controlling input validation, changing reporter
 * parameters or seeding random generation.
 * 
 * @author Tyson Jones
 *
 * @defgroup debug Debug
 * @ingroup api
 * @brief Utilities for controlling QuEST behaviour such as seeding, input validation and printing.
 * @{
 */

#ifndef DEBUG_H
#define DEBUG_H

#include "quest/include/types.h"



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif



/** 
 * @defgroup debug_seed Seeding
 * @brief Functions for seeding QuEST's random generators.
 * @details Re-seeding with identical seeds will determine all of QuEST's subsequent 
 *          random outputs (such as measurement and random state preparation), and can
 *          be done at any stage of execution. When seeding is not explicitly performed,
 *          QuEST will attempt to use a cryptographically secure pseudorandom number generator
 *          (CSPRNG) if locally available, else fall back to a standard PRNG, via using
 *          the standard C++ `random_device` class.
 * @{
 */


/** Sets the seeds used by QuEST's random number generation.
 * 
 * This affects and determines all pseudorandom decisions made the QuEST library,
 * such as qubit measurement outcomes and random state preparation. The seeds are
 * passed without modification to QuEST's 19,937-bit
 * [Mersenne Twister](https://cplusplus.com/reference/random/mt19937_64/) generator.
 * As such, fully specifying the generator's initial state requires 19,937 seed bits, 
 * or 623 unsigned integers, although this is rarely necessary.
 * 
 * Seeding can be performed at any time before a random decision (such as measurement);
 * no random outcomes/state is stored on any object beforehand.
 *
 * > [!NOTE]
 * > In distributed simulation, only the root node's seeds are consulted (and broadcast
 * > to all nodes) while those passed on all other nodes are ignored.
 * 
 * @myexample
 * ```
    // make report() print only one trailing newline
    setQuESTNumReportedNewlines(1);

    // seed
    unsigned seeds[] = {123456789u, 987654321u, 546372819u};
    setQuESTSeeds(seeds, 3);

    // randomly collapse the plus state
    initPlusState(qureg);
    reportScalar("outcome 0", applyQubitMeasurement(qureg, 0));
    reportScalar("outcome 1", applyQubitMeasurement(qureg, 1));
    reportScalar("outcome 2", applyQubitMeasurement(qureg, 2));

    // restore RNG to state prior to collapse
    setQuESTSeeds(seeds, 3);

    // randomly collapse the plus state again, obtaining same outcomes
    reportStr("");
    initPlusState(qureg);
    reportScalar("outcome 0", applyQubitMeasurement(qureg, 0));
    reportScalar("outcome 1", applyQubitMeasurement(qureg, 1));
    reportScalar("outcome 2", applyQubitMeasurement(qureg, 2));
 * ```
 * Example output:
 * ```
    outcome 0: 0
    outcome 1: 1
    outcome 2: 1

    outcome 0: 0
    outcome 1: 1
    outcome 2: 1
 * ```
 * 
 * @param[in] seeds    a list of seeds.
 * @param[in] numSeeds the length of @p seeds.
 * @throws @validationerror
 * - if @p seeds is a null pointer.
 * - if @p numSeeds is less than one.
 * @see
 * - getQuESTSeeds()
 * - getQuESTNumSeeds()
 * - setQuESTSeedsToDefault()
 * @author Tyson Jones
 */
void setQuESTSeeds(unsigned* seeds, int numSeeds);


/** Re-randomizes QuEST, seeding its random number generator with new seeds as pseudorandomly
 * produced by the @c C++ [@c std::random_device](https://en.cppreference.com/cpp/numeric/random/random_device)
 * function, potentially using a hardware RNG.
 * 
 * Similar to setQuESTSeeds(), this affects all pseudorandom decisions made the QuEST library,
 * such as qubit measurement outcomes and random state preparation. Unlike setQuESTSeeds()
 * however, repeated calls will not guarantee identical random generation; the internally
 * processed seeds may differ at each invocation. This is, in fact, the function first
 * internally called during QuEST environment initialisation.
 * 
 * Seeding can be performed at any time after QuEST library initialisation, before a
 * random decision (such as measurement); no random outcomes/state is stored on any object beforehand.
 * 
 * This function consults @c std::random_device to produce @c DEFAULT_NUM_RNG_SEEDS=4
 * seeds, which is infact insufficient to fully specify an initial state of QuEST's
 * [19,937-bit  Mersenne Twister](https://cplusplus.com/reference/random/mt19937_64/)
 * generator. Presently, @c DEFAULT_NUM_RNG_SEEDS cannot be changed except through
 * manual modification and recompilation.
 *
 * > [!NOTE]
 * > In distributed simulation, only the root node's seeds are consulted (and broadcast
 * > to all nodes) while those passed on all other nodes are ignored.
 * 
 * Note that many data structures (e.g. CompMatr, DiagMatr, KrausMap) will assess epsilon-dependent
 * validation properties such as unitarity _once_, recording the result in a persistent heap field
 * (like KrausMap.isApproxCPTP) to avoid superfluous re-calculation.
 * Updating the global validation epsilon via this function will update all persistent heap fields,
 * marking epsilon-dependent properties as "unknown", which will be lazily re-evaluated when validation
 * is next performed. Ergo, calling this function can cause later additional function overheads.
 * 
 * @myexample
 * ```
    // make report() print only one trailing newline
    setQuESTNumReportedNewlines(1);

    // randomise the RNG state
    setQuESTSeedsToDefault();

    // randomly collapse the plus state
    initPlusState(qureg);
    reportScalar("outcome 0", applyQubitMeasurement(qureg, 0));
    reportScalar("outcome 1", applyQubitMeasurement(qureg, 1));
    reportScalar("outcome 2", applyQubitMeasurement(qureg, 2));

    // re-randomise the RNG state
    setQuESTSeedsToDefault();

    // randomly collapse the plus state again, obtaining new outcomes
    reportStr("");
    initPlusState(qureg);
    reportScalar("outcome 0", applyQubitMeasurement(qureg, 0));
    reportScalar("outcome 1", applyQubitMeasurement(qureg, 1));
    reportScalar("outcome 2", applyQubitMeasurement(qureg, 2));
 * ```
 * Example output:
 * ```
    outcome 0: 0
    outcome 1: 1
    outcome 2: 0

    outcome 0: 1
    outcome 1: 1
    outcome 2: 1
 * ```
 * We can view the random seeds chosen by QuEST using getQuESTSeeds():
 * ```
    int numSeeds = getQuESTNumSeeds(); // fixed=4
    unsigned seeds[1000]; // lazily/dangerously assume num<=1000
    
    setQuESTSeedsToDefault();
    getQuESTSeeds(seeds);

    for (int i=0; i<numSeeds; i++)
        printf("%u ", seeds[i]);
    printf("\n\n");

    setQuESTSeedsToDefault();
    getQuESTSeeds(seeds);
    
    for (int i=0; i<numSeeds; i++)
        printf("%u ", seeds[i]);
    printf("\n");
 * ```
 * Example output:
 * ```
    766738893 3449946700 2005286535 3261547242 

    3996404805 474611433 1261615115 3403708295 
 * ```
 * 
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * @see
 * - setQuESTSeeds()
 * - getQuESTSeeds()
 * - getQuESTNumSeeds()
 * @author Tyson Jones
 */
void setQuESTSeedsToDefault();


/** Populates @p seeds with those which last seeded QuEST's random number generation.
 * 
 * This will return the last seeds passed to setQuESTSeeds(), or those internally
 * chosen within setQuESTSeedsToDefault() (and equivalently, those initially chosen
 * during the QuEST environment initialisation). The number of seeds, which informs
 * the necessary capacity of @p seeds, is obtained by getQuESTNumSeeds().
 * 
 * @myexample
 * 
 * ```
    int numSeeds = getQuESTNumSeeds();
    unsigned *seeds = malloc(numSeeds * sizeof *seeds);

    getQuESTSeeds(seeds);
    printf("seeds[%d] = { ", numSeeds);

    for (int i=0; i<numSeeds; i++)
        printf("%u ", seeds[i]);
    printf("} \n");

    free(seeds);
 * ```
 * may output
 * ```
    seeds[4] = { 3674477761 2499034318 2614242128 475980445 } 
 * ```
 * 
 * See setQuESTSeedsToDefault() for an example of how getQuESTSeeds() interacts with
 * other seeding calls.
 * 
 * @param[out] seeds the list of seeds
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * - if @p seeds is a nullptr.
 * @throws seg-fault
 * - if @p seeds does not have capacity to write as many @c unsigned as getQuESTNumSeeds() returns.
 * @see
 * - getQuESTNumSeeds()
 * - setQuESTSeeds()
 * @author Tyson Jones
 */
void getQuESTSeeds(unsigned* seeds);


/** Returns the number of seeds used by QuEST in its latest round of random number generator seeding.
 * 
 * This is the length of the list output by getQuESTSeeds().
 * 
 * @returns the number of seeds.
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * @see
 * - getQuESTSeeds()
 * - setQuESTSeeds()
 * @author Tyson Jones
 */
int getQuESTNumSeeds();


/** @} */



/** 
 * @defgroup debug_validation Validation
 * @brief Functions to control QuEST's user-input validation.
 * @details These can be used to adjust the precision with which properties like unitarity 
 *          are checked/enforced, or otherwise disable all input validation (e.g. is the
 *          given qubit index valid?). Note passing erroneous input while validation is 
 *          disabled can result in runtime errors like segmentation faults. 
 * @{
 */


/** Sets the function which QuEST will call when encountering an invalid input.
 * 
 * By default, when a user passes an invalid input to QuEST (such as a negative qubit index),
 * an internal function @c default_inputErrorHandler() is called which prints a message to
 * @c stdout, attempts to gracefully clean up communication in distributed settings, then
 * exits execution with @c exit(EXIT_FAILURE). If this is undesired, setQuESTInputErrorHandler()
 * allows the user to substitute @p callback for the default handler, which will receive the
 * throwing API function name @p func, and the error message string @p msg.
 * 
 * QuEST endeavours to perform input validation upfront before proceeding to any mutation
 * of passed objects (like a Qureg). This permits gracefully catching validation errors through
 * custom handling without corrupting QuEST's state. For example, @c C++ users may wish 
 * to throw an exception within @p callback, or MPI superusers may wish to perform custom 
 * communicator cleanup before exiting.
 * 
 * > [!IMPORTANT]
 * > It is crucial that @p callback does not return execution back to the throwing
 * > QuEST function, which is likely to cause a segmentation fault or other internal
 * > error. Instead, @p callback should exit or throw an exception, caught by the
 * > user's control flow.
 * 
 * Note validation can be changed or disabled with setQuESTValidationEpsilon() and
 * setQuESTValidationOff(), which affects when @p callback will be called. This function
 * can be called at any time to update the error handler.
 * 
 * @myexample
 * 
 * ```
    void myErrorHandler(const char* errFunc, const char* errMsg) {
        printf("Ruh-roh, Raggy! Function '%s' has reported '%s'.\n", errFunc, errMsg);
        printf("We will now be very good children and exit immediately!\n");
        exit(0);
    }

    int main() {
        initQuESTEnv();
        setQuESTInputErrorHandler(myErrorHandler);
        createQureg(9999); // invokes myErrorHandler
        ...
    }
 * ```
 *
 * @param[in] callback a pointer to a function which accepts two `const char*` arguments.
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * @throws seg-fault
 * - if @p callback is a null-ptr and an invalid input is later encountered.
 * @see
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/setting_errorhandler.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/setting_errorhandler.cpp) examples
 * - setQuESTValidationOff()
 * @author Tyson Jones
 */
void setQuESTInputErrorHandler(void (*callback)(const char* func, const char* msg));


/** Restores QuEST's input validation.
 * 
 * This means that invalid inputs to other QuEST API functions will call the error handler
 * (the default, or one passed to setQuESTInputErrorHandler()), rather than be silently ignored.
 * 
 * This function only has an affect if setQuESTValidationOff() was prior called. This function
 * does not affect the validation epsilon controlled with setQuESTValidationEpsilon(), which when
 * zero, will still see the skipping of numerically-approximated validations (such as unitarity checks).
 * 
 * @see
 * - setQuESTValidationOff()
 * @author Tyson Jones
 */
void setQuESTValidationOn();


/** Disables all of QuEST's input validation.
 * 
 * When a QuEST API function encounters an invalid input, the error will be ignored and execution of the
 * function will proceed. This is useful in order to call a QuEST function in a manner which is ordinarily
 * forbidden to avoid user mistakes. 
 * 
 * > [!IMPORTANT]
 * > This function disables all of QuEST's runtime input validation, meaning invalid inputs such as
 * > negative qubit indices will be accepted and trusted, likely causing internal errors and segmentation faults.
 * 
 * Users wishing only to adjust numerical validation tolerances, or disable numerical validations such as
 * matrix unitarity checks, should instead use the safer setQuESTValidationEpsilon(). If it is essential to
 * disable all validation, it should be later restored with setQuESTValidationOn().
 * 
 * @myexample
 *
 * ```
    Qureg qureg = createDensityQureg(3); // ~ 6-qubit statevector
    CompMatr1 matr = getInlineCompMatr1({{1,2},{3,4}});

    int target = 4; // 0-2 valid, 3-5 hacky, 6+ seg-fault

    setQuESTValidationOff();
    leftapplyCompMatr1(qureg, target, matr);
    setQuESTValidationOn();
 * ```
 * 
 * @see
 * - setQuESTValidationOn()
 * - setQuESTValidationEpsilon()
 * - setQuESTInputErrorHandler()
 * @author Tyson Jones
 */
void setQuESTValidationOff();


/** Restores QuEST's validation threshold for testing numerical or approximate quantities, to its default value.
 * 
 * The default value is informed by the environment variable `QUEST_DEFAULT_VALIDATION_EPSILON`. If the
 * environment variable was not specified during the launch of the QuEST executable, then the default
 * validation epsilon is specific to the precision of `qreal`, as controlled by `QUEST_FLOAT_PRECISION`.
 * These are:
 * | @c QUEST_FLOAT_PRECISION | @c qreal  | default epsilon |
 * |--------------------------|-----------|-----------------|
 * | 1                        | @c float  | @c 1E-5         |
 * | 2                        | @c double | @c 1E-12        |
 * | 4                        | @c long @c double | @c 1E-15 |
 *
 * @see
 * - setQuESTValidationEpsilon()
 * @author Tyson Jones
 */
void setQuESTValidationEpsilonToDefault();


/** Modifies QuEST's validation threshold for testing numerical or approximate quantities, to @p eps.
 * 
 * Many of QuEST's API functions validate that expected numerical properties of the input are satisfied.
 * For example, that the matrix passed to applyCompMatr1() is unitary, and ergo that the
 * product of the matrix with its own adjoint produces the identity matrix. Due to floating-point error,
 * such properties cannot be evaluated exactly, and small disagreement between the expected and given
 * property is tolerated. This difference is the validation epsilon, as overridden by this function.
 * Precisely how the validation epsilon is used by numerical validation is function specific, and
 * individually documented.
 * 
 * > [!TIP]
 * > Passing @p eps=0 effectively encodes @p eps=infinity, and *disables* all numerical validation.
 * 
 * The validation epsilon has no effect on non-numerical validation, such as whether qubit incices
 * are valid. In general, it is therefore safe to modify and disable numerical validation via this
 * function.
 * 
 * The default validation epsilon, which can itself be controlled by the `QUEST_DEFAULT_VALIDATION_EPSILON`
 * environment variable, is restored via setQuESTValidationEpsilonToDefault().
 * 
 * Note that many data structures (e.g. CompMatr, DiagMatr, KrausMap) will assess epsilon-dependent
 * validation properties such as unitarity _once_, recording the result in a persistent heap field
 * (like KrausMap.isApproxCPTP) to avoid superfluous re-calculation.
 * Updating the global validation epsilon via this function will update all persistent heap fields,
 * marking epsilon-dependent properties as "unknown", which will be lazily re-evaluated when validation
 * is next performed. Ergo, calling this function can cause later additional function overheads.
 * 
 * @myexample
 * 
 * ```
    // | max [matr . adj(matr) - identity] |^2 = 576
    CompMatr1 matr = getInlineCompMatr1({{1,2},{3,4}}); // non-unitary
    setQuESTValidationEpsilon(576.0);
    applyCompMatr1(qureg, 0, matr); // no error

    matr.elems[0][0] = 999;
    setQuESTValidationEpsilon(0); // disable all numerical validation
    applyCompMatr1(qureg, 0, matr); // no error

    // target=-1 would still trigger an error
    // applyCompMatr1(qureg, -1, matr);
 * ```
 * 
 * @param[in] eps the new validation epsilon.
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * - if @p eps is negative.
 * @see
 * - setQuESTValidationEpsilonToDefault()
 * - getQuESTValidationEpsilon()
 * - setQuESTValidationOff()
 * @author Tyson Jones
 */
void setQuESTValidationEpsilon(qreal eps);


/** Returns the threshold used by QuEST's numerical validation.
 * 
 * This is the value last passed to setQuESTValidationEpsilon(), unless overridden by
 * setQuESTValidationEpsilonToDefault(), or similarly if never called. It indicates the
 * precision and correctness demanded of input numerical quantities to the QuEST API. A larger
 * epsilon corresponds to more permissive validation, while a smaller epsilon means validation
 * is harder to pass and input quantities must be more carefully prepared. 
 * 
 * > [!NOTE]
 * > A validation epsilon of @c 0 indicates numerical validation is disabled. 
 * 
 * The exact usage of the validation epsilon is function specific. For an example, see applyCompMatr1().
 * The validation has no effect when validation has been disabled entirely via setQuESTValidationOff().
 * 
 * @returns The validation epsilon.
 * @see
 * - setQuESTValidationEpsilon()
 * - setQuESTValidationEpsilonToDefault()
 * - setQuESTValidationOff()
 * @author Tyson Jones
 */
qreal getQuESTValidationEpsilon();


/** @} */



/** 
 * @defgroup debug_reporting Reporting
 * @brief Functions to control how QuEST's reporters display and truncate information.
 * @{
 */


/** Sets the maximum number of rows and columns of data structures subsequently printed by
 * QuEST's report functions, with remaining elements ellipted.
 * 
 * This function is useful for keeping @c stdout concise and readable when reporting
 * large data structures, such as matrices and Qureg. The ellipted elements are in
 * the _center_ of the reported data structure. The specified maximums persist until
 * changed again with this function.
 * 
 * > [!TIP]
 * > Specifying @p numRows=0 or @p numCols=0 respectively disables ellipsis across
 * > rows and columns respectively. Beware this means that functions like reportQureg()
 * > will display the entirety of their data, which can be very large.
 * 
 * This function presently affects the output of functions:
 * - reportQureg()
 * - reportPauliStrSum()
 * - reportCompMatr1()
 * - reportCompMatr2()
 * - reportCompMatr()
 * - reportDiagMatr1()
 * - reportDiagMatr2()
 * - reportDiagMatr()
 * - reportFullStateDiagMatr()
 * - reportKrausMap()
 * - reportSuperOp()
 * 
 * In contrast, it has no effect on the output of functions:
 * - reportPauliStr()
 * - reportScalar()
 * - reportStr()
 * 
 * > When this function is not called, the QuEST defaults are adopted, which are presently:
 * > - @c numRows=32
 * > - @c numCols=4
 * 
 * @myexample
 * 
 * ```cpp
    Qureg qureg = createDensityQureg(6);
    initRandomPureState(qureg);
    setQuESTMaxNumReportedItems(7,4);
    reportQureg(qureg);
 * ```
 * will report 7 rows and 4 columns of @c qureg.
 * ```text
    Qureg (6 qubit density matrix, 64x64 qcomps, 64.1 KiB):
        0.012273+(3.2329e-19)i  -0.0056906+0.0072884i   …  -0.0052104-0.0070839i    0.00092891-0.0086299i    
        -0.0056906-0.0072884i   0.0069669-(1.363e-20)i  …  -0.0017909+0.0063788i    -0.0055556+0.0034498i    
        0.014757-0.0083751i     -0.0018689+0.012647i    …  -0.011099-0.0049622i     -0.0047721-0.011011i     
        0.00078401+0.0014319i   -0.0012139-0.00019834i  …  0.00049365-0.0010604i    0.0010662-0.00044291i    
                  ⋮
        -0.010456+0.012047i     -0.0023056-0.011795i    …  0.011393+0.00092115i     0.0076794+0.0082644i     
        -0.0052104+0.0070839i   -0.0017909-0.0063788i   …  0.0063009-(3.3799e-20)i  0.0045868+0.0041999i     
        0.00092891+0.0086299i   -0.0055556-0.0034498i   …  0.0045868-0.0041999i     0.0061386-(1.4941e-19)i  
 * ```
 *
 * @param[in] numRows the max number of rows to report (all if @c =0).
 * @param[in] numCols the max number of columns to report (all if @c =0).
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * - if either @p numRows or @p numCols is negative.
 * @see
 * - [C](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/reporting_matrices.c) and 
 *   [C++](https://github.com/QuEST-Kit/QuEST/blob/main/examples/isolated/reporting_matrices.cpp) examples
 * - setQuESTMaxNumReportedSigFigs()
 * - setQuESTNumReportedNewlines()
 */
void setQuESTMaxNumReportedItems(qindex numRows, qindex numCols);


/** Sets the maximum number of significant figures in floating-point numbers printed by
 * QuEST's reporter functions.
 * 
 * This function is useful for keeping @c stdout concise and readable when reporting
 * numerical quantites, such as matrices and Qureg amplitudes. The specified number of
 * significant figures persists until changed again with this function.
 * 
 * > [!TIP]
 * > Numbers with _fewer_ non-zero significant figures will _not_ be padded with zeros,
 * > and so will print fewer digits than @p numSigFigs, keeping the output concise.
 * 
 * > [!NOTE]
 * > Complex quantities will have their real and imaginary components separately printed,
 * > each with the specified number of significant figures.
 * 
 * > [!IMPORTANT]
 * > This function does not affect the significant figures in printed memory sizes
 * > (e.g. `5.32 KiB`) which is always shown with three significant figures 
 * > (or four when in bytes, e.g. `1023 bytes`).
 * 
 * @myexample
 * 
 * ```cpp
    Qureg qureg = createQureg(3);
    initRandomPureState(qureg);

    setQuESTMaxNumReportedSigFigs(2);
    reportQureg(qureg);

    setQuESTMaxNumReportedSigFigs(10);
    reportQureg(qureg);
 * ```
 * may output
 * ```text
    Qureg (3 qubit statevector, 8 qcomps, 232 bytes):
        0.49+0.35i    |0⟩
        -0.085+0.17i  |1⟩
        0.37+0.099i   |2⟩
        -0.14+0.088i  |3⟩
        -0.29+0.078i  |4⟩
        0.19+0.019i   |5⟩
        0.02-0.5i     |6⟩
        -0.032+0.22i  |7⟩

    Qureg (3 qubit statevector, 8 qcomps, 232 bytes):
        0.4915502074+0.3471270204i    |0⟩
        -0.08473589829+0.1685819182i  |1⟩
        0.3702183582+0.0991921232i    |2⟩
        -0.142193399+0.08829704821i   |3⟩
        -0.2881909331+0.07795511511i  |4⟩
        0.1868915394+0.01946703883i   |5⟩
        0.01956017542-0.5029788111i   |6⟩
        -0.03171267079+0.220342331i   |7⟩
 * ```
 * Meanwhile,
 * ```cpp
    CompMatr1 matr = getInlineCompMatr1({{1,2},{3,4.123456789}});
    setQuESTMaxNumReportedSigFigs(3);
    reportCompMatr1(matr);
 * ```
 * will output
 * ```text
    CompMatr1 (1 qubit, 2x2 qcomps, 80 bytes):
        1  2     
        3  4.12
 * ```
 *
 * @param[in] numSigFigs the max number of significant figures to print in subsequent report functions.
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * - if @p numSigFigs is negative.
 * @see
 * - setQuESTMaxNumReportedItems()
 * - setQuESTNumReportedNewlines()
 * @author Tyson Jones
 */
void setQuESTMaxNumReportedSigFigs(int numSigFigs);


/** Sets the number of trailing newlines printed at the end of QuEST's report functions.
 * 
 * These newlines are merely a convenience so that users do not have to manually intersperse
 * newlines in @c stdout between functions like reportScalar() and reportPauliStr(),
 * which becomes a greater pain in distributed settings when avoiding duplicated output
 * across processes. The specified number of newlines persists until changed again with
 * this function.
 * 
 * @myexample
 * 
 * By default, the sequence
 * ```cpp
    reportCompMatr1(getInlineCompMatr1({{1,2},{3,4}}));
    reportPauliStr(getInlinePauliStr("XYZ", {0,2,4}));
    reportScalar("x", 5);
    reportStr("hello world!");
 * ```
 * will print with the internal default of @p numNewLines=2
 * ```text
    CompMatr1 (1 qubit, 2x2 qcomps, 80 bytes):
        1  2  
        3  4  

    ZIYIX

    x: 5

    hello world!
 * ```
 * but if called after @c setQuESTNumReportedNewlines(1), will output
 * ```text
    CompMatr1 (1 qubit, 2x2 qcomps, 80 bytes):
        1  2  
        3  4  
    ZIYIX
    x: 5
    hello world!
 * ```
 * It is possible to forego all trailing newlines, and also to write to @c stdout between
 * report functions.
 * ```cpp
    setQuESTNumReportedNewlines(0);
    reportPauliStr(getInlinePauliStr("XYZ", {0,2,4}));
    printf(" * ");
    reportPauliStr(getInlinePauliStr("ZZZ", {0,1,2}));
    printf(" = ");
    reportPauliStr(getInlinePauliStr("YZXZ", {0,1,2,4}));
    printf("\n");
 * ```
 * ```text
    ZIYIX * ZZZ = ZIXZY
 * ```
 *
 * @param[in] numNewlines the new number of trailing newlines.
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * - if @p numNewlines is negative.
 * @author Tyson Jones
 */
void setQuESTNumReportedNewlines(int numNewlines);


/** Sets the characters used by reportPauliStr() and reportPauliStrSum() to indicate the
 * @c I, @c X, @c Y and @c Z Pauli operators.
 * 
 * @myexample
 * ```cpp
   PauliStr str = getInlinePauliStr("XYZZ", {0,10,13,20});
   reportPauliStr(str);

   setQuESTReportedPauliChars(".xyz");
   reportPauliStr(str);

   setQuESTReportedPauliChars(" !!!");
   reportPauliStr(str);
 * ```
 * ```text
    ZIIIIIIZIIYIIIIIIIIIX

    z......z..y.........x

    !      !  !         !
 * ```
 * These symbols are used across all styles accepted by setQuESTReportedPauliStrStyle().
 * ```cpp
   setQuESTReportedPauliStrStyle(1); // style=1
   setQuESTReportedPauliChars(".XyS");
   reportPauliStr(str);
 * ```
 * ```text
   X0 y10 S13 S20 
 * ```
 *
 * @param[in] paulis four characters to indicate @c I, @c X, @c Y and @c Z respectively.
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * - if @p paulis is not length @c 4.
 * @throws seg-fault
 * - if @p paulis is a null pointer.
 * - if @p paulis does not contain a terminal character and is smaller than 4 bytes in size.
 * @see
 * - setQuESTReportedPauliStrStyle()
 * @author Tyson Jones
 */
void setQuESTReportedPauliChars(const char* paulis);


/** Sets the visual style of Pauli strings printed by reportPauliStr() and reportPauliStrSum().
 * 
 * - When @p style=0 (default), every Pauli operator in the string is printed, _until_ the
 *   largest-index non-identity operator. The rightmost printed operator is the least
 *   significant (i.e. of qubit index @c 0).
 * 
 *   For example:
 *   ```cpp
     setQuESTReportedPauliStrStyle(0); // default
     reportPauliStr(getInlinePauliStr("XYZZ", {0,10,13,20}));
 *   ```
 *   ```text
     ZIIIIIIZIIYIIIIIIIIIX
 *   ```
 * 
 * - When @p style=1, only non-identity Pauli operators are printed, followed by their indices.
 *   
 *   For example:
 *   ```cpp
     setQuESTReportedPauliStrStyle(1);
     reportPauliStr(getInlinePauliStr("XYZZ", {0,10,13,20}));
 *   ```
 *   ```text
     X0 Y10 Z13 Z20
 *   ```
 *
 * The symbols for @c I, @c X, @c Y and @c Z can be overridden with setQuESTReportedPauliChars().
 * 
 * @param[in] style either @c 0 or @c 1 to respectively indicate the above styles.
 * @throws @validationerror
 * - if the QuEST environment has not been initialised via initQuESTEnv().
 * - if @p style is not @c 0 or @c 1.
 * @see
 * - setQuESTReportedPauliChars()
 * @author Tyson Jones
 */
void setQuESTReportedPauliStrStyle(int style);


/** @} */



/** 
 * @defgroup debug_cache Caching
 * @brief Functions to control temporary memory used by the QuEST process.
 * @{
 */


/** Returns the current size (in bytes) of the persistent GPU cache used to accelerate
 * QuEST.
 * 
 * Some QuEST functions involve allocating non-trivial device data, which persists between
 * invocations for effiency. For example, applyCompMatr() when given a matrix of more than
 * five qubits, will allocate data proportional in size to the matrix dimension and the
 * simulated Qureg size. This allocation happens on-the-fly.
 * 
 * This function returns the size of that cache. It does not include device memory consumed
 * by persistent QuEST objects, such as Qureg and CompMatr.
 * When QuEST is not running in GPU-accelerated mode, this function always returns @c 0.
 * 
 * @returns The current GPU cache size in bytes.
 * @see
 * - clearQuESTGpuCache()
 * @author Tyson Jones
 */
qindex getQuESTGpuCacheSize();


/** Clears QuEST's GPU cache.
 * 
 * Some QuEST functions involve allocating non-trivial device data, which persists between
 * invocations for effiency. For example, applyCompMatr() when given a matrix of more than
 * five qubits, will allocate data proportional in size to the matrix dimension and the
 * simulated Qureg size. This allocation happens on-the-fly.
 * 
 * This function clears this cache, freeing memory at the cost of increased runtime when
 * calling functions like applyCompMatr(), due to the overhead of reallocation therein.
 * This function has no effect on the persistent device memory owned by QuEST objects like
 * Qureg and CompMatr, nor does it have any effect when QuEST is not using GPU acceleration.
 * 
 * @see
 * - getQuESTGpuCacheSize()
 * @author Tyson Jones
 */
void clearQuESTGpuCache();


/** @} */



/** 
 * @defgroup debug_info Info
 * @brief Functions for getting debugging information.
 * @{
 */


/** Writes a compact description of the active QuEST environment into @p str.
 *
 * The caller must provide a character buffer of length at least @c 200. The
 * resulting string is null-terminated and summarises the current deployment.
 *
 * @param[out] str  a character buffer of length @c 200 to overwrite.
 * @throws @validationerror
 * - if the QuEST environment is not initialised.
 * @notyettested
 * @author Tyson Jones
 */
void getQuESTEnvironmentString(char str[200]);


/** @} */



// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). We
 * manually add these to their respective Doxygen doc groups.
 */

#ifdef __cplusplus

#include <vector>


/// @ingroup debug_seed
/// @notyettested
/// @notyetdoced
/// @cppvectoroverload
/// @see setQuESTSeeds()
void setQuESTSeeds(std::vector<unsigned> seeds);


/// @ingroup debug_seed
/// @notyettested
/// @notyetdoced
/// @cpponly
/// @see getQuESTSeeds()
std::vector<unsigned> getQuESTSeeds();


#endif // __cplusplus



#endif // DEBUG_H

/** @} */ // (end file-wide doxygen defgroup)
