/** @file
 * Backwards-compatible definitions of deprecated v3 functions which
 * work as expected, but issue a user-toggleable warning during 
 * compilation. The deprecated functions are necessarily instantiated
 * here as macros so that they are only resolved (and the compiler
 * warnings therein triggered) when a user actually calls them. Note
 * that uses should NOT dely upon these, but instead use them to ease
 * the process of refactoring v3-compatible QuEST code to v4. Further,
 * they do not enable every v3 function signatures to compile - some 
 * manual modification of the v3 QuEST program will be required, as 
 * attemptedly reported by this header during compilation.
 * 
 * @author Tyson Jones
 * 
 * (excluded from doxygen doc)
 */

#ifndef DEPRECATED_H
#define DEPRECATED_H

#include "quest.h"

#include "stdlib.h"



/*
 * INITIAL WARNING
 */

#if !defined(DISABLE_DEPRECATION_WARNINGS) || DISABLE_DEPRECATION_WARNINGS == 0

    // #warning command is always recognised (deprecated API is not MSVC-compatible)
    #warning "\
Deprecated functions have been included in compilation. The QuEST v3 API will be attemptedly \
automatically substituted for the v4 API, although some uses of the old API will still fail to \
compile. For example, access to v3 struct fields (e.g. 'ComplexMatrixN.real[0][0]') must be \
manually replaced with v4 struct fields (e.g. 'CompMatr.cpuElems[0][0]'), though new v4 functions \
make struct field access mostly redundant (e.g. via 'setCompMatr()'). Even when successfully \
compiling, use of the deprecated v3 functions is dangerous; these are not unit-tested, and the \
auto-porting to v4 may introduce new bugs. As such, you should only use this facility to help \
refactor your code to v4, and should absolutely not continue to use the old v3 API for simulation."

#endif



/*
 * TOGGLEABLE WARNING MESSAGES
 *
 * users can define precompiler constant DISABLE_DEPRECATION_WARNINGS=1
 * in order to disable compile-time deprecation warnings. This will
 * make most of the QuEST v3 API silently work by casting to the 
 * v4 API at compile-time. Note that _Pragma() are resolved at
 * compile-time, AFTER pre-processing; some compilers (like gcc) will
 * ergo treat them like statements which can fill an 'if() (token)', so
 * that (e.g.) 'if(..) oldFunc()' will issue a warning but execute
 * 'oldFunc()' outside of the conditional. That won't compile if the
 * function depends on control-flow params like a for-loop index. Yuck!
 */

#define _EFFECT_PRAGMA(cmd) _Pragma(#cmd)

#if DISABLE_DEPRECATION_WARNINGS

    #define _WARN_TYPE_RENAMED(oldname, newname)

    #define _WARN_FUNC_RENAMED(oldname, newname)

    #define _WARN_FUNC_NOW_HAS_FEWER_ARGS(oldsig, newsig)

    #define _WARN_GENERAL_MSG(msg)

#else

    #define _WARN_TYPE_RENAMED(oldname, newname) \
        _EFFECT_PRAGMA(message( \
            "The QuEST type '" oldname "' is deprecated. " \
            "Please instead use '" newname "' which has been automatically invoked."))
    
    #define _WARN_FUNC_RENAMED(oldname, newname) \
        _EFFECT_PRAGMA(message( \
            "The QuEST function '" oldname "' is deprecated. " \
            "Please instead use '" newname "' which has been automatically invoked."))

    #define _WARN_FUNC_NOW_HAS_FEWER_ARGS(oldsig, newsig) \
        _EFFECT_PRAGMA(message( \
            "The QuEST function '" oldsig "' now accepts fewer arguments, and has " \
            "new signature '" newsig "' which has been automatically invoked."))

    #define _WARN_GENERAL_MSG(msg) \
        _EFFECT_PRAGMA(message(msg))

#endif



/*
 * NON-TOGGLEABLE WARNING
 *
 * which cannot be suppressed because it sometimes precedes an error
 * which is obfuscated without prior warning
 */

#define _WARN_UNSUPPRESSABLE_MSG(msg) \
    _EFFECT_PRAGMA(message(msg))



/*
 * NON-TOGGLEABLE ERRORS
 *
 * which cannot be disabled, and after which compilation cannot continue.
 */


/// @todo: 
/// we use placeholder message(...) below which only issues a warning,
/// and does not kill compilation; there alas seems to be no cross-platform
/// method of aborting compilation with a _Pragma (and aborting during
/// preprocessing via #error is too soon). We currently work around this by
/// referring to an undefined symbol. Ew!


#define _FORCE_COMPILATION_TO_FAIL() \
    _FORCING_UNKNOWN_SYMBOL_ERROR_TO_STOP_COMPILATION


#define _ERROR_GENERAL_MSG(msg) \
    _EFFECT_PRAGMA(message(msg)) \
    _FORCE_COMPILATION_TO_FAIL()


#define _ERROR_FUNC_RENAMED(oldname, newfunc) \
    _ERROR_GENERAL_MSG( \
        "The QuEST function '" oldname "' is deprecated. " \
        "Please instead use '" newfunc "' which could not here be automatically invoked.") \
    _FORCE_COMPILATION_TO_FAIL()
    

#define _ERROR_FUNC_REMOVED(oldname) \
    _ERROR_GENERAL_MSG( \
        "The QuEST function '" oldname "' is deprecated, and has no replacement.") \
    _FORCE_COMPILATION_TO_FAIL()


#define _ERROR_PREPROCESSOR_RENAMED(oldname, newname) \
    _ERROR_GENERAL_MSG( \
        "The QuEST preprocessor '" oldname "' is deprecated. " \
        "Please instead use '" newname "' which could not here be automatically invoked.") \
    _FORCE_COMPILATION_TO_FAIL()


#define _ERROR_PREPROCESSOR_REMOVED(oldname) \
    _ERROR_GENERAL_MSG( \
        "The QuEST preprocessor '" oldname "' is deprecated, and has no replacement.") \
    _FORCE_COMPILATION_TO_FAIL()


#define _ERROR_PHASE_FUNC_REMOVED(oldname) \
    _ERROR_GENERAL_MSG( \
        "The QuEST function '" oldname "' is deprecated. Please instead create a 'DiagMatr' or 'FullStateDiagMatr', initialise it " \
        "via functions like 'setDiagMatrFromMultiVarFunc()' or 'setDiagMatrFromMultiDimLists()', and apply it via 'applyDiagMatr() or " \
        "'applyFullStateDiagMatr()'. This procedure cannot be automatically performed here." ) \
    _FORCE_COMPILATION_TO_FAIL()



/*
 * PREPROCESSOR CONSTANTS
 *
 * We cannot simply replace format specifiers like 'REAL_STRING_FORMAT' because our
 * warning code interrupts the string joining done by the preprocessor. For example,
 * '"x = " REAL_STRING_FORMAT' would not compile due to the warning inside REAL_*.
 */

#define REAL_STRING_FORMAT \
    _ERROR_PREPROCESSOR_RENAMED("REAL_STRING_FORMAT",  "QREAL_FORMAT_SPECIFIER")

#define REAL_QASM_FORMAT \
    _ERROR_PREPROCESSOR_REMOVED("REAL_QASM_FORMAT")

#define MPI_QuEST_REAL \
    _ERROR_PREPROCESSOR_REMOVED("MPI_QuEST_REAL")

#define MPI_MAX_AMPS_IN_MSG \
    _ERROR_PREPROCESSOR_REMOVED("MPI_MAX_AMPS_IN_MSG")

#define REAL_EPS \
    _ERROR_PREPROCESSOR_REMOVED("REAL_EPS")

#define REAL_SPECIFIER \
    _ERROR_PREPROCESSOR_REMOVED("REAL_SPECIFIER")

#define absReal(...) \
    _ERROR_PREPROCESSOR_REMOVED("absReal()")



/*
 * STACK-MEMORY STRUCT TYPES
 *
 * which involves us defining entirely the old v3 structs, and
 * warning users to migrate to the new v4 versions. We avoided
 * just replacing the v3 structs with the v4 equivalents 
 * because then direct access of v3 fields (.real[0][0]) would
 * cause a compiler error with no error messages. We accept
 * that only for heap-memory structs which always have a
 * proceeding create() function which will dispatch the warning.
 */


typedef struct ComplexMatrix2
{
    qreal real[2][2];
    qreal imag[2][2];
} ComplexMatrix2;

typedef struct ComplexMatrix4
{
    qreal real[4][4];
    qreal imag[4][4];
} ComplexMatrix4;


// enable referencess to ComplexMatrix2/ 4 inside this header without
// causing a warning unrelated to the user's code (must define these
// before below macro which would then seek to replace the values)

typedef ComplexMatrix2 _NoWarnComplexMatrix2;
typedef ComplexMatrix4 _NoWarnComplexMatrix4;


// warn about use of ComplexMatrix1/2, but allow it

#define ComplexMatrix2 \
    _WARN_GENERAL_MSG( \
        "The QuEST type 'ComplexMatrix2' is deprecated in favour of 'CompMatr1' which has a different memory layout. We will attempt to " \
        "automatically replace your 'ComplexMatrix2' with a 'CompMatr1' instance in subsequently invoked functions.") \
    ComplexMatrix2

#define ComplexMatrix4 \
    _WARN_GENERAL_MSG( \
        "The QuEST type 'ComplexMatrix4' is deprecated in favour of 'CompMatr2' which has a different memory layout. We will attempt to " \
        "automatically replace your 'ComplexMatrix4' with a 'CompMatr2' instance in subsequently invoked functions.") \
    ComplexMatrix4


// automatically convert ComplexMatrix2/4 to CompMatr1/2 when
// passed to a v3 function which has a v4 port

#define _GET_COMP_MATR_1_FROM_COMPLEX_MATRIX_2(m) \
    (CompMatr1) { \
        .numQubits = 1, \
        .numRows = 2, \
        .elems = { \
            {getQcomp(m.real[0][0], m.imag[0][0]), getQcomp(m.real[0][1], m.imag[0][1])}, \
            {getQcomp(m.real[1][0], m.imag[1][0]), getQcomp(m.real[1][1], m.imag[1][1])}}}

#define _GET_COMP_MATR_2_FROM_COMPLEX_MATRIX_4(m) \
    (CompMatr2) { \
        .numQubits = 2, \
        .numRows = 4, \
        .elems = { \
            {getQcomp(m.real[0][0], m.imag[0][0]), getQcomp(m.real[0][1], m.imag[0][1]), getQcomp(m.real[0][2], m.imag[0][2]), getQcomp(m.real[0][3], m.imag[0][3])}, \
            {getQcomp(m.real[1][0], m.imag[1][0]), getQcomp(m.real[1][1], m.imag[1][1]), getQcomp(m.real[1][2], m.imag[1][2]), getQcomp(m.real[1][3], m.imag[1][3])}, \
            {getQcomp(m.real[2][0], m.imag[2][0]), getQcomp(m.real[2][1], m.imag[2][1]), getQcomp(m.real[2][2], m.imag[2][2]), getQcomp(m.real[2][3], m.imag[2][3])}, \
            {getQcomp(m.real[3][0], m.imag[3][0]), getQcomp(m.real[3][1], m.imag[3][1]), getQcomp(m.real[3][2], m.imag[3][2]), getQcomp(m.real[3][3], m.imag[3][3])}}}



/*
 * PAULI OPERATOR TYPES
 */

enum pauliOpType {PAULI_I=0, PAULI_X=1, PAULI_Y=2, PAULI_Z=3};

#define _WARN_DEPRECATED_PAULI_ENUM(enum, intcode, char1, char2, char3) \
    _WARN_GENERAL_MSG( \
        "The QuEST enum '" #enum "' is deprecated, although has been defined in the deprecation header for convenience. " \
        "Please instead use integer " #intcode " or characters '" #char1 "', '" #char2 "', '" #char3 "'.") \

typedef enum pauliOpType _NoWarnPauliOpType;

#define pauliOpType \
    _WARN_GENERAL_MSG("The QuEST type 'enum pauliOpType' is deprecated, although it is still defined in the deprecation header for convenience.") \
    pauliOpType

#define PAULI_I \
    _WARN_DEPRECATED_PAULI_ENUM(PAULI_I, 0, I, i, 0) \
    PAULI_I

#define PAULI_X \
    _WARN_DEPRECATED_PAULI_ENUM(PAULI_X, 1, X, x, 1) \
    PAULI_X

#define PAULI_Y \
    _WARN_DEPRECATED_PAULI_ENUM(PAULI_Y, 2, Y, y, 2) \
    PAULI_Y

#define PAULI_Z \
    _WARN_DEPRECATED_PAULI_ENUM(PAULI_Z, 3, Z, z, 3) \
    PAULI_Z



/*
 * REMAINING TYPES
 */


#define _GET_QCOMP_FROM_COMPLEX_STRUCT(x) \
    getQcomp((x).real, (x).imag)

#define Complex \
    _WARN_GENERAL_MSG( \
        "The QuEST type 'Complex' is deprecated, and replaced with complex scalar primitive 'qcomp' which can be instantiated " \
        "with literals and modified with overloaded arithmetic operators. References to 'Complex' will be automatically replaced " \
        "with an anonymous inline struct, though the QuEST will only ever return 'qcomp' instances.") \
    struct { qreal real; qreal imag; }


#define ComplexMatrixN \
    _WARN_TYPE_RENAMED("ComplexMatrixN", "CompMatr") \
    CompMatr


#define PauliHamil \
    _WARN_TYPE_RENAMED("PauliHamil", "PauliStrSum") \
    PauliStrSum


#define DiagonalOp \
    _WARN_TYPE_RENAMED("DiagonalOp", "FullStateDiagMatr") \
    FullStateDiagMatr


#define SubDiagonalOp \
    _WARN_TYPE_RENAMED("SubDiagonalOp", "DiagMatr") \
    DiagMatr


#define Vector \
    _WARN_GENERAL_MSG("The QuEST type 'Vector' is deprecated, and will be automatically replaced with an inline struct.") \
    struct { qreal x; qreal y; qreal z; }


#define phaseFunc \
    _ERROR_GENERAL_MSG( \
        "The QuEST type 'enum phaseFunc' is deprecated. The functions which replace the 'applyPhaseFunc' family, such as " \
        "'setDiagMatrFromMultiVarFunc()' and 'setDiagMatrFromMultiDimArray()', do not use pre-defined enums." )
    
#define bitEncoding \
    _ERROR_GENERAL_MSG( \
        "The QuEST type 'enum bitEncoding' is deprecated. The new v4 function 'setDiagMatrFromMultiVarFunc()' instead accepts " \
        "an int flag to indicate whether (1) or not (0) to interpret the basis state bits under two's complement singed encoding." )



/*
 * REMOVED FUNCTIONS WITH NO REPLACEMENT
 */


#define toComplex(...) \
    _ERROR_FUNC_REMOVED("toComplex()")

#define fromComplex(...) \
    _ERROR_FUNC_REMOVED("fromComplex()")


#define applyMultiControlledMatrixN(...) \
    _ERROR_FUNC_REMOVED("applyMultiControlledMatrixN()") // our new multiplyCompMatr doesn't accept controls


#define syncQuESTSuccess(...) \
    _ERROR_FUNC_REMOVED("syncQuESTSuccess()")


#define startRecordingQASM(...) \
    _ERROR_FUNC_REMOVED("startRecordingQASM()")

#define stopRecordingQASM(...) \
    _ERROR_FUNC_REMOVED("stopRecordingQASM()")

#define clearRecordedQASM(...) \
    _ERROR_FUNC_REMOVED("clearRecordedQASM()")

#define printRecordedQASM(...) \
    _ERROR_FUNC_REMOVED("printRecordedQASM()")

#define writeRecordedQASMToFile(...) \
    _ERROR_FUNC_REMOVED("writeRecordedQASMToFile()")


#define bindArraysToStackComplexMatrixN(...) \
    _ERROR_FUNC_REMOVED("bindArraysToStackComplexMatrixN()")

#define getStaticComplexMatrixN(...) \
    _ERROR_FUNC_REMOVED("getStaticComplexMatrixN()")


#define reportState(qureg) \
    _ERROR_FUNC_REMOVED("reportState(qureg)", "reportQuregToFile(qureg, char* fn)")



/*
 * FUNCTIONS WITH REPLACEMENTS WHICH CANNOT BE AUTOMATICALLY CALLED
 */


#define initComplexMatrixN(...) \
    _ERROR_FUNC_RENAMED("initComplexMatrixN(ComplexMatrixN, qreal[][], qreal[][])", "setCompMatr(CompMatr, qcomp[][])")


#define createPauliHamil(...) \
    _ERROR_FUNC_RENAMED( \
        "createPauliHamil(int numQubits, int numTerms)", \
        "createPauliStrSum(PauliStr* strings, qcomp* coeffs, qindex numTerms)")

#define initPauliHamil(...) \
    _ERROR_FUNC_RENAMED( \
        "initPauliHamil(PauliHamil hamil, qreal* coeffs, enum pauliOpType* codes)", \
        "createPauliStrSum(PauliStr* strings, qcomp* coeffs, qindex numTerms)")


#define initDiagonalOp(...) \
    _ERROR_FUNC_RENAMED( \
        "initDiagonalOp(DiagonalOp obj, qreal* allElemsRe, qreal* allElemsIm)", \
        "setFullStateDiagMatr(FullStateDiagMatr obj, qindex startInd, qcomp* someElems, qindex numElems)")

#define setDiagonalOpElems(...) \
    _ERROR_FUNC_RENAMED( \
        "setDiagonalOpElems(DiagonalOp, qindex, qreal*, qreal*, qindex)", \
        "setFullStateDiagMatr(FullStateDiagMatr, qindex, qcomp*, qindex)")


#define initStateFromAmps(...) \
    _ERROR_FUNC_RENAMED("initStateFromAmps(Qureg, qreal*, qreal*)", "initArbitraryPureState(Qureg, qcomp*)")

#define setAmps(...) \
    _ERROR_FUNC_RENAMED("setAmps(Qureg, qindex, qreal*, qreal*, qindex)", "setQuregAmps(Qureg, qindex, qcomp*, qindex)")

#define setDensityAmps(...) \
    _ERROR_FUNC_RENAMED( \
        "setDensityAmps(Qureg, qindex startRow, qindex startCol, qreal* flatAmpsRe, qreal* flatAmpsIm, qindex numFlatAmps)", \
        "setDensityQuregAmps(Qureg, qindex startRow, qindex startCol, qcomp** amps, qindex numRows, qindex numCols)")


#define getQuESTSeeds(...) \
    _ERROR_GENERAL_MSG( \
        "The QuEST function 'getQuESTSeeds(QuESTEnv env, unsigned long int* out, int numOut)' has been deprecated. " \
        "Please instead use 'getSeeds(unsigned* out)' which accepts a pointer to pre-allocated memory of length " \
        "equal to that returned by 'getNumSeeds()'. We cannot automatically invoke this replacement routine." )


#define applyPhaseFunc(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyPhaseFunc")

#define applyPhaseFuncOverrides(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyPhaseFuncOverrides")

#define applyMultiVarPhaseFunc(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyMultiVarPhaseFunc")

#define applyMultiVarPhaseFuncOverrides(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyMultiVarPhaseFuncOverrides")

#define applyNamedPhaseFunc(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyNamedPhaseFunc")

#define applyNamedPhaseFuncOverrides(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyNamedPhaseFuncOverrides")

#define applyParamNamedPhaseFunc(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyParamNamedPhaseFunc")

#define applyParamNamedPhaseFuncOverrides(...) \
    _ERROR_PHASE_FUNC_REMOVED("applyParamNamedPhaseFuncOverrides")



/*
 * FUNCTIONS WITH THE SAME NAME BUT 1 INSTEAD OF 2 ARGS
 *
 * can be redirected to the new single-arg form, issuing a warning
 * when given the second, superfluous parameter. Note that deprecated 
 * functions which can accept C99 inline temporary arrays, e.g. 
 * f((qcomp[]) {1,2,3}), can not be overloaded in this way, because 
 * the array elements confuse variadic macros
 */


#define _GET_MACRO_WITH_1_OR_2_ARGS(_1, _2, macroname, ...) macroname

#define _CALL_MACRO_WITH_1_OR_2_ARGS(prefix, ...) \
    _GET_MACRO_WITH_1_OR_2_ARGS(__VA_ARGS__, prefix##_2, prefix##_1)(__VA_ARGS__)


#define _CREATE_QUREG_1(n) \
    createQureg(n)

#define _CREATE_QUREG_2(n, env) \
    _WARN_FUNC_NOW_HAS_FEWER_ARGS("createQureg(int, QuESTEnv)", "createQureg(int)") \
    _CREATE_QUREG_1(n)

#define createQureg(...) \
    _CALL_MACRO_WITH_1_OR_2_ARGS(_CREATE_QUREG, __VA_ARGS__)



#define _CREATE_DENSITY_QUREG_1(n) \
    createDensityQureg(n)

#define _CREATE_DENSITY_QUREG_2(n, env) \
    _WARN_FUNC_NOW_HAS_FEWER_ARGS("createDensityQureg(int, QuESTEnv)", "createDensityQureg(int)") \
    _CREATE_DENSITY_QUREG_1(n)

#define createDensityQureg(...) \
    _CALL_MACRO_WITH_1_OR_2_ARGS(_CREATE_DENSITY_QUREG, __VA_ARGS__)



#define _CREATE_CLONE_QUREG_1(n) \
    createCloneQureg(n)

#define _CREATE_CLONE_QUREG_2(n, env) \
    _WARN_FUNC_NOW_HAS_FEWER_ARGS("createCloneQureg(Qureg, QuESTEnv)", "createCloneQureg(Qureg)") \
    _CREATE_CLONE_QUREG_1(n)

#define createCloneQureg(...) \
    _CALL_MACRO_WITH_1_OR_2_ARGS(_CREATE_CLONE_QUREG, __VA_ARGS__)



#define _DESTROY_QUREG_1(qureg) \
    destroyQureg(qureg)

#define _DESTROY_QUREG_2(n, env) \
    _WARN_FUNC_NOW_HAS_FEWER_ARGS("destroyQureg(Qureg, QuESTEnv)", "destroyQureg(Qureg)") \
    _DESTROY_QUREG_1(n)

#define destroyQureg(...) \
    _CALL_MACRO_WITH_1_OR_2_ARGS(_DESTROY_QUREG, __VA_ARGS__)



#define _GET_ENVIRONMENT_STRING_1(str) \
    getEnvironmentString(str)

#define _GET_ENVIRONMENT_STRING_2(str) \
    _WARN_FUNC_NOW_HAS_FEWER_ARGS("getEnvironmentString(QuESTEnv, char[200])", "getEnvironmentString(char[200])") \
    _GET_ENVIRONMENT_STRING_1(str)

#define getEnvironmentString(...) \
    _CALL_MACRO_WITH_1_OR_2_ARGS(_GET_ENVIRONMENT_STRING, __VA_ARGS__)



/*
 * FUNCTIONS WITH THE SAME NAME BUT 0 INSTEAD OF 1 ARGS
 *
 * which are handled similar to above, but require more
 * pre-processor trickery to handle the no-arg case
 */


#define _COUNT_ZERO_OR_ONE_ARGS(...) _COUNT_ZERO_OR_ONE_ARGS_INNER(__VA_ARGS__, 1, 0,)
#define _COUNT_ZERO_OR_ONE_ARGS_INNER(_0, _1, X, ...) X
#define _CONCAT_SYMBOLS(A, B) _CONCAT_SYMBOLS_INNER(A, B)
#define _CONCAT_SYMBOLS_INNER(A, B) A ## B
#define _COMMA_SEPERATOR ,



#define _SYNC_QUEST_ENV_0(_) \
    syncQuESTEnv()

#define _SYNC_QUEST_ENV_1(env) \
    _WARN_FUNC_NOW_HAS_FEWER_ARGS("syncQuESTEnv(QuESTEnv)", "syncQuESTEnv()") \
    syncQuESTEnv()

#define syncQuESTEnv(X) \
    _CONCAT_SYMBOLS(_SYNC_QUEST_ENV_, _COUNT_ZERO_OR_ONE_ARGS(_COMMA_SEPERATOR ## X))(X)    



#define _REPORT_QUEST_ENV_0(_) \
    reportQuESTEnv()

#define _REPORT_QUEST_ENV_1(env) \
    _WARN_FUNC_NOW_HAS_FEWER_ARGS("reportQuESTEnv(QuESTEnv)", "reportQuESTEnv()") \
    reportQuESTEnv()

#define reportQuESTEnv(X) \
    _CONCAT_SYMBOLS(_REPORT_QUEST_ENV_, _COUNT_ZERO_OR_ONE_ARGS(_COMMA_SEPERATOR ## X))(X)   



/*
 * KRAUS MAPS
 *
 * which we handle explicitly, because the name 'mixKrausMap' is retained
 * (in v3, it is a 1-qubit Kraus map, but in v4, it has any-qubit) with
 * the same number of arguments, albeit with different types. This requires
 * us to use a C11 _Generic which must dispatch to compile-time (NOT
 * preprocessor) functions. We re-use an inner macro for all Kraus map funcs,
 * which accepts the variable number of target qubits at the end.
 */


#define _MIX_KRAUS_MAP_INNER(qureg, ops, numOps, targs, numTargs) \
    int dim = 1 << numTargs; \
    qcomp*** ptrs = (qcomp***) malloc(numOps * sizeof *ptrs); \
    for (int n=0; n<numOps; n++) { \
        ptrs[n] = (qcomp**) malloc(dim * sizeof **ptrs); \
        for (int r=0; r<dim; r++) { \
            ptrs[n][r] = (qcomp*) malloc(dim * sizeof ***ptrs); \
            for (int c=0; c<dim; c++) \
                ptrs[n][r][c] = getQcomp(ops[n].real[r][c], ops[n].imag[r][c]); \
        } \
    } \
    \
    KrausMap map = createKrausMap(numTargs, numOps); \
    setKrausMap(map, ptrs); \
    (mixKrausMap)(qureg, (targs), numTargs, map); /* calls below macro, wrapped to avoid warning */ \
    destroyKrausMap(map); \
    \
    for (int n=0; n<numOps; n++) { \
        for (int r=0; r<dim; r++) \
            free(ptrs[n][r]); \
        free(ptrs[n]); \
    } \
    free(ptrs);

static inline void v3_mixKrausMap(Qureg qureg, int targ, _NoWarnComplexMatrix2 *ops, int numOps) {
    _MIX_KRAUS_MAP_INNER(qureg, ops, numOps, &targ, 1);
}

#define mixKrausMap(...) \
    _WARN_UNSUPPRESSABLE_MSG( \
        "The QuEST function 'mixKrausMap()' has changed from v3 to v4. The old signature was " \
        "'mixKrausMap(Qureg, int targ, ComplexMatrix2* ops, int numOps)' while the new signature is " \
        "'mixKrausMap(Qureg, int* targs, int numTargs, KrausMap)'. Because the function name and its " \
        "number of arguments have not changed, we are unable to automatically map a v3 invocation to " \
        "the v4 API, nor can we detect when the v4 API has been used in order to avoid this message. " \
        "Hence, this warning will either be obsolete (because you have already migrated to the v4 API), " \
        "or precedes a compilation failure below due to continued use of the v3 API. To continue using " \
        "the v3 API, replace 'mixKrausMap' with 'v3_mixKrausMap' which will issue no warning nor error. " \
        "This warning cannot be suppressed.") \
    mixKrausMap(__VA_ARGS__)
    




static inline void _mixNonTPKrausMap(Qureg qureg, int targ, _NoWarnComplexMatrix2 *ops, int numOps) {
    qreal eps = getValidationEpsilon();
    setValidationEpsilon(0);
    _MIX_KRAUS_MAP_INNER(qureg, ops, numOps, &targ, 1);
    setValidationEpsilon(eps);
}

#define mixNonTPKrausMap(...) \
    _WARN_FUNC_RENAMED( \
        "mixNonTPKrausMap(Qureg, int targ, ComplexMatrix2* ops, int numOps)", \
        "mixUnnormalizedKrausMap(Qureg, int* targs, int numTargs, KrausMap)") \
    _mixNonTPKrausMap(__VA_ARGS__)



static inline void _mixTwoQubitKrausMap(Qureg qureg, int targ1, int targ2, _NoWarnComplexMatrix4 *ops, int numOps, int isNonCPTP) {
    int targs[] = {targ1, targ2};
    qreal eps = getValidationEpsilon();
    if (isNonCPTP) setValidationEpsilon(0);
    _MIX_KRAUS_MAP_INNER(qureg, ops, numOps, targs, 2);
    setValidationEpsilon(eps);
}

#define mixTwoQubitKrausMap(...) \
    _WARN_FUNC_RENAMED( \
        "mixTwoQubitKrausMap(Qureg, int targ1, int targ2, ComplexMatrix4* ops, int numOps)", \
        "mixKrausMap(Qureg, int* targs, int numTargs, KrausMap)") \
    _mixTwoQubitKrausMap(__VA_ARGS__, 0)

#define mixNonTPTwoQubitKrausMap(...) \
    _WARN_FUNC_RENAMED( \
        "mixNonTPTwoQubitKrausMap(Qureg, int targ1, int targ2, ComplexMatrix4* ops, int numOps)", \
        "mixUnnormalizedKrausMap(Qureg, int* targs, int numTargs, KrausMap)") \
    _mixTwoQubitKrausMap(__VA_ARGS__, 1)



static inline void _mixMultiQubitKrausMap(Qureg qureg, int* targs, int numTargs, CompMatr *ops, int numOps, int isNonCPTP) {

    qcomp*** ptrs = (qcomp***) malloc(numOps * sizeof *ptrs);
    for (int n=0; n<numOps; n++)
        ptrs[n] = ops[n].cpuElems;

    KrausMap map = createKrausMap(numTargs, numOps);
    setKrausMap(map, ptrs);
    free(ptrs);

    qreal eps = getValidationEpsilon();
    if (isNonCPTP) setValidationEpsilon(0);
    (mixKrausMap)(qureg, targs, numTargs, map); // calls above macro, wrapped to avoid warning */
    destroyKrausMap(map);
    setValidationEpsilon(eps);
}

#define mixMultiQubitKrausMap(...) \
    _WARN_FUNC_RENAMED( \
        "mixMultiQubitKrausMap(Qureg, int* targs, int numTargs, ComplexMatrixN* ops, int numOps)", \
        "mixKrausMap(Qureg, int* targs, int numTargs, KrausMap)") \
    _mixMultiQubitKrausMap(__VA_ARGS__, 0)

#define mixNonTPMultiQubitKrausMap(...) \
    _WARN_FUNC_RENAMED( \
        "mixNonTPMultiQubitKrausMap(Qureg, int* targs, int numTargs, ComplexMatrixN* ops, int numOps)", \
        "mixUnnormalisedKrausMap(Qureg, int* targs, int numTargs, KrausMap)") \
    _mixMultiQubitKrausMap(__VA_ARGS__, 1)



/*
 * FUNCTIONS WITH CHANGED NAMES (AND POSSIBLY NUM ARGS)
 *
 * We use variadic args wherever possible so that users passing 
 * the wrong number of args will trigger a function-related
 * compiler error, rather than the less-readable macro one. This
 * also enables users to pass C99 inline compound literal arrays,
 * which otherwise confuse the macros.
 */


static inline QuESTEnv _createQuESTEnv() {
    initQuESTEnv();
    return getQuESTEnv();
}

#define createQuESTEnv() \
    _WARN_GENERAL_MSG( \
        "The QuEST function 'createQuESTEnv()' is deprecated in favour of " \
        "'initQuESTEnv()' which returns nothing, and 'getQuESTEnv()' which " \
        "returns an immutable copy of the internally-managed environment. " \
        "This second function is for convenience; the API no longer needs " \
        "to receive the environment instance as a parameter. The above " \
        "mentioned functions have been automatically invoked.") \
    _createQuESTEnv()

#define destroyQuESTEnv(...) \
    _WARN_FUNC_RENAMED("destroyQuESTEnv(QuESTEnv)", "finalizeQuESTEnv()") \
    finalizeQuESTEnv()



#define createComplexMatrixN(...) \
    _WARN_FUNC_RENAMED("createComplexMatrixN()", "createCompMatr()") \
    createCompMatr(__VA_ARGS__)

#define destroyComplexMatrixN(...) \
    _WARN_FUNC_RENAMED("destroyComplexMatrixN()", "destroyCompMatr()") \
    destroyCompMatr(__VA_ARGS__)



#define destroyPauliHamil(...) \
    _WARN_FUNC_RENAMED("destroyPauliHamil()", "destroyPauliStrSum()") \
    destroyPauliStrSum(__VA_ARGS__)

#define createPauliHamilFromFile(...) \
    _WARN_FUNC_RENAMED("createPauliHamilFromFile()", "createPauliStrSumFromReversedFile()") \
    createPauliStrSumFromReversedFile(__VA_ARGS__)

#define reportPauliHamil(...) \
    _WARN_FUNC_RENAMED("reportPauliHamil()", "reportPauliStrSum()") \
    reportPauliStrSum(__VA_ARGS__)



#define createDiagonalOp(numQb, env) \
    _WARN_FUNC_RENAMED("createDiagonalOp(int, QuESTEnv)", "createFullStateDiagMatr(int)") \
    createFullStateDiagMatr(numQb)

#define destroyDiagonalOp(diagOp, env) \
    _WARN_FUNC_RENAMED("destroyDiagonalOp(DiagonalOp, QuESTEnv)", "destroyFullStateDiagMatr(FullStateDiagMatr)") \
    destroyFullStateDiagMatr(diagOp)

#define syncDiagonalOp(...) \
    _WARN_FUNC_RENAMED("syncDiagonalOp()", "syncFullStateDiagMatr()") \
    syncFullStateDiagMatr(__VA_ARGS__)

#define initDiagonalOpFromPauliHamil(...) \
    _WARN_FUNC_RENAMED("initDiagonalOpFromPauliHamil()", "setFullStateDiagMatrFromPauliStrSum()") \
    setFullStateDiagMatrFromPauliStrSum(__VA_ARGS__)

#define createDiagonalOpFromPauliHamilFile(fn, env) \
    _WARN_FUNC_RENAMED("createDiagonalOpFromPauliHamilFile(char*, QuESTEnv)", "createFullStateDiagMatrFromPauliStrSumFile(char*)") \
    createFullStateDiagMatrFromPauliStrSumFile(fn)

#define applyDiagonalOp(...) \
    _WARN_FUNC_RENAMED("applyDiagonalOp()", "multiplyFullStateDiagMatr()") \
    multiplyFullStateDiagMatr(__VA_ARGS__)

#define calcExpecDiagonalOp(...) \
    _WARN_FUNC_RENAMED("calcExpecDiagonalOp()", "calcExpecNonHermitianFullStateDiagMatr()") \
    calcExpecNonHermitianFullStateDiagMatr(__VA_ARGS__)



#define createSubDiagonalOp(...) \
    _WARN_FUNC_RENAMED("createSubDiagonalOp()", "createDiagMatr()") \
    createDiagMatr(__VA_ARGS__)

#define destroySubDiagonalOp(...) \
    _WARN_FUNC_RENAMED("destroySubDiagonalOp()", "destroyDiagMatr()") \
    destroyDiagMatr(__VA_ARGS__)

#define diagonalUnitary(...) \
    _WARN_FUNC_RENAMED("diagonalUnitary()", "applyDiagMatr()") \
    applyDiagMatr(__VA_ARGS__)

#define applySubDiagonalOp(...) \
    _WARN_FUNC_RENAMED("applySubDiagonalOp()", "multiplyDiagMatr()") \
    multiplyDiagMatr(__VA_ARGS__)

static inline void _applyGateSubDiagonalOp(Qureg qureg, int* targets, int numTargets, DiagMatr op) {
    qreal eps = getValidationEpsilon();
    setValidationEpsilon(0);
    applyDiagMatr(qureg, targets, numTargets, op);
    setValidationEpsilon(eps);
}
#define applyGateSubDiagonalOp(...) \
    _WARN_GENERAL_MSG( \
        "The QuEST function 'applyGateSubDiagonalOp()' is deprecated. To achieve the same thing, disable " \
        "numerical validation via 'setValidationEpsilon(0)' before calling 'applyDiagMatr()'. You can " \
        "save the existing epsilon via 'getValidationEpsilon()' to thereafter restore. This procedure " \
        "has been performed here automatically.") \
    _applyGateSubDiagonalOp(__VA_ARGS__)



#define reportStateToScreen(qureg, env, rank) \
    _WARN_FUNC_RENAMED("reportStateToScreen(qureg, env, rank)", "reportQureg(qureg)") \
    reportQureg(qureg)



#define getNumQubits(qureg) \
    _WARN_FUNC_RENAMED("getNumQubits(qureg)", "qureg.numQubits") \
    (qureg).numQubits 

 #define getNumAmps(qureg) \
    _WARN_FUNC_RENAMED("getNumAmps(qureg)", "qureg.numAmps") \
    (qureg).numAmps 



#define setQuregToPauliHamil(...) \
    _WARN_FUNC_RENAMED("setQuregToPauliHamil()", "setQuregToPauliStrSum()") \
    setQuregToPauliStrSum(__VA_ARGS__)

#define cloneQureg(...) \
    _WARN_FUNC_RENAMED("cloneQureg()", "setQuregToClone()") \
    setQuregToClone(__VA_ARGS__)

#define setWeightedQureg(f1, q1, f2, q2, fOut, qOut) \
    _WARN_GENERAL_MSG( \
        "The QuEST function 'setWeightedQureg(f1,q1, f2,q2, fOut,qOut)' is deprecated, and replaced with " \
        "'setQuregToSuperposition(fOut,qOut, f1,q1, f2,q2)' which has been automatically invoked. The new " \
        "fucntion however accepts only statevectors, not density matrices, so may error at runtime. Beware " \
        "that the order of the arguments has changed, so that the first supplied Qureg is modified." ) \
    setQuregToSuperposition( \
        getQcomp(fOut.real, fOut.imag), qOut, \
        getQcomp(f1.real, f1.imag), q1, \
        getQcomp(f2.real, f2.imag), q2)



#define phaseShift(...) \
    _WARN_FUNC_RENAMED("phaseShift()", "applyPhaseShift()") \
    applyPhaseShift(__VA_ARGS__)

#define controlledPhaseShift(...) \
    _WARN_FUNC_RENAMED("controlledPhaseShift()", "applyTwoQubitPhaseShift()") \
    applyTwoQubitPhaseShift(__VA_ARGS__)

#define multiControlledPhaseShift(...) \
    _WARN_FUNC_RENAMED("multiControlledPhaseShift()", "applyMultiQubitPhaseShift()") \
    applyMultiQubitPhaseShift(__VA_ARGS__)

#define controlledPhaseFlip(...) \
    _WARN_FUNC_RENAMED("controlledPhaseFlip()", "applyTwoQubitPhaseFlip()") \
    applyTwoQubitPhaseFlip(__VA_ARGS__)

#define multiControlledPhaseFlip(...) \
    _WARN_FUNC_RENAMED("multiControlledPhaseFlip()", "applyMultiQubitPhaseFlip()") \
    applyMultiQubitPhaseFlip(__VA_ARGS__)



#define sGate(...) \
    _WARN_FUNC_RENAMED("sGate()", "applyS()") \
    applyS(__VA_ARGS__)

#define tGate(...) \
    _WARN_FUNC_RENAMED("tGate()", "applyT()") \
    applyT(__VA_ARGS__)



#define copyStateToGPU(...) \
    _WARN_FUNC_RENAMED("copyStateToGPU()", "syncQuregToGpu()") \
    syncQuregToGpu(__VA_ARGS__)

#define copyStateFromGPU(...) \
    _WARN_FUNC_RENAMED("copyStateFromGPU()", "syncQuregFromGpu()") \
    syncQuregFromGpu(__VA_ARGS__)

#define copySubstateToGPU(...) \
    _WARN_FUNC_RENAMED("copySubstateToGPU()", "syncSubQuregToGpu()") \
    syncSubQuregToGpu(__VA_ARGS__)

#define copySubstateFromGPU(...) \
    _WARN_FUNC_RENAMED("copySubstateFromGPU()", "syncSubQuregFromGpu()") \
    syncSubQuregFromGpu(__VA_ARGS__)



#define getAmp(...) \
    _WARN_FUNC_RENAMED("getAmp()", "getQuregAmp()") \
    getQuregAmp(__VA_ARGS__)

#define getDensityAmp(...) \
    _WARN_FUNC_RENAMED("getDensityAmp()", "getDensityQuregAmp()") \
    getDensityQuregAmp(__VA_ARGS__)

#define getRealAmp(...) \
    _WARN_FUNC_RENAMED("getRealAmp()", "real(getQuregAmp())") \
    real(getQuregAmp(__VA_ARGS__))

#define getImagAmp(...) \
    _WARN_FUNC_RENAMED("getImagAmp()", "imag(getImagAmp())") \
    imag(getQuregAmp(__VA_ARGS__))

#define getProbAmp(...) \
    _WARN_FUNC_RENAMED("getProbAmp()", "calcProbOfBasisState()") \
    calcProbOfBasisState(__VA_ARGS__)

#define calcProbOfOutcome(...) \
    _WARN_FUNC_RENAMED("calcProbOfOutcome()", "calcProbOfQubitOutcome()") \
    calcProbOfQubitOutcome(__VA_ARGS__)

#define calcProbOfAllOutcomes(...) \
    _WARN_FUNC_RENAMED("calcProbOfAllOutcomes()", "calcProbsOfAllMultiQubitOutcomes()") \
    calcProbsOfAllMultiQubitOutcomes(__VA_ARGS__)

#define calcDensityInnerProduct(...) \
    _WARN_FUNC_RENAMED("calcDensityInnerProduct()", "calcInnerProduct()") \
    calcInnerProduct(__VA_ARGS__)

#define calcHilbertSchmidtDistance(...) \
    _WARN_FUNC_RENAMED("calcHilbertSchmidtDistance()", "calcDistance()") \
    calcDistance(__VA_ARGS__)

#define calcExpecPauliHamil(qureg, hamil, workspace) \
    _WARN_FUNC_RENAMED("calcExpecPauliHamil(Qureg, PauliHamil, Qureg)", "calcExpecPauliStrSum(Qureg, PauliStrSum)") \
    calcExpecPauliStrSum(qureg, hamil)



static inline qreal _calcExpecPauliStr(Qureg qureg, int* targs, _NoWarnPauliOpType* enums, int numTargs) {
    
    int codes[100];
    for (int i=0; i<numTargs && i<100; i++)
        codes[i] = enums[i];
    
    return calcExpecPauliStr(qureg, getPauliStr(codes, targs, numTargs));
}

#define calcExpecPauliProd(qureg, targs, paulis, numTargs, workspace) \
    _WARN_FUNC_RENAMED( \
        "calcExpecPauliProd(qureg, targs, paulis, numTargs, workspace)", \
        "calcExpecPauliStr(qureg, getPauliStr(paulis, targs, numTargs))") \
    _calcExpecPauliStr(qureg, targs, paulis, numTargs)



static inline PauliStrSum _createPauliStrSumFromCodes(int numQubits, _NoWarnPauliOpType* allPauliCodes, qreal* termCoeffs, int numTerms) {

    int* targs = (int*) malloc(numQubits * sizeof *targs);
    for (int i=0; i<numQubits; i++)
        targs[i] = i;

    PauliStr* strings = (PauliStr*) malloc(numTerms * sizeof *strings);
    for (int i=0; i<numTerms; i++) {

        int codes[100]; // assumes numQubits<=100
        for (int j=0; j<numQubits && j<100; j++)
            codes[j] = (int) allPauliCodes[i*numQubits+j];

        strings[i] = getPauliStr(codes, targs, numQubits);
    }

    qcomp* coeffs = (qcomp*) malloc(numTerms * sizeof *coeffs);
    for (int i=0; i<numTerms; i++)
        coeffs[i] = termCoeffs[i];

    PauliStrSum sum = createPauliStrSum(strings, coeffs, numTerms);
    free(targs);
    free(strings);
    free(coeffs);

    return sum;
}

static inline qreal _calcExpecPauliSum(Qureg qureg, _NoWarnPauliOpType* allPauliCodes, qreal* termCoeffs, int numTerms) {
    PauliStrSum sum = _createPauliStrSumFromCodes(qureg.numQubits, allPauliCodes, termCoeffs, numTerms);
    qreal out = calcExpecPauliStrSum(qureg, sum);
    destroyPauliStrSum(sum);
    return out;
}

#define calcExpecPauliSum(qureg, paulis, coeffs, numTerms, workspace) \
    _WARN_FUNC_RENAMED("calcExpecPauliSum(Qureg, ...)", "calcExpecPauliStrSum(Qureg, PauliStrSum)") \
    _calcExpecPauliSum(qureg, paulis, coeffs, numTerms)

static inline void _applyPauliSum(Qureg inQureg, _NoWarnPauliOpType* allPauliCodes, qreal* termCoeffs, int numSumTerms, Qureg outQureg) {
    PauliStrSum sum = _createPauliStrSumFromCodes(inQureg.numQubits, allPauliCodes, termCoeffs, numSumTerms);
    setQuregToClone(outQureg, inQureg); 
    multiplyPauliStrSum(outQureg, sum, inQureg);
    destroyPauliStrSum(sum);
}

#define applyPauliSum(...) \
    _WARN_FUNC_RENAMED("applyPauliSum(inQureg, ..., outQureg)", "multiplyPauliStrSum(outQureg, PauliStrSum)") \
    _applyPauliSum(__VA_ARGS__)

static inline void _applyPauliHamil(Qureg inQureg, PauliStrSum hamil, Qureg outQureg) {
    setQuregToClone(outQureg, inQureg); 
    multiplyPauliStrSum(outQureg, hamil, inQureg);
}

#define applyPauliHamil(...) \
    _WARN_FUNC_RENAMED("applyPauliHamil(inQureg, PauliHamil, outQureg)", "multiplyPauliStrSum(qureg, PauliStrSum, workspace)") \
    _applyPauliHamil(__VA_ARGS__)



#define compactUnitary(q, t, a, b) \
    _WARN_FUNC_RENAMED( \
        "compactUnitary(qureg, t, a, b)", \
        "applyCompMatr1(qureg, t, getInlineCompMatr1({{a,-conj(b)},{b,conj(a)}}))") \
    applyCompMatr1(q, t, getInlineCompMatr1({ \
        {_GET_QCOMP_FROM_COMPLEX_STRUCT(a), - conj(_GET_QCOMP_FROM_COMPLEX_STRUCT(b))}, \
        {_GET_QCOMP_FROM_COMPLEX_STRUCT(b),   conj(_GET_QCOMP_FROM_COMPLEX_STRUCT(a))}} ))

#define controlledCompactUnitary(q, c, t, a, b) \
    _WARN_FUNC_RENAMED( \
        "controlledCompactUnitary(qureg, c, t, a, b)", \
        "applyControlledCompMatr1(qureg, c, t, getInlineCompMatr1({{a,-conj(b)},{b,conj(a)}}))") \
    applyControlledCompMatr1(q, c, t, getInlineCompMatr1({ \
        {_GET_QCOMP_FROM_COMPLEX_STRUCT(a), - conj(_GET_QCOMP_FROM_COMPLEX_STRUCT(b))}, \
        {_GET_QCOMP_FROM_COMPLEX_STRUCT(b),   conj(_GET_QCOMP_FROM_COMPLEX_STRUCT(a))}} ))



#define unitary(qureg, targ, ...) \
    _WARN_FUNC_RENAMED("unitary()", "applyCompMatr1()") \
    applyCompMatr1(qureg, targ, _GET_COMP_MATR_1_FROM_COMPLEX_MATRIX_2(__VA_ARGS__))

#define controlledUnitary(qureg, ctrl, targ, ...) \
    _WARN_FUNC_RENAMED("controlledUnitary()", "applyControlledCompMatr1()") \
    applyControlledCompMatr1(qureg, ctrl, targ, _GET_COMP_MATR_1_FROM_COMPLEX_MATRIX_2(__VA_ARGS__))

#define multiControlledUnitary(qureg, ctrls, numctrls, targ, ...) \
    _WARN_FUNC_RENAMED("multiControlledUnitary()", "applyMultiControlledCompMatr1()") \
    applyMultiControlledCompMatr1(qureg, ctrls, numctrls, targ, _GET_COMP_MATR_1_FROM_COMPLEX_MATRIX_2(__VA_ARGS__))

#define multiStateControlledUnitary(qureg, ctrls, states, numctrls, targ, ...) \
    _WARN_FUNC_RENAMED("multiStateControlledUnitary()", "applyMultiStateControlledCompMatr1()") \
    applyMultiStateControlledCompMatr1(qureg, ctrls, states, numctrls, targ, _GET_COMP_MATR_1_FROM_COMPLEX_MATRIX_2(__VA_ARGS__))

#define twoQubitUnitary(qureg, t1, t2, ...) \
    _WARN_FUNC_RENAMED("twoQubitUnitary()", "applyCompMatr2()") \
    applyCompMatr2(qureg, t1, t2, _GET_COMP_MATR_2_FROM_COMPLEX_MATRIX_4(__VA_ARGS__))

#define controlledTwoQubitUnitary(qureg, c, t1, t2, ...) \
    _WARN_FUNC_RENAMED("controlledTwoQubitUnitary()", "applyControlledCompMatr2()") \
    applyControlledCompMatr2(qureg, c, t1, t2, _GET_COMP_MATR_2_FROM_COMPLEX_MATRIX_4(__VA_ARGS__))

#define multiControlledTwoQubitUnitary(qureg, ctrls, nctrls, t1, t2, ...) \
    _WARN_FUNC_RENAMED("multiControlledTwoQubitUnitary()", "applyMultiControlledCompMatr2()") \
    applyMultiControlledCompMatr2(qureg, ctrls, nctrls, t1, t2, _GET_COMP_MATR_2_FROM_COMPLEX_MATRIX_4(__VA_ARGS__))

#define multiQubitUnitary(...) \
    _WARN_FUNC_RENAMED("multiQubitUnitary()", "applyCompMatr()") \
    applyCompMatr(__VA_ARGS__)

#define controlledMultiQubitUnitary(...) \
    _WARN_FUNC_RENAMED("controlledMultiQubitUnitary()", "applyControlledCompMatr()") \
    applyControlledCompMatr(__VA_ARGS__)

#define multiControlledMultiQubitUnitary(...) \
    _WARN_FUNC_RENAMED("multiControlledMultiQubitUnitary()", "applyMultiControlledCompMatr()") \
    applyMultiControlledCompMatr(__VA_ARGS__)



#define applyMatrix2(qureg, targ, ...) \
    _WARN_FUNC_RENAMED("applyMatrix2()", "multiplyCompMatr1()") \
    multiplyCompMatr1(qureg, targ, _GET_COMP_MATR_1_FROM_COMPLEX_MATRIX_2(__VA_ARGS__))

#define applyMatrix4(qureg, targ1, targ2, ...) \
    _WARN_FUNC_RENAMED("applyMatrix4()", "multiplyCompMatr2()") \
    multiplyCompMatr2(qureg, targ1, targ2, _GET_COMP_MATR_2_FROM_COMPLEX_MATRIX_4(__VA_ARGS__))

#define applyMatrixN(...) \
    _WARN_FUNC_RENAMED("applyMatrixN()", "multiplyCompMatr()") \
    multiplyCompMatr(__VA_ARGS__)



static inline void _applyGateMatrixN(Qureg qureg, int* targs, int numTargs, CompMatr u) {
    qreal eps = getValidationEpsilon();
    setValidationEpsilon(0);
    applyCompMatr(qureg, targs, numTargs, u);
    setValidationEpsilon(eps);
}

#define applyGateMatrixN(...) \
    _WARN_GENERAL_MSG( \
        "The QuEST function 'applyGateMatrixN()' is deprecated. To achieve the same thing, disable " \
        "numerical validation via 'setValidationEpsilon(0)' before calling 'applyCompMatr()'. You can " \
        "save the existing epsilon via 'getValidationEpsilon()' to thereafter restore. This procedure " \
        "has been performed here automatically.") \
    _applyGateMatrixN(__VA_ARGS__)

static inline void _applyMultiControlledGateMatrixN(Qureg qureg, int* ctrls, int numCtrls, int* targs, int numTargs, CompMatr u) {
    qreal eps = getValidationEpsilon();
    setValidationEpsilon(0);
    applyMultiControlledCompMatr(qureg, ctrls, numCtrls, targs, numTargs, u);
    setValidationEpsilon(eps);
}

#define applyMultiControlledGateMatrixN(...) \
    _WARN_GENERAL_MSG( \
        "The QuEST function 'applyMultiControlledGateMatrixN()' is deprecated. To achieve the same thing, disable " \
        "numerical validation via 'setValidationEpsilon(0)' before calling 'applyMultiControlledCompMatr()'. You can " \
        "save the existing epsilon via 'getValidationEpsilon()' to thereafter restore. This procedure has been " \
        "performed here automatically.") \
    _applyMultiControlledGateMatrixN(__VA_ARGS__)



#define pauliX(...) \
    _WARN_FUNC_RENAMED("pauliX()", "applyPauliX()") \
    applyPauliX(__VA_ARGS__)

#define pauliY(...) \
    _WARN_FUNC_RENAMED("pauliY()", "applyPauliY()") \
    applyPauliY(__VA_ARGS__)

#define pauliZ(...) \
    _WARN_FUNC_RENAMED("pauliZ()", "applyPauliZ()") \
    applyPauliZ(__VA_ARGS__)

#define controlledPauliX(...) \
    _WARN_FUNC_RENAMED("controlledPauliX()", "applyControlledPauliX()") \
    applyControlledPauliX(__VA_ARGS__)

#define controlledPauliY(...) \
    _WARN_FUNC_RENAMED("controlledPauliY()", "applyControlledPauliY()") \
    applyControlledPauliY(__VA_ARGS__)

#define controlledPauliZ(...) \
    _WARN_FUNC_RENAMED("controlledPauliZ()", "applyControlledPauliZ()") \
    applyControlledPauliZ(__VA_ARGS__)



#define rotateX(...) \
    _WARN_FUNC_RENAMED("rotateX()", "applyRotateX()") \
    applyRotateX(__VA_ARGS__)

#define rotateY(...) \
    _WARN_FUNC_RENAMED("rotateY()", "applyRotateY()") \
    applyRotateY(__VA_ARGS__)

#define rotateZ(...) \
    _WARN_FUNC_RENAMED("rotateZ()", "applyRotateZ()") \
    applyRotateZ(__VA_ARGS__)

#define rotateAroundAxis(q, t, a, v) \
    _WARN_FUNC_RENAMED( \
        "rotateAroundAxis(Qureg, int, qreal, Vector)", \
        "applyRotateAroundAxis(Qureg, int, qreal, qreal vx, qreal vy, qreal vz)") \
    applyRotateAroundAxis(q, t, a, v.x, v.y, v.z)

#define controlledRotateX(...) \
    _WARN_FUNC_RENAMED("controlledRotateX()", "applyControlledRotateX()") \
    applyControlledRotateX(__VA_ARGS__)

#define controlledRotateY(...) \
    _WARN_FUNC_RENAMED("controlledRotateY()", "applyControlledRotateY()") \
    applyControlledRotateY(__VA_ARGS__)

#define controlledRotateZ(...) \
    _WARN_FUNC_RENAMED("controlledRotateZ()", "applyControlledRotateZ()") \
    applyControlledRotateZ(__VA_ARGS__)

#define controlledRotateAroundAxis(q, c, t, a, v) \
    _WARN_FUNC_RENAMED( \
        "controlledRotateAroundAxis(Qureg, int, int, qreal, Vector)", \
        "applyControlledRotateAroundAxis(Qureg, int, int, qreal, qreal vx, qreal vy, qreal vz)") \
    applyControlledRotateAroundAxis(q, c, t, a, v.x, v.y, v.z)



#define hadamard(...) \
    _WARN_FUNC_RENAMED("hadamard()", "applyHadamard()") \
    applyHadamard(__VA_ARGS__) 

#define controlledNot(...) \
    _WARN_FUNC_RENAMED("controlledNot()", "applyControlledPauliX()") \
    applyControlledPauliX(__VA_ARGS__)



#define multiRotateZ(...) \
    _WARN_FUNC_RENAMED("multiRotateZ()", "applyPhaseGadget()") \
    applyPhaseGadget(__VA_ARGS__)

#define multiControlledMultiRotateZ(...) \
    _WARN_FUNC_RENAMED("multiControlledMultiRotateZ()", "applyMultiControlledPhaseGadget()") \
    applyMultiControlledPhaseGadget(__VA_ARGS__)

static inline void _multiRotatePauli(Qureg qureg, int* targs, _NoWarnPauliOpType* enums, int numTargs, qreal angle) {
    int codes[100];
    for (int i=0; i<numTargs && i<100; i++)
        codes[i] = enums[i];
    applyPauliGadget(qureg, getPauliStr(codes, targs, numTargs), angle);
}
#define multiRotatePauli(...) \
    _WARN_FUNC_RENAMED( \
        "multiRotatePauli(qureg, targs, paulis, numTargs, angle)", \
        "applyPauliGadget(qureg, getPauliStr(paulis, targs, numTargs), angle)") \
    _multiRotatePauli(__VA_ARGS__)

static inline void _multiControlledMultiRotatePauli(Qureg qureg, int* ctrls, int numCtrls, int* targs, _NoWarnPauliOpType* enums, int numTargs, qreal angle) {
    int codes[100];
    for (int i=0; i<numTargs && i<100; i++)
        codes[i] = enums[i];
    applyMultiControlledPauliGadget(qureg, ctrls, numCtrls, getPauliStr(codes, targs, numTargs), angle);
}
#define multiControlledMultiRotatePauli(...) \
    _WARN_FUNC_RENAMED( \
        "multiControlledMultiRotatePauli(qureg, ctrls, numCtrls, targs, paulis, numTargs, angle)", \
        "applyMultiControlledPauliGadget(qureg, ctrls, numCtrls, getPauliStr(paulis, targs, numTargs), angle)") \
    _multiControlledMultiRotatePauli(__VA_ARGS__)

#define multiQubitNot(...) \
    _WARN_FUNC_RENAMED("multiQubitNot()", "applyMultiQubitNot()") \
    applyMultiQubitNot(__VA_ARGS__)

#define multiControlledMultiQubitNot(...) \
    _WARN_FUNC_RENAMED("multiControlledMultiQubitNot()", "applyMultiControlledMultiQubitNot()") \
    applyMultiControlledMultiQubitNot(__VA_ARGS__)



#define collapseToOutcome(...) \
    _WARN_FUNC_RENAMED("collapseToOutcome()", "applyForcedQubitMeasurement()") \
    applyForcedQubitMeasurement(__VA_ARGS__)

#define measure(...) \
    _WARN_FUNC_RENAMED("measure()", "applyQubitMeasurement()") \
    applyQubitMeasurement(__VA_ARGS__)

#define measureWithStats(...) \
    _WARN_FUNC_RENAMED("measureWithStats()", "applyQubitMeasurementAndGetProb()") \
    applyQubitMeasurementAndGetProb(__VA_ARGS__)

#define applyProjector(...) \
    _WARN_FUNC_RENAMED("applyProjector()", "applyQubitProjector()") \
    applyQubitProjector(__VA_ARGS__)



#define mixPauli(...) \
    _WARN_FUNC_RENAMED("mixPauli()", "mixPaulis()") \
    mixPaulis(__VA_ARGS__)

#define mixDensityMatrix(outQureg, prob, inQureg) \
    _WARN_FUNC_RENAMED("mixDensityMatrix(outQureg, prob, inQureg)", "mixQureg(outQureg, inQureg, prob)") \
    mixQureg(outQureg, inQureg, prob)



#define swapGate(...) \
    _WARN_FUNC_RENAMED("swapGate()", "applySwap()") \
    applySwap(__VA_ARGS__)

#define sqrtSwapGate(...) \
    _WARN_FUNC_RENAMED("sqrtSwapGate()", "applySqrtSwap()") \
    applySqrtSwap(__VA_ARGS__) 



#define applyTrotterCircuit(...) \
    _WARN_FUNC_RENAMED("applyTrotterCircuit(..., PauliHamil, ...)", "applyTrotterizedPauliStrSumGadget(..., PauliStrSum, ...)") \
    applyTrotterizedPauliStrSumGadget(__VA_ARGS__)

#define applyFullQFT(...) \
    _WARN_FUNC_RENAMED("applyFullQFT()", "applyFullQuantumFourierTransform()") \
    applyFullQuantumFourierTransform(__VA_ARGS__)

#define applyQFT(...) \
    _WARN_FUNC_RENAMED("applyQFT()", "applyQuantumFourierTransform()") \
    applyQuantumFourierTransform(__VA_ARGS__)



#define seedQuESTDefault(...) \
    _WARN_FUNC_RENAMED("seedQuESTDefault(QuESTEnv)", "setSeedsToDefault()") \
    setSeedsToDefault()

#define seedQuEST(env, seeds, numSeeds) \
    _WARN_FUNC_RENAMED("seedQuEST(QuESTEnv, unsigned long int*, int)", "setSeeds(unsigned*, int)") \
    setSeeds(seeds, numSeeds)



#endif // DEPRECATED_H
