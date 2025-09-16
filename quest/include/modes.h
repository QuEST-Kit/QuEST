/** @file
 * Constants related to configuring QuEST runtime modes,
 * and documentation of environment variables
 * 
 * @author Tyson Jones
 * 
 * @defgroup modes Modes
 * @ingroup api
 * @brief Constants and environment variables for controlling QuEST execution.
 * @{
 */

#ifndef MODES_H
#define MODES_H



// document environment variables

// spoof env-vars as consts to doc (hackily and hopefully temporarily)
#if 0


    /** @envvardoc
     * 
     * Specifies whether to permit multiple MPI processes to deploy to the same GPU.
     * 
     * @attention 
     * This environment variable has no effect when either (or both) of distribution or 
     * GPU-acceleration are disabled.
     * 
     * In multi-GPU execution, which combines distribution with GPU-acceleration, it is 
     * prudent to assign each GPU to at most one MPI process in order to avoid superfluous 
     * slowdown. Hence by default, initQuESTEnv() will forbid assigning multiple MPI processes 
     * to the same GPU. This environment variable can be set to `1` to disable this validation, 
     * permitting sharing of a single GPU, as is often useful for debugging or unit testing 
     * (for example, testing multi-GPU execution when only a single GPU is available).
     * 
     * @warning
     * Permitting GPU sharing may cause unintended behaviour when additionally using cuQuantum.
     * 
     * @envvarvalues
     *  - forbid sharing: @p 0, @p '0', @p '', @p , (unspecified)
     *  - permit sharing: @p 1, @p '1'
     * 
     * @author Tyson Jones
     */
    const int PERMIT_NODES_TO_SHARE_GPU = 0;


    /** @envvardoc
     * 
     * Specifies the default validation epsilon. 
     * 
     * Specifying `DEFAULT_VALIDATION_EPSILON` to a positive, real number overrides the 
     * precision-specific default (`1E-5`, `1E-12`, `1E-13` for single, double and quadruple 
     * precision respectively). The specified epsilon is used by QuEST for numerical validation
     * unless overriden at runtime via setValidationEpsilon(), in which case it can be
     * restored to that specified by this environment variable using setValidationEpsilonToDefault().
     * 
     * @envvarvalues
     *  - setting @p DEFAULT_VALIDATION_EPSILON=0 disables numerical validation, as if the value
     *    were instead infinity.
     *  - setting @p DEFAULT_VALIDATION_EPSILON='' is equivalent to _not_ specifying the variable,
     *    adopting instead the precision-specific default above.
     *  - setting @p DEFAULT_VALIDATION_EPSILON=x where `x` is a positive, valid `qreal` in any
     *    format accepted by `C` or `C++` (e.g. `0.01`, `1E-2`, `+1e-2`) will use `x` as the
     *    default validation epsilon.
     * 
     * @constraints
     * The function initQuESTEnv() will throw a validation error if:
     *   - The specified epsilon must be `0` or positive.
     *   - The specified epsilon must not exceed that maximum or minimum value which can be stored
     *     in a `qreal`, which is specific to its precision.
     * 
     * @author Tyson Jones
     */
    const qreal DEFAULT_VALIDATION_EPSILON = 0;


#endif



// user flags for choosing automatic deployment; only accessible by C++ 
// backend and C++ users; C users must hardcode -1 

#ifdef __cplusplus

namespace modeflag { 

    extern int USE_AUTO;
}

#endif // __cplusplus



#endif // MODES_H

/** @} */ // (end file-wide doxygen defgroup)
