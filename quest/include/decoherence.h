/** @file
 * API signatures for effecting decohering channels upon Quregs
 * which are instantiated as density matrices.
 * 
 * @author Tyson Jones
 * 
 * @defgroup decoherence Decoherence
 * @ingroup api
 * @brief Functions for effecting decoherence channels upon density matrices.
 * @{
 */

#ifndef DECOHERENCE_H
#define DECOHERENCE_H

#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/channels.h"



/*
 * C AND C++ AGNOSTIC FUNCTIONS
 */

// enable invocation by both C and C++ binaries
#ifdef __cplusplus
extern "C" {
#endif


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob and @f$ t = @f$ @p target. 
 * 
 * This function effects
 * @f[ 
        \dmrho 
            \;\rightarrow\;
        (1 - p) \, \dmrho 
            \,+\, 
        p \, \hat{Z}_t \,\dmrho\, \hat{Z}_t.
 * @f]
 * 
 * @equivalences
 * This function is equivalent to (but much faster than):
 * - mixPaulis() with a zero probability for the @f$\hat{X}@f$ and @f$\hat{Y}@f$ components.
 *   ```
    mixPaulis(qureg, target, 0, 0, prob);
 *   ```
 * - mixKrausMap() with (scaled) @f$\hat{\id}@f$ and @f$\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-prob);
    qreal b = sqrt(prob);

    KrausMap map = createInlineKrausMap(1, 2, {
        {{a,0},{0, a}}, // a * I
        {{b,0},{0,-b}}  // b * Z
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 * - mixQureg() with a duplicated Qureg modified under applyPauliZ().
 *   ```
    Qureg clone = createCloneQureg(qureg);
    applyPauliZ(clone);
    mixQureg(qureg, other, prob);
 *   ```
 *
 * @notyetvalidated
 */
void mixDephasing(Qureg qureg, int target, qreal prob);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob, @f$ t_1 = @f$ @p target1 and @f$ t_2 = @f$ @p target2.
 * 
 * This function effects
 * @f[ 
        \dmrho 
            \;\rightarrow\;
            (1 - p) \, \dmrho 
                \,+\,
            \frac{p}{3} \left(
            \hat{Z}_{t_1} \dmrho \hat{Z}_{t_1}
                \,+\,
            \hat{Z}_{t_1} \dmrho \hat{Z}_{t_1}
                \,+\,
            \hat{Z}_{t_1}  \hat{Z}_{t_2} \dmrho \hat{Z}_{t_1} \hat{Z}_{t_2}
        \right).
 * @f]
 * 
 * @equivalences
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with (scaled) @f$\hat{\id}\otimes\hat{\id}@f$, @f$\hat{\id}\otimes\hat{Z}@f$, 
 *   @f$\hat{Z}\otimes\hat{\id}@f$ and @f$\hat{Z}\otimes\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-prob);
    qreal b = sqrt(prob/3);

    KrausMap map = createInlineKrausMap(2, 4, {
        {{a,0,0,0},{0, a,0,0},{0,0, a,0},{0,0,0, a}}, // a * II
        {{b,0,0,0},{0,-b,0,0},{0,0, b,0},{0,0,0,-b}}, // b * IZ
        {{b,0,0,0},{0, b,0,0},{0,0,-b,0},{0,0,0,-b}}, // b * ZI
        {{b,0,0,0},{0,-b,0,0},{0,0,-b,0},{0,0,0, b}}  // b * ZZ
    });
    
    int targets[] = {target1, target2};
    mixKrausMap(qureg, targets, 2, map);
 *   ```
 *
 * @notyetvalidated
 */
void mixTwoQubitDephasing(Qureg qureg, int target1, int target2, qreal prob);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob and @f$ t = @f$ @p target. 
 * 
 * This function effects
 * @f[ 
        \dmrho \;\rightarrow\; 
            (1 - p) \, \dmrho \,+\, \frac{p}{3} \left( 
            \hat{X}_t \dmrho \hat{X}_t \,+\,
            \hat{Y}_t \dmrho \hat{Y}_t \,+\,
            \hat{Z}_t \dmrho \hat{Z}_t
        \right).
 * @f]
 * 
 * @equivalences
 * This function is equivalent to (but much faster than):
 * - mixPaulis() with a uniform probability.
 *   ```
    mixPaulis(qureg, target, prob/3, prob/3, prob/3);
 *   ```
 * - mixKrausMap() with (scaled) @f$\hat{\id}@f$, @f$\hat{X}@f$, @f$\hat{Y}@f$ and @f$\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-prob);
    qreal b = sqrt(prob/3);

    KrausMap map = createInlineKrausMap(1, 4, {
        {{a,0},{0, a}}, // a * I
        {{0,b},{b, 0}}, // b * X
        {{b,0},{0,-b}}  // b * Z
        {{0,-1i*b},{1i*b,0}}, // b * Y
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 *
 * @notyetvalidated
 */
void mixDepolarising(Qureg qureg, int target, qreal prob);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob, @f$ t_1 = @f$ @p target1 and @f$ t_2 = @f$ @p target2.
 * 
 * This function effects:
 * @f[
        \dmrho \; \rightarrow \;
            (1 - p) \dmrho 
                +
            \frac{p}{15} \left(
                \sum_{\hat{\sigma} \in \{\hat{\id},\hat{X},\hat{Y},\hat{Z}\}}
                \sum_{\hat{\sigma}' \in \{\hat{\id},\hat{X},\hat{Y},\hat{Z}\}}
                \hat{\sigma}_{t_1} \hat{\sigma}_{t_2}'
                \; \dmrho \;
                \hat{\sigma}_{t_1} \hat{\sigma}_{t_2}'
            \right)
            - \frac{p}{15} \hat{\id}_{t_1} \hat{\id}_{t_2} \dmrho \hat{\id}_{t_1} \hat{\id}_{t_2},
 * @f]
 *
 * or verbosely:
 * 
 * @f[
    \dmrho \; \rightarrow \;
    (1 - p) \, \rho + \frac{p}{15} \; 
    \left( 
    \begin{gathered}
        \hat{X}_{t_1} \, \rho \, \hat{X}_{t_1} + 
        \hat{Y}_{t_1} \, \rho \, \hat{Y}_{t_1} + 
        \hat{Z}_{t_1} \, \rho \, \hat{Z}_{t_1} + 
            \\
        \hat{X}_{t_2} \, \rho \, \hat{X}_{t_2} + 
        \hat{Y}_{t_2} \, \rho \, \hat{Y}_{t_2} + 
        \hat{Z}_{t_2} \, \rho \, \hat{Z}_{t_2} + 
            \\
        \hat{X}_{t_1} \hat{X}_{t_2} \, \rho \, \hat{X}_{t_1} \hat{X}_{t_2} + 
        \hat{Y}_{t_1} \hat{Y}_{t_2} \, \rho \, \hat{Y}_{t_1} \hat{Y}_{t_2} + 
        \hat{Z}_{t_1} \hat{Z}_{t_2} \, \rho \, \hat{Z}_{t_1} \hat{Z}_{t_2} +
            \\
        \hat{X}_{t_1} \hat{Y}_{t_2} \, \rho \, \hat{X}_{t_1} \hat{Y}_{t_2} + 
        \hat{Y}_{t_1} \hat{Z}_{t_2} \, \rho \, \hat{Y}_{t_1} \hat{Z}_{t_2} + 
        \hat{Z}_{t_1} \hat{X}_{t_2} \, \rho \, \hat{Z}_{t_1} \hat{X}_{t_2} +
            \\
        \hat{X}_{t_1} \hat{Z}_{t_2} \, \rho \, \hat{X}_{t_1} \hat{Z}_{t_2} + 
        \hat{Y}_{t_1} \hat{X}_{t_2} \, \rho \, \hat{Y}_{t_1} \hat{Z}_{t_2} + 
        \hat{Z}_{t_1} \hat{Y}_{t_2} \, \rho \, \hat{Z}_{t_1} \hat{Y}_{t_2}
    \end{gathered}
    \right).
 * @f]
 *
 * @equivalences
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with Kraus operators containing every possible tensor product 
 *   of two Pauli matrices, all scaled by @f$ (p/15)^{1/2} @f$, _except_ for
 *   @f$ \hat{\id} \otimes \hat{\id} @f$ which is scaled by @f$ (1-16p/15)^{1/2} @f$.
 *
 * @notyetvalidated
 */
void mixTwoQubitDepolarising(Qureg qureg, int target1, int target2, qreal prob);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ p = @f$ @p prob and @f$ t = @f$ @p target.
 * 
 * This function effects
 * @f[ 
        \dmrho \; \rightarrow \; 
        \hat{K}_t^{(1)} \dmrho \, {\hat{K}_t^{(2)}}^\dagger 
            \,+\,
        \hat{K}_t^{(2)} \dmrho \, {\hat{K}_t^{(2)}}^\dagger
 * @f]
 * where @f$ \hat{K}^{(1)} @f$ and @f$ \hat{K}^{(2)} @f$ are one-qubit Kraus operators
 * @f[
    \hat{K}^{(1)} = \begin{pmatrix} 1 & 0 \\ 0 & \sqrt{1-p} \end{pmatrix},
    \;\;
    \hat{K}^{(1)} = \begin{pmatrix} 0 & \sqrt{p} \\ 0 & 0 \end{pmatrix}.
 * @f]
 * 
 * @equivalences
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with the above Kraus operators.
 *   ```
    KrausMap map = createInlineKrausMap(1, 2, {
        {{1,0},{0,sqrt(1-prob)}}, // K1
        {{0,sqrt(p)},{0,0}}       // K2
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 *
 * @notyetvalidated
 */
void mixDamping(Qureg qureg, int target, qreal prob);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ t = @f$ @p target, and
 * @f$ p_x = @f$ @p probX, @f$ p_y = @f$ @p probY, @f$ p_z = @f$ @p probZ.
 * 
 * This function effects
 * @f[ 
        \dmrho \;\rightarrow\;
             (1 - p) \, \dmrho 
             \,+\,
            p_x \, \hat{X}_t \dmrho \hat{X}_t
             \,+\,
            p_y \, \hat{Y}_t \dmrho \hat{Y}_t
             \,+\,
            p_z \, \hat{Z}_t \dmrho \hat{Z}_t.
 * @f]
 * 
 * @equivalences
 * This function is equivalent to (but much faster than):
 * - mixKrausMap() with (scaled) @f$\hat{\id}@f$, @f$\hat{X}@f$, @f$\hat{Y}@f$ and @f$\hat{Z}@f$ Kraus operators.
 *   ```
    qreal a = sqrt(1-probX-probY-probZ);
    qreal b = sqrt(probX);
    qreal c = sqrt(probY);
    qreal d = sqrt(probZ);

    KrausMap map = createInlineKrausMap(1, 4, {
        {{a,0},{0, a}}, // a * I
        {{0,b},{b, 0}}, // b * X
        {{d,0},{0,-d}}  // d * Z
        {{0,-1i*c},{1i*c,0}}, // c * Y
    });
    
    mixKrausMap(qureg, &target, 1, map);
 *   ```
 *
 * @notyetvalidated
 */
void mixPaulis(Qureg qureg, int target, qreal probX, qreal probY, qreal probZ);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho_1 = @f$ @p qureg, @f$ \dmrho_2 = @f$ @p other and @f$ p = @f$ @p prob.
 * 
 * This function effects
 * @f[ 
        \dmrho_1 \;\rightarrow \;
            (1 - p) \, \dmrho_1 
                \,+\,
            p \, \dmrho_2.
 * @f]
 * 
 * @notyetvalidated
 */
void mixQureg(Qureg qureg, Qureg other, qreal prob);


/** @notyetdoced
 * 
 * @formulae
 * Let @f$ \dmrho = @f$ @p qureg, @f$ \vec{t} = @f$ @p targets and @f$ \hat{K}^{(i)} @f$
 * denote the @f$i@f$-th Kraus operator in @p map.
 * 
 * This function effects
 * @f[ 
        \dmrho \; \rightarrow \; 
        \sum\limits_i
        \hat{K}_{\vec{t}}^{(i)} \dmrho \, {\hat{K}_{\vec{t}}^{(i)}}^\dagger
 * @f]
 *
 * @notyetvalidated
 */
void mixKrausMap(Qureg qureg, int* targets, int numTargets, KrausMap map);


/// @notyetdoced
/// @notyetvalidated
void mixSuperOp(Qureg qureg, int* targets, int numTargets, SuperOp superop);


// end de-mangler
#ifdef __cplusplus
}
#endif



/*
 * C++ OVERLOADS
 *
 * which are only accessible to C++ binaries, and accept
 * arguments more natural to C++ (e.g. std::vector). These
 * are included in the file-wide doxygen group (no subgroups).
 */

#ifdef __cplusplus

#include <vector>

/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cppvectoroverload
/// @see mixKrausMap()
void mixKrausMap(Qureg qureg, std::vector<int> targets, KrausMap map);

/// @notyettested
/// @notyetdoced
/// @notyetvalidated
/// @cppvectoroverload
/// @see mixSuperOp()
void mixSuperOp(Qureg qureg, std::vector<int> targets, SuperOp superop);

#endif // __cplusplus



#endif // DECOHERENCE_H

/** @} */ // (end file-wide doxygen defgroup)
