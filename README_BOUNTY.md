# unitaryHACK26 Submission: QuEST Core Optimizations

This pull request resolves two open bounties for the QuEST (Quantum Exact Simulation Toolkit) framework: the implementation of an $\mathcal{O}(k)$ graph routing algorithm for SWAP fusion, and the integration of Neumaier compensated summation for numerically stable dense matrix operations.

---

## 1. SWAP Fusion: $\mathcal{O}(1)$ Fused State Vector Permutation

### The Bottleneck
Applying $k$ sequential SWAP operations to a quantum state requires $k$ independent passes over the $\mathcal{O}(2^n)$ state vector array. For large systems, this severely bottlenecks on memory bandwidth and destroys CPU cache locality.

### The Implementation
We introduce a new public API function: `applyMultiSwap(Qureg qureg, const int* targs1, const int* targs2, int numSwaps)`. 

Instead of executing the SWAPs sequentially, the algorithm processes the execution queue to form a bipartite graph of logical-to-physical memory maps. The bit-permutation is calculated in $\mathcal{O}(k)$ time. The final amplitudes are then moved to their precise memory addresses in a single, cache-friendly $\mathcal{O}(2^n)$ pass.

### Rigorous Correctness Proof
**Theorem:** *Given a set of $k$ disjoint qubit-index transpositions, the induced state vector permutation can be fully resolved in-place in a single $\mathcal{O}(2^n)$ iteration using the guard $\pi(i) > i$.*

**Proof:**
Let the state vector amplitudes be indexed by $i \in \{0, 1, \dots, 2^n - 1\}$. Let $T = \{ \tau_1, \tau_2, \dots, \tau_k \}$ be a set of disjoint transpositions acting on the $n$ qubit indices, where $\tau_m = (a_m, b_m)$ and $\{a_m, b_m\} \cap \{a_{m'}, b_{m'}\} = \emptyset$ for $m \neq m'$. 

The induced permutation $\pi$ on the amplitude index $i$ flips the $a_m$-th and $b_m$-th bits of $i$ if and only if those bits differ. Because the transpositions are disjoint, their bit-flips are completely orthogonal. An index $i$ may have its bit-pairs flipped by multiple independent transpositions simultaneously, but because these flips commute, applying $\pi$ twice restores all original bit values:

$$\pi(\pi(i)) = i \quad \forall i \in \{0, \dots, 2^n - 1\}$$

Because $\pi^2 = \text{id}$, the permutation $\pi$ is an involution. The disjoint cycle decomposition of an involution consists exclusively of fixed points and 2-cycles.
* **Fixed Points ($\pi(i) = i$):** The guard $\pi(i) > i$ evaluates to `false`. No memory swap occurs.
* **2-Cycles ($\pi(i) = j, j \neq i$):** For every cycle $(i, j)$, exactly one element satisfies $i < j$. When the iterator reaches the smaller index, the guard $\pi(i) > i$ is `true`, and the amplitudes are swapped. When the iterator reaches the larger index $j$, because $\pi$ is an involution, $\pi(j) = i < j$. The guard $\pi(j) > j$ is `false`, strictly preventing the reversion of the swap.

Every 2-cycle is processed exactly once. The complete multi-qubit permutation is applied accurately in-place, requiring only $\mathcal{O}(1)$ auxiliary memory. $\blacksquare$

---

## 2. Numerical Stability: Neumaier Compensated Summation

### The Bottleneck
In `cpu_statevec_anyCtrlAnyTargDenseMatr_sub`, the application of an $N$-target dense complex matrix requires iterating $2^N$ times, linearly combining dynamic amplitudes via the standard `+=` accumulation. In IEEE 754 arithmetic, this standard addition suffers from catastrophic cancellation when summing highly oscillatory complex amplitudes. The arithmetic error scales as $\mathcal{O}(n \varepsilon)$, destroying the strict unitarity ($\text{Tr}(\rho) = 1$) of the density matrix for massive state vectors.

### The Implementation
We replaced the naive accumulation inner loops with a **Neumaier Summation** algorithm (an improvement over Kahan summation that covers instances where the next term is larger than the running sum). This isolates the low-order bits lost during floating-point rounding into an independent error accumulator, reducing the overall algorithmic error bound to effectively $\mathcal{O}(\varepsilon)$.

### Compiler Flag Override
Compiling QuEST with `-Ofast` or `-ffast-math` forces the compiler to reassociate floating-point operations, which mathematically annihilates the Neumaier error compensation logic. To prevent this without disabling global optimization, the specific inner loop function is safeguarded using a function-specific compiler attribute:

```cpp
__attribute__((optimize("no-fast-math")))
void cpu_statevec_anyCtrlAnyTargDenseMatr_sub(...) {
    // Neumaier logic implemented via kahan.hpp
}