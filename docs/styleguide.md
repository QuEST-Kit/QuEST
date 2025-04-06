# 🎨  Style guide

<!--
  A style guide for QuEST contributors
  (this comment must be under the title for valid doxygen rendering)
  
  @author Tyson Jones
-->



Don't agonise about style - write your code as you see fit and we can address major issues in review/PR.
Some encouraged conventions include:

- use `camelCase` for everything except:
  - constants which use `CAPITALS_AND_UNDERSCORES`
  - related function prefixes, like `prefix_someFunction()`
- favour clarity over concision, for example
  ```cpp
  qcomp elem = state[ind][ind];
  qreal prob = std::real(elem);
  return prob;
  ```
  over
  ```cpp
  return std::real(state[ind][ind]);
  ```
- never ever do:
  ```cpp
  using namespace std;
  ```
  but _do_ shorten common containers like `vector`:
  ```cpp
  using std::vector;

  vector<int> mylist;
  ```
- whitespace is free; use it wherever it can improve clarity, like to separate subroutines.
  ```cpp
  // i000 = nth local index where all suffix bits are 0
  qindex i000 = insertThreeZeroBits(n, braQb1, ketQb2, ketQb1);
  qindex i0b0 = setBit(i000, ketQb2, braBit2);
  qindex i1b1 = flipTwoBits(i0b0, braQb1, ketQb1);

  // j = nth received amp in buffer
  qindex j = n + offset;

  // mix pair of amps using buffer
  qcomp amp0b0 = qureg.cpuAmps[i0b0];
  qcomp amp1b1 = qureg.cpuAmps[i1b1];

  qureg.cpuAmps[i0b0] = c1*amp0b0 + c2*(amp1b1 + qureg.cpuCommBuffer[j]);
  qureg.cpuAmps[i1b1] = c1*amp1b1 + c2*(amp0b0 + qureg.cpuCommBuffer[j]);
  ```
- use `auto` where it improves readability, discretionarily. Obviously it is better than massive, unimportant types of objects or heavily templated collections, but sometimes knowing the precise type of a primitive is helpful
- It is permissable to avoid superfluous braces around single-line branches:
  ```cpp
  if (cond)
      return x;
  ```
- always prefix calls to mathematical functions like `abs()` with the `std` namespace, i.e. `std::abs()`. This avoids ambiguity with `C` overloads like `abs(int)` which can cause insidious bugs! The full list of functions to prefix are:
  - `abs`
  - `real`
  - `imag`
  - `conj`
  - `norm`
  - `sin`
  - `cos`
  - `log`
  - `log2`
  - `exp`
  - `pow`
  - `sqrt`
  - `floor`
  - `ceil`
  - `atan2`
  - `min`
  - `max`
