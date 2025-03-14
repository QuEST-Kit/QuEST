<!--
  A style guide for QuEST contributors
  
  @author Tyson Jones
-->

>TODO!

Don't agonise about style - write your code as you see fit and we can address major issues in review/PR.

Some stylistic choices Tyson personally makes:
- use `auto` where it improves readability, discretionarily. Obviously it is better than massive, unimportant types of objects or heavily templated collections, but sometimes knowing the precise type of a primitive is helpful
- I avoid superfluous braces (even despite being a major proponent of defensive design):
  ```C++
  if (cond)
      return x;
  ```
  Sue me!


- _never_ use `abs()`, always use `std::abs()`. This is because `abs` may actually call `C`'s `abs(int)` function which is utilised by the v3 deprecated unit tests. As such, `abs(0.999)` will return `0` causing insidious bugs. Use of `real` and `imag` (in lieu of `std::real` and `std::imag`) are okay since always defined in terms of floats.