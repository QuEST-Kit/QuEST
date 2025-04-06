<!--
  Tests
  
  @author Tyson Jones
-->

# 🧪  Tests

This folder contains QuEST's extensive tests. See [`compile.md`](/docs/compile.md#tests) and [`launch.md`](/docs/launch.md#tests) to get them running.

The subdirectories are:
- [`utils`](utils/) containing test utilities, including non-optimised functions against which QuEST output is compared.
- [`unit`](unit) containing [unit tests](https://en.wikipedia.org/wiki/Unit_testing) which test individual QuEST functions in isolation, under their entire input domains (where feasible).
- [`integration`](integration/) containing [integration tests](https://en.wikipedia.org/wiki/Integration_testing) which test multiple QuEST functions working at scale.
- [`deprecated`](deprecated/) containing `v3`'s tests and utilities, only used when explicitly [activated](/docs/compile.md#v3).

The tests use [Catch2](https://github.com/catchorg/Catch2) and are generally structured as
```cpp
TEST_CASE( "someApiFunc", "[funcs]" ) {

    PREPARE_TEST(...)

    SECTION( "correctness" ) {

        SECTION( "statevector" ) {

            auto a = getApiResult();
            auto b = getReferenceResult();
            REQUIRE_AGREE(a, b);
        }

        SECTION( "density matrix" ) {

            auto a = getApiResult();
            auto b = getReferenceResult();s
        }
    }

    SECTION( "validation" ) {

        SECTION( "some way to mess it up" ) {

            REQUIRE_THROWS( someApiFunc(badArgs) );
        }
    }
}
```
