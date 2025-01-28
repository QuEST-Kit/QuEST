#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include "macros.hpp"

#include <vector>
using std::vector;

// TODO: 
// use CMake to obtain catch2 and replace #include "catch.hpp" with:
// #include <catch2/generators/catch_generators.hpp>
// #include <catch2/generators/catch_generators_adapters.hpp>



/*
 * RANGE
 */


vector<int> getRange(int start, int endExcl) {
    DEMAND( endExcl >= start );

    vector<int> out(endExcl - start);

    for (size_t i=0; i<out.size(); i++)
        out[i] = start + i;

    return out;
}


vector<int> getRange(int endExcl) {
    return getRange(0, endExcl);
}



/*
 * GENERATORS
 */


// TODO:
// rework v3 generators to use vectors
