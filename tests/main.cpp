#include <catch2/catch_session.hpp>

// TODO:
// when we switch to CMake-supplied Catch2,
// we must replace the above include with:
// #include <catch2/catch_session.hpp>


#include "quest.h"
#include <stdexcept>


// TODO:
// implement a custom reporter in order to avoid output
// duplication when running tests distributed


// recast QuEST errors into exceptions which Catch can catch  
extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {

    throw std::runtime_error(errMsg);
}


// custom catch2 main so that we can prepare QuEST (needed due to MPI)
int main(int argc, char* argv[]) {

    initQuESTEnv();

    int result = Catch::Session().run( argc, argv );

    finalizeQuESTEnv();
    return result;
}
