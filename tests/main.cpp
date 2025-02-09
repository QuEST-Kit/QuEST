#include <catch2/catch_session.hpp>
#include <stdexcept>

#include "quest.h"
#include "tests/utils/cache.hpp"


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
    createCachedQuregs();

    // TODO:
    // is there some way for us to announce what deployment modes will be run?!?!

    int result = Catch::Session().run( argc, argv );

    destroyCachedQuregs();
    finalizeQuESTEnv();
    return result;
}
