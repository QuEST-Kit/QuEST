#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

#include <stdexcept>
#include <iostream>

#include "quest.h"
#include "tests/utils/cache.hpp"


/*
 * recast QuEST errors into exceptions which Catch can intercept
 */

extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {

    throw std::runtime_error(errMsg);
}


/*
 * report QuEST env when tests start
 */

class startListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;
    void testRunStarting(Catch::TestRunInfo const&) override {

        // a full report is too verbose...
        // reportQuESTEnv();

        // so we summarise the important info
        QuESTEnv env = getQuESTEnv();
        std::cout << std::endl;
        std::cout << "QuEST testing execution environment:" << std::endl;
        std::cout << "  precision:       " << FLOAT_PRECISION << std::endl;
        std::cout << "  multithreaded:   " << getQuESTEnv().isMultithreaded << std::endl;
        std::cout << "  distributed:     " << getQuESTEnv().isMultithreaded << std::endl;
        std::cout << "  GPU-accelerated: " << getQuESTEnv().isMultithreaded << std::endl;
        std::cout << "  cuQuantum:       " << COMPILE_CUQUANTUM << std::endl;
        std::cout << std::endl;
    }
};

CATCH_REGISTER_LISTENER(startListener)


/*
 * setup QuEST before Catch2 session
 */

int main(int argc, char* argv[]) {

    initQuESTEnv();
    createCachedQuregs();

    // disable Catch2 output on non-root nodes
    if (getQuESTEnv().rank != 0)
        std::cout.rdbuf(NULL);

    // launch Catch2, triggering above event listener
    int result = Catch::Session().run( argc, argv );

    destroyCachedQuregs();
    finalizeQuESTEnv();
    return result;
}
