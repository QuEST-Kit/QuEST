/** @file
 * Entry-point for all unit and integration tests.
 *
 * @author Tyson Jones
 * 
 * @defgroup tests Tests
 * 
 * @defgroup testutils Utilities
 * @ingroup tests
 * 
 * @defgroup unittests Unit tests
 * @ingroup tests
 * @todo
 *   we cannot really document the v4 unit tests this way because
 *   many are produced by macros rather than individual TEST_CASE, hmm
 * 
 * @defgroup integrationtests Integration tests
 * @ingroup tests
 * 
 * @defgroup deprecatedtests Deprecated tests
 * @ingroup tests
 * 
 * @defgroup deprecatedutils Deprecated utilities
 * @ingroup tests
 */

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>

#include <stdexcept>
#include <iostream>
#include <string>

#include "quest/include/quest.h"
#include "tests/utils/cache.hpp"
#include "tests/utils/macros.hpp"


/*
 * recast QuEST errors into exceptions which Catch can intercept
 */

extern "C" void invalidQuESTInputError(const char* errMsg, const char* errFunc) {

    throw std::runtime_error(std::string(errFunc) + ": " + std::string(errMsg));
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
        std::cout << "QuEST execution environment:" << std::endl;
        std::cout << "  precision:       " << FLOAT_PRECISION      << std::endl;
        std::cout << "  multithreaded:   " << env.isMultithreaded  << std::endl;
        std::cout << "  distributed:     " << env.isDistributed    << std::endl;
        std::cout << "  GPU-accelerated: " << env.isGpuAccelerated << std::endl;
        std::cout << "  cuQuantum:       " << (env.isGpuAccelerated && COMPILE_CUQUANTUM) << std::endl;
        std::cout << "  node count:      " << env.numNodes         << std::endl;
        std::cout << "  unit Qureg size: " << getNumCachedQubits() << std::endl;
        std::cout << std::endl;

        std::cout << "Tested Qureg deployments:" << std::endl;
        for (auto& [label, qureg]: getCachedStatevecs())
            std::cout << "  " << label << std::endl;
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
