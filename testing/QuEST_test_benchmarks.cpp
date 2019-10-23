/** @file
 * @author Tyson Jones
 */

#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "QuEST.h"
#include "QuEST_test_utils.hpp"

TEST_CASE( "benchmark compactUnitary", "[!benchmark]" ) {
    
    int numQb = 10;
    QuESTEnv env = createQuESTEnv();
    Qureg stat = createQureg(numQb, env);
    Qureg dens = createDensityQureg(numQb, env);
    
    qcomp a = .3 * exp(2i);
    qcomp b = sqrt(1-abs(a)*abs(a)) * exp(-3i);
    Complex alpha = toComplex( a );
    Complex beta = toComplex( b );

    BENCHMARK( "state vector", t) {
        compactUnitary(stat, t%numQb, alpha, beta);
    };
    BENCHMARK( "density matrix", t ) {
        compactUnitary(dens, t%numQb, alpha, beta);
    };

    destroyQureg(stat, env);
    destroyQureg(dens, env);
    destroyQuESTEnv(env);
}