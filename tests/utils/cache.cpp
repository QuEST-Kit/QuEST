/** @file
 * Testing utilities which create Quregs across all
 * available hardware deployments
 *
 * @author Tyson Jones
 */

#include "quest.h"

#include "macros.hpp"
#include "qvector.hpp"
#include "qmatrix.hpp"
#include "macros.hpp"
#include "linalg.hpp"
#include "config.hpp"
#include "cache.hpp"

#include <unordered_map>
#include <tuple>
#include <vector>
#include <string>

using std::tuple;
using std::vector;
using std::string;



/*
 * caches of Quregs and FullStateDiagMatr
 * which persist between calls of get-cache()
 */

quregCache statevecs1;
quregCache statevecs2;
quregCache densmatrs1;
quregCache densmatrs2;
matrixCache matrices;

int getNumCachedQubits() {

    // we are merely aliasing the below env-var fetching function
    // to minimise a diff since pre-runtime controlling the tested
    // Qureg sizes is experimental and not fully designed (we may
    // eventually wish to specify different sizes for statevectors
    // vs density matrices, or control integration test Qureg sizes
    // also through environment variables, etc)
    return getNumQubitsInUnitTestedQuregs();
}



/*
 * deployments in cache
 */

deployInfo getSupportedDeployments() {

    deployInfo out;

    // determine which Qureg deployments are supported
    QuESTEnv env = getQuESTEnv();
    bool omp = env.isMultithreaded;
    bool mpi = env.isDistributed;
    bool gpu = env.isGpuAccelerated;

    // return only the "most-accelerated" deployment, unless all are desired
    bool one = ! getWhetherToTestAllDeployments();

    // add only those supported to the output list, in order of preference.
    // flag order is (MPI, GPU, OMP), matching createCustomQureg
    if (gpu && omp && mpi) { out.push_back({"GPU + OMP + MPI", 1, 1, 1}); if (one) return out; }
    if (gpu && mpi)        { out.push_back({"GPU + MPI",       1, 1, 0}); if (one) return out; }
    if (gpu && omp)        { out.push_back({"GPU + OMP",       0, 1, 1}); if (one) return out; }
    if (gpu)               { out.push_back({"GPU",             0, 1, 0}); if (one) return out; }
    if (mpi && omp)        { out.push_back({"CPU + OMP + MPI", 1, 0, 1}); if (one) return out; }
    if (mpi)               { out.push_back({"CPU + MPI",       1, 0, 0}); if (one) return out; }
    if (omp)               { out.push_back({"CPU + OMP",       0, 0, 1}); if (one) return out; }
    if (true)              { out.push_back({"CPU",             0, 0, 0}); if (one) return out; }

    // always contains CPU obviously, but this makes it explicit
    DEMAND( !out.empty() );
    
    // return all supported deployments
    return out;
}



/*
 * manage cached quregs
 */

quregCache createCachedStatevecsOrDensmatrs(bool isDensMatr) {

    quregCache out;

    // only add supported-deployment quregs to the cache
    for (auto [label, mpi, gpu, omp] : getSupportedDeployments())
        out[label] = createCustomQureg(getNumCachedQubits(), isDensMatr, mpi, gpu, omp);

    return out;
}

void createCachedQuregs() {

    // must not be called twice nor pre-creation
    DEMAND( statevecs1.empty() );
    DEMAND( statevecs2.empty() );
    DEMAND( densmatrs1.empty() );
    DEMAND( densmatrs2.empty() );

    statevecs1 = createCachedStatevecsOrDensmatrs(false);
    statevecs2 = createCachedStatevecsOrDensmatrs(false);
    densmatrs1 = createCachedStatevecsOrDensmatrs(true);
    densmatrs2 = createCachedStatevecsOrDensmatrs(true);
}

void destroyCachedQuregs() {

    // must not be called twice nor pre-creation
    DEMAND( ! statevecs1.empty() );
    DEMAND( ! statevecs2.empty() );
    DEMAND( ! densmatrs1.empty() );
    DEMAND( ! densmatrs2.empty() );

    auto caches = {
        statevecs1, statevecs2, 
        densmatrs1, densmatrs2};

    for (auto& cache : caches)
        for (auto& [label, qureg]: cache)
            destroyQureg(qureg);

    statevecs1.clear();
    statevecs2.clear();
    densmatrs1.clear();
    densmatrs2.clear();
}

quregCache getCachedStatevecs() {

    // must not be called pre-creation nor post-destruction
    DEMAND( !statevecs1.empty() );

    return statevecs1;
}
quregCache getCachedDensmatrs() {

    // must not be called pre-creation nor post-destruction
    DEMAND( !densmatrs1.empty() );

    return densmatrs1;
}

quregCache getAltCachedStatevecs() {

    // must not be called pre-creation nor post-destruction
    DEMAND( !statevecs2.empty() );

    return statevecs2;
}
quregCache getAltCachedDensmatrs() {

    // must not be called pre-creation nor post-destruction
    DEMAND( !densmatrs2.empty() );

    return densmatrs2;
}



/*
 * manage cached FullStateDiagMatr
 */

void createCachedFullStateDiagMatrs() {

    // must not be called twice
    DEMAND( matrices.empty() );

    // only add supported-deployment matrices to the cache
    for (auto [label, mpi, gpu, omp] : getSupportedDeployments())
        matrices[label] = createCustomFullStateDiagMatr(getNumCachedQubits(), mpi, gpu, omp);
}

void destroyCachedFullStateDiagMatrs() {

    // must not be called twice
    DEMAND( !matrices.empty() );

    for (auto& [label, matrix]: matrices)
        destroyFullStateDiagMatr(matrix);

    matrices.clear();
}

matrixCache getCachedFullStateDiagMatrs() {

    // must not be called pre-creation nor post-destruction
    DEMAND( !matrices.empty() );

    return matrices;
}


/*
 * reference states of equivalent 
 * dimension to the cached quregs
 */

qvector getRefStatevec() {
    return getZeroVector(getPow2(getNumCachedQubits()));
}
qmatrix getRefDensmatr() {
    return getZeroMatrix(getPow2(getNumCachedQubits()));
}
