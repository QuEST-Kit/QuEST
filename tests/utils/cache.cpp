/** @file
 * Testing utilities which create Quregs across all
 * available hardware deployments
 *
 * @author Tyson Jones
 */

#include "quest/include/quest.h"

#include "macros.hpp"
#include "qvector.hpp"
#include "qmatrix.hpp"
#include "linalg.hpp"
#include "cache.hpp"

#include <unordered_map>
#include <tuple>
#include <vector>
#include <string>

using std::tuple;
using std::vector;
using std::string;



/*
 * cached quregs which persist between
 * calls of getCachedStatevecs() and
 * getCachedDensmatrs()
 */

quregCache statevecs1;
quregCache statevecs2;
quregCache densmatrs1;
quregCache densmatrs2;



/*
 * while the number of qubits in the unit-test Quregs
 * is fixed, it is defined privately here (with internal
 * linkage) so that it can be changed between compilations
 * without having to recompiling the entire test suite
 */

static constexpr int NUM_QUBITS_IN_QUREGS = 6;

int getNumCachedQubits() {
    return NUM_QUBITS_IN_QUREGS;
}



/*
 * manage cached quregs
 */

deployInfo getSupportedDeployments() {

    deployInfo out;

    // determine which Qureg deployments are supported
    QuESTEnv env = getQuESTEnv();
    bool omp = env.isMultithreaded;
    bool mpi = env.isDistributed;
    bool gpu = env.isGpuAccelerated;

    // return only the "most-accelerated" deployment, unless all are desired
    bool one = ! TEST_ALL_DEPLOYMENTS;

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
    
    // return all supported deployments
    return out;
}

quregCache createCachedStatevecsOrDensmatrs(bool isDensMatr) {

    quregCache out;

    // only add supported-deployment quregs to the cache
    for (auto [label, mpi, gpu, omp] : getSupportedDeployments())
        out[label] = createCustomQureg(NUM_QUBITS_IN_QUREGS, isDensMatr, mpi, gpu, omp);

    return out;
}

void createCachedQuregs() {
    statevecs1 = createCachedStatevecsOrDensmatrs(false);
    statevecs2 = createCachedStatevecsOrDensmatrs(false);
    densmatrs1 = createCachedStatevecsOrDensmatrs(true);
    densmatrs2 = createCachedStatevecsOrDensmatrs(true);
}

void destroyCachedQuregs() {

    auto caches = {
        statevecs1, statevecs2, 
        densmatrs1, densmatrs2};

    for (auto& cache : caches)
        for (auto& [label, qureg]: cache)
            destroyQureg(qureg);
}

quregCache getCachedStatevecs() {
    return statevecs1;
}
quregCache getCachedDensmatrs() {
    return densmatrs1;
}

quregCache getAltCachedStatevecs() {
    return statevecs2;
}
quregCache getAltCachedDensmatrs() {
    return densmatrs2;
}



/*
 * reference states of equivalent 
 * dimension to the cached quregs
 */

qvector getRefStatevec() {
    return getZeroVector(getPow2(NUM_QUBITS_IN_QUREGS));
}
qmatrix getRefDensmatr() {
    return getZeroMatrix(getPow2(NUM_QUBITS_IN_QUREGS));
}
