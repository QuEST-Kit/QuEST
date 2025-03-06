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

quregCache statevecs;
quregCache densmatrs;

using std::tuple;
using std::vector;
using std::string;



deployInfo getSupportedDeployments() {

    deployInfo out;

    // determine which Qureg deployments are supported
    QuESTEnv env = getQuESTEnv();
    bool omp = env.isMultithreaded;
    bool mpi = env.isDistributed;
    bool gpu = env.isGpuAccelerated;

    // add only those supported to the output list
    // order is (MPI, GPU, OMP), matching createCustomQureg
    if (true)              out.push_back({"CPU",             0, 0, 0});
    if (mpi)               out.push_back({"CPU + MPI",       1, 0, 0});
    if (omp)               out.push_back({"CPU + OMP",       0, 0, 1});
    if (mpi && omp)        out.push_back({"CPU + OMP + MPI", 1, 0, 1});
    if (gpu)               out.push_back({"GPU",             0, 1, 0});
    if (gpu && omp)        out.push_back({"GPU + OMP",       0, 1, 1});
    if (gpu && mpi)        out.push_back({"GPU + MPI",       1, 1, 0});
    if (gpu && omp && mpi) out.push_back({"GPU + OMP + MPI", 1, 1, 1});

    return out;
}


quregCache createCachedStatevecsOrDensmatrs(bool isDensMatr) {

    quregCache out;

    // only add supported-deployment quregs to the cache
    for (auto [label, mpi, gpu, omp] : getSupportedDeployments())
        out[label] = createCustomQureg(NUM_UNIT_QUREG_QUBITS, isDensMatr, mpi, gpu, omp);

    return out;
}


void createCachedQuregs() {
    statevecs = createCachedStatevecsOrDensmatrs(false);
    densmatrs = createCachedStatevecsOrDensmatrs(true);
}


void destroyCachedQuregs() {

    for (auto& [label, qureg]: statevecs)
        destroyQureg(qureg);

    for (auto& [label, qureg]: densmatrs)
        destroyQureg(qureg);
}


quregCache getCachedStatevecs() {
    return statevecs;
}
quregCache getCachedDensmatrs() {
    return densmatrs;
}


qvector getRefStatevec() {
    return getZeroVector(getPow2(NUM_UNIT_QUREG_QUBITS));
}
qmatrix getRefDensmatr() {
    return getZeroMatrix(getPow2(NUM_UNIT_QUREG_QUBITS));
}
