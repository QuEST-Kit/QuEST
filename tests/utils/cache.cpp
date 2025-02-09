#include "quest.h"
#include "macros.hpp"
#include "cache.hpp"

#include <unordered_map>

quregCache statevecs;
quregCache densmatrs;


quregCache createCachedStatevecsOrDensmatrs(bool isDensMatr) {

    // determine which Qureg deployments are supported
    QuESTEnv env = getQuESTEnv();
    bool omp = env.isMultithreaded;
    bool mpi = env.isDistributed;
    bool gpu = env.isGpuAccelerated;

    // only add supported-deployment quregs to the cache
    quregCache cache;

    // signature is createCustomQureg(..., MPI, GPU, OMP)
    if (true)              cache["CPU"]             = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 0, 0, 0);
    if (mpi)               cache["CPU + MPI"]       = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 1, 0, 0);
    if (omp)               cache["CPU + MPI"]       = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 0, 0, 1);
    if (mpi && omp)        cache["CPU + OMP + MPI"] = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 1, 0, 1);
    if (gpu)               cache["GPU"]             = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 0, 1, 0);
    if (gpu && omp)        cache["GPU + OMP"]       = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 0, 1, 1);
    if (gpu && mpi)        cache["GPU + OMP + MPI"] = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 1, 1, 0);
    if (gpu && omp && mpi) cache["GPU + OMP + MPI"] = createCustomQureg(NUM_QUREG_QUBITS, isDensMatr, 1, 1, 1);

    return cache;
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
