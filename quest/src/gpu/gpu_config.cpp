/** @file
 * Utility functions for querying GPU hardware,
 * and allocating and copying GPU VRAM data.
 * 
 * @author Tyson Jones
 */

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/matrices.h"
#include "quest/include/channels.h"
#include "quest/include/environment.h"

#include "quest/src/core/errors.hpp"
#include "quest/src/core/memory.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"
#include "quest/src/gpu/gpu_config.hpp"

#include <set>
#include <array>
#include <string>
#include <algorithm>


#if COMPILE_CUDA && ! (defined(__NVCC__) || defined(__HIP__))
    #error \
        "Attempted to compile gpu_config.cpp in GPU-accelerated mode with a non-GPU compiler. "\
        "Please compile this file with a CUDA (nvcc) or ROCm (hipcc) compiler."
#endif


#if COMPILE_CUDA && defined(__NVCC__)
    #include <cuda.h>
    #include <cuda_runtime.h>
#endif
#if COMPILE_CUDA && defined(__HIP__)
    #include "quest/src/gpu/cuda_to_hip.hpp"
#endif



/*
 * CUDA ERROR HANDLING
 *
 * which is only defined when CUDA-compiling, since it is invoked only a macro (defined
 * in gpu_config.hpp) which wraps CUDA API calls
 */

#if COMPILE_CUDA

void assertCudaCallSucceeded(int result, const char* call, const char* caller, const char* file, int line) {

    // result (int) is actually type cudaError_t but we cannot use this CUDA-defined type
    // in gpu_config.hpp (since it's included by non-CUDA-compiled files), and we wish to keep
    // the signature consistent.
    cudaError_t code = (cudaError_t) result;

    if (result != cudaSuccess)
        error_cudaCallFailed(cudaGetErrorString(code), call, caller, file, line);
}

#endif



/*
 * CUQUANTUM MANAGEMENT
 *
 * these functions are defined in gpu_cuquantum.hpp when
 * COMPILE_CUQUANTUM is 1, but are otherwise defaulted to
 * the internal errors below. This slight inelegance
 * enables us to keep gpu_cuquantum.hpp as a single header
 * file, without exposing it to code beyond gpu/
 */


#if ! COMPILE_CUQUANTUM

void gpu_initCuQuantum() {
    error_cuQuantumInitOrFinalizedButNotCompiled();
}

void gpu_finalizeCuQuantum() {
    error_cuQuantumInitOrFinalizedButNotCompiled();
}

#endif



/*
 * HARDWARE AVAILABILITY
 */


// many of the below functions must be assessed
// per-device, when MPI + CUDA hybridisation means
// MPI processes are bound to unique GPUs. Such
// functions must not be called before we have
// explicitly bound GPUs to MPI processes, which
// we use this flag to defensively ensure
bool hasGpuBeenBound = false;


int getBoundGpuId() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    int id;
    CUDA_CHECK( cudaGetDevice(&id) );
    return id;
    
#else
    error_gpuQueriedButGpuNotCompiled();
    return -1;
#endif
}


int gpu_getComputeCapability() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    cudaDeviceProp props;
    CUDA_CHECK( cudaGetDeviceProperties(&props, getBoundGpuId()) );
    return props.major * 10 + props.minor;

#else
    error_gpuQueriedButGpuNotCompiled();
    return -1;
#endif
}


bool gpu_isGpuCompiled() {
    return (bool) COMPILE_CUDA;
}


bool gpu_isCuQuantumCompiled() {
    return (bool) COMPILE_CUQUANTUM;
}


int gpu_getNumberOfLocalGpus() {
#if COMPILE_CUDA

    // HIP throws an error when a CUDA API function
    // is called but no devices exist, which we handle
    int num;
    auto status = cudaGetDeviceCount(&num);
    return (status == cudaSuccess)? num : 0;

#else
    error_gpuQueriedButGpuNotCompiled();
    return -1;
#endif
}


bool gpu_isGpuAvailable() {
#if COMPILE_CUDA

    int numDevices = gpu_getNumberOfLocalGpus();
    if (numDevices == 0)
        return false;

    // check if any reported device is a valid GPU
    for (int deviceInd=0; deviceInd < numDevices; deviceInd++) {

        // by checking the properties of each device
        struct cudaDeviceProp props;
        auto status = cudaGetDeviceProperties(&props, deviceInd);

        // if the query failed, device is anyway unusable
        if (status != cudaSuccess) 
            continue;

        // if the device is a real GPU, it's 'major' compute capability is != 9999 (meaning emulation)
        if (props.major != 9999)
            return true;
    }

    // no non-emulation devices were found
    return false;

#else
    error_gpuQueriedButGpuNotCompiled();
    return false;
#endif
}


bool gpu_isDirectGpuCommPossible() {
#if COMPILE_CUDA

    if (!comm_isMpiGpuAware())
        return false;

    if (!gpu_isGpuAvailable())
        return false;

    /// @todo
    /// and are GPUs compatible?
    /// (the above might need to call a GPU-compiled func)

    return true;

#else
    error_gpuQueriedButGpuNotCompiled();
    return false;
#endif
}


size_t gpu_getCurrentAvailableMemoryInBytes() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    // note that in distributed settings, all GPUs
    // are being simultaneously queried, and it is
    // possible their values differ per-node

    size_t free, total;
    CUDA_CHECK( cudaMemGetInfo(&free, &total) );
    return free;

#else
    error_gpuQueriedButGpuNotCompiled();
    return 0;
#endif
}


size_t gpu_getTotalMemoryInBytes() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    size_t free, total;
    CUDA_CHECK( cudaMemGetInfo(&free, &total) );
    return total;

#else
    error_gpuQueriedButGpuNotCompiled();
    return 0;
#endif
}


bool gpu_doesGpuSupportMemPools() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    int supports;
    CUDA_CHECK( cudaDeviceGetAttribute(&supports, cudaDevAttrMemoryPoolsSupported, getBoundGpuId()) );
    return (bool) supports;

#else
    error_gpuQueriedButGpuNotCompiled();
    return false;
#endif
}


qindex gpu_getMaxNumConcurrentThreads() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    int deviceId = getBoundGpuId();

    // this may differ between nodes (which have different GPUs), which is fine
    int maxThreadsPerBlock;
    int maxNumBlocks;
    CUDA_CHECK( cudaDeviceGetAttribute(&maxThreadsPerBlock, cudaDevAttrMaxThreadsPerBlock,  deviceId) );
    CUDA_CHECK( cudaDeviceGetAttribute(&maxNumBlocks,       cudaDevAttrMultiProcessorCount, deviceId) );

    return maxThreadsPerBlock * static_cast<qindex>(maxNumBlocks); // avoid overflow

#else
    error_gpuQueriedButGpuNotCompiled();
    return -1;
#endif
}



/*
 * ENVIRONMENT MANAGEMENT
 */


std::array<char,17> getBoundGpuUuid() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    constexpr int numUuidChars = 16;
    constexpr int numOutChars = numUuidChars + 1; // terminal char

    std::array<char,numOutChars> out;
    cudaUUID_t uuid;

    // ROCm v5's cudaDeviceProp doesn't have a uuid field
    #if defined(__HIP__)
        hipDevice_t device;
        CUDA_CHECK( hipDeviceGet(&device, getBoundGpuId()) );
        CUDA_CHECK( hipDeviceGetUuid(&uuid, device) );
    #else
        cudaDeviceProp prop;
        CUDA_CHECK( cudaGetDeviceProperties(&prop, getBoundGpuId()) );
        uuid = prop.uuid;
    #endif

    // copy char[16] to out[17] (not human readable)
    std::copy_n(uuid.bytes, numUuidChars, out.begin());

    // include terminal char so that subsequent to-string
    // operations will succeed without knowing string length
    out[numOutChars-1] = '\0';
    return out;

#else
    error_gpuQueriedButGpuNotCompiled();
    return {};
#endif
}


void gpu_bindLocalGPUsToNodes() {
#if COMPILE_CUDA

    // distribute local MPI processes across local GPUs;
    int numLocalGpus = gpu_getNumberOfLocalGpus();
    int localGpuInd = comm_getRank() % numLocalGpus;
    CUDA_CHECK( cudaSetDevice(localGpuInd) );

    // note it is possible for multiple MPI processes
    // to bind to the same local GPU (as can be assessed
    // with gpu_doAnyMpiProcessesShareLocalGpu()), but
    // this will incur slowdowns due to context-switching
    // and has no benefit - is if further illegal when
    // using cuStateVec, as the caller will validate

    // indicate that other GPU-queries are now legal
    hasGpuBeenBound = true;

#else
    error_gpuQueriedButGpuNotCompiled();
#endif 
}


bool gpu_areAnyNodesBoundToSameGpu() {
#if COMPILE_CUDA
    assert_gpuHasBeenBound(hasGpuBeenBound);

    if (!comm_isInit())
        return false;

    // obtain bound GPU's UUID; a unique identifier 16-char identifier
    auto uuidStr = getBoundGpuUuid();

    // we can repurpose string-to-root sending to collect all uuids
    auto allUuids = comm_gatherStringsToRoot(uuidStr.data(), uuidStr.size());
    auto uniqueUuids = std::set<std::string>(allUuids.begin(), allUuids.end());

    // and assess whether they're all unique (non-root's bools are overwritten)
    bool localGpusAreUnique = allUuids.size() == uniqueUuids.size();
    bool globalGpusAreUnique = comm_isTrueOnRootNode(localGpusAreUnique);
    return ! globalGpusAreUnique;

#else
    error_gpuQueriedButGpuNotCompiled();
    return false;
#endif 
}


void gpu_sync() {
#if COMPILE_CUDA

    CUDA_CHECK( cudaDeviceSynchronize() );

#else
    error_gpuSyncedButGpuNotCompiled();
#endif
}



/*
 * MEMORY ALLOCATION
 */


qcomp* gpu_allocArray(qindex length) {
#if COMPILE_CUDA

    size_t numBytes = mem_getLocalQuregMemoryRequired(length);

    // attempt to malloc
    qcomp* ptr;
    cudaError_t errCode = cudaMalloc(&ptr, numBytes);

    // intercept memory-alloc error and merely return nullptr pointer (to be handled by validation)
    if (errCode == cudaErrorMemoryAllocation)
        return nullptr;

    // pass all other unexpected errors to internal error handling
    CUDA_CHECK(errCode);

    return ptr;

#else
    error_gpuAllocButGpuNotCompiled();
    return nullptr;
#endif
}


void gpu_deallocArray(qcomp* amps) {
#if COMPILE_CUDA

    // cudaFree on nullptr is fine
    CUDA_CHECK( cudaFree(amps) );

#else
    error_gpuDeallocButGpuNotCompiled();
#endif
}



/*
 * MEMORY MOVEMENT
 */


// flags to make the memory transfer direction visually distinct
enum CopyDirection {
    TO_HOST,
    TO_DEVICE
};


void copyArrayIfGpuCompiled(qcomp* cpuArr, qcomp* gpuArr, qindex numElems, enum CopyDirection direction) {
#if COMPILE_CUDA

    // must ensure gpu amps are up to date
    gpu_sync();

    auto flag = (direction == TO_HOST)? 
        cudaMemcpyDeviceToHost:
        cudaMemcpyHostToDevice;

    auto src = (direction == TO_HOST)? gpuArr : cpuArr;
    auto dst = (direction == TO_HOST)? cpuArr : gpuArr;

    // synchronous memory copy
    size_t numBytes = numElems * sizeof(qcomp);
    CUDA_CHECK( cudaMemcpy(dst, src, numBytes, flag) );

#else
    error_gpuCopyButGpuNotCompiled();
#endif
}


void copyMatrixIfGpuCompiled(qcomp** cpuMatr, qcomp* gpuArr, qindex matrDim, enum CopyDirection direction) {
#if COMPILE_CUDA

    // NOTE:
    // this function copies a 2D CPU matrix into a 1D row-major GPU array,
    // although this is not actually needed by the QuEST backend which
    // maintains 1D row-major CPU memories merely aliased by 2D structures
    // for the user's benefit. As such, this is dead code, but preserved in
    // case it is ever needed (like if custom user 2D data was needed in GPU).
    error_gpuDeadCopyMatrixFunctionCalled();

    // for completeness, we permit copying from the 1D GPU memory to the 2D CPU memory,
    // although we never actually have the need to do this!
    auto flag = (direction == TO_HOST)? 
        cudaMemcpyDeviceToHost:
        cudaMemcpyHostToDevice;

    // copy each CPU row into flattened GPU memory. we make each memcpy asynch,
    // but it's unclear it helps, nor whether single-stream sync is necessary
    size_t numBytesPerRow = matrDim * sizeof(**cpuMatr);
    gpu_sync();

    for (qindex r=0; r<matrDim; r++) {
        qcomp* cpuRow = cpuMatr[r];
        qcomp* gpuSlice = &gpuArr[r*matrDim];

        auto src = (direction == TO_HOST)? gpuSlice : cpuRow;
        auto dst = (direction == TO_HOST)? cpuRow   : gpuSlice;

        CUDA_CHECK( cudaMemcpyAsync(dst, src, numBytesPerRow, flag) );
    }

    // wait for async copies to complete
    gpu_sync();

#else
    error_gpuCopyButGpuNotCompiled();
#endif
}


template <typename T>
void assertHeapObjectGpuMemIsAllocated(T obj) {

    if (! mem_isAllocated(util_getGpuMemPtr(obj)) || ! getQuESTEnv().isGpuAccelerated)
        error_gpuCopyButMatrixNotGpuAccelerated();
}


void gpu_copyArray(qcomp* dest, qcomp* src, qindex dim) {
#if COMPILE_CUDA

    // ensure src and dest aren't being modified
    gpu_sync();

    CUDA_CHECK( cudaMemcpy(dest, src, dim * sizeof(qcomp), cudaMemcpyDeviceToDevice) );

#else
    error_gpuCopyButGpuNotCompiled();
#endif
}


void gpu_copyCpuToGpu(qcomp* cpuArr, qcomp* gpuArr, qindex numElems) {
    copyArrayIfGpuCompiled(cpuArr, gpuArr, numElems, TO_DEVICE);
}
void gpu_copyGpuToCpu(qcomp* gpuArr, qcomp* cpuArr, qindex numElems) {
    copyArrayIfGpuCompiled(cpuArr, gpuArr, numElems, TO_HOST);
}


void gpu_copyCpuToGpu(Qureg qureg, qcomp* cpuArr, qcomp* gpuArr, qindex numElems) {
    // used for moving memory within the same Qureg, hence the Qureg arg only for assertion
    assert_quregIsGpuAccelerated(qureg);
    copyArrayIfGpuCompiled(cpuArr, gpuArr, numElems, TO_DEVICE);
}
void gpu_copyGpuToCpu(Qureg qureg, qcomp* gpuArr, qcomp* cpuArr, qindex numElems) {
    assert_quregIsGpuAccelerated(qureg);
    copyArrayIfGpuCompiled(cpuArr, gpuArr, numElems, TO_HOST);
}


void gpu_copyCpuToGpu(Qureg qureg) {
    gpu_copyCpuToGpu(qureg, qureg.cpuAmps, qureg.gpuAmps, qureg.numAmpsPerNode);
}
void gpu_copyGpuToCpu(Qureg qureg) {
    gpu_copyGpuToCpu(qureg, qureg.gpuAmps, qureg.cpuAmps, qureg.numAmpsPerNode);
}


void gpu_copyCpuToGpu(CompMatr matr) {
    assertHeapObjectGpuMemIsAllocated(matr);

    // note matr.cpuElems is merely a 2D alias for matr.cpuElemsFlat, which
    // matches the format of matr.gpuElemsFlat. Ergo, we do not invoke 
    // copyMatrixIfGpuCompiled(), and instead more efficiently overwrite
    // the contiguous memory, which retains any user changes to .cpuElems

    qindex numElems = matr.numRows * matr.numRows;
    copyArrayIfGpuCompiled(matr.cpuElemsFlat, util_getGpuMemPtr(matr), numElems, TO_DEVICE);
}
void gpu_copyGpuToCpu(CompMatr matr) {
    assertHeapObjectGpuMemIsAllocated(matr);

    // note matr.cpuElems is merely a 2D alias for matr.cpuElemsFlat, which
    // matches the format of matr.gpuElemsFlat. Ergo, we do not invoke 
    // copyMatrixIfGpuCompiled(), and instead more efficiently overwrite
    // the contiguous matr.cpuElemsFlat, which users can access via .cpuElems

    qindex numElems = matr.numRows * matr.numRows;
    copyArrayIfGpuCompiled(matr.cpuElemsFlat, util_getGpuMemPtr(matr), numElems, TO_HOST);
}


void gpu_copyCpuToGpu(DiagMatr matr) {
    assertHeapObjectGpuMemIsAllocated(matr);
    copyArrayIfGpuCompiled(matr.cpuElems, util_getGpuMemPtr(matr), matr.numElems, TO_DEVICE);
}
void gpu_copyGpuToCpu(DiagMatr matr) {
    assertHeapObjectGpuMemIsAllocated(matr);
    copyArrayIfGpuCompiled(matr.cpuElems, util_getGpuMemPtr(matr), matr.numElems, TO_HOST);
}


void gpu_copyCpuToGpu(SuperOp op) {
    assertHeapObjectGpuMemIsAllocated(op);

    // note op.cpuElems is merely a 2D alias for op.cpuElemsFlat, which
    // matches the format of op.gpuElemsFlat. Ergo, we do not invoke 
    // copyMatrixIfGpuCompiled(), and instead more efficiently overwrite
    // the contiguous memory, which retains any user changes to .cpuElems

    qindex numElems = op.numRows * op.numRows;
    copyArrayIfGpuCompiled(op.cpuElemsFlat, util_getGpuMemPtr(op), numElems, TO_DEVICE);
}
void gpu_copyGpuToCpu(SuperOp op) {
    assertHeapObjectGpuMemIsAllocated(op);

    // note op.cpuElems is merely a 2D alias for op.cpuElemsFlat, which
    // matches the format of op.gpuElemsFlat. Ergo, we do not invoke 
    // copyMatrixIfGpuCompiled(), and instead more efficiently overwrite
    // the contiguous op.cpuElemsFlat, which users can access via .cpuElems

    qindex numElems = op.numRows * op.numRows;
    copyArrayIfGpuCompiled(op.cpuElemsFlat, util_getGpuMemPtr(op), numElems, TO_HOST);
}


void gpu_copyCpuToGpu(FullStateDiagMatr matr) {
    assertHeapObjectGpuMemIsAllocated(matr);
    copyArrayIfGpuCompiled(matr.cpuElems, util_getGpuMemPtr(matr), matr.numElemsPerNode, TO_DEVICE);
}



/*
 * CACHE MANAGEMENT
 */


// persistent but variably-sized cache memory used by the any-targ dense
// matrix kernel as working memory, which is lazily runtime expanded when
// necessary, and only ever cleared when triggered by the user
qcomp* gpuCache = nullptr;
qindex gpuCacheLen = 0;


qcomp* gpu_getCacheOfSize(qindex numElemsPerThread, qindex numThreads) {
#if COMPILE_CUDA

    // do not interfere with existing kernels using the cache
    gpu_sync();

    qindex numNewElems = numElemsPerThread * numThreads;

    // return existing cache if it's already sufficiently big
    if (numNewElems <= gpuCacheLen)
        return gpuCache;

    // otherwise, resize the cache
    gpuCacheLen = numNewElems;
    CUDA_CHECK( cudaFree(gpuCache) ); // nullptr fine to free
    CUDA_CHECK( cudaMalloc(&gpuCache, gpuCacheLen * sizeof *gpuCache) );

    return gpuCache;

#else
    error_gpuCacheModifiedButGpuNotCompiled();
    return nullptr;
#endif
}


void gpu_clearCache() {
#if COMPILE_CUDA

    // do not interfere with existing kernels using the cache
    gpu_sync();

    // cudaFree on nullptr is fine
    CUDA_CHECK( cudaFree(gpuCache) );

    gpuCache = nullptr;
    gpuCacheLen = 0;

#else
    error_gpuCacheModifiedButGpuNotCompiled();
#endif
}


size_t gpu_getCacheMemoryInBytes() {

    // query permitted even when not GPU accelerated
    return gpuCacheLen * sizeof *gpuCache;
}
