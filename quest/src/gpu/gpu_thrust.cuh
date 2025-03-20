/** @file
 * Subroutines which invoke Thrust. This file is only ever included
 * when COMPILE_CUDA=1 so it can safely invoke CUDA signatures without 
 * guards. Further, as it is entirely a header, it can declare templated
 * times without explicitly instantiating them across all parameter values.
 * 
 * Where possible, we use Thrust vector iterators in lieu of explicitly
 * iterating indices (and consulting a vector bound to a Thrust functor),
 * since the latter may lead to sub-optimal memory access (e.g. non-coalesced).
 * Alas, some functors are too complicated (e.g. they access multiple non-contiguous
 * amps per natural iteration) and fall-back to binding the GPU memory.
 * 
 * This file is a (CUDA) header since only ever included by gpu_subroutines.cpp.
 * 
 * When compiling for AMD GPUs, the Thrust and CUDA symbols invoked herein are
 * mapped to RocmThrust and HIP symbols by cuda_to_hip.h.
 * 
 * @author Tyson Jones
 */

#ifndef GPU_THRUST_HPP
#define GPU_THRUST_HPP

#if ! COMPILE_CUDA
    #error "A file being compiled somehow included gpu_thrust.hpp despite QuEST not being compiled in GPU-accelerated mode."
#endif

#include "quest/include/modes.h"
#include "quest/include/types.h"
#include "quest/include/qureg.h"
#include "quest/include/paulis.h"
#include "quest/include/matrices.h"

#include "quest/src/gpu/gpu_types.cuh"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/constants.hpp"
#include "quest/src/core/utilities.hpp"
#include "quest/src/core/randomiser.hpp"
#include "quest/src/comm/comm_config.hpp"

// kernels/thrust must use cu_qcomp, never qcomp
#define USE_CU_QCOMP
#include "quest/src/core/fastmath.hpp"
#undef USE_CU_QCOMP

#include <thrust/random.h>
#include <thrust/complex.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/for_each.h>
#include <thrust/inner_product.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/system/system_error.h>



/*
 * LISTS
 *
 * are copied to device memory using thrust's device_vector's 
 * copy constructor (devicevec d_vec = hostvec). The pointer 
 * to the data (d_vec.data()) can be cast into a raw pointer
 * and passed directly to CUDA kernels
 */


using devints = thrust::device_vector<int>;

int* getPtr(devints& qubits) {

    return thrust::raw_pointer_cast(qubits.data());
}


using devreals = thrust::device_vector<qreal>;

qreal* getPtr(devreals& reals) {

    return thrust::raw_pointer_cast(reals.data());
}

void copyFromDeviceVec(devreals& reals, qreal* out) {

    thrust::copy(reals.begin(), reals.end(), out);
}

devreals getDeviceRealsVec(qindex dim) {

    devreals out;

    try  {
        out.resize(dim);
        thrust::fill(out.begin(), out.end(), 0.);

    } catch (thrust::system_error &e) {
        error_thrustTempGpuAllocFailed();
    }

    return out;
}



/*
 * AMP POINTERS
 *
 * used to enumerate GPU amps of matrices, quregs and
 * full-state diagonal matrices inside thrust functions
 */


thrust::device_ptr<cu_qcomp> getStartPtr(cu_qcomp* amps) {

    return thrust::device_pointer_cast(amps);
}
auto getStartPtr(qcomp* amps) {

    return getStartPtr(toCuQcomps(amps));
}


auto getStartPtr(Qureg qureg) {

    return getStartPtr(qureg.gpuAmps);
}
auto getEndPtr(Qureg qureg) {

    return getStartPtr(qureg) + qureg.numAmpsPerNode;
}


auto getStartPtr(FullStateDiagMatr matr) {

    return getStartPtr(matr.gpuElems);
}
auto getEndPtr(FullStateDiagMatr matr) {

    return getStartPtr(matr) + matr.numElemsPerNode;
}



/*
 * CUSTOM FUNCTORS
 *
 * used to effect custom transformations upon GPU
 * amps using thrust functions
 */


struct functor_getAmpConj : public thrust::unary_function<cu_qcomp,cu_qcomp> {

    __host__ __device__ cu_qcomp operator()(cu_qcomp amp) {
        return getCompConj(amp);
    }
};

struct functor_getAmpNorm : public thrust::unary_function<cu_qcomp,qreal> {

    __host__ __device__ qreal operator()(cu_qcomp amp) {
        return getCompNorm(amp);
    }
};

struct functor_getAmpReal : public thrust::unary_function<cu_qcomp,qreal> {

    __host__ __device__ qreal operator()(cu_qcomp amp) {
        return getCompReal(amp);
    }
};


struct functor_getAmpConjProd : public thrust::binary_function<cu_qcomp,cu_qcomp,cu_qcomp> {

    __host__ __device__ cu_qcomp operator()(cu_qcomp braAmp, cu_qcomp ketAmp) { 
        return getCompConj(braAmp) * ketAmp;
    }
};

struct functor_getNormOfAmpDif : public thrust::binary_function<cu_qcomp,cu_qcomp,qreal> {

    __host__ __device__ qreal operator()(cu_qcomp amp1, cu_qcomp amp2) { 
        return getCompNorm(amp1 - amp2);
    }
};


struct functor_getExpecStateVecZTerm : public thrust::binary_function<qindex,cu_qcomp,qreal> {

    // this functor computes a single term from the sum
    // in the expectation value of Z of a statevector

    qindex targMask;
    functor_getExpecStateVecZTerm(qindex mask) : targMask(mask) {}

    __device__ qreal operator()(qindex ind, cu_qcomp amp) {
        
        int par = cudaGetBitMaskParity(ind & targMask); // device-only
        int sign = fast_getPlusOrMinusOne(par);
        return sign * getCompNorm(amp);
    }
};


struct functor_getExpecDensMatrZTerm : public thrust::unary_function<qindex,cu_qcomp> {

    // this functor computes a single term from the sum
    // in the expectation value of Z of a density matrix

    qindex numAmpsPerCol, firstDiagInd, targMask;
    cu_qcomp* amps;

    functor_getExpecDensMatrZTerm(qindex dim, qindex diagInd, qindex mask, cu_qcomp* _amps) : 
        numAmpsPerCol(dim), firstDiagInd(diagInd), targMask(mask), amps(_amps) {}

    __device__ cu_qcomp operator()(qindex n) {

        qindex i = fast_getQuregLocalIndexOfDiagonalAmp(n, firstDiagInd, numAmpsPerCol);
        qindex r = n + firstDiagInd;

        int par = cudaGetBitMaskParity(r & targMask); // device-only
        int sign = fast_getPlusOrMinusOne(par);
        return sign * amps[i];
    }
};


struct functor_getExpecStateVecPauliTerm : public thrust::unary_function<qindex,cu_qcomp> {

    // this functor computes a single term from the sum in the
    // expectation value of a Pauli str (which necessarily contains 
    // at least one X or Y) of a statevector

    qindex maskXY, maskYZ;
    cu_qcomp *amps, *pairAmps;

    functor_getExpecStateVecPauliTerm(qindex _maskXY, qindex _maskYZ, cu_qcomp* _amps, cu_qcomp* _pairAmps) : 
        maskXY(_maskXY), maskYZ(_maskYZ), amps(_amps), pairAmps(_pairAmps) {}

    __device__ cu_qcomp operator()(qindex n) {

        qindex j = flipBits(n, maskXY);
        int par = cudaGetBitMaskParity(j & maskYZ); // device-only
        int sign = fast_getPlusOrMinusOne(par);

        // sign excludes i^numY contribution
        return sign * getCompConj(amps[n]) * pairAmps[j]; // pairAmps may be amps or buffer
    }
};


struct functor_getExpecDensMatrPauliTerm : public thrust::unary_function<qindex,cu_qcomp> {

    // this functor computes a single term from the sum in the
    // expectation value of a Pauli str (which necessarily contains 
    // at least one X or Y) of a density matrix

    qindex maskXY, maskYZ;
    qindex numAmpsPerCol, firstDiagInd;
    cu_qcomp *amps;

    functor_getExpecDensMatrPauliTerm(qindex _maskXY, qindex _maskYZ, qindex _numAmpsPerCol, qindex _firstDiagInd, cu_qcomp* _amps) :
        maskXY(_maskXY), maskYZ(_maskYZ), numAmpsPerCol(_numAmpsPerCol), firstDiagInd(_firstDiagInd), amps(_amps) {}

    __device__ cu_qcomp operator()(qindex n) {

        qindex r = n + firstDiagInd;
        qindex i = flipBits(r, maskXY);
        qindex m = fast_getQuregLocalFlatIndex(i, n, numAmpsPerCol);

        // sign excludes i^numY contribution
        int par = cudaGetBitMaskParity(i & maskYZ); // device-only
        int sign = fast_getPlusOrMinusOne(par);
        return sign * amps[m];
    }
};


template <bool HasPower, bool UseRealPow>
struct functor_getExpecDensMatrDiagMatrTerm : public thrust::unary_function<qindex,cu_qcomp> {

    // this functor computes a single term from the sum in the expectation 
    // value of a FullStateDiagMatr upon a density matrix

    qindex numAmpsPerCol, firstDiagInd;
    cu_qcomp *amps, *elems, expo;

    functor_getExpecDensMatrDiagMatrTerm(qindex dim, qindex diagInd, cu_qcomp* _amps, cu_qcomp* _elems, cu_qcomp _expo) : 
        numAmpsPerCol(dim), firstDiagInd(diagInd), amps(_amps), elems(_elems), expo(_expo) {}

    __device__ cu_qcomp operator()(qindex n) {

        cu_qcomp elem = elems[n];

        if constexpr (HasPower && ! UseRealPow)
            elem = getCompPower(elem, expo);
        if constexpr (HasPower &&   UseRealPow)
            elem = getCuQcomp(pow(getCompReal(elem), getCompReal(expo)),0); // CUDA pow(qreal,qreal)

        qindex i = fast_getQuregLocalIndexOfDiagonalAmp(n, firstDiagInd, numAmpsPerCol);

        return amps[i] * elem;
    }
};


template <bool IsDiag>
struct functor_setAmpToPauliStrSumElem {

    int rank;
    qindex dim;
    qindex suffixLen;
    qindex numTerms;

    cu_qcomp* amps;
    cu_qcomp* coeffs;
    PauliStr* strings;
    
    functor_setAmpToPauliStrSumElem(
        int rank, qindex dim, qindex suffixLen, qindex numTerms, 
        cu_qcomp* amps, cu_qcomp* coeffs, PauliStr* strings
    ) :
        rank(rank), dim(dim), suffixLen(suffixLen), numTerms(numTerms), 
        amps(amps), coeffs(coeffs), strings(strings)
    {}

    __device__ void operator()(qindex n) {

        qindex i = concatenateBits(rank, n, suffixLen);

        // repurpose this functor for populating both density matrices
        // and full-state diagonal operators (for which strings=I,Z)
        qindex r, c;

        if constexpr (IsDiag) {
            r = i;
            c = i;
        } else {
            r = fast_getQuregGlobalRowFromFlatIndex(i, dim);
            c = fast_getQuregGlobalColFromFlatIndex(i, dim);
        }

        amps[n] = fast_getPauliStrSumElem(coeffs, strings, numTerms, r, c);
    }
};


struct functor_mixAmps : public thrust::binary_function<cu_qcomp,cu_qcomp,cu_qcomp> {

    // this functor linearly combines the given pair
    // of amplitudes, weighted by the fixed qreals,
    // and is used by mixQureg upon density matrices

    qreal outProb, inProb;
    functor_mixAmps(qreal out, qreal in) : outProb(out), inProb(in) {}

    __host__ __device__ cu_qcomp operator()(cu_qcomp outAmp, cu_qcomp inAmp) {
        
        return (outProb * outAmp) + (inProb * inAmp);
    }
};


struct functor_superposeAmps {

    // this functor linearly combines the given trio
    // of amplitudes, weighted by fixed qcomps, and is
    // used by setQuregToSuperposition

    cu_qcomp fac0, fac1, fac2;
    functor_superposeAmps(cu_qcomp f0, cu_qcomp f1, cu_qcomp f2) : fac0(f0), fac1(f1), fac2(f2) {}

    template <typename Tuple> __host__ __device__ void operator()(Tuple t) {
        thrust::get<0>(t) = fac0*thrust::get<0>(t) + fac1*thrust::get<1>(t) + fac2*thrust::get<2>(t);
    }
};


template <bool HasPower, bool UseRealPow, bool Norm>
struct functor_multiplyElemPowerWithAmpOrNorm : public thrust::binary_function<cu_qcomp,cu_qcomp,cu_qcomp> {

    // this functor multiplies a diagonal matrix element 
    // raised to a power (templated to optimise away the 
    // exponentiation at compile-time when power==1) upon
    // a statevector amp (used when modifying the state)
    // or its norm (used when calculating expected values)

    cu_qcomp exponent;
    functor_multiplyElemPowerWithAmpOrNorm(cu_qcomp power) : exponent(power) {}

    __host__ __device__ cu_qcomp operator()(cu_qcomp quregAmp, cu_qcomp matrElem) {

        if constexpr (HasPower && ! UseRealPow)
            matrElem = getCompPower(matrElem, exponent);
        if constexpr (HasPower &&   UseRealPow)
            matrElem = getCuQcomp(pow(getCompReal(matrElem), getCompReal(exponent)),0); // CUDA pow(qreal,qreal)

        if constexpr (Norm)
            quregAmp = getCuQcomp(getCompNorm(quregAmp), 0);

        return matrElem * quregAmp;
    }
};


struct functor_getDiagInd : public thrust::unary_function<qindex,qindex> {

    // this functor accepts the index of a statevector 
    // basis-state and produces the index of a density 
    // matrix's corresponding diagonal basis-state

    qindex numAmpsPerCol, firstDiagInd;

    functor_getDiagInd(Qureg qureg) {
        firstDiagInd = util_getLocalIndexOfFirstDiagonalAmp(qureg);
        numAmpsPerCol = powerOf2(qureg.numQubits);
    }

    __host__ __device__ qindex operator()(qindex i) {
        return fast_getQuregLocalIndexOfDiagonalAmp(i, firstDiagInd, numAmpsPerCol);
    }
};


template <int NumBits>
struct functor_insertBits : public thrust::unary_function<qindex,qindex> {

    // this functor inserts bits into a qindex value, and
    // is used to enumerate specific basis-state indices
    // with qubits in the specified bit values

    int* sortedIndsPtr;
    qindex valueMask;
    int numBits;

    functor_insertBits(int* ptr, qindex mask, int nBits) :
        sortedIndsPtr(ptr), valueMask(mask), numBits(nBits)
    {
        assert_numTargsMatchesTemplateParam(nBits, NumBits);
    }

    __host__ __device__ qindex operator()(qindex i) {

        // use the compile-time value if possible, to auto-unroll the insertBits loop
        SET_VAR_AT_COMPILE_TIME(int, nbits, NumBits, numBits);

        // return ith local index where bits have the specified values at the specified indices
        return insertBitsWithMaskedValues(i, sortedIndsPtr, nbits, valueMask);
    }
};


template <bool Conj>
struct functor_getFidelityTerm : public thrust::unary_function<qindex,cu_qcomp>
{
    int rank, numQubits;
    qindex logNumAmpsPerNode, numAmpsPerCol;
    cu_qcomp *rho, *psi;

    functor_getFidelityTerm(
        int _rank, int _numQubits, qindex _logNumAmpsPerNode, qindex _numAmpsPerCol, cu_qcomp* _rho, cu_qcomp* _psi
    ) :
        rank(_rank), numQubits(_numQubits), logNumAmpsPerNode(_logNumAmpsPerNode), numAmpsPerCol(_numAmpsPerCol), rho(_rho), psi(_psi)
    {}

    __host__ __device__ cu_qcomp operator()(qindex n) {

        // i = global index of nth local amp of rho
        qindex i = concatenateBits(rank, n, logNumAmpsPerNode);

        // r, c = global row and column indices corresponding to i
        qindex r = getBitsRightOfIndex(i, numQubits);
        qindex c = getBitsLeftOfIndex(i, numQubits-1);

        // collect amps involved in this term
        cu_qcomp rhoAmp = rho[n];
        cu_qcomp rowAmp = psi[r];
        cu_qcomp colAmp = psi[c];

        // compute term of <psi|rho^dagger|psi> or <psi|rho|psi>
        if constexpr (Conj) {
            rhoAmp = getCompConj(rhoAmp);
            colAmp = getCompConj(colAmp);
        } else
            rowAmp = getCompConj(rowAmp);

        cu_qcomp fid = rhoAmp * rowAmp * colAmp;
        return fid;
    }
};


template <int NumTargets>
struct functor_projectStateVec : public thrust::binary_function<qindex,cu_qcomp,cu_qcomp> {

    // this functor multiplies an amp with zero or a 
    // renormalisation codfficient, depending on whether
    // the basis state of the amp has qubits in a particular
    // configuration. This is used to project statevector
    // qubits into a particular measurement outcome

    int* targetsPtr;
    int numTargets, rank;
    qindex retainValue;
    qreal renorm;

    functor_projectStateVec(
        int* targetsPtr, int numTargets, 
        qindex retainValue, qreal renorm
    ) :
        targetsPtr(targetsPtr), numTargets(numTargets),
        retainValue(retainValue), renorm(renorm)
    { 
        assert_numTargsMatchesTemplateParam(numTargets, NumTargets);
    }

    __host__ __device__ cu_qcomp operator()(qindex n, cu_qcomp amp) {

        // use the compile-time value if possible, to auto-unroll the getValueOfBits() loop below
        SET_VAR_AT_COMPILE_TIME(int, numBits, NumTargets, numTargets);

        // return amp scaled by zero or renorm, depending on whether n has projected substate
        qindex val = getValueOfBits(n, targetsPtr, numBits);
        qreal fac = renorm * (val == retainValue);
        return fac * amp;
    }
};


template <int NumTargets>
struct functor_projectDensMatr : public thrust::binary_function<qindex,cu_qcomp,cu_qcomp> {

    // this functor multiplies an amp with zero or a 
    // renormalisation coefficient, depending on whether
    // the basis state of the amp has qubits in a particular
    // configuration. This is used to project density matrix
    // qubits into a particular measurement outcome

    int* targetsPtr;
    int numTargets, rank, numQuregQubits;
    qindex logNumAmpsPerNode, retainValue;
    qreal renorm;

    functor_projectDensMatr(
        int* targetsPtr, int numTargets, int rank, int numQuregQubits,
        qindex logNumAmpsPerNode, qindex retainValue, qreal renorm
    ) :
        targetsPtr(targetsPtr), numTargets(numTargets), rank(rank), numQuregQubits(numQuregQubits),
        logNumAmpsPerNode(logNumAmpsPerNode), retainValue(retainValue), renorm(renorm)
    { 
        assert_numTargsMatchesTemplateParam(numTargets, NumTargets);
    }

    __host__ __device__ cu_qcomp operator()(qindex n, cu_qcomp amp) {

        // use the compile-time value if possible, to auto-unroll the getValueOfBits() loop below
        SET_VAR_AT_COMPILE_TIME(int, numBits, NumTargets, numTargets);

        // i = global index of nth local amp
        qindex i = concatenateBits(rank, n, logNumAmpsPerNode);

        // r, c = global row and column indices of nth local amp
        qindex r = getBitsRightOfIndex(i, numQuregQubits);
        qindex c = getBitsLeftOfIndex(i, numQuregQubits-1);

        qindex v1 = getValueOfBits(r, targetsPtr, numBits);
        qindex v2 = getValueOfBits(c, targetsPtr, numBits);

        // multiply amp with renorm or zero if values disagree with given outcomes
        qreal fac = renorm * (v1 == v2) * (retainValue == v1);
        return fac * amp;
    }
};


struct functor_setRandomStateVecAmp : public thrust::unary_function<qindex,cu_qcomp> {

    // this functor generates a random, unnormalised
    // statevector amplitude which, after normalisation
    // of all amps, produces uniformly-random pure states

    /// @todo
    /// this method of parallel RNG is slow, since every
    /// amplitude uses an independent freshly-created 
    /// generator, as per the limitations of Thrust. A
    /// lower-level method like use of cuRAND may prove
    /// faster, though we should ensure continued compatibility
    /// with AMD GPUs via HIP/ROCm and RocThrust. We should
    /// first quantify the speed of this function in comparison to
    /// a single-qubit gate; being slower than <5 gates is acceptable

    unsigned baseSeed;
    functor_setRandomStateVecAmp(unsigned seed) : baseSeed(seed) {}

    __host__ __device__ cu_qcomp operator()(qindex ampInd) {

        // wastefully create new distributions for every amp
        thrust::random::normal_distribution<qreal> normDist(0, 1); // mean=0, var=1
        thrust::random::uniform_real_distribution<qreal> phaseDist(0, 2*const_PI); // ~ [0, 2pi]

        // wastefully initialise a new generator for every amp...
        thrust::random::default_random_engine gen;

        // which we uniquely seed, as opposed to commonly-seeding and advancing
        // each thread by a different amount; this avoids gen.discard() having 
        // to serially perform 2^N advances in the final thread. Alas, we have to
        // prepare the unique seed (combining baseSeed and ampInd) ourselves, being 
        // unable to use std::seed_seq; we'll use SplitMix64, adapted from:
        // https://xoshiro.di.unimi.it/splitmix64.c

        unsigned long long uniqueSeed = ampInd + baseSeed;
        uniqueSeed = (uniqueSeed ^ (uniqueSeed >> 30)) * 0xbf58476d1ce4e5b9ULL;
        uniqueSeed = (uniqueSeed ^ (uniqueSeed >> 27)) * 0x94d049bb133111ebULL;
        uniqueSeed = (uniqueSeed ^ (uniqueSeed >> 31));
        gen.seed(uniqueSeed);

        // using the uniquely-seeded generator, produce a few variates to inform one amp;
        // see an explanation of the maths in randomiser.cpp's rand_getThreadPrivateRandomAmp()
        qreal n1 = normDist(gen);
        qreal n2 = normDist(gen);
        qreal prob = n1*n1 + n2*n2;
        qreal phase = phaseDist(gen);
        auto iphase = thrust::complex<qreal>(0, phase);
        auto amp = sqrt(prob) * thrust::exp(iphase); // CUDA sqrt

        // cast thrust::complex to cu_qcomp
        return getCuQcomp(amp.real(), amp.imag());
    }
};



/*
 * MATRIX INITIALISATION
 */


void thrust_fullstatediagmatr_setElemsToPauliStrSum(FullStateDiagMatr out, PauliStrSum in) {

    // copy 'in' lists into GPU memory, which is not a big deal even when 'in'
    // is very large, because we only do this during FullStateDiagMatr initialisation
    thrust::device_vector<qcomp> devCoeffs(in.coeffs, in.coeffs + in.numTerms);
    thrust::device_vector<PauliStr> devStrings(in.strings, in.strings + in.numTerms);
    
    // obtain raw pointers which can be passed to fastmath.hpp routines
    cu_qcomp* devCoeffsPtr = toCuQcomps(thrust::raw_pointer_cast(devCoeffs.data()));
    PauliStr* devStringsPtr = thrust::raw_pointer_cast(devStrings.data());

    int rank = out.isDistributed? comm_getRank() : 0;
    qindex logNumElemsPerNode = logBase2(out.numElemsPerNode);

    // <true> indicates the PauliStrSum is diagonal (contains only I or Z)
    auto functor = functor_setAmpToPauliStrSumElem<true>(
        rank, out.numElems, logNumElemsPerNode,
        in.numTerms, toCuQcomps(out.gpuElems), devCoeffsPtr, devStringsPtr);

    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + out.numElemsPerNode;
    thrust::for_each(indIter, endIter, functor);
}



/*
 * STATE MODIFICATION 
 */


void thrust_setElemsToConjugate(cu_qcomp* matrElemsPtr, qindex matrElemsLen) {

    auto ptr = getStartPtr(matrElemsPtr);
    thrust::transform(ptr, ptr + matrElemsLen, ptr, functor_getAmpConj());
}


void thrust_densmatr_setAmpsToPauliStrSum_sub(Qureg qureg, PauliStrSum sum) {

    // this assertion exists because fast_getPauliStrElem() (invoked in functor)
    // previously assumed str.highPaulis=0 for all str in sum (for a speedup)
    // which is gauranteed satisfied for all sum compatible with a density-matrix.
    // This is no longer essential, since fast_getPauliStrElem() has relaxed this
    // requirement and foregone the optimisation, but we retain this check in
    // case a similar optimisation is restored in the future
    assert_highPauliStrSumMaskIsZero(sum);

    // copy sum lists into GPU memory, which is not a big deal even when sum
    // is very large, because we only do this during Qureg initialisation (infrequent)
    thrust::device_vector<qcomp> devCoeffs(sum.coeffs, sum.coeffs + sum.numTerms);
    thrust::device_vector<PauliStr> devStrings(sum.strings, sum.strings + sum.numTerms);
    
    // obtain raw pointers which can be passed to fastmath.hpp routines
    cu_qcomp* devCoeffsPtr = toCuQcomps(thrust::raw_pointer_cast(devCoeffs.data()));
    PauliStr* devStringsPtr = thrust::raw_pointer_cast(devStrings.data());

    // <false> indicates the PauliStrSum is not diagonal (contains X or Y)
    auto functor = functor_setAmpToPauliStrSumElem<false>(
        qureg.rank, powerOf2(qureg.numQubits), qureg.logNumAmpsPerNode,
        sum.numTerms, toCuQcomps(qureg.gpuAmps), devCoeffsPtr, devStringsPtr);

    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + qureg.numAmpsPerNode;
    thrust::for_each(indIter, endIter, functor);
}


void thrust_densmatr_mixQureg_subA(qreal outProb, Qureg outQureg, qreal inProb, Qureg inQureg) {

    thrust::transform(
        getStartPtr(outQureg), getEndPtr(outQureg), 
        getStartPtr(inQureg),  getStartPtr(outQureg), // 4th arg is output pointer
        functor_mixAmps(outProb, inProb));
}


void thrust_statevec_setQuregToSuperposition_sub(cu_qcomp facOut, Qureg outQureg, cu_qcomp fac1, Qureg inQureg1, cu_qcomp fac2, Qureg inQureg2) {

    thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(getStartPtr(outQureg), getStartPtr(inQureg1), getStartPtr(inQureg2))),
        thrust::make_zip_iterator(thrust::make_tuple(getEndPtr(outQureg),   getEndPtr(inQureg1),   getEndPtr(inQureg2))),
        functor_superposeAmps(facOut, fac1, fac2));
}


template <bool HasPower>
void thrust_statevec_allTargDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, cu_qcomp exponent) {

    thrust::transform(
        getStartPtr(qureg), getEndPtr(qureg), 
        getStartPtr(matr),  getStartPtr(qureg), // 4th arg is output pointer
        functor_multiplyElemPowerWithAmpOrNorm<HasPower,false,false>(exponent));
}



/*
 * PROBABILITIES
 */


qreal thrust_statevec_calcTotalProb_sub(Qureg qureg) {

    /// @todo
    /// this function should be more accurate than other
    /// functions (like calcProbOfQubitOutcome) because it
    /// is primarily used for checking normalisation and ergo
    /// numerical accuracy. We such, naively trusting thrust's
    /// reduction may be suboptimal; it may be necessary to
    /// implement a custom numerically-stable CUDA reduction.

    // this being a 0 integral literal instead of a 0. float 
    // literal causes a silent Thrust error. Grr...
    qreal init = 0.0;

    qreal prob = thrust::transform_reduce(
        getStartPtr(qureg), getEndPtr(qureg), 
        functor_getAmpNorm(), init, thrust::plus<qreal>());

    return prob;
}


qreal thrust_densmatr_calcTotalProb_sub(Qureg qureg) {

    /// @todo
    /// this function should be more accurate than other
    /// functions (like calcProbOfQubitOutcome) because it
    /// is primarily used for checking normalisation and ergo
    /// numerical accuracy. We such, naively trusting thrust's
    /// reduction may be suboptimal; it may be necessary to
    /// implement a custom numerically-stable CUDA reduction.

    auto rawIter = thrust::make_counting_iterator(0);
    auto indIter = thrust::make_transform_iterator(rawIter, functor_getDiagInd(qureg));
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), indIter);
    auto probIter= thrust::make_transform_iterator(ampIter, functor_getAmpReal());

    qindex numIts = powerOf2(qureg.logNumColsPerNode);
    qreal prob = thrust::reduce(probIter, probIter + numIts);
    return prob;
}


template <int NumQubits>
qreal thrust_statevec_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    devints sortedQubits = util_getSorted(qubits);
    qindex valueMask = util_getBitMask(qubits, outcomes);

    auto indFunctor = functor_insertBits<NumQubits>(getPtr(sortedQubits), valueMask, qubits.size());
    auto probFunctor = functor_getAmpNorm();

    auto rawIter = thrust::make_counting_iterator(0);
    auto indIter = thrust::make_transform_iterator(rawIter, indFunctor);
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), indIter);
    auto probIter = thrust::make_transform_iterator(ampIter, probFunctor);

    qindex numIts = qureg.numAmpsPerNode / powerOf2(qubits.size());
    qreal prob = thrust::reduce(probIter, probIter + numIts);
    return prob;
}


template <int NumQubits>
qreal thrust_densmatr_calcProbOfMultiQubitOutcome_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes) {

    // cannot move these into functor_insertBits constructor, since the memory
    // would dangle - and we cannot bind deviceints as an attribute - it's host-only!
    devints sortedQubits = util_getSorted(qubits);
    qindex valueMask = util_getBitMask(qubits, outcomes);

    auto basisIndFunctor = functor_insertBits<NumQubits>(getPtr(sortedQubits), valueMask, qubits.size());
    auto diagIndFunctor = functor_getDiagInd(qureg);
    auto probFunctor = functor_getAmpReal();

    auto rawIter = thrust::make_counting_iterator(0);
    auto indIter = thrust::make_transform_iterator(rawIter, basisIndFunctor);
    auto diagIter= thrust::make_transform_iterator(indIter, diagIndFunctor);
    auto ampIter = thrust::make_permutation_iterator(getStartPtr(qureg), diagIter);
    auto probIter= thrust::make_transform_iterator(ampIter, probFunctor);

    qindex numIts = powerOf2(qureg.logNumColsPerNode - qubits.size());
    qreal prob = thrust::reduce(probIter, probIter + numIts);
    return prob;
}



/*
 * INNER PRODUCTS AND MEASURES
 */


cu_qcomp thrust_statevec_calcInnerProduct_sub(Qureg quregA, Qureg quregB) {

    cu_qcomp init = getCuQcomp(0, 0);

    cu_qcomp prod = thrust::inner_product(
        getStartPtr(quregA), getEndPtr(quregA), getStartPtr(quregB), 
        init, thrust::plus<cu_qcomp>(), functor_getAmpConjProd());

    return prod;
}


qreal thrust_densmatr_calcHilbertSchmidtDistance_sub(Qureg quregA, Qureg quregB) {

    qreal init = 0;

    qreal dist = thrust::inner_product(
        getStartPtr(quregA), getEndPtr(quregA), getStartPtr(quregB), 
        init, thrust::plus<qreal>(), functor_getNormOfAmpDif());

    return dist;
}


template <bool Conj>
cu_qcomp thrust_densmatr_calcFidelityWithPureState_sub(Qureg rho, Qureg psi) {

    // functor accepts an index and produces a cu_qcomp
    auto functor = functor_getFidelityTerm<Conj>(
        rho.rank, rho.numQubits, rho.logNumAmpsPerNode, 
        psi.numAmps, toCuQcomps(rho.gpuAmps), toCuQcomps(psi.gpuAmps));

    auto indIter = thrust::make_counting_iterator(0);
    qindex numIts = rho.numAmpsPerNode;

    cu_qcomp init = getCuQcomp(0, 0);
    cu_qcomp fid = thrust::transform_reduce(
        indIter, indIter + numIts, 
        functor, init, thrust::plus<cu_qcomp>());

    return fid;
}



/*
 * PAULI EXPECTATION VALUES
 */


qreal thrust_statevec_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs) {

    qindex mask = util_getBitMask(targs);
    auto functor = functor_getExpecStateVecZTerm(mask);

    qreal init = 0;
    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + qureg.numAmpsPerNode;

    return thrust::inner_product(
        indIter, endIter, getStartPtr(qureg), 
        init, thrust::plus<qreal>(), functor);
}


cu_qcomp thrust_densmatr_calcExpecAnyTargZ_sub(Qureg qureg, vector<int> targs) {

    qindex dim = powerOf2(qureg.numQubits);
    qindex ind = util_getLocalIndexOfFirstDiagonalAmp(qureg);
    qindex mask = util_getBitMask(targs);
    auto functor = functor_getExpecDensMatrZTerm(dim, ind, mask, toCuQcomps(qureg.gpuAmps));

    cu_qcomp init = getCuQcomp(0, 0);
    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + powerOf2(qureg.logNumColsPerNode);

    return thrust::transform_reduce(indIter, endIter, functor, init, thrust::plus<cu_qcomp>());
}


cu_qcomp thrust_statevec_calcExpecPauliStr_subA(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

    qindex maskXY = util_getBitMask(util_getConcatenated(x, y));
    qindex maskYZ = util_getBitMask(util_getConcatenated(y, z));
    auto ampsPtr = toCuQcomps(qureg.gpuAmps);
    auto functor = functor_getExpecStateVecPauliTerm(maskXY, maskYZ, ampsPtr, ampsPtr); // amps=pairAmps

    cu_qcomp init = getCuQcomp(0, 0);
    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + qureg.numAmpsPerNode;

    cu_qcomp value = thrust::transform_reduce(indIter, endIter, functor, init, thrust::plus<cu_qcomp>());

    return value * toCuQcomp(util_getPowerOfI(y.size()));
}


cu_qcomp thrust_statevec_calcExpecPauliStr_subB(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

    qindex maskXY = util_getBitMask(util_getConcatenated(x, y));
    qindex maskYZ = util_getBitMask(util_getConcatenated(y, z));
    auto ampsPtr = toCuQcomps(qureg.gpuAmps);
    auto buffPtr = toCuQcomps(qureg.gpuCommBuffer);
    auto functor = functor_getExpecStateVecPauliTerm(maskXY, maskYZ, ampsPtr, buffPtr);

    cu_qcomp init = getCuQcomp(0, 0);
    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + qureg.numAmpsPerNode;

    cu_qcomp value = thrust::transform_reduce(indIter, endIter, functor, init, thrust::plus<cu_qcomp>());

    return value * toCuQcomp(util_getPowerOfI(y.size()));
}


cu_qcomp thrust_densmatr_calcExpecPauliStr_sub(Qureg qureg, vector<int> x, vector<int> y, vector<int> z) {

    qindex mXY = util_getBitMask(util_getConcatenated(x, y));
    qindex mYZ = util_getBitMask(util_getConcatenated(y, z));
    qindex dim = powerOf2(qureg.numQubits);
    qindex ind = util_getLocalIndexOfFirstDiagonalAmp(qureg);
    auto functor = functor_getExpecDensMatrPauliTerm(mXY, mYZ, dim, ind, toCuQcomps(qureg.gpuAmps));

    cu_qcomp init = getCuQcomp(0, 0);
    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + powerOf2(qureg.logNumColsPerNode);

    cu_qcomp value = thrust::transform_reduce(indIter, endIter, functor, init, thrust::plus<cu_qcomp>());

    return value * toCuQcomp(util_getPowerOfI(y.size()));
}



/*
 * DIAGONAL MATRIX EXPECTATION VALUES
 */


template <bool HasPower, bool UseRealPow> 
cu_qcomp thrust_statevec_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, cu_qcomp expo) {

    cu_qcomp init = getCuQcomp(0, 0);
    auto functor = functor_multiplyElemPowerWithAmpOrNorm<HasPower,UseRealPow,true>(expo);

    cu_qcomp value = thrust::inner_product(
        getStartPtr(qureg), getEndPtr(qureg), getStartPtr(matr), 
        init, thrust::plus<cu_qcomp>(), functor);

    return value;
}


template <bool HasPower, bool UseRealPow> 
cu_qcomp thrust_densmatr_calcExpecFullStateDiagMatr_sub(Qureg qureg, FullStateDiagMatr matr, cu_qcomp expo) {

    qindex dim = powerOf2(qureg.numQubits);
    qindex ind = util_getLocalIndexOfFirstDiagonalAmp(qureg);
    auto ampsPtr = toCuQcomps(qureg.gpuAmps);
    auto elemsPtr = toCuQcomps(matr.gpuElems);
    auto functor = functor_getExpecDensMatrDiagMatrTerm<HasPower,UseRealPow>(dim, ind, ampsPtr, elemsPtr, expo);

    cu_qcomp init = getCuQcomp(0, 0);
    auto indIter = thrust::make_counting_iterator(0);
    auto endIter = indIter + powerOf2(qureg.logNumColsPerNode);

    return thrust::transform_reduce(indIter, endIter, functor, init, thrust::plus<cu_qcomp>());
}



/*
 * PROJECTORS
 */


template <int NumQubits>
void thrust_statevec_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal renorm) {

    devints devQubits = qubits;
    qindex retainValue = getIntegerFromBits(outcomes.data(), outcomes.size());
    auto projFunctor = functor_projectStateVec<NumQubits>(
        getPtr(devQubits), qubits.size(), retainValue, renorm);

    auto indIter = thrust::make_counting_iterator(0);
    auto ampIter = getStartPtr(qureg);

    qindex numIts = qureg.numAmpsPerNode;
    thrust::transform(indIter, indIter + numIts, ampIter, ampIter, projFunctor); // 4th arg gets modified
}


template <int NumQubits>
void thrust_densmatr_multiQubitProjector_sub(Qureg qureg, vector<int> qubits, vector<int> outcomes, qreal renorm) {

    devints devQubits = qubits;
    qindex retainValue = getIntegerFromBits(outcomes.data(), outcomes.size());
    auto projFunctor = functor_projectDensMatr<NumQubits>(
        getPtr(devQubits), qubits.size(), qureg.rank, qureg.numQubits,
        qureg.logNumAmpsPerNode, retainValue, renorm);

    auto indIter = thrust::make_counting_iterator(0);
    auto ampIter = getStartPtr(qureg);

    qindex numIts = qureg.numAmpsPerNode;
    thrust::transform(indIter, indIter + numIts, ampIter, ampIter, projFunctor); // 4th arg gets modified
}



/*
 * STATE INITIALISATION
 */


void thrust_statevec_initUniformState(Qureg qureg, cu_qcomp amp) {

    thrust::fill(getStartPtr(qureg), getEndPtr(qureg), amp);
}


void thrust_statevec_initDebugState_sub(Qureg qureg) {

    // globally, |n> gains coefficient 2n/10 + i(2n+1)/10,
    // which is a step-size of 2/10 + i(2/10)...
    cu_qcomp step = getCuQcomp(2/10., 2/10.);

    // and each node begins from a unique n (if distributed)
    qindex n = util_getGlobalIndexOfFirstLocalAmp(qureg);
    cu_qcomp init = getCuQcomp(2*n/10., (2*n+1)/10.);

    thrust::sequence(getStartPtr(qureg), getEndPtr(qureg), init, step);
}


void thrust_statevec_initUnnormalisedUniformlyRandomPureStateAmps_sub(Qureg qureg) {

    // thread amp generators uniquely perturb a common base seed
    unsigned seed = rand_getThreadSharedRandomSeed(qureg.isDistributed);
    auto functor = functor_setRandomStateVecAmp(seed);

    auto indIter = thrust::make_counting_iterator(0);
    auto ampIter = getStartPtr(qureg);

    qindex numIts = qureg.numAmpsPerNode;
    thrust::transform(indIter, indIter + numIts, ampIter, functor); // 3rd arg gets modified
}



#endif // GPU_THRUST_HPP