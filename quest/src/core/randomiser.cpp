/** @file
 * Functions for generating random numbers consistent
 * between distributed nodes, used for emulating
 * quantum measurements and randomly initialising Quregs.
 * 
 * @author Tyson Jones
 * @author Balint Koczor (patched v3 MSVC seeding)
 */

#include "quest/include/types.h"

#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/core/constants.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <cmath>
#include <complex>
#include <random>
#include <limits>
#include <vector>

using std::vector;



/*
 * RNG HYPERPARAMTERS
 */

#define DEFAULT_NUM_RNG_SEEDS 4



/*
 * PRIVATE RNG SINGLETONS
 *
 * The use of static ensures we never accidentally expose the
 * singletons to other files. We must always sample/advance the 
 * RNGm on every node at the same time to avoid divergence of
 * outcomes (producing non-L2 states), though this does not
 * require communication; we just ensure rand_() funcs are 
 * never elsewhere invoked inside rank-dependent control flow. 
 */

static vector<unsigned> currentSeeds;

static std::mt19937_64 mainGenerator;



/*
 * SEEDING
 */


void rand_setSeeds(vector<unsigned> seeds) {

    // this function consults only root-node seeds, broadcasting 
    // to other nodes which might not have the existing space to
    // receive them, so we explicitly reserve the needed space

    // all nodes learn root node's #seeds
    unsigned numRootSeeds = seeds.size();
    if (comm_isInit())
        comm_broadcastUnsignedsFromRoot(&numRootSeeds, 1);

    // all nodes ensure they have space to receive root node's seeds
    seeds.resize(numRootSeeds);
    
    // all nodes receive root seeds
    if (comm_isInit())
        comm_broadcastUnsignedsFromRoot(seeds.data(), seeds.size());

    // all nodes remember seeds (in case user wishes to later recall them)
    currentSeeds = seeds;

    // all nodes pass all seeds to RNG
    std::seed_seq seq(seeds.begin(), seeds.end());
    mainGenerator.seed(seq);
}


void rand_setSeedsToDefault() {

    // seeds will attempt to use device's CSPRNG
    std::random_device device;

    // use CSPRNG to generate well-spread seeds
    vector<unsigned> seeds(DEFAULT_NUM_RNG_SEEDS);

    // root node produces seeds (non-root seeds would be ignored)
    if (comm_getRank() == 0)
        for (auto &seed : seeds)
            seed = device();

    rand_setSeeds(seeds);
}


int rand_getNumSeeds() {

    return currentSeeds.size();
}


vector<unsigned> rand_getSeeds() {

    // returns a copy, so currentSeeds stays unmutated
    return currentSeeds;
}



/*
 * MEASUREMENT SAMPLING
 *
 * which is always consistent between nodes, as
 * maintained by identical seeding and RNG advancing
 */


int rand_getRandomSingleQubitOutcome(qreal probOfZero) {

    // assumes 0 <= probOfZero <= 1, and produces a Bernoulli variate,
    // which is always consistent between distributed nodes

    std::uniform_real_distribution<qreal> distrib(0, 1); // ~[0,1]

    // advances generator on every node, retaining consensus
    qreal sample = distrib(mainGenerator);

    // produces 0 with probOfZero probability
    return sample > probOfZero;
}


qindex rand_getRandomMultiQubitOutcome(vector<qreal> probs) {

    // assumes sum(probs) = 1, and produces a multinomial variate
    // which is always consistent between distributed nodes

    std::uniform_real_distribution<qreal> distrib(0, 1); // ~[0,1]

    // advances generator on every node, retaining consensus
    qreal sample = distrib(mainGenerator);

    // map sample to an element of probs
    qreal cumProb = 0;
    for (size_t i=0; i<probs.size(); i++) {
        cumProb += probs[i];

        if (sample < cumProb)
            return i;
    }

    // it is principally possible that cumProb == 1 - eps, and that
    // sample lies within [1-eps, 1], and should take final elem
    if (1 - cumProb <= validateconfig_getEpsilon())
        return probs.size() - 1;

    error_randomiserGivenNonNormalisedProbList();
    return probs.size();
}



/*
 * STATE AMPLITUDE SAMPLING
 */


unsigned rand_getThreadSharedRandomSeed(bool distinctPerNode) {

    // (callable by multiple nodes, but NOT by multiple threads!)

    // this function produces a seed which can be subsequently
    // used by many threads (each perturbing their seed via their
    // thread rank) to generate random amplitudes in parallel.
    // The dispatched seed is informed by the mainGenerator,
    // such that the user's API seeding determines both all
    // measurement outcomes and state amps, whilst retaining
    // independence amps for each subsequent random state

    // a single uniform seed suffices to seed local random amp generation
    unsigned maxSeed = std::numeric_limits<unsigned>::max();
    std::uniform_int_distribution<unsigned> distrib(0, maxSeed); // ~[0 .. max]

    // if we want the same seed on every node, advance all generators
    if (!distinctPerNode)
        return distrib(mainGenerator);

    // otherwise, if we wish for distinct per-node seeds (e.g. a
    // Qureg is distributed in a distributed environment)...
    unsigned mySeed = 0;
    int myRank = comm_getRank();
    int numNodes = comm_getNumNodes();

    // then we must advance all generators #ranks times...
    for (int rank=0; rank<numNodes; rank++) {

        unsigned globalSeed = distrib(mainGenerator);

        // but keep a distinct seed per-node
        if (rank == myRank)
            mySeed = globalSeed;
    }

    return mySeed;
}


std::mt19937_64 rand_getThreadPrivateGenerator(unsigned sharedSeed, int threadId) {

    // combine the shared seed (same per-thread, potentially per-node)
    // with the thread ID so that thread RNG outcomes differ
    std::seed_seq seeds { sharedSeed, static_cast<unsigned>(threadId) };

    // give each thread an independent generator, for thread safety
    std::mt19937_64 threadGenerator;
    threadGenerator.seed(seeds);
    return threadGenerator;
}


std::normal_distribution<qreal> rand_getThreadPrivateAmpAbsDistribution() {

    // there's nothing special about our instantiation that
    // makes it thread-specific; the name of this function
    // reminds the user to instantiate one distribution
    // per-thread just to avoid inefficient re-instantiation
    // per iteration of a hot-loop

    return std::normal_distribution<qreal>(0, 1); // mean=0, var=1
}


std::uniform_real_distribution<qreal> rand_getThreadPrivateAmpPhaseDistribution() {

    // there's nothing special about our instantiation that
    // makes it thread-specific; the name of this function
    // reminds the user to instantiate one distribution
    // per-thread just to avoid inefficient re-instantiation
    // per iteration of a hot-loop

    return std::uniform_real_distribution<qreal>(0, 2*const_PI);
}


qcomp rand_getThreadPrivateRandomAmp(std::mt19937_64 &gen, std::normal_distribution<qreal> &normDist, std::uniform_real_distribution<qreal> &phaseDist) {

    // generate the amp's probability (unnormalised)
    qreal n1 = normDist(gen);
    qreal n2 = normDist(gen);
    qreal prob = n1*n1 + n2*n2;

    // generate the amp's complex phase
    qreal phase = phaseDist(gen);

    // produce an unnormalised amp; subsequent normalisation
    // (divising by all probs between all nodes) will make
    // prob a chi-squared variate, and amp uniformly random.
    // https://sumeetkhatri.com/wp-content/uploads/2020/05/random_pure_states.pdf
    qcomp amp = std::sqrt(prob) * std::exp(phase * 1_i);
    return amp;
}
