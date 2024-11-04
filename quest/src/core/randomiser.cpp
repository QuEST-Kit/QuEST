/** @file
 * Functions for generating random numbers consistent
 * between distributed nodes, used for emulating
 * quantum measurements
 */

#include "quest/include/types.h"

#include "quest/src/core/bitwise.hpp"
#include "quest/src/core/errors.hpp"
#include "quest/src/core/validation.hpp"
#include "quest/src/comm/comm_config.hpp"
#include "quest/src/comm/comm_routines.hpp"

#include <random>
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
static std::mt19937_64 generator;
static std::uniform_real_distribution<qreal> distribution(0, 1);



/*
 * SEEDING
 */


void rand_setSeeds(vector<unsigned> seeds) {

    // remember seeds (in case user wishes to later recall them)
    currentSeeds = seeds;

    // pass all seeds to RNG; we just seeds is consistent between nodes
    std::seed_seq seq(seeds.begin(), seeds.end());
    generator.seed(seq);
}


void rand_setSeedsToDefault() {

    // seeds will attempt to use device's CSPRNG
    std::random_device device;

    // use CSPRNG to generate well-spread seeds
    vector<unsigned> seeds(DEFAULT_NUM_RNG_SEEDS);

    // root node produces seeds
    if (comm_getRank() == 0)
        for (auto &seed : seeds)
            seed = device();
    
    // and broadcasts them to all other nodes...
    if (comm_isInit())
        comm_broadcastUnsignedsFromRoot(seeds.data(), seeds.size());

    // so that all nodes set the same seeds
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
 * SAMPLING
 */


int rand_getRandomSingleQubitOutcome(qreal probOfZero) {

    // assumes 0 <= probOfZero <= 1

    // advances generator on every node, retaining consensus
    qreal sample = distribution(generator);

    // produces 0 with probOfZero probability
    return sample > probOfZero;
}


qindex rand_getRandomMultiQubitOutcome(vector<qreal> probs) {

    // assumes sum(probs) = 1

    // advances generator on every node, retaining consensus
    qreal sample = distribution(generator);

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
