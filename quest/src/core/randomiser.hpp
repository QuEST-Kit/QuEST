/** @file
 * Functions for generating random numbers consistent
 * between distributed nodes, used for emulating
 * quantum measurements and randomly initialising Quregs.
 * 
 * @author Tyson Jones
 */

#ifndef RANDOMISER_HPP
#define RANDOMISER_HPP

#include "quest/include/types.h"

#include <vector>
#include <random>

using std::vector;



/*
 * SEEDING
 */


void rand_setSeeds(vector<unsigned> seeds);

void rand_setSeedsToDefault();

int rand_getNumSeeds();

vector<unsigned> rand_getSeeds();



/*
 * SAMPLING
 */


int rand_getRandomSingleQubitOutcome(qreal probOfZero);

qindex rand_getRandomMultiQubitOutcome(vector<qreal> probs);



/*
 * STATE AMPLITUDE SAMPLING
 */


unsigned rand_getThreadSharedRandomSeed(bool distinctPerNode);

std::mt19937_64 rand_getThreadPrivateGenerator(unsigned sharedSeed, int threadId);

std::normal_distribution<qreal> rand_getThreadPrivateAmpAbsDistribution();

std::uniform_real_distribution<qreal> rand_getThreadPrivateAmpPhaseDistribution();

qcomp rand_getThreadPrivateRandomAmp(std::mt19937_64 &gen, std::normal_distribution<qreal> &normDist, std::uniform_real_distribution<qreal> &phaseDist);



#endif // RANDOMISER_HPP