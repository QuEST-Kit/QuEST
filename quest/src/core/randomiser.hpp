/** @file
 * Functions for generating random numbers consistent
 * between distributed nodes, used for emulating
 * quantum measurements
 */

#ifndef RANDOMISER_HPP
#define RANDOMISER_HPP

#include "quest/include/types.h"

#include <vector>

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



#endif // RANDOMISER_HPP