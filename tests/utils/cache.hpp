/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilscache Cache
 * @ingroup testutils
 * @brief 
 * Testing utilities which create Quregs across all
 * available hardware deployments
 * @{
 */

#ifndef CACHE_HPP
#define CACHE_HPP

#include "quest.h"

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <unordered_map>
#include <vector>
#include <string>
#include <tuple>

using quregCache  = std::unordered_map<std::string,Qureg>;
using matrixCache = std::unordered_map<std::string,FullStateDiagMatr>;
using deployInfo  = std::vector<std::tuple<std::string,int,int,int>>;


/*
 * CUSTOM CACHE
 *
 * for obtaining and manually maintaining Quregs of all possible
 * deployment types, but with the specified dimensions
 */

quregCache createCustomCachedQuregs(int numQubits, bool isDensityMatrix);
void destroyCustomCachedQuregs(quregCache& cache);


/*
 * MAIN CACHE
 *
 * for obtaining Quregs and FullStateDiagMatrs of all possible
 * deployment types, managed by the test utils, and with 
 * fixed dimensions specific to the test config
 */

int getNumCachedQubits();
deployInfo getSupportedDeployments();

void createCachedFullStateDiagMatrs();
void destroyCachedFullStateDiagMatrs();
matrixCache getCachedFullStateDiagMatrs();

void createCachedQuregs();
void destroyCachedQuregs();
quregCache getCachedStatevecs();
quregCache getCachedDensmatrs();
quregCache getAltCachedStatevecs();
quregCache getAltCachedDensmatrs();

Qureg getArbitraryCachedStatevec();
Qureg getArbitraryCachedDensmatr();


/*
 * REFERENCE STATES
 */

qvector getRefStatevec();
qmatrix getRefDensmatr();


#endif // CACHE_HPP

/** @} (end defgroup) */
