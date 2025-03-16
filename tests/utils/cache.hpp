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

#include "quest/include/quest.h"

#include "qvector.hpp"
#include "qmatrix.hpp"

#include <unordered_map>
#include <vector>
#include <string>
#include <tuple>

using quregCache  = std::unordered_map<std::string,Qureg>;
using matrixCache = std::unordered_map<std::string,FullStateDiagMatr>;
using deployInfo  = std::vector<std::tuple<std::string,int,int,int>>;

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

qvector getRefStatevec();
qmatrix getRefDensmatr();


#endif // CACHE_HPP

/** @} (end defgroup) */
