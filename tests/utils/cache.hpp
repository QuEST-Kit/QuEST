/** @file
 * @author Tyson Jones
 * @defgroup testutilscache Cache
 * @brief 
 *   Testing utilities which create Quregs across all
 *   available hardware deployments
 * @ingroup testutils
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

using quregCache = std::unordered_map<std::string,Qureg>;
using deployInfo = std::vector<std::tuple<std::string,int,int,int>>;

int getNumCachedQubits();
deployInfo getSupportedDeployments();

void createCachedQuregs();
void destroyCachedQuregs();

quregCache getCachedStatevecs();
quregCache getCachedDensmatrs();

qvector getRefStatevec();
qmatrix getRefDensmatr();


#endif // CACHE_HPP

/** @} (end defgroup) */
