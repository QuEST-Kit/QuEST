#ifndef CACHE_HPP
#define CACHE_HPP

#include "quest.h"
#include <unordered_map>
#include <vector>
#include <string>
#include <tuple>

using quregCache = std::unordered_map<std::string,Qureg>;
using deployInfo = std::vector<std::tuple<std::string,int,int,int>>;

deployInfo getSupportedDeployments();

void createCachedQuregs();
void destroyCachedQuregs();

quregCache getCachedStatevecs();
quregCache getCachedDensmatrs();


#endif // CACHE_HPP