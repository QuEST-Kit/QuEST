#ifndef CACHE_HPP
#define CACHE_HPP

#include "quest.h"
#include <unordered_map>

using quregCache = std::unordered_map<std::string,Qureg>;

void createCachedQuregs();
void destroyCachedQuregs();

quregCache getCachedStatevecs();
quregCache getCachedDensmatrs();


#endif // CACHE_HPP