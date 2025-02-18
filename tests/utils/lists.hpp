#ifndef LISTS_HPP
#define LISTS_HPP

#include <catch2/generators/catch_generators_adapters.hpp>

#include <vector>
#include <tuple>


using std::vector;

vector<int> getRange(int start, int endExcl);
vector<int> getRange(int endExcl);

vector<int> getSublist(vector<int> list, int start, int len);


template<class T> 
using CatchGen = Catch::Generators::GeneratorWrapper<T>;
using listpair = std::tuple<vector<int>,vector<int>>;

CatchGen<vector<int>> sublists(CatchGen<int>&& gen, int sublen);
CatchGen<listpair> disjointsublists(CatchGen<int>&& gen, int sublen1, int sublen2);


#endif // LISTS_HPP