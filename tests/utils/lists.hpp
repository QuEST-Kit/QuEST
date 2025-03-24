/** @file
 * @author Tyson Jones
 * 
 * @defgroup testutilslists Lists
 * @ingroup testutils
 * @brief
 * Testing utilities which generate lists of integers.
 * @{
 */

#ifndef LISTS_HPP
#define LISTS_HPP

#include <catch2/generators/catch_generators_adapters.hpp>

#include "quest/include/quest.h"

#include <vector>
#include <tuple>

using std::vector;

using listpair = std::tuple<vector<int>,vector<int>>;
using listtrio = std::tuple<vector<int>,vector<int>,vector<int>>;


vector<int> getRange(int start, int endExcl);
vector<int> getRange(int endExcl);
vector<int> getComplement(vector<int> listA, vector<int> listB);
vector<int>   getSublist(vector<int>   list, int start, int len);
vector<qcomp> getSublist(vector<qcomp> list, int start, int len);


template<class T> using CatchGen = Catch::Generators::GeneratorWrapper<T>;
CatchGen<vector<int>> sublists(CatchGen<int>&& gen, int sublen);
CatchGen<listpair> disjointsublists(CatchGen<int>&& gen, int sublen1, int sublen2);


vector<int> GENERATE_TARGS(int numQubits, int numTargs);
listpair GENERATE_CTRLS_AND_TARGS(int numQubits, int numCtrls, int numTargs);


#endif // LISTS_HPP

/** @} (end defgroup) */
