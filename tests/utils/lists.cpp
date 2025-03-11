/** @file
 * Testing utilities which generate lists of integers.
 *
 * @author Tyson Jones
 */

#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include "lists.hpp"
#include "macros.hpp"
#include "random.hpp"
#include "linalg.hpp"

#include <tuple>
#include <vector>
#include <algorithm>

using std::tuple;
using std::vector;

using Catch::Generators::IGenerator;
using Catch::Generators::GeneratorWrapper;



/*
 * PRIVATE
 */


vector<int> getGeneratorElems(GeneratorWrapper<int>& gen) {

    vector<int> list;
    do { list.push_back(gen.get()); } while (gen.next());

    return list;
}



/*
 * SUBLISTS
 */


vector<int> getSublist(vector<int> list, int start, int len) {

    return vector<int>(list.begin() + start, list.begin() + start + len);
}



/*
 * RANGE
 */


vector<int> getRange(int start, int endExcl) {
    DEMAND( endExcl >= start );

    vector<int> out(endExcl - start);

    for (size_t i=0; i<out.size(); i++)
        out[i] = start + i;

    return out;
}


vector<int> getRange(int endExcl) {
    return getRange(0, endExcl);
}



/*
 * COMPLEMENTS
 */


vector<int> getComplement(vector<int> listA, vector<int> listB) {

    std::sort(listA.begin(), listA.end());
    std::sort(listB.begin(), listB.end());

    vector<int> out;
    std::set_difference(
        listA.begin(), listA.end(), 
        listB.begin(), listB.end(), 
        std::back_inserter(out));

    return out;
}



/*
 * SUBLIST GENERATORs
 */


class SublistGenerator final : public IGenerator<vector<int>> {
    
    // list of all possible (unique) elements
    vector<int> list;

    // output sublist of user-specified length
    int sublen;
    vector<int> sublist;

    // indicates which elements of list are in sublist
    vector<bool> featured;

private:

    void prepareSublist() {
        
        // choose the next combination
        int j=0;
        for (size_t i=0; i<list.size(); i++)
            if (featured[i])
                sublist[j++] = list[i];
                
        // prepare for permuting
        std::sort(sublist.begin(), sublist.end());
    }

public:

    SublistGenerator(vector<int> list, int sublen): 
        list(list),
        sublen(sublen)
    {
        // ensure sublist would not be longer than list
        DEMAND( sublen <= list.size() );

        // populate sublist with first elements
        sublist.resize(sublen);
        featured.resize(list.size());
        fill(featured.end() - sublen, featured.end(), true);

        prepareSublist();
    }

    vector<int> const& get() const override {

        // return a copy of sublist
        return sublist;
    }
    
    bool next() override {
        
        // offer next permutation of the current combination
        if (std::next_permutation(sublist.begin(), sublist.end()))
            return true;

        // else generate the next combination
        if (std::next_permutation(featured.begin(), featured.end())) {

            prepareSublist();
            return true;
        }
        
        // else indicate generator is exhausted
        return false;
    }
};


GeneratorWrapper<vector<int>> sublists(GeneratorWrapper<int>&& gen, int sublen) {

    vector<int> list;
    do { list.push_back(gen.get()); } while (gen.next());

    return GeneratorWrapper<vector<int>>(
        Catch::Detail::make_unique<SublistGenerator>(
            list, sublen));
}



/*
 * DISJOINT SUBLISTS GENERATOR
 */


class DisjointSublistsGenerator : public IGenerator<listpair> {

    // list of all possible (unique) elements
    vector<int> list;

    // lengths of each sublist
    int sublen1;
    int sublen2;
    listpair sublists;

    // underlying generator
    SublistGenerator subGen;

public:

    DisjointSublistsGenerator(vector<int> list, int sublen1, int sublen2):
        list(list),
        sublen1(sublen1),
        sublen2(sublen2),
        subGen(list, sublen1 + sublen2)
    {
        static_cast<void>(next()); // ignore bool return
    }

    listpair const& get() const override {
        return sublists;
    }
    
    bool next() override {

        // stop when sublist generator is exhausted
        if (!subGen.next())
            return false;

        // else, obtain the next (sublist) sequence
        vector<int> combined = subGen.get();
        vector<int> sublist1 = getSublist(combined, 0, sublen1);
        vector<int> sublist2 = getSublist(combined, sublen1, sublen2);

        sublists = {sublist1, sublist2};
        return true;
    }
};


GeneratorWrapper<listpair> disjointsublists(GeneratorWrapper<int>&& gen, int sublen1, int sublen2) {

    vector<int> list;
    do { list.push_back(gen.get()); } while (gen.next());

    return GeneratorWrapper<listpair>(
        Catch::Detail::make_unique<DisjointSublistsGenerator>(
            list, sublen1, sublen2));
}



/*
 * GENERATOR WRAPPERS
 */


listpair GENERATE_CTRLS_AND_TARGS(int numQubits, int numCtrls, int numTargs) {
    DEMAND( numQubits >= numCtrls + numTargs );

#if FAST_UNIT_TESTS

    // avoid testing all combinations (for speed), and instead repeatedly
    // generate a random combination of disjoint control and target qubits,
    // avoiding superfluously repeating permutations when numRepeats is small

    int numRepeats = getNumPermutations(numQubits, numCtrls + numTargs);
    if (numRepeats > MAX_NUM_FAST_QUBIT_PERMUTATIONS)
        numRepeats = MAX_NUM_FAST_QUBIT_PERMUTATIONS;
    
    GENERATE_COPY( range(0,numRepeats) );

    return getRandomFixedNumCtrlsTargs(numQubits, numCtrls, numTargs);

#else

    // generate every possible combination of disjoint control and target lists,
    // upon all possible qubits, in every possible order. Note this even generates
    // every ordering of the control qubits (e.g. {1,2} and {2,1}) which should
    // have no effect on the resulting operation, but we check it to be super
    // careful anyway!

    return GENERATE_COPY( disjointsublists(range(0,numQubits), numCtrls, numTargs) );

#endif
}


vector<int> GENERATE_TARGS(int numQubits, int numTargs) {
    DEMAND( numQubits >= numTargs );

    auto [ctrls, targs] = GENERATE_CTRLS_AND_TARGS(numQubits, 0, numTargs);
    return targs;
}
