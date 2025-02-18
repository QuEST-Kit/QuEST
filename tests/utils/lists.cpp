#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include "macros.hpp"
#include "random.hpp"

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


using listpair = tuple<
    vector<int>,
    vector<int>>;


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
