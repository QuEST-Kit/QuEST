/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#include "quest/include/types.h"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <vector>
#include <tuple>
#include <string>


// I like to live dangerously; and I am willing to debate with the
// most dogmatic C++ enthusiast about the unreadability of namespaces
using namespace std;



/*
 * TYPE STRINGS
 */


// try safely check for de-mangler library
#if defined __has_include
    #if __has_include (<cxxabi.h>)
        #include <cxxabi.h>
        #define DEMANGLE_TYPE true
    #else
        #define DEMANGLE_TYPE false
    #endif

// fall-back to dangerous compiler-specific check (we know MSVC doesn't mangle)
#elif ! defined _MSC_VER

    #warning "Attempting to include compiler-specific library <cxxabi.h>"

    #include <cxxabi.h>
    #define DEMANGLE_TYPE true
#else
    #define DEMANGLE_TYPE false
#endif


// macros for printing multi-word macros (necessary indirection)
#define GET_STR_INTERNAL(x) #x
#define GET_STR(x) GET_STR_INTERNAL(x)


template<typename T> string getTypeName(T) { 

    // get possibly-mangled type name
    string name = typeid(T).name();

    // non-MSVC compilers can de-mangle
    #ifdef DEMANGLE_TYPE
        name = abi::__cxa_demangle(name.c_str(), nullptr, nullptr, nullptr);
    #endif

    return name;
}


string form_getQcompType() {

    qcomp x;
    return getTypeName(x);
}

string form_getQrealType() {

    // more portable than getTypeName
    return GET_STR( FLOAT_TYPE );
}

string form_getQindexType() {

    // more portable than getTypeName
    return GET_STR( INDEX_TYPE );
}

string form_getFloatPrecisionFlag() {

    return GET_STR( FLOAT_PRECISION );
}



/*
 * STRING CASTS
 */

// all types not explicitly overloaded below use template in header

string form_str(qreal num) {

	// it's unbelievable this was the best C++17 compatible method
	char repr[64];
	sprintf(repr, "%g", num);
	return string(repr);
}



/*
 * TABLE PRINTING
 */


void form_printTable(string title, vector<tuple<string, string>> rows, string indent) {

	int minColumnGap = 5;

	// find max-width of left column
	int maxWidth = 0;
	for (auto const& [key, value] : rows)
		if (key.length() > maxWidth)
			maxWidth = key.length();

	// print table title (indented)
	cout << indent << "[" << title << "]" << endl;

	// pad left column (which is doubly indented) to align right column
	for (auto const& [key, value] : rows)
		cout << indent << indent << left << setw(maxWidth + minColumnGap) << setfill('.') << key << value << endl;

}


void form_printTable(string title, vector<tuple<string, long long int>> rows, string indent) {

	// convert all values to strings
	vector<tuple<string, string>> casted;
	for (auto const& [key, value] : rows)
		casted.push_back({key, to_string(value)});

	form_printTable(title, casted, indent);
}
