/** @file
 * Unit tests of the paulis module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitpaulis Paulis
 * @ingroup unittests
 */

#include "quest/include/quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/qvector.hpp"
#include "tests/utils/qmatrix.hpp"
#include "tests/utils/macros.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/linalg.hpp"
#include "tests/utils/random.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

using std::vector;
using Catch::Matchers::ContainsSubstring;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[paulis]"



/** 
 * TESTS
 * 
 * @ingroup unitpaulis
 * @{
 */


TEST_CASE( "getPauliStr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        GENERATE( range(0,10) );

        std::string charSet = GENERATE( "IXYZ", "ixyz", "0123", "iX2z" );

        int numPaulis = GENERATE( 1, 2, 5, 10, 20, 30, 31, 32, 33, 60, 64 );
        auto targets = getRandomSubRange(0, 64, numPaulis);

        // not necessary to add a terminal char when we specify numPaulis
        vector<char> pauliChars(numPaulis, 'I');
        vector<int>  pauliInts(numPaulis, 0);
        long long unsigned lowValue = 0;
        long long unsigned highValue = 0;

        for (int i=0; i<numPaulis; i++) {
            int v = getRandomInt(1, 3+1);
            pauliInts[i] = v;
            pauliChars[i] = charSet[v];

            int t = targets[i];
            (t < 32)? 
                (lowValue  += v * getPow2(2 * t)): // pow4
                (highValue += v * getPow2(2 * (t-32))); // pow4
        }

        SECTION( LABEL_C_INTERFACE ) {
            
            SECTION( "from chars" ) {

                CAPTURE( targets, pauliChars );

                PauliStr str = getPauliStr(pauliChars.data(), targets.data(), numPaulis);

                REQUIRE( str.lowPaulis  == lowValue );
                REQUIRE( str.highPaulis == highValue );
            }

            SECTION( "from ints" ) {

                CAPTURE( targets, pauliChars );

                PauliStr str = getPauliStr(pauliInts.data(), targets.data(), numPaulis);

                REQUIRE( str.lowPaulis  == lowValue );
                REQUIRE( str.highPaulis == highValue );
            }

            SECTION( "from literal" ) {

                // lazily ignores some above prepared vars

                int targ = targets[0];
                CAPTURE( targ );
                
                const char* in = "X";
                PauliStr str = getPauliStr(in, &targ, 1);
                REQUIRE( str.lowPaulis  == ((targ <  32)? getPow2(2*targ)      : 0) );
                REQUIRE( str.highPaulis == ((targ >= 32)? getPow2(2*(targ-32)) : 0) );
            }
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            // std::string() requires pauliChars includes a terminal char
            pauliChars.push_back('\0');
            std::string in = std::string(pauliChars.data());

            SECTION( "from string" ) {

                CAPTURE( targets, pauliChars );

                PauliStr str = getPauliStr(in, targets.data(), numPaulis);

                REQUIRE( str.lowPaulis  == lowValue );
                REQUIRE( str.highPaulis == highValue );
            }

            SECTION( "from vector" ) {

                CAPTURE( targets, pauliChars );

                PauliStr str = getPauliStr(in, targets);

                REQUIRE( str.lowPaulis  == lowValue );
                REQUIRE( str.highPaulis == highValue );
            }

            SECTION( "from literal" ) {

                // lazily ignores some above prepared vars

                int targ = targets[0];
                CAPTURE( targ );
                
                PauliStr str = getPauliStr("X", &targ, 1);
                REQUIRE( str.lowPaulis  == ((targ <  32)? getPow2(2*targ)      : 0) );
                REQUIRE( str.highPaulis == ((targ >= 32)? getPow2(2*(targ-32)) : 0) );
            }

            SECTION( "from only string" ) {

                CAPTURE( targets, pauliChars );

                char chars[65]; // 64 + terminal char
                chars[64] = '\0';
                for (int i=0; i<64; i++)
                    chars[i] = 'I';

                // string is DECREASING significance
                for (int i=0; i<numPaulis; i++)
                    chars[64-targets[i]-1] = pauliChars[i];

                std::string all = std::string(chars);
                PauliStr str = getPauliStr(all);

                CAPTURE( all );
                REQUIRE( str.lowPaulis  == lowValue );
                REQUIRE( str.highPaulis == highValue );
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "invalid target" ) {

            int target = GENERATE( -1, 64, 65, 9999 );

            REQUIRE_THROWS_WITH( getPauliStr("X", {target}), ContainsSubstring("Invalid index") );
        }

        SECTION( "duplicated target" ) {

            REQUIRE_THROWS_WITH( getPauliStr("XY", {0,0}), ContainsSubstring("duplicate") );
        }

        SECTION( "invalid number of paulis" ) {

            int numPaulis = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( getPauliStr("X", nullptr, numPaulis), ContainsSubstring("must contain at least one Pauli operator") );
        }

        SECTION( "string terminated early" ) {

            REQUIRE_THROWS_WITH( getPauliStr("X", {1,2}), ContainsSubstring("different number of Pauli operators") && ContainsSubstring("qubit indices") );
        }

        SECTION( "unrecognised char" ) {

            REQUIRE_THROWS_WITH( getPauliStr("hi", {1,2}), ContainsSubstring("unrecognised Pauli character") );
        }
    }
}


TEST_CASE( "getInlinePauliStr", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        PauliStr str;
        long long unsigned val;
        
        val = 1 + 2*4 + 3*4*4 + 0*4*4*4;
        str = getInlinePauliStr("XYZI", {0,1,2,3});
        REQUIRE( str.lowPaulis == val );
        str = getInlinePauliStr("xyzi", {0,1,2,3});
        REQUIRE( str.lowPaulis == val );
        str = getInlinePauliStr("1230", {0,1,2,3});
        REQUIRE( str.lowPaulis == val );

        val = (1*4*4*4 + 2*4*4 + 3*4 + 0*4);
        str = getInlinePauliStr("XYZI",  {3,2,1,0});
        REQUIRE( str.lowPaulis == val );
        str = getInlinePauliStr("xyzi",  {3,2,1,0});
        REQUIRE( str.lowPaulis == val );
        str = getInlinePauliStr("1230",  {3,2,1,0});
        REQUIRE( str.lowPaulis == val );

        val *= (1ULL << (2 * (63 - 3 - 32)));
        str = getInlinePauliStr("XYZ", {63,62,61});
        REQUIRE( str.highPaulis == val );
        str = getInlinePauliStr("xyz", {63,62,61});
        REQUIRE( str.highPaulis == val );
        str = getInlinePauliStr("123", {63,62,61});
        REQUIRE( str.highPaulis == val );
    }

    SECTION( LABEL_VALIDATION ) {

        // here, only the C++ interface can be tested

        SECTION( "invalid target" ) {

            int target = GENERATE( -1, 64, 65, 9999 );

            REQUIRE_THROWS_WITH( getInlinePauliStr("X", {target}), ContainsSubstring("Invalid index") );
        }

        SECTION( "duplicated target" ) {

            REQUIRE_THROWS_WITH( getInlinePauliStr("XY", {0,0}), ContainsSubstring("duplicate") );
        }

        SECTION( "string terminated early" ) {

            REQUIRE_THROWS_WITH( getInlinePauliStr("X", {1,2}), ContainsSubstring("different number of Pauli operators") && ContainsSubstring("qubit indices") );
        }

        SECTION( "unrecognised char" ) {

            REQUIRE_THROWS_WITH( getInlinePauliStr("ABC", {1,2,3}), ContainsSubstring("unrecognised Pauli character") );
        }
    }
}


TEST_CASE( "createPauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // the utilities for automatically checking
        // correctness of this function would "beg the
        // question", so we instead compare to a
        // hardcoded problem solved externally

        int numQubits = 2;

        vector<PauliStr> strings = {
            getPauliStr("XY", {0,1}),
            getPauliStr("ZX", {0,1}),
            getPauliStr("YZ", {0,1})
        };
        vector<qcomp> coeffs = {1, 2, 3_i};

        qmatrix ref = {
            {  0,   3,    2, -1_i},
            { -3,   0, -1_i,   -2},
            {  2, 1_i,    0,   -3},
            {1_i,  -2,    3,    0}};

        SECTION( LABEL_C_INTERFACE ) {

            PauliStrSum sum = createPauliStrSum(strings.data(), coeffs.data(), coeffs.size());

            REQUIRE( sum.numTerms == (qindex) strings.size() );
            REQUIRE( *(sum.isApproxHermitian) == -1 );

            REQUIRE_AGREE( getMatrix(sum,numQubits), ref );

            destroyPauliStrSum(sum);
        }

        SECTION( LABEL_CPP_INTERFACE ) {

            PauliStrSum sum = createPauliStrSum(strings, coeffs);

            REQUIRE( sum.numTerms == (qindex) strings.size() );
            REQUIRE( *(sum.isApproxHermitian) == -1 );

            REQUIRE_AGREE( getMatrix(sum,numQubits), ref );

            destroyPauliStrSum(sum);
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "number of terms" ) {

            int numTerms = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( createPauliStrSum(nullptr, nullptr, numTerms), ContainsSubstring("number of terms must be a positive integer") );
        }

        SECTION( "mismatching lengths" ) {

            // specific to the C++ interface

            REQUIRE_THROWS_WITH( createPauliStrSum({}, {.1}), ContainsSubstring("different number of Pauli strings") && ContainsSubstring("coefficients") );
        }
    }
}


TEST_CASE( "createInlinePauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // C++ interface does not mandate literal

        SECTION( "pauli parsing" ) {

            // ZYXI
            unsigned ref = 0 + 1*4 + 2*4*4 + 3*4*4*4;
            auto str = GENERATE(
                ".1i ZYXI",  ".1I Z Y X I",  ".1j   ZY   X I  ", 
                ".1i zyxi",  ".1I z y x i",  ".1j   zy   x i  ", 
                ".1i 3210",  ".1I 3 2 1 0",  ".1J   32   1 0  ",
                ".1i Zy1i",  ".1I 3 Y X 0",  ".1J   Zy   1 I  "
            );

            PauliStrSum sum = createInlinePauliStrSum(str);
            REQUIRE( sum.strings[0].lowPaulis == ref );
            destroyPauliStrSum(sum);
        }

        SECTION( "coefficient parsing" ) {

            vector<std::string> strs = {"1 X", "0 X", "0.1 X", "5E2-1i X", "-1E-50i X",  "1 - 6E-5i X", "-1.5E-15  -   5.123E-30i  0"};
            vector<qcomp> coeffs     = { 1,     0,     0.1,     5E2-1_i,   -(1E-50)*1_i,  1 -(6E-5)*1_i, qcomp(-1.5E-15, -5.123E-30) };

            size_t i = GENERATE_REF( range(0, (int) strs.size()) );
            CAPTURE( strs[i], coeffs[i] );

            PauliStrSum sum = createInlinePauliStrSum(strs[i]);
            REQUIRE_AGREE( sum.coeffs[0], coeffs[i] ); // should be strict
            destroyPauliStrSum(sum);
        }
        
        SECTION( "newlines" ) {

            PauliStrSum sum = createInlinePauliStrSum(R"(
                + 5E2-1i     XYZ 
                - 1E-50i     IXY 
                + 1 - 6E-5i  IIX 
                  0          III 
                  5.         XXX 
                  .5         ZZZ 
            )");

            REQUIRE( sum.numTerms == 6 );
            REQUIRE( sum.strings[3].lowPaulis == 0ULL );
            REQUIRE_AGREE( sum.coeffs[5], qcomp(.5,0) ); // should be strict

            destroyPauliStrSum(sum);
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not init" ) {

            // no way to test this
            SUCCEED( );
        }

        SECTION( "empty" ) {

            auto str = GENERATE( "", " ", "\n", " \n " );

            REQUIRE_THROWS_WITH( createInlinePauliStrSum(str), ContainsSubstring("empty") );
        }

        SECTION( "uninterpretable" ) {

            auto str = GENERATE( "X", "1", "a X", "-1 H", "0 .3", "1 23456" );

            REQUIRE_THROWS_WITH( createInlinePauliStrSum(str), ContainsSubstring("Could not interpret") );

            REQUIRE_NOTHROW( createInlinePauliStrSum("1 2 3") ); // = 1 * YZ and is legal
        }

        SECTION( "inconsistent number of qubits" ) {

            REQUIRE_THROWS_WITH( createInlinePauliStrSum("3 XYZ \n 2 YX"), ContainsSubstring("inconsistent") );
        }

        SECTION( "too many qubits" ) {

            // C++ interface permits both literal AND passing existing string
            std::string str = "1 XXXXXYYYYYZZZZZIIIIIXXXXXYYYYYZZZZZIIIIIXXXXXYYYYYZZZZZIIIIIXXXXX"; // 65 paulis

            REQUIRE_THROWS_WITH( createInlinePauliStrSum(str), ContainsSubstring("exceeds the maximum of 64") );
        }
    }
}


TEST_CASE( "createPauliStrSumFromFile", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        std::string fn = "test.txt";
    
        // file contents can be identical to createInlinePauliStrSum input above
        if (getQuESTEnv().rank == 0) {
            std::ofstream file;
            file.open(fn);
            file << R"(
                + 5E2-1i     XYZ 
                - 1E-50i     IXY 
                + 1 - 6E-5i  IIX 
                0            III 
                5.           IXX
                .5           ZYX 
            )";
            file.close();
        }

        // all nodes must wait for root to finish writing
        syncQuESTEnv();
    
        PauliStrSum sum = createPauliStrSumFromFile(fn);

        REQUIRE( sum.strings[0].lowPaulis == 3 + 2*4 + 1*4*4 );

        REQUIRE( sum.numTerms == 6 );
        REQUIRE( sum.coeffs[0] == qcomp(500, -1) );
        REQUIRE( sum.strings[3].lowPaulis == 0ULL );
        REQUIRE( sum.coeffs[5] == qcomp(.5,0) );

        destroyPauliStrSum(sum);
    }

    SECTION( LABEL_VALIDATION ) {

        // we skip all tests which overlap createInlinePauliStrSum
        // above, since the function simplify parses the file then
        // calls createInlinePauliStrSum(), and writing everything
        // to file to subsequently test this function is a chore

        SECTION( "bad file name" ) {

            auto fn = GENERATE( "", " ", "\n", "nonexistentfile.txt" );

            REQUIRE_THROWS_WITH( createPauliStrSumFromFile(fn), ContainsSubstring("Could not load and read the given file") );
        }
    }
}


TEST_CASE( "createPauliStrSumFromReversedFile", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        std::string fn = "test.txt";
    
        // file contents can be identical to createInlinePauliStrSum input above
        if (getQuESTEnv().rank == 0) {
            std::ofstream file;
            file.open(fn);
            file << R"(
                + 5E2-1i     XYZ 
                - 1E-50i     IXY 
                + 1 - 6E-5i  IIX 
                0            III 
                5.           IXX
                .5           ZYX 
            )";
            file.close();
        }

        // all nodes must wait for root to finish writing
        syncQuESTEnv();
    
        PauliStrSum sum = createPauliStrSumFromReversedFile(fn);

        // reversed order from createPauliStrSumFromFile() above
        REQUIRE( sum.strings[0].lowPaulis == 1 + 2*4 + 3*4*4 );

        REQUIRE( sum.numTerms == 6 );
        REQUIRE( sum.coeffs[0] == qcomp(500, -1) );
        REQUIRE( sum.strings[3].lowPaulis == 0ULL );
        REQUIRE( sum.coeffs[5] == qcomp(.5,0) );

        destroyPauliStrSum(sum);
    }

    SECTION( LABEL_VALIDATION ) {

        // we skip all tests which overlap createInlinePauliStrSum
        // above, since the function simplify parses the file then
        // calls createInlinePauliStrSum(), and writing everything
        // to file to subsequently test this function is a chore

        SECTION( "bad file name" ) {

            auto fn = GENERATE( "", " ", "\n", "nonexistentfile.txt" );

            REQUIRE_THROWS_WITH( createPauliStrSumFromReversedFile(fn), ContainsSubstring("Could not load and read the given file") );
        }
    }
}


TEST_CASE( "destroyPauliStrSum", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        PauliStrSum sum = createInlinePauliStrSum("1 X");
        REQUIRE_NOTHROW( destroyPauliStrSum(sum) );
    }

    SECTION( LABEL_VALIDATION ) {

        /// @todo fails in MSVC for unknown reason
        #ifndef _MSC_VER
        // sanitizer messes with default initialisation
        #ifndef SANITIZER_IS_ACTIVE
        SECTION( "not created" ) {

            PauliStrSum sum;

            // uninitialised sum fields can be coincidentally
            // valid on some platforms (Github Actions linux
            // gcc), so we force invalidity
            sum.numTerms = -1;
            sum.isApproxHermitian = nullptr; // else seg-faults if field accessed

            REQUIRE_THROWS_WITH( destroyPauliStrSum(sum), 
                ContainsSubstring("invalid fields") || 
                ContainsSubstring("heap pointers was unexpectedly NULL") ||
                ContainsSubstring("It is likely the structure was not created by its proper function")
            );
        }
        #endif
        #endif
    }
}


/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */


void reportPauliStr(PauliStr str);

void reportPauliStrSum(PauliStrSum str);
