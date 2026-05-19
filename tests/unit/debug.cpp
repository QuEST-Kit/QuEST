/** @file
 * Unit tests of the debug module.
 *
 * @author Tyson Jones
 * 
 * @defgroup unitdebug Debug
 * @ingroup unittests
 */

#include "quest.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "tests/utils/macros.hpp"
#include "tests/utils/cache.hpp"
#include "tests/utils/convert.hpp"
#include "tests/utils/compare.hpp"
#include "tests/utils/random.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

using Catch::Matchers::ContainsSubstring;
using std::vector;



/*
 * UTILITIES
 */

#define TEST_CATEGORY \
    LABEL_UNIT_TAG "[debug]"



/** 
 * TESTS
 * 
 * @ingroup unitdebug
 * @{
 */


TEST_CASE( "setQuESTInputErrorHandler", TEST_CATEGORY ) {

    /// @todo
    /// We can test this by saving the current handler,
    /// overwriting it to local handlers of our choosing
    /// and verifying e.g. their error messages are detected,
    /// then restoring the original handler. Alas, this
    /// requires exposing the handler in main, else adding
    /// a quest function for obtaining the user-set handler.
    /// Both are annoying - I forego them presently! This
    /// function has anyway been manually tested via examples.
    
    SUCCEED( );
}


TEST_CASE( "setQuESTMaxNumReportedSigFigs", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        qcomp scalar = getQcomp(0.12345, 0.12345);

        vector<std::string> refs = {
            "0.1+0.1i",
            "0.12+0.12i",
            "0.123+0.123i",
            "0.1235+0.1235i", // rounded
            "0.12345+0.12345i"
        };

        // disable auto \n after lines
        setQuESTNumReportedNewlines(0);

        for (size_t numSigFigs=1; numSigFigs<=refs.size(); numSigFigs++) {

            setQuESTMaxNumReportedSigFigs(numSigFigs);

            // redirect stdout to buffer
            std::stringstream buffer;
            std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
            reportScalar("", scalar);
            std::cout.rdbuf(old);
            std::string out = buffer.str();

            std::string ref = refs[numSigFigs-1];

            CAPTURE( numSigFigs, ref );
            REQUIRE( out == ref );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "number" ) {

            int num = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( setQuESTMaxNumReportedSigFigs(num), ContainsSubstring("Cannot be less than one") );
        }
    }

    // restore to QuEST default for future tests
    setQuESTMaxNumReportedSigFigs(5);
}


TEST_CASE( "setQuESTNumReportedNewlines", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        for (int numNewlines=0; numNewlines<3; numNewlines++) {

            setQuESTNumReportedNewlines(numNewlines);

            // redirect stdout to buffer
            std::stringstream buffer;
            std::streambuf * old = std::cout.rdbuf(buffer.rdbuf());
            reportStr("x");
            std::cout.rdbuf(old);
            std::string out = buffer.str();

            std::string ref = "x" + std::string(numNewlines, '\n');
            
            CAPTURE( numNewlines, ref );
            REQUIRE( out == ref );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "number" ) {

            REQUIRE_THROWS_WITH( setQuESTNumReportedNewlines(-1), ContainsSubstring("Cannot generally be less than zero") );
        }

        SECTION( "multine number" ) {

            setQuESTNumReportedNewlines(0);

            REQUIRE_THROWS_WITH( reportQuESTEnv(), ContainsSubstring("zero") && ContainsSubstring("not permitted when calling multi-line") );
        }
    }

    // restore to QuEST default for future tests
    setQuESTNumReportedNewlines(2);
}


TEST_CASE( "setQuESTSeeds", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // perform the below test for every combination
        // of deployments, since their RNGs differ

        for (auto& [label, qureg]: getCachedDensmatrs()) {

            DYNAMIC_SECTION( label ) {

                SECTION( "same seed consistency" ) {

                    unsigned seeds[] = {123, 543, 755};
                    const int numSeeds = 3;
                    const int numMixedStates = 10;
                    const int numReps = 5;

                    // set an arbitrary fixed seed...
                    setQuESTSeeds(seeds, numSeeds);

                    // generate and remember a random state
                    initRandomMixedState(qureg, numMixedStates);
                    qmatrix ref = getMatrix(qureg);

                    // get a set of random measurement outcomes
                    vector<int> outcomes(qureg.numQubits);
                    for (int i=0; i<qureg.numQubits; i++)
                        outcomes[i] = applyQubitMeasurement(qureg, i);

                    // repeatedly...
                    for (int r=0; r<numReps; r++) {

                        // reset the seed
                        setQuESTSeeds(seeds, numSeeds);

                        // and confirm all random states are re-produced
                        initRandomMixedState(qureg, numMixedStates);
                        REQUIRE_AGREE( qureg, ref);

                        // as are all measurement outcomes
                        for (int i=0; i<qureg.numQubits; i++)
                            REQUIRE( outcomes[i] == applyQubitMeasurement(qureg,i) );
                    }
                }

                SECTION( "different key inconsistency" ) {

                    unsigned seeds[] = {123, 543, 755};
                    const int numSeeds = 3;
                    const int ampInd = 0;

                    // set arbitrary seed and collect random-state amp
                    setQuESTSeeds(seeds, numSeeds);
                    initRandomPureState(qureg);
                    qcomp amp1 = getDensityQuregAmp(qureg, ampInd, ampInd);

                    // change one passed seed and re-collect random-state amp
                    int i = GENERATE_COPY( range(0,numSeeds) );
                    seeds[i] = 987654321;
                    setQuESTSeeds(seeds, numSeeds);
                    initRandomPureState(qureg);
                    qcomp amp2 = getDensityQuregAmp(qureg, ampInd, ampInd);

                    // amps differ if seed successfully impacts outcomes
                    REQUIRE( amp1 != amp2 );
                }
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // no way to test this
            SUCCEED( );
        }

        SECTION( "number of seeds" ) {

            int numSeeds = GENERATE( -1, 0 );

            REQUIRE_THROWS_WITH( setQuESTSeeds(nullptr, numSeeds), ContainsSubstring("Invalid number of random seeds") );
        }

        // inconsistency between nodes is permitted
    }

    // re-randomise seeds for remaining tests
    setQuESTSeedsToDefault();
}


TEST_CASE( "setQuESTSeedsToDefault", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // perform the below test for every combination
        // of deployments, since their RNGs differ

        for (auto& [label, qureg]: getCachedDensmatrs()) {

            DYNAMIC_SECTION( label ) {

                SECTION( "different key inconsistency" ) {

                    const int ampInd = 0;

                    // randomise seed and collect random-state amp
                    setQuESTSeedsToDefault();
                    initRandomPureState(qureg);
                    qcomp amp1 = getDensityQuregAmp(qureg, ampInd, ampInd);

                    // re-randomise seed and collect new random-state amp
                    setQuESTSeedsToDefault();
                    initRandomPureState(qureg);
                    qcomp amp2 = getDensityQuregAmp(qureg, ampInd, ampInd);

                    // amps differ if seed re-randomisation impacts outcomes
                    REQUIRE( amp1 != amp2 );
                }
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // no way to test this
            SUCCEED( );
        }
    }

    // re-randomise seeds for remaining tests
    setQuESTSeedsToDefault();
}


TEST_CASE( "getQuESTSeeds", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( "can be called immediately" ) {

            REQUIRE_NOTHROW( getQuESTNumSeeds() );

            int numSeeds = getQuESTNumSeeds();
            vector<unsigned> out(numSeeds);

            REQUIRE_NOTHROW( getQuESTSeeds(out.data()) );
        }

        SECTION( "correct output" ) {

            GENERATE( range(0,10) );

            // prepare random seeds (using test utils RNG, not QuEST's)
            int numSeeds = getRandomInt(1, 10);
            vector<unsigned> in(numSeeds);
            for (int i=0; i<numSeeds; i++)
                in[i] = static_cast<unsigned>(getRandomInt(0, 99999));

            // pass seeds to QuEST
            setQuESTSeeds(in.data(), numSeeds);

            // check we get them back
            vector<unsigned> out(numSeeds);
            getQuESTSeeds(out.data());
            for (int i=0; i<numSeeds; i++)
                REQUIRE( in[i] == out[i] );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // no way to test this
            SUCCEED( );
        }
    }

    // re-randomise seeds for remaining tests
    setQuESTSeedsToDefault();
}


TEST_CASE( "getQuESTNumSeeds", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( "can be called immediately" ) {

            REQUIRE_NOTHROW( getQuESTNumSeeds() );
        }

        SECTION( "correct output" ) {

            GENERATE( range(0,10) );

            // prepare random seeds (using test utils RNG, not QuEST's)
            int numSeeds = getRandomInt(1, 10);
            vector<unsigned> in(numSeeds);
            for (int i=0; i<numSeeds; i++)
                in[i] = static_cast<unsigned>(getRandomInt(0, 99999));

            // pass seeds to QuEST
            setQuESTSeeds(in.data(), numSeeds);

            // confirm we get out correct number
            REQUIRE( getQuESTNumSeeds() == numSeeds );
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "env not initialised" ) {

            // no way to test this
            SUCCEED( );
        }
    }

    // re-randomise seeds for remaining tests
    setQuESTSeedsToDefault();
}


TEST_CASE( "setQuESTValidationOn", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // always safe to call
        for (int i=0; i<3; i++)
            REQUIRE_NOTHROW( setQuESTValidationOn() );

        // illegal and caught
        REQUIRE_THROWS( setQuESTSeeds(nullptr, -99) );
    }

    SECTION( LABEL_VALIDATION ) {

        // has no validation
        SUCCEED( );
    }
}


TEST_CASE( "setQuESTValidationOff", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // confirm always safe to call
        for (int i=0; i<3; i++)
            REQUIRE_NOTHROW( setQuESTValidationOff() );

        // prepare non-unitary matrix
        CompMatr1 m = getCompMatr1({{1,2},{3,4}});
        Qureg qureg = createQureg(1);

        // confirm not-unitary error suppressed...
        REQUIRE_NOTHROW( applyCompMatr1(qureg, 0, m) );

        // which otherwise triggers
        setQuESTValidationOn();
        REQUIRE_THROWS( applyCompMatr1(qureg, 0, m) );

        destroyQureg(qureg);
    }

    SECTION( LABEL_VALIDATION ) {

        // has no validation
        SUCCEED( );
    }

    // ensure validation is on for remaining tests
    setQuESTValidationOn();
}


TEST_CASE( "setQuESTValidationEpsilon", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( "affects validation" ) {

            Qureg qureg = createQureg(5);

            SECTION( "unitarity" ) {

                // prepare non-unitary matrix
                CompMatr1 m = getCompMatr1({{1,2},{3,4}});
                
                // confirm it throws non-unitary error
                REQUIRE_THROWS( applyCompMatr1(qureg, 0, m) );

                // confirm setting = 0 disables epsilon errors...
                setQuESTValidationEpsilon(0);
                REQUIRE_NOTHROW( applyCompMatr1(qureg, 0, m) );

                // but does not disable absolute errors
                REQUIRE_THROWS( applyCompMatr1(qureg, -1, m) );

                // confirm non-zero (forgive all) works
                setQuESTValidationEpsilon(9999); // bigger than dist of m*conj(m) from identity squared
                REQUIRE_NOTHROW( applyCompMatr1(qureg, 0, m) );
            }

            /// @todo
            /// to be completely rigorous, we should test
            /// that unitarity, hermiticity, CPTPness and 
            /// non-zero-ness of all of CompMatr, DiagMatr 
            /// and FullStateDiagMatr are affected! This
            /// is quite a chore of course

            destroyQureg(qureg);
        }
        
        SECTION( "affects struct fields" ) {

            SECTION( "CompMatr" ) {

                CompMatr m = createCompMatr(1);
                *(m.isApproxUnitary)   = 1;
                *(m.isApproxHermitian) = 1;

                setQuESTValidationEpsilon(.1);
                REQUIRE( *(m.isApproxUnitary)   == -1 );
                REQUIRE( *(m.isApproxHermitian) == -1 );

                destroyCompMatr(m);
            }

            SECTION( "DiagMatr" ) {

                DiagMatr m = createDiagMatr(1);
                *(m.isApproxUnitary)   = 1;
                *(m.isApproxHermitian) = 0;
                *(m.isApproxNonZero)   = 1;

                setQuESTValidationEpsilon(.1);
                REQUIRE( *(m.isApproxUnitary)   == -1 );
                REQUIRE( *(m.isApproxHermitian) == -1 );
                REQUIRE( *(m.isApproxNonZero)   == -1 );

                destroyDiagMatr(m);
            }

            SECTION( "FullStateDiagMatr" ) {

                FullStateDiagMatr m = createFullStateDiagMatr(1);
                *(m.isApproxUnitary)   = 1;
                *(m.isApproxHermitian) = 0;
                *(m.isApproxNonZero)   = 1;

                setQuESTValidationEpsilon(.1);
                REQUIRE( *(m.isApproxUnitary)   == -1 );
                REQUIRE( *(m.isApproxHermitian) == -1 );
                REQUIRE( *(m.isApproxNonZero)   == -1 );

                destroyFullStateDiagMatr(m);
            }

            SECTION( "KrausMap" ) {

                KrausMap k = createKrausMap(1, 3);
                *(k.isApproxCPTP) = 1;

                setQuESTValidationEpsilon(.1);
                REQUIRE( *(k.isApproxCPTP) == -1 );

                destroyKrausMap(k);
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( "negative epsilon" ) {

            qreal eps = GENERATE( -0.5, -1, -100 );

            REQUIRE_THROWS_WITH( setQuESTValidationEpsilon(eps), ContainsSubstring("positive number") );
        }
    }

    // ensure validation epsilon is default for remaining tests
    setQuESTValidationEpsilonToDefault();
}


TEST_CASE( "getQuESTValidationEpsilon", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // confirm always safe to call
        for (int i=0; i<3; i++)
            REQUIRE_NOTHROW( getQuESTValidationEpsilon() ); // ignores output

        GENERATE( range(0,10) );

        // confirm set correctly
        qreal eps = getRandomReal(0, 99999);
        setQuESTValidationEpsilon(eps);

        REQUIRE( getQuESTValidationEpsilon() == eps );
    }

    SECTION( LABEL_VALIDATION ) {

        // has no validation
        SUCCEED( );
    }

    // ensure validation epsilon is default for remaining tests
    setQuESTValidationEpsilonToDefault();
}


TEST_CASE( "setQuESTValidationEpsilonToDefault", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        SECTION( "always safe to call" ) {

            for (int i=0; i<3; i++)
                REQUIRE_NOTHROW( setQuESTValidationEpsilonToDefault() );
        }

        SECTION( "affects validation" ) {

            // prepare non-unitary matrix
            CompMatr1 m = getCompMatr1({{1,2},{3,4}});
            Qureg qureg = createQureg(1);

            // confirm it throws non-unitary error
            REQUIRE_THROWS( applyCompMatr1(qureg, 0, m) );

            // confirm setting = 0 disables epsilon errors...
            setQuESTValidationEpsilon(0);
            REQUIRE_NOTHROW( applyCompMatr1(qureg, 0, m) );

            // which returns when stored to default
            setQuESTValidationEpsilonToDefault();
            REQUIRE_THROWS( applyCompMatr1(qureg, 0, m) );

            destroyQureg(qureg);
        }

        SECTION( "affects struct fields" ) {

            SECTION( "CompMatr" ) {

                CompMatr m = createCompMatr(1);
                *(m.isApproxUnitary)   = 1;
                *(m.isApproxHermitian) = 1;

                setQuESTValidationEpsilonToDefault();
                REQUIRE( *(m.isApproxUnitary)   == -1 );
                REQUIRE( *(m.isApproxHermitian) == -1 );

                destroyCompMatr(m);
            }

            SECTION( "DiagMatr" ) {

                DiagMatr m = createDiagMatr(1);
                *(m.isApproxUnitary)   = 1;
                *(m.isApproxHermitian) = 0;
                *(m.isApproxNonZero)   = 1;

                setQuESTValidationEpsilonToDefault();
                REQUIRE( *(m.isApproxUnitary)   == -1 );
                REQUIRE( *(m.isApproxHermitian) == -1 );
                REQUIRE( *(m.isApproxNonZero)   == -1 );

                destroyDiagMatr(m);
            }

            SECTION( "FullStateDiagMatr" ) {

                FullStateDiagMatr m = createFullStateDiagMatr(1);
                *(m.isApproxUnitary)   = 1;
                *(m.isApproxHermitian) = 0;
                *(m.isApproxNonZero)   = 1;

                setQuESTValidationEpsilonToDefault();
                REQUIRE( *(m.isApproxUnitary)   == -1 );
                REQUIRE( *(m.isApproxHermitian) == -1 );
                REQUIRE( *(m.isApproxNonZero)   == -1 );

                destroyFullStateDiagMatr(m);
            }

            SECTION( "KrausMap" ) {

                KrausMap k = createKrausMap(1, 3);
                *(k.isApproxCPTP) = 1;

                setQuESTValidationEpsilonToDefault();
                REQUIRE( *(k.isApproxCPTP) == -1 );

                destroyKrausMap(k);
            }
        }
    }

    SECTION( LABEL_VALIDATION ) {

        SECTION( LABEL_VALIDATION ) {

            // performs no validation
            SUCCEED( );
        }
    }

    // ensure validation epsilon is default for remaining tests
    setQuESTValidationEpsilonToDefault();
}


TEST_CASE( "getQuESTGpuCacheSize", TEST_CATEGORY ) {

    SECTION( LABEL_CORRECTNESS ) {

        // confirm cache begins empty
        clearQuESTGpuCache();
        REQUIRE( getQuESTGpuCacheSize() == 0 );

        // hackily detect cuQuantum
        char envStr[200];
        getQuESTEnvironmentString(envStr);
        bool usingCuQuantum = std::string(envStr).find("cuQuantum=0") == std::string::npos;

        // proceed only if we're ever using our own GPU cache
        if (getQuESTEnv().isGpuAccelerated && !usingCuQuantum) {

            // GPU cache created for >= 6 qubit matrices
            Qureg qureg = createCustomQureg(10, 0, 0,1,0); // gpu only
            CompMatr matr = createCompMatr(6); // always in gpu
            for (qindex i=0; i<matr.numRows; i++)
                matr.cpuElems[i][i] = 1;
            syncCompMatr(matr);

            // each control qubit halves the necessary cache size;
            // so we start with many controls and remove each in-turn,
            // expanding the GPU cache in every call
            int targs[] = {0,1,2,3,4,5};
            int ctrls[] = {6,7,8,9};
            qindex cacheSize = 0;

            for (int numCtrls=4; numCtrls>=0; numCtrls--) {

                // expand the cache
                applyMultiControlledCompMatr(qureg, ctrls, numCtrls, targs, 6, matr);

                // confirm it expanded, OR stayed the same, which happens when
                // the total number of simultaneous threads needed hits/exceeds
                // the number available in the hardware
                qindex newSize = getQuESTGpuCacheSize();
                CAPTURE( cacheSize, newSize );
                REQUIRE( newSize >= cacheSize );

                cacheSize = newSize;
            }

            destroyQureg(qureg);
            destroyCompMatr(matr);
        }
    }

    SECTION( LABEL_VALIDATION ) {

        // performs no validation
        SUCCEED( );
    }
}
 

/** @} (end defgroup) */



/**
 * @todo
 * UNTESTED FUNCTIONS
 */


void setQuESTMaxNumReportedItems(qindex numRows, qindex numCols);

void getQuESTEnvironmentString(char str[200]);

void setQuESTReportedPauliChars(const char* paulis);

void setQuESTReportedPauliStrStyle(int style);
