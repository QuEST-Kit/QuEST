#include "quest.h"

#include <stdio.h>
#include <iostream>
#include <complex>
#include <chrono>

using namespace std;
using namespace std::chrono;




#define TIME_FUNC(VAR, CALL) \
    start = high_resolution_clock::now(); \
    for (int n=0; n<NUM_REPS; n++) { \
        for (int t=0; t<numQubits; t++) { \
            for (int c=0; c<numQubits; c++) { \
                if (t==c) \
                    continue; \
                int s = (n+t+c)%2; \
                CALL; \
            } \
        } \
    } \
    stop = high_resolution_clock::now(); \
    auto VAR = duration_cast<microseconds>(stop - start).count();



const int NUM_REPS = 1000;


int main(int argc, char* argv[]) {

    if (argc != 2)
        printf("numqb\n");

    int numQubits = atoi(argv[1]);

    initQuESTEnv();
    Qureg qureg = createQureg(numQubits);
    CompMatr1 m = getCompMatr1({{.1i,.2i},{-4+.2i,8}});
    setValidationOff();
    
    //reportQureg(qureg);



    // try to eliminate warm-up effects
    for (int n=0; n<NUM_REPS; n++) {
        for (int t=0; t<numQubits; t++) {
            for (int c=0; c<numQubits; c++) {
                if (t==c)
                    continue;


                int s = (n+t+c)%2;
                oneStateCtrlGate(qureg, c, s, t, m);
                oneCtrlGate(qureg, c, t, m);
                noCtrlGate(qureg, t, m);


                // NEW_oneStateCtrlGate(qureg, c, s, t, m);
                // NEW_oneCtrlGate(qureg, c, t, m);
                // NEW_noCtrlGate(qureg, t, m);
            }
        }
    }


    auto start = high_resolution_clock::now();
    auto stop   = high_resolution_clock::now();



    // TIME_FUNC(dur_OLD_noCtrlGate, noCtrlGate(qureg, t, m); )
    // TIME_FUNC(dur_NEW_noCtrlGate, NEW_noCtrlGate(qureg, t, m); )

    // TIME_FUNC(dur_OLD_oneCtrlGate, oneCtrlGate(qureg, c, t, m); )
    // TIME_FUNC(dur_NEW_oneCtrlGate, NEW_oneCtrlGate(qureg, c, t, m); )

    // TIME_FUNC(dur_OLD_oneStateCtrlGate, oneStateCtrlGate(qureg, c, s, t, m); )
    // TIME_FUNC(dur_NEW_oneStateCtrlGate, NEW_oneStateCtrlGate(qureg, c, s, t, m); )


    // cout << "dur_OLD_noCtrlGate:\t" << dur_OLD_noCtrlGate << std::endl;
    // cout << "dur_NEW_noCtrlGate:\t" << dur_NEW_noCtrlGate << std::endl;
    // cout << std::endl;

    // cout << "dur_OLD_oneCtrlGate:\t" << dur_OLD_oneCtrlGate << std::endl;
    // cout << "dur_NEW_oneCtrlGate:\t" << dur_NEW_oneCtrlGate << std::endl;
    // cout << std::endl;

    // cout << "dur_OLD_oneStateCtrlGate:\t" << dur_OLD_oneStateCtrlGate << std::endl;
    // cout << "dur_NEW_oneStateCtrlGate:\t" << dur_NEW_oneStateCtrlGate << std::endl;
    // cout << std::endl;


    // start = high_resolution_clock::now();

    // for (int n=0; n<NUM_REPS; n++) {
    //     for (int t=0; t<numQubits; t++) {
    //         for (int c=0; c<numQubits; c++) {
    //             if (t==c)
    //                 continue;

    //             //int s = (n+t+c)%2;
    //             const int s = 1;
    //             oneStateCtrlGate(qureg, c, s, t, m);
    //         }
    //     }
    // }

    // stop = high_resolution_clock::now();
    // auto durB = duration_cast<microseconds>(stop - start).count();

    // cout << "durA (oneCtrlGate):      " << durA << endl;
    // cout << "durB (oneStateCtrlGate): " << durB << endl;


    destroyQureg(qureg);
    finalizeQuESTEnv();
    return 0;
}
