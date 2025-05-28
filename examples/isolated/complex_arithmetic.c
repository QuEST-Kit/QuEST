/** @file
 * Examples of using QuEST's precision-agnostic
 * overloaded arithmetic operators between qcomp 
 * and other types like ints and floats, in C11.
 * MSVC does not support C complex arithmetic.
 * 
 * @author Tyson Jones
*/

#include "quest.h"


int main() {

    initQuESTEnv();

#if !defined(_MSC_VER)

    qcomp x = 1.2 + 3.4i;

    // C + R
    x = x + (int) 3;
    x = x + (qindex) 3;
    x = x + (float) 2;
    x = x + (double) 2;
    x = x + (long double) 2;

    // R + C
    x = (int) 3         + x;
    x = (qindex) 3      + x;
    x = (float) 2       + x;
    x = (double) 2      + x;
    x = (long double) 2 + x;

    // C - R
    x = x - (int) 3;
    x = x - (qindex) 3;
    x = x - (float) 2;
    x = x - (double) 2;
    x = x - (long double) 2;

    // R - C
    x = (int) 3         - x;
    x = (qindex) 3      - x;
    x = (float) 2       - x;
    x = (double) 2      - x;
    x = (long double) 2 - x;

    // C * R
    x = x * (int) 3;
    x = x * (qindex) 3;
    x = x * (float) 2;
    x = x * (double) 2;
    x = x * (long double) 2;

    // R * C
    x = (int) 3         * x;    
    x = (qindex) 3      * x;
    x = (float) 2       * x;
    x = (double) 2      * x;
    x = (long double) 2 * x;

    // C / R
    x = x / (int) 3;
    x = x / (qindex) 3;
    x = x / (float) 2;
    x = x / (double) 2;
    x = x / (long double) 2;

    // R / C
    x = (int) 3         / x;
    x = (qindex) 3      / x;
    x = (float) 2       / x;
    x = (double) 2      / x;
    x = (long double) 2 / x;

    // C += R
    x += (int) 3;
    x += (qindex) 3;
    x += (float) 2;
    x += (double) 2;
    x += (long double) 2;

    // C -= R
    x -= (int) 3;
    x -= (qindex) 3;
    x -= (float) 2;
    x -= (double) 2;
    x -= (long double) 2;

    // C *= R
    x *= (int) 3;
    x *= (qindex) 3;
    x *= (float) 2;
    x *= (double) 2;
    x *= (long double) 2;

    // C /= R
    x /= (int) 3;
    x /= (qindex) 3;
    x /= (float) 2;
    x /= (double) 2;
    x /= (long double) 2;

    reportScalar("x", x);
#endif

    finalizeQuESTEnv();
    return 0;
}