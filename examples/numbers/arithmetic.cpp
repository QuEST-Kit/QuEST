#include "quest.h"

int main() {

    initQuESTEnv();

    // use 1_i instead of 1i to make the complex literal precision agnostic
    qcomp x = 1.2 + 3.4_i;

    // C with non-same-prec C
    x *= 1.5i - 2i - x + 3.5*x / 10.5_i;

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

    reportQcomp(x);
    
    finalizeQuESTEnv();
}
