/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#ifndef FORMATTER_HPP
#define FORMATTER_HPP

#include "quest/include/structures.h"

#include <vector>
#include <tuple>
#include <complex>
#include <string>
#include <sstream>
#include <type_traits>



/*
 * TYPE NAME STRINGS
 */

std::string form_getQrealType();

std::string form_getQcompType();

std::string form_getQindexType();

std::string form_getFloatPrecisionFlag();



/*
 * STRING CASTS
 *
 * which are aggravatingly necessarily defined here in the header
 * because of their use of templating.
 */

template<typename T> std::string form_str(T expr) {

    // note that passing qcomp (and complex<T> generally) will see them
    // formatted as tuples (re, im). they should instead be passed to
    // compToStr() which is currently a private member of formatter.cpp.
    // I tried and failed to conjure the templating magic necessary to
    // automatically invoke it when T is a complex<K>.

    // write to buffer (rather than use to_string()) so that floating-point numbers
    // are automatically converted to scientific notation when necessary
    std::ostringstream buffer;
    buffer << expr;
    return buffer.str();
}



/*
 * SUBSTRINGS USED BY REPORTERS
 */

namespace form_substrings { 
    extern std::string eq;
    extern std::string mu;
    extern std::string by;
    extern std::string pn;
    extern std::string pg;
    extern std::string pm;
    extern std::string bt;
    extern std::string na;
    extern std::string un;
}



/*
 * MATRIX PRINTING
 */

void form_printMatrix(CompMatr1 matr, std::string indent="    ");

void form_printMatrix(CompMatr2 matr, std::string indent="    ");

void form_printMatrix(CompMatrN matr, std::string indent="    ");



/*
 * TABLE PRINTING
 */

void form_printTable(std::string title, std::vector<std::tuple<std::string, std::string>> rows, std::string indent="  ");

void form_printTable(std::string title, std::vector<std::tuple<std::string, long long int>> rows, std::string indent="  ");



#endif // FORMATTER_HPP