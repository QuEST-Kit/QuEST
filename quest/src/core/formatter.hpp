/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#ifndef FORMATTER_HPP
#define FORMATTER_HPP

#include <vector>
#include <tuple>
#include <string>
#include <sstream>
#include <type_traits>



/*
 * TYPE STRINGS
 */

std::string form_getQrealType();

std::string form_getQcompType();

std::string form_getQindexType();

std::string form_getFloatPrecisionFlag();



/*
 * STRING CASTS
 */

template<typename T> std::string form_str(T expr) {
    std::ostringstream buffer;

    // floats are printed in scientific notation, to avoid precision loss
    if constexpr (std::is_floating_point_v<T>)
        buffer << std::scientific;

    buffer << expr;
    return buffer.str();
}
}



/*
 * TABLE PRINTING
 */

void form_printTable(std::string title, std::vector<std::tuple<std::string, std::string>> rows, std::string indent="  ");

void form_printTable(std::string title, std::vector<std::tuple<std::string, long long int>> rows, std::string indent="  ");



#endif // FORMATTER_HPP