/** @file
 * String formatting functions, primarily used by reportQureg and reportQuESTEnv()
 */

#ifndef FORMATTER_HPP
#define FORMATTER_HPP

#include <vector>
#include <tuple>
#include <stdlib.h>



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

std::string form_str(qreal num);

template<typename T>
std::string form_str(T expr) {
	return std::to_string(expr);
}



/*
 * TABLE PRINTING
 */

void form_printTable(std::string title, std::vector<std::tuple<std::string, std::string>> rows, std::string indent="  ");

void form_printTable(std::string title, std::vector<std::tuple<std::string, long long int>> rows, std::string indent="  ");



#endif // FORMATTER_HPP