#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <map>
#include "define.hpp"

/* Print out vector of string in a line */
Str printStrVec(const StrVec& str_vec);

/* Get residue ranges from string */
Int2d retrieveRange(Str& inp_str);

/* Get selection context from string */
SeleContext retrieveSeleContext(StrQue inp_str);

/* Print out error message then return error code */
template<typename T>
T error_exit(const Str& inp_str, const T err_code) {
    std::cout << "ERROR> " << inp_str << std::endl;
    return err_code;
}

/* Print DEBUG info for simple data types */
template<typename T>
void DEBUG(const T& debug_info) {
    std::cout << "DEBUG> " << debug_info << std::endl;
}

/* Print program info */
void printProgName();

#endif