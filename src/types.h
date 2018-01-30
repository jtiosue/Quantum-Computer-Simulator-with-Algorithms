#ifndef TYPES_INCLUDE
#define TYPES_INCLUDE

#include <map>
#include <string>
#include <complex>
#include <vector>

typedef std::complex<double> amp;
typedef std::map<std::string, amp> state_map;
typedef std::vector<unsigned int> vec_int;
typedef std::vector<std::string> vec_states;

#endif