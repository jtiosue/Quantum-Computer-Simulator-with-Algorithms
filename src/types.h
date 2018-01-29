#ifndef TYPES_INCLUDE
#define TYPES_INCLUDE

#include <map>
#include <string>
#include <complex>
#include <vector>

using namespace std;

typedef complex<double> amp;
typedef map<string, amp> state_map;
typedef vector<int> vec_int;
typedef vector<string> vec_states;

#endif