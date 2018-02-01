#include "methods.h"

// Useful functions

unsigned int char_to_int(char c) { return c - '0'; }

unsigned int gcd(unsigned int a, unsigned int b) {
	// Find the greatest common divisor between a and b.
	unsigned int t;
	while (b != 0) { t = b; b = a % b; a = t; }
	return a;
}

unsigned int binary_to_base10(std::string s) {
	/*
	Convert binary number to base 10 with bit shifts.

	This function also does the following:
	Find row of matrix corresponding to state
	i.e. for 3 qubits, the unitary matrix would be 8x8.
	The basis is, in order,
	{|000>, |001>, |010>, |011>, |100>, |101>, |110>, |111>}
	Thus, if we input the state "000", this function returns 0.
	if we input the state "101", this function returns 5.
	*/


	/* RECURSIVE WAY
	if (!state.length()) return 0;
	int i = 0;
	if (state[0] == '1') i = pow(2, state.length()) / 2;
	return i + binary_to_base10(state.substr(1, state.length() - 1));
	*/

	// BITSHIFT WAY
	unsigned int result = 0;
	for (unsigned int i = 0; i < s.length(); i++) {
		result ^= (char_to_int(s[i]) << (s.length() - i - 1));
	}
	return result;
}

std::string base10_to_binary(unsigned int x) {
	if (x == 0) return "0";
	std::string s = "";
	while (x > 0) {
		if (x % 2 == 0) s = "0" + s;
		else s = "1" + s;
		x /= 2;
	}
	return s;
}

unsigned int mod_power(unsigned int a, unsigned int x, unsigned int N) {
	/*
	For large x, doing pow(a, x) % N causes overflow. So instead we
	use this function to compute a^x mod N for large x.

	This function operates on the fact that x^y mod n = (x mod n)^y mod n
		(see the paper for proof).
	*/
	
	unsigned int res = 1;

	a = a % N;  // Update a if it is more than or equal to N

	while (x > 0) {
		// If x is odd, multiply a with result
		if (x & 1) res = (res*a) % N;

		// x must be even now
		x = x >> 1; // x = floor(x / 2)
		a = (a*a) % N;
	}
	return res;
}
