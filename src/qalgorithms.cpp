
#include <cmath>
#include <algorithm>
#include <functional>
#include "quantum.h"
#include "qalgorithms.h"
#include "rand.h"


const double pi = acos(-1.0);

// Useful functions

int char_to_int(char c) {return c - '0';}

int gcd(unsigned int a, unsigned int b) {
	// Find the greatest common divisor between a and b.
	int t;
	while (b != 0) { t = b; b = a % b; a = t; }
	return a;
}

int binary_to_base10(string s) {
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
	int result = 0;
	for (int i = 0; i < s.length(); i++) {
		result ^= (char_to_int(s[i]) << (s.length() - i - 1));
	}
	return result;
}

// QUANTUM FOURIER TRANSFORM

void QFT(Register *reg, int start, int end) {
	/*
	Apply the quantum fourier transform to qubits
	start through end - 1. int end is noninclusive.

	start and end are set by default be zero in the header file.
	*/

	if (end == 0) end = reg->num_qubits;

	// WITH SINGLE UNITARY TRANSFORMATION
	/*
	int n = end - start; vec_int v;
	for (int i = start; i < end; i++) v.push_back(i);
	reg->apply_gate(Unitary::QFT(n), v);
	*/

	// REALISTICALLY, WITH SINGLE AND DOUBLE QUBIT GATES
	for (int j = start; j < end; j++) {
		reg->Hadamard(j);
		for (int k = 1; k < end - j; k++) {
			reg->ControlledPhaseShift(j+k, j, pi/pow(2, k));
		}
	}
	for (int i = 0; i < floor((end-start) / 2.0); i++) reg->Swap(start+i, end-i-1);
}

// INVERSE QUANTUM FOURIER TRANSFORM

void IQFT(Register *reg, int start, int end) {
	// Just run QFT backwords for qubits start through end - 1.
	// start and end are set by default be zero in the header file.

	if (end == 0) end = reg->num_qubits;

	// WITH SINGLE UNITARY TRANSFORMATION
	/*
	int n = end - start; vec_int v;
	for (int i = start; i < end; i++) v.push_back(i);
	reg->apply_gate(Unitary::IQFT(n), v);
	*/

	// REALISTCALLY, WITH SINGLE AND DOUBLE QUBIT GATES

	for (int i = 0; i < floor((end-start) / 2.0); i++) reg->Swap(start+i, end-i-1);

	for (int j = end-1; j >= start; j--) {
		for (int k = end-j-1; k >= 1; k--) {
			reg->ControlledPhaseShift(j + k, j, -pi / pow(2, k));
		}
		reg->Hadamard(j);
	}
}


// FACTORIZTION

int find_period_classical(std::function<int(int)> f) {
	/*
	classical algorithm (Brent's algorithm) to find the period r of the function
		int f(int x)
	*/

	int r = 1; int tortoise = 0; int hare = 1; int power = 1;

	while (tortoise != hare) {
		if (power == r) {tortoise = hare; power *= 2; r = 0;}
		hare = f(hare);
		r++;
	}
	return r;
}

int find_period_quantum(std::function<int(int)> f) {
	/*
	Quantum algorithm to find the period r of the function
		int f(int x)
	*/
	return 0;
}

int find_period(int a, int N) {
	/*
	algorithm to find the period r of the function
		f(x) = a^x mod N
	*/

	auto f = [a, N](int x) -> int { return int(pow(a, x)) % N; };
	return find_period_classical(f);
	// return find_period_quantum(f);
}

int Shor(unsigned int N, unsigned int depth_limit) { 
	/*
	Find a single factor of N. N must not be an integer
	power of a prime. See the notes for an explaination of
	how Shor's algorithm works. The quantum computation only
	occurs in the find_period function.

	depth_limit is set by default in header. This limits the recursive depth.

	set_srand() must be called before using this function.
	*/
	if (depth_limit <= 0) return 1;
	int a = int(floor(get_rand()*N)); int g = gcd(a, N);
	if (g != 1) return g; // we're done.
	int r = find_period(a, N); int n = pow(a, r / 2);
	if (r % 2 == 1 || n % N == N-1) return Shor(N, depth_limit-1); // start over
	return gcd(n - 1, N); // v = gcd(n+1, N); // we're done.
}

// RIPPLE CARRY ADDITION

void HalfAdder(int qubit1, int qubit2, int carry, Register *s) {
	s->Toffoli(qubit1, qubit2, carry); // Carry (AND)
	s->ControlledNot(qubit1, qubit2); // Sum (XOR)
	// Sum is stored in qubit2.
}

void FullAdder(int qubit1, int qubit2, int carryin, int carryout, Register *s) {
	s->Toffoli(qubit1, qubit2, carryout); // Carry (AND)
	s->ControlledNot(qubit1, qubit2); // Sum (XOR)
	s->Toffoli(qubit2, carryin, carryout);	// Carry (AND)
	s->ControlledNot(carryin, qubit2); // Sum (XOR)
	// Sum is stored in b
}

void Add(Register *reg) {
	/*
	Semi-classical ripple-carry adder.
	Uses only Controlled-NOT and Toffoli (CC-NOT) gates to compute the sum.

	If the register has 3*n qubits, this function computes the sum of the
	number respresented by the first n qubits (0 to n-1) and the number represented by
	the second n qubits (n to 2n-1). The extra n qubits (2n to 3n-1) are used to aid
	in the computation AND SHOULD BE SET TO ZERO. The last qubit (3*n-1) is the overflow digit.

	This is deterministic if the register is in a particular state,

		i.e. |101011000> is equivalent to 101 + 011 in binary or 5 + 3 
		in base 10 and will always evaluate to 8.

	but if the register is in a superposition of states, the result is probabilistic.

	NOTE that this method affects the inputed register, but the register will not be measured.
	*/

	int num_bits = reg->num_qubits / 3; // number of bits per number.
	
	// Half adder on least significant qubit
	HalfAdder(num_bits - 1, 2 * num_bits - 1, 2 * num_bits, reg);

	for (int i = 2; i <= num_bits; i++)
		FullAdder(num_bits-i, 2*num_bits-i, 2*num_bits+i-2, 2*num_bits+i-1, reg);
}

int Add(unsigned int a, unsigned int b) {
	/*
	Semi-classical ripple-carry adder.
	Uses only Controlled-NOT and Toffoli (CC-NOT) gates to compute the sum.

	Sets up a register so that when calling Add(register) (see the function above)
	it deterministically evaluates a+b correctly.

	Example:
		If a is 5 == 101 in binary
		and b is 12 == 1100 in binary
		
		then we set up a Register such that it is in the state
		|0101 1100 0000> with probability one. The last four qubits
		are used to aid in the computation that occurs in the Add(register)
		function above.
	*/
	int num_bits = int(log2(max(a, b))) + 1; // num of bits per number.

	Register reg(num_bits * 3);

	// Set up the first section of qubits for the integer a
	// and the second section of qubits for the integer b.

	int ta, tb;
	for (int i = 1; i <= num_bits; i++) {

		// Get binary one or zero for this digit of a.
		if (pow(2, num_bits - i) > a) ta = 0;
		else {ta = 1; a = a % int(pow(2, num_bits - i));}

		// Get binary one or zero for this digit of b.
		if (pow(2, num_bits - i) > b) tb = 0;
		else { tb = 1; b = b % int(pow(2, num_bits - i)); }

		// Set them up in the register.
		if (ta) reg.PauliX(i-1);
		if (tb) reg.PauliX(num_bits+i-1);
	}

	Add(&reg);

	string s = reg.measure();
	// our result is the qubits in the place of number b; i.e. qubits num_bits
	// through 2*num_bits - 1. But we also have an overflow at qubit 3*num_bits - 1.

	return binary_to_base10(s.substr(3 * num_bits - 1, 1) + s.substr(num_bits, num_bits));
}

void ModAdd(Register *reg) { // DOESN'T WORK. THE QFT AND IQFT ARE CORRECT.
	// https://arxiv.org/pdf/quant-ph/0008033.pdf
	// must be an even number of bits

	int n = reg->num_qubits / 2;

	QFT(reg, 0, n);
	
	int phase;
	
	for (int i = n - 1; i >= 0; i--) { // targets are qubit i
		phase = 0;
		for (int j = n + i; j >= n; j--) { // controls are j
			reg->ControlledPhaseShift(j, i, pi/pow(2, phase));
			phase++;
		}
	}

	IQFT(reg, 0, n);
}

int ModAdd(unsigned int a, unsigned int b, unsigned int num_bits) {
	// compute a + b mod pow(2, num_bits)
	// num_bits is per number.

	if (num_bits < int(log2(max(a, b))) + 1) {
		printf("Not enough bits to compute %d + %d mod %d\n", int(a), int(b), (int)pow(2, num_bits)); 
		return 0;
	}


	Register reg(num_bits * 2);

	// Set up the first section of qubits for the integer a
	// and the second section of qubits for the integer b.

	int ta, tb;
	for (int i = 1; i <= num_bits; i++) {

		// Get binary one or zero for this digit of a.
		if (pow(2, num_bits - i) > a) ta = 0;
		else { ta = 1; a = a % int(pow(2, num_bits - i)); }

		// Get binary one or zero for this digit of b.
		if (pow(2, num_bits - i) > b) tb = 0;
		else { tb = 1; b = b % int(pow(2, num_bits - i)); }

		// Set them up in the register.
		if (ta) reg.PauliX(i - 1);
		if (tb) reg.PauliX(num_bits + i - 1);
	}

	ModAdd(&reg);

	// reg.print_states();

	string s = reg.measure();

	return binary_to_base10(s.substr(0, num_bits));
}

// GROVER'S SEARCH

int Grover(unsigned int omega, unsigned int num_bits, bool verbose) {
	/*
	Perform a Grover search to find what omega is given the black box
	operator Uomega. Of course, here we know omega and make Uomega from
	that, but the principle is that we are given Uomega such that
		Uomega |x> = |x>   if x != omega
		Uomega |x> = -|x>  if x == omega
	and we want to find omega. This is the simplest appication of the
	Grover search.

	omega must be in the range 0 <= omega < pow(2, num_bits).
	If verbose is true, the register will be printed prior to the measurement.
	set_srand() must be called before calling this function.
	*/

	int N = pow(2, num_bits);
	if (omega >= N) { printf("%d is not enough bits for omega = %d\n", num_bits, omega); return 0; }

	// Make Uomega, our black box operator.
	Unitary Uomega = Unitary::Identity(N); Uomega[omega][omega] = -1.0;

	// The Grover diffusion operator, 2|s><s| - I, where |s>
	// is equal superposition of all states.
	Unitary D(N);
	for (int i = 0; i < D.dimension; i++) {
		for (int j = 0; j < D.dimension; j++) {
			D[i][j] = 2.0 / double(N);
			if (i == j) D[i][j] -= 1;
		}
	}

	// // Here I define Us such that Hadamard^n * Us * Hadamard^n is the
	// // Grover diffusion operator, where Hadamard^n means the Hadamard
	// // gate applied to each qubit.
	// Unitary Us = Unitary::Identity(N) * -1.0; Us[0][0] = 1.0;

	Register r(num_bits); vec_int v;

	// Begin with equal superposition.
	for (int i = 0; i < r.num_qubits; i++) { r.Hadamard(i); v.push_back(i); }

	// iterate O(sqrt(N)) times. where we stop is important!
	for (int _ = 0; _ < round(pi / (4.0*asin(1/sqrt(N)))-1.0/2.0); _++) {
		// Uomega operator is applied to the whole system
		r.apply_gate(Uomega, v);

		// Apply the Grover diffusion operator.
		/*
		Instead of r.apply_gate(D, v), could do the following

		for (int i = 0; i < r.num_qubits; i++) r.Hadamard(i);
		r.apply_gate(Us, v);
		for (int i = 0; i < r.num_qubits; i++) r.Hadamard(i);
		*/

		r.apply_gate(D, v);
	}

	// When printing states, you can see that the basis element
	// corresponding to omega has a much higher amplitude than
	// the rest of the basis elements.
	if(verbose) r.print_states();

	// Collapse the system. There is a high probability that we get
	// the basis element corresponding to omega.
	string s = r.measure(); // result

	return binary_to_base10(s);
}
