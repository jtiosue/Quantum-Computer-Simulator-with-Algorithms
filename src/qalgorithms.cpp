#include <cmath>
#include <algorithm>
#include "qalgorithms.h"
#include "rand.h"
#include "methods.h"

#include <iostream>

// QUANTUM FOURIER TRANSFORM

void QFT(Register *reg, unsigned int start, unsigned int end) {
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
	for (unsigned int j = start; j < end; j++) {
		reg->Hadamard(j);
		for (unsigned int k = 1; k < end - j; k++) {
			reg->ControlledPhaseShift(j + k, j, pi/double(1 << k)); // 1 << k is pow(2, k)
		}
	}
	for (unsigned int i = 0; i < floor((end-start) / 2.0); i++) reg->Swap(start+i, end-i-1);
}

// INVERSE QUANTUM FOURIER TRANSFORM

void IQFT(Register *reg, unsigned int start, unsigned int end) {
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
	///*
	for (unsigned int i = 0; i < floor((end-start) / 2.0); i++) reg->Swap(start+i, end-i-1);

	// note: can't use unsigned int's here because unsigned j=-1
	// is always greater than zero. NEED to convert end and start to int
	// in order to properly compare.
	for (int j = int(end) - 1; j >= int(start); j--) {
		for (int k = int(end)-j-1; k >= 1; k--) {
			// don't need to explicilty convert to unsigned int here, but might as well.
			reg->ControlledPhaseShift((unsigned int)(j + k), (unsigned int)j, -pi / double(1 << k)); // 1 << k is pow(2, k)
		}
		reg->Hadamard((unsigned int)j);
	}
	//*/
}


// FINDING THE PERIOD OF f(x) = a^x mod N USING THE QUANTUM ALGORITHM

unsigned int find_Shor_period(unsigned int a, unsigned int N, unsigned int depth_limit) {
	/*
	Find the period r of the function
		f(x) = a^x mod N
	using the quantum algorithm
	*/

	if (depth_limit <= 0) {
		printf("Reached maximum depth limit in period find. Returning 1.\n");
		return 1;
	}
		
	// q is the number of numbers register 1 can hold.
	// q must be at least 2*N so that even if the period r = N-1, the first register
	// can still hold enough numbers to have at least two x such that f(x0) = f(x1)
	// because x0 + r = x1. Most literature I have found says that we should initialize
	// q such that N^2 <= q <= 2N^2 so that q / r > N even when r = N - 1, thus
	// there will be at least N x such that f(x0) = f(x1) = ... = f(xN) where
	// xi = x0 + i*r. But for the code, I've found that q = 2*N works fine. With a smaller
	// q, we have a smaller probability of measuring the correct r, so we may have to recurse
	// through this function more. But simulating a quanutum register on a classical computer is
	// exponential, so recursing through is still more faster than simulating more qubits.
	unsigned int q = 2 * N;

	unsigned int L1 = floor(log2(q)) + 1; // number of qubits in register 1
	unsigned int L2 = floor(log2(N)) + 1; // number of qubits in register 2
	// printf("Initialized register 1 with %d qubits and register 2 with %d\n", L1, L2);

	// This is very important! I messed this up for a while. q is 2^L1.
	// What we set earlier was just to help pick L1.
	q = (unsigned int)(1 << L1); // equiv to q = (unsigned int)pow(2, L1);


	Register reg(L1 + L2);

	// Make equal superposition in register 1.
	for (unsigned int i = 0; i < L1; i++) reg.Hadamard(i);
	// Could have also just QFTed the first L1 qubits. Has same effect.


	auto f = [a, N, L1, L2](string s) {
		/*
		Given a state of the two registers |x>|0>, this function
		will make it |x>|a^x mod N>. The first register has L1 
		qubits, the second has L2 qubits.
		*/
		string t = s.substr(0, L1);
		unsigned int x = binary_to_base10(t);
		// use my defined mod_power function so that we don't overflow
		string res = base10_to_binary(mod_power(a, x, N));
		while (res.size() < L2) res = "0" + res;
		return t + res;
	};
	// This entangles the registers. Sends |x> to |f(x)>.
	// In our case, we define f so that this sends |x, y> to |x, f(x)>
	reg.apply_function(f);

	// Don't technically need to measure yet, I don't think.
	// Could just wait until the end, but this reduces the number
	// of states, and so the QFT will perform faster (on the computer,
	// in real life, I don't think this affects the speed).
	for (unsigned int i = L1; i < L1 + L2; i++) reg.measure(i); // measure register 2.

	// Quantum fourier transform the first register.
	QFT(&reg, 0, L1);

	// m will be an integer multiple of q / r with high prbability
	// where r is the period.
	unsigned int m = binary_to_base10(reg.measure().substr(0, L1)); // Measurement of register 1.

	if (m == 0) {
		// printf("Quantum period find failed; trying again\n");
		return find_Shor_period(a, N, depth_limit - 1);
	}

	// with high probability, m = lambda * q / r for some
	// integer lambda where r is the period. Find the period.
	// Let c = m / q. We know c = lambda / r. Expand c with
	// continuous fractions (see the paper) to find the 
	// denominator r. This is O(log N) time, I think, so not bad.
	
	unsigned int r = 0; // If r is still zero afte r the while loop, we've failed.
	double c = double(m) / double(q); // this equals lambda / r for some lambda. Find the denomitor = r.
	// printf("Beginning continuous fractions expansion to find denominator for %g\n", c);
	unsigned int a1, d0, d1, dtemp; // a1 will hold cf expansion coefficients. The d's are denominators
	a1 = floor(1 / c);
	d0 = 1; d1 = a1;

	// we know the denominator is not greater than q, because lam / r = m / q.
	// if c == 0 then we've completed the cf expansion.
	// if c > 1, then it means that at some point c was very small, basically zero
	// and somehow this caused c = 1.0 / c - floor(1 / c) to not be less than 1.
	// Basically, it means we finished the cf expansion and didn't find r.
	// This almost never happens when q >= N^2, but happens often when Q = 2*N.
	while (c != 0.0 && c <= 1.0 && d1 < q) {
		if (!(d1 & 1) && mod_power(a, d1, N) == 1) { r = d1; break; } // Make sure r is even. Found it.
		else if (2 * d1 < q && mod_power(a, 2 * d1, N) == 1) { r = 2 * d1; break; } // Try two times, which will be even.
		c = 1.0 / c - double(a1);
		a1 = floor(1.0 / c);
		dtemp = d1; d1 = d1 * a1 + d0; d0 = dtemp;
	}

	if (r == 0) {
		// printf("Quantum period find failed; trying again\n");
		return find_Shor_period(a, N, depth_limit - 1);
	}

	// We already made sure that r was even, so we're good.

	printf("Quantum period find found the period of %d^x mod %d to be %d\n", a, N, r);
	return r;
}

// FACTORIZTION

unsigned int Shor(unsigned int N, unsigned int depth_limit) { 
	/*
	Find a single factor of N. N must not be an integer
	power of a prime. See the notes for an explaination of
	how Shor's algorithm works. The quantum computation only
	occurs in the find_period function.

	depth_limit is set by default in header. This limits the recursive depth.

	set_srand() must be called before using this function.
	*/
	if (depth_limit <= 0) {
		printf("Reached maximum depth limit in Shor. Try again, or increase depth limit\n");
		return 1;
	}

	if (N % 2 == 0) return 2;
	unsigned int a = (unsigned int)(floor(get_rand()*(N-1)+1)); if (a == 1) a++; 
	unsigned int g = gcd(a, N);
	if (g != 1 && g != N) {
		printf("Completed Shor's algorithm classically. Found a factor of %d to be %d\n", N, g); 
        printf("But we want to solve it quantumly! So starting over...\n");
		// return g;
        return Shor(N, depth_limit);
	} 

	// printf("Using quantum period finding algorithm to find the period of %d ^ x mod %d\n", a, N);
	unsigned int r = find_Shor_period(a, N); unsigned int n = mod_power(a, r / 2, N);
	// if (r % 2 == 1 || n % N == N-1) return Shor(N, depth_limit-1); // start over
    
	unsigned int res = gcd(n - 1, N);
	if (res != 1 && res != N) {
		printf("Shor found a factor of %d to be %d\n", N, res);
		return res;
	}
	res = gcd(n + 1, N);
	if (res != 1 && res != N) {
		printf("Shor found a factor of %d to be %d\n", N, res);
		return res;
	}

	// printf("Shor did not find a factor; trying again\n");
	return Shor(N, depth_limit - 1);
}

// RIPPLE CARRY ADDITION

void HalfAdder(unsigned int qubit1, unsigned int qubit2, unsigned int carry, Register *s) {
	s->Toffoli(qubit1, qubit2, carry); // Carry (AND)
	s->ControlledNot(qubit1, qubit2); // Sum (XOR)
	// Sum is stored in qubit2.
}

void FullAdder(unsigned int qubit1, unsigned int qubit2, unsigned int carryin, unsigned int carryout, Register *s) {
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

	unsigned int num_bits = reg->num_qubits / 3; // number of bits per number.
	
	// Half adder on least significant qubit
	HalfAdder(num_bits - 1, 2 * num_bits - 1, 2 * num_bits, reg);

	for (unsigned int i = 2; i <= num_bits; i++)
		FullAdder(num_bits-i, 2*num_bits-i, 2*num_bits+i-2, 2*num_bits+i-1, reg);
}

unsigned int Add(unsigned int a, unsigned int b) {
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
	unsigned int num_bits = (unsigned int)(log2(max(a, b))) + 1; // num of bits per number.

	Register reg(num_bits * 3);

	// Set up the first section of qubits for the integer a
	// and the second section of qubits for the integer b.

	unsigned int ta, tb;
	for (unsigned int i = 1; i <= num_bits; i++) {

		// Get binary one or zero for this digit of a.
		if ((unsigned int)(1 << (num_bits - i)) > a) ta = 0; // 1 << x is pow(2, x)
		else {ta = 1; a = a % (unsigned int)(1 << (num_bits - i));}

		// Get binary one or zero for this digit of b.
		if ((unsigned int)(1 << (num_bits - i)) > b) tb = 0;
		else { tb = 1; b = b % (unsigned int)(1 << (num_bits - i)); }

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

// MODULAR ADDITION USING QFT AND IQFT

void ModAdd(Register *reg) {
	/*
	For each state in the register, transforms a state
		|x, y>
	into the state
		|(x+y) mod pow(2, reg->num_bits), z>
	where we don't care about z.

	x and y must be represented by the same number of bits, so
	reg must have an even number of bits.
	*/

	unsigned int n = reg->num_qubits / 2;

	QFT(reg, 0, n); // Apply to just the bits representing x.

	// Bits representing y act as controls.
	int power;
	for (unsigned int control = n; control < reg->num_qubits; control++) {
		power = control - n;
		for (int target = int(n) - 1; target >= int(2 * n - control) - 1; target--) {
			reg->ControlledPhaseShift(control, (unsigned int)target, pi / double(1 << power)); // 1 << power is pow(2, power)
			power--;
		}
	}

	IQFT(reg, 0, n); // Apply to just the bits representing x.
}

unsigned int ModAdd(unsigned int a, unsigned int b, unsigned int num_bits) {
	/*
	Compute (a + b) mod pow(2, num_bits) where num_bits is per number.
	We set up a Register to be in the state |a b>
	where a is reprented in binary with num_bits and likewise
	for b. The result is stored on the bits the represent a.
	*/

	if (num_bits < (unsigned int)(log2(max(a, b))) + 1) {
		printf("Not enough bits to compute %d + %d mod %d\n", a, b, (int)(1 << num_bits)); // 1 << num_bits is pow(2, num_bits)
		return 0;
	}

	Register reg(num_bits * 2);

	// Set up the first section of qubits for the integer a
	// and the second section of qubits for the integer b.

	unsigned int ta, tb;
	for (unsigned int i = 1; i <= num_bits; i++) {

		// Get binary one or zero for this digit of a.
		if ((unsigned int)(1 << (num_bits - i)) > a) ta = 0; // 1 << x is pow(2, x)
		else { ta = 1; a = a % (unsigned int)(1 << (num_bits - i)); }

		// Get binary one or zero for this digit of b.
		if ((unsigned int)(1 << (num_bits - i)) > b) tb = 0;
		else { tb = 1; b = b % (unsigned int)(1 << (num_bits - i)); }

		// Set them up in the register.
		if (ta) reg.PauliX(i - 1);
		if (tb) reg.PauliX(num_bits + i - 1);
	}

	ModAdd(&reg);

	// reg will now be in a single state.

	string s = reg.measure();

	return binary_to_base10(s.substr(0, num_bits));
}

// GROVER'S SEARCH

unsigned int Grover(unsigned int omega, unsigned int num_bits, bool verbose) {
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
		Instead of r.apply_gate(D, v), could do the following, which is more physically realistic I think.

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
