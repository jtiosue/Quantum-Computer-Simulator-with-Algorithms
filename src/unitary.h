#ifndef UNITARY_INCLUDE
#define UNITARY_INCLUDE

#include "types.h"

class Unitary {

	private:
		amp **matrix;

	public:
		int dimension;
		Unitary(int dimension);
		// ~Unitary();
		amp * operator[](int i);
		Unitary operator*(amp x);
		friend std::ostream &operator<<(std::ostream &os, Unitary &U);

		static Unitary Identity(int dimension);
		static Unitary Hadamard();
		static Unitary PauliX();
		static Unitary PauliY();
		static Unitary PauliZ();
		static Unitary PhaseShift(double theta);
		static Unitary PiOverEight();
		static Unitary ControlledNot();
		static Unitary Toffoli();
		static Unitary ControlledPhaseShift(double theta);
		static Unitary Swap();
		static Unitary Ising(double theta);
		static Unitary QFT(int num_qubits);
		static Unitary IQFT(int num_qubits);

		// Matrices to be created upon initialization.
		// const static Unitary MHadamard;

};

// For left multiplication.
Unitary operator*(amp x, Unitary &U);

#endif