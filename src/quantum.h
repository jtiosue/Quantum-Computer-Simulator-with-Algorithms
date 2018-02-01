#ifndef QUANTUM_INCLUDE
#define QUANTUM_INCLUDE

#include <cstdio>
// #include "types.h" defined from unitary
#include "unitary.h"
#include <functional>

using namespace std;

class Register {

	private:
		state_map states;
		bool check_state(string state);

	public:
		unsigned int num_qubits;
		static vec_states all_states(unsigned int n);
		Register(unsigned int num_qubits);
		void set_nonzero_states(state_map &s);
		amp amplitude(string state);
		double probability(string state);
		string measure();
		char measure(unsigned int qubit);
		void print_states();
		friend std::ostream &operator<<(std::ostream &os, Register &reg);
		static state_map copy_map(state_map &s);
		vec_states nonzero_states();
		amp & operator[](string state);

		// Gates
		void apply_gate(Unitary u, vec_int qubits);

		void Hadamard(unsigned int qubit);
		void PhaseShift(unsigned int qubit, double theta);
		void PiOverEight(unsigned int qubit);
		void PauliX(unsigned int qubit);
		void PauliY(unsigned int qubit);
		void PauliZ(unsigned int qubit);
		void ControlledNot(unsigned int control_qubit, unsigned int target_qubit);
		void Toffoli(unsigned int control_qubit1, unsigned int control_qubit2, unsigned int target_qubit);
		void ControlledPhaseShift(unsigned int control_qubit, unsigned int target_qubit, double theta);
		void Swap(unsigned int qubit1, unsigned int qubit2);
		void Ising(unsigned int qubit1, unsigned int qubit2, double theta);

		// Sort of a gate
		void apply_function(function<string(string)> f);
};

#endif