#ifndef QUANTUM_INCLUDE
#define QUANTUM_INCLUDE

#include <cstdio>
#include <string>
//#include "types.h"
#include "unitary.h"

using namespace std;

class Register {

	private:
		state_map states;
		bool check_state(string state);

	public:
		int num_qubits;
		static vec_states all_states(int n);
		Register(int num_qubits);
		// ~Register();
		void set_nonzero_states(state_map &s);
		amp amplitude(string state);
		double probability(string state);
		string measure();
		char measure(int qubit);
		void print_states();
		friend std::ostream &operator<<(std::ostream &os, Register &reg);
		static state_map copy_map(state_map &s);

		// Gates
		void apply_gate(Unitary u, vec_int qubits);

		void Hadamard(int qubit);
		void PhaseShift(int qubit, double theta);
		void PiOverEight(int qubit);
		void PauliX(int qubit);
		void PauliY(int qubit);
		void PauliZ(int qubit);
		void ControlledNot(int control_qubit, int target_qubit);
		void Toffoli(int control_qubit1, int control_qubit2, int target_qubit);
		void ControlledPhaseShift(int control_qubit, int target_qubit, double theta);
		void Swap(int qubit1, int qubit2);
		void Ising(int qubit1, int qubit2, double theta);
};

#endif