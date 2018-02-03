#include "unitary.h"
#include <iostream>


Unitary::Unitary(unsigned int dimension) {
	// Initialize square matrix to zero.
	this->dimension = dimension;

	matrix = new amp *[dimension];
	for (unsigned int i = 0; i < dimension; i++) matrix[i] = new amp[dimension];
	/*
	for (unsigned int i = 0; i < dimension; i++) {
		for (unsigned int j = 0; j < dimension; j++) {
			matrix[i][j] = 0.0;
		}
	}
	*/
}

amp *Unitary::operator[](unsigned int i) {
	return matrix[i];
}

Unitary Unitary::operator*(amp x) {
	Unitary n(dimension);
	for (unsigned int r = 0; r < dimension; r++) {
		for (unsigned int c = 0; c < dimension; c++) {
			n[r][c] = matrix[r][c] * x;
		}
	}
	return n;
}

// Left scalar multiplication, nonmember function
Unitary operator*(amp x, Unitary &U) {
	return U * x;
}

Unitary Unitary::operator*(Unitary &U) {
	Unitary f(dimension);
	if (U.dimension != dimension) {
		printf("Matrices cannot be multiplied; different dimensions\n");
		return f;
	}
	for (unsigned int row = 0; row < dimension; row++) {
		for (unsigned int column = 0; column < dimension; column++) {
			for (unsigned int t = 0; t < dimension; t++) {
				f[row][column] += matrix[row][t] * U[t][column];
			}
		}
	}
	return f;
}

std::ostream &operator<<(std::ostream &os, Unitary &U) {
	// Show all nonzero states
	for (unsigned int r = 0; r < U.dimension; r++) {
		for (unsigned int c = 0; c < U.dimension; c++) {
			os << U[r][c] << "\t";
		}
		os << "\n";
	}
	return os;
}

void Unitary::print_matrix() {
	std::cout << *this;
}

Unitary Unitary::Identity(unsigned int dimension) {
	// Unitary is initialized to all zeros
	// I don't think this is compiler specific.
	Unitary u(dimension);
	for (unsigned int i = 0; i < dimension; i++) u[i][i] = 1;
	return u;
}

Unitary Unitary::Hadamard() {
	// H = ((|0> + |1>) <0 | + (|0> - |1>) <1|) / sqrt(2)
	Unitary u(2); amp c = 1 / sqrt(2);
	u[0][0] = c; u[0][1] = c; u[1][0] = c; u[1][1] = -c;
	return u;
}

Unitary Unitary::PauliX() {
	// PX = |1><0| + |0><1|
	Unitary u(2);
	u[0][1] = 1.0; u[1][0] = 1.0;
	return u;
}

Unitary Unitary::PauliY() {
	// PY = i(| 1><0 | -| 0><1 | )
	Unitary u(2); amp c(0.0, 1.0); 
	u[0][1] = -c; u[1][0] = c;
	return u;
}

Unitary Unitary::PauliZ() {
	// PZ = |1><0| - |0><1|
	Unitary u(2);
	u[0][0] = 1.0; u[1][1] = -1.0;
	return u;
}

Unitary Unitary::PhaseShift(double theta) {
	// P = |0><0| +exp(i theta) |1><1|
	Unitary u(2);
	u[0][0] = 1; u[1][1] = exp(amp(0, theta));
	return u;
}

Unitary Unitary::PiOverEight() {
	// PhaseShift(Pi/4)
	Unitary u(2);
	u[0][0] = 1; u[1][1] = amp(1, 1) / sqrt(2);
	return u;
}

Unitary Unitary::ControlledNot() {
	Unitary u(4); 
	u[0][0] = 1.0; u[1][1] = 1.0; u[2][3] = 1.0; u[3][2] = 1.0;
	return u;
}

Unitary Unitary::Toffoli() {
	Unitary u(8);
	for (unsigned int i = 0; i < 6; i++) u[i][i] = 1;
	u[6][7] = 1; u[7][6] = 1;
	return u;
}

Unitary Unitary::ControlledPhaseShift(double theta) {
	// |00> + |01> + |10> + exp(i theta) |11>
	Unitary u(4);
	u[0][0] = 1.0; u[1][1] = 1.0; u[2][2] = 1.0; u[3][3] = exp(amp(0, theta));
	return u;
}

Unitary Unitary::Swap() {
	/*
	Swapping q1 and q2 is equivalent to 
	ControlledNot(q1, q2), ControlledNot(q2, q1), ControlledNot(q1, q2)
	*/
	Unitary u(4);
	u[0][0] = 1.0; u[3][3] = 1.0; u[1][2] = 1.0; u[2][1] = 1.0;
	return u;
}

Unitary Unitary::Ising(double theta) {
	Unitary u(4); amp s = 1.0 / sqrt(2);
	for (unsigned int i = 0; i < 4; i++) u[i][i] = s;
	u[1][2] = s * amp(0, -1); u[2][1] = s * amp(0, -1);
	// acos(-1) is pi.
	u[0][3] = s * exp(amp(0, 1)*(theta - acos(-1.0) / 2.0));
	u[3][0] = s * exp(amp(0, 1)*(-theta - acos(-1.0) / 2.0));
	return u;
}

Unitary Unitary::QFT(unsigned int num_qubits) {
	// Makes the matrix that transforms the system.
	unsigned int N = 1 << num_qubits; Unitary u(N); // 1 << num_qubits is pow(2, qubits)
	amp omega = exp(2*acos(-1.)*amp(0, 1) / double(N)); // acos(-1) is pi.
	amp c = 1 / sqrt(N);
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			u[i][j] = pow(omega, i*j) * c;
		}
	}
	return u;
}

Unitary Unitary::IQFT(unsigned int num_qubits) {
	// Complex transpose of the QFT.
	unsigned int N = 1 << num_qubits; Unitary u(N);
	amp omega = exp(2 * acos(-1.)*amp(0, -1) / double(N)); // acos(-1) is pi.
	amp c = 1 / sqrt(N);
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			u[i][j] = pow(omega, i*j) * c;
		}
	}
	return u;
}
