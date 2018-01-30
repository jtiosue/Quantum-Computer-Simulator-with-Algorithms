#include "test.h"

using namespace std;

void _test_unitary() {
	cout << "Testing unitary...\n" << endl;
	Unitary *u;

	cout << "Check default initialization to zero:" << endl;
	u = new Unitary(2);
	cout << *u << endl;

	cout << "Check right scalar multiplication: the matrix below should be negative two times the following matrix" << endl;
	u = new Unitary(3); for (int r = 0; r < u->dimension; r++) { for (int c = 0; c < u->dimension; c++) { (*u)[r][c] = get_rand(); } }
	cout << *u << endl;
	*u = (*u)*(-2.0);
	cout << *u << endl;

	cout << "Check left scalar multiplication: the matrix below should be three times the following matrix" << endl;
	u = new Unitary(4); for (int r = 0; r < u->dimension; r++) { for (int c = 0; c < u->dimension; c++) { (*u)[r][c] = get_rand(); } }
	cout << *u << endl;
	*u = 3.0 * (*u);
	cout << *u << endl;
}

void _test_quantum_add() {
	cout << "Testing ripple carry quantum addition..." << endl;
	unsigned int a, b;

	a = 7; b = 25;
	cout << "\t" << a << "+" << b << " = " << a + b << ", Quantum add give: " << Add(a, b) << endl;

	a = 543; b = 7;
	cout << "\t" << a << "+" << b << " = " << a + b << ", Quantum add give: " << Add(a, b) << endl;

	a = 457; b = 323;
	cout << "\t" << a << "+" << b << " = " << a + b << ", Quantum add give: " << Add(a, b) << endl;

	cout << "Testing quantum fourier transform modular addition..." << endl;
	unsigned int N, n;

	a = 6; b = 8; n = 4; N = pow(2, n);
	cout << "\t" << a << "+" << b << " mod " << N << " = " << (a + b) % N << ", Quantum mod add gave: " << ModAdd(a, b, n) << endl;

	a = 50; b = 101; n = 7; N = pow(2, n);
	cout << "\t" << a << "+" << b << " mod " << N << " = " << (a + b) % N << ", Quantum mod add gave: " << ModAdd(a, b, n) << endl;

	a = 5; b = 3; n = 7; N = pow(2, n);
	cout << "\t" << a << "+" << b << " mod " << N << " = " << (a + b) % N << ", Quantum mod add gave: " << ModAdd(a, b, n) << endl;

	a = 5; b = 3; n = 3; N = pow(2, n);
	cout << "\t" << a << "+" << b << " mod " << N << " = " << (a + b) % N << ", Quantum mod add gave: " << ModAdd(a, b, n) << endl;
}

void _test_collapse() {
	cout << "Testing collapsing a register...\n" << endl;
	int n;

	n = 3; Register *r = new Register(n);
	for (int i = 0; i < r->num_qubits; i++) {
		r->Hadamard(i);
		r->PauliY(i);
		r->PiOverEight(i);
		r->PhaseShift(i, 0.2);
		r->ControlledNot(i, (i + 1) % r->num_qubits);
		r->Toffoli(i, (i + 1) % r->num_qubits, (i + 2) % r->num_qubits);
	}
	//r->print_states();
	cout << *r << endl;
	cout << "\nMeasure all qubits, got state " << r->measure() << "\n" << endl;
	r->print_states();
	cout << endl;

	n = 2; r = new Register(n);
	for (int i = 0; i < r->num_qubits; i++) {
		r->Hadamard(i);
		r->PauliX(i);
		r->PiOverEight(i);
		r->PhaseShift(i, 0.7);
		r->ControlledNot((i + 1) % r->num_qubits, i);
		r->PauliZ(i);
	}
	r->print_states();
	cout << "\nMeasure 0th qubit, got " << r->measure(0) << "\n" << endl;
	r->print_states();
	cout << endl;

	n = 4; r = new Register(n);
	for (int i = 0; i < r->num_qubits - 1; i++) {
		r->Hadamard(i);
		r->PauliX(i);
		r->PiOverEight(i);
		r->PhaseShift(i, 0.7);
		r->ControlledNot((i + 1) % r->num_qubits, i);
		r->PauliZ(i);
	}
	r->print_states();
	cout << "\nMeasure 2nd qubit, got " << r->measure(2) << "\n" << endl;
	r->print_states();
	cout << endl;
}

void _test_Grover() {
	cout << "Testing Grover search algorithm..." << endl;
	int omega; int num_bits; int N; int result;

	num_bits = 3; N = pow(2, num_bits); omega = int(get_rand()*N); result = Grover(omega, num_bits);
	cout << "\t" << "want " << omega << ", got " << result << "\n" << endl;

	num_bits = 4; N = pow(2, num_bits); omega = int(get_rand()*N); result = Grover(omega, num_bits);
	cout << "\t" << "want " << omega << ", got " << result << "\n" << endl;

	num_bits = 5; N = pow(2, num_bits); omega = int(get_rand()*N); result = Grover(omega, num_bits);
	cout << "\t" << "want " << omega << ", got " << result << "\n" << endl;

	num_bits = 6; N = pow(2, num_bits); omega = int(get_rand()*N); result = Grover(omega, num_bits, false);
	cout << "\t" << "want " << omega << ", got " << result << "\n" << endl;
}

void _test_QFT() {
	cout << "Testing quantum Fourier transform..." << endl;
	Register *r;

	r = new Register(3);
	cout << *r << endl;
	cout << "\tNow apply QFT..." << endl;
	QFT(r);
	cout << *r << endl;
	cout << "\tNow apply inverse QFT..." << endl;
	IQFT(r);
	cout << *r << endl;

	cout << "\tNew register" << endl;

	r = new Register(4); r->Toffoli(0, 2, 3); r->ControlledNot(2, 0); r->Hadamard(1);
	cout << *r << endl;
	cout << "\tNow apply QFT..." << endl;
	QFT(r);
	cout << *r << endl;
	cout << "\tNow apply inverse QFT..." << endl;
	IQFT(r);
	cout << *r << endl;

	cout << "\tNew register" << endl;

	r = new Register(4); r->Hadamard(0); r->Hadamard(1); r->Hadamard(2); r->Hadamard(3);
	r->Toffoli(0, 2, 3); r->ControlledNot(2, 0);
	r->PhaseShift(2, 2.1); r->ControlledPhaseShift(2, 1, 4.568); r->Ising(0, 2, 1.05548);
	cout << *r << endl;
	cout << "\tNow apply inverse QFT..." << endl;
	IQFT(r);
	cout << *r << endl;
	cout << "\tNow apply QFT..." << endl;
	QFT(r);
	cout << *r << endl;

	cout << "\tNew register" << endl;

	r = new Register(4);
	cout << *r << endl;
	cout << "\tNow apply QFT to just qubits 1 and 2..." << endl;
	QFT(r, 1, 3);
	cout << *r << endl;
	cout << "\tNow apply inverse QFT to just 1 and 2..." << endl;
	IQFT(r, 1, 3);
	cout << *r << endl;
}