#include "test.h"

int main() {
	set_srand(); // always need this at the start. initializes random.
	
	_test_unitary();
	_test_quantum_add();
	_test_Shor_factorization();
	_test_collapse();
	_test_Grover();
	_test_QFT();

	return 0;
}