#include "test.h"

int main() {
    set_srand(); // Always need this at the start. Initializes random.
    
    _test_unitary();
    _test_QFT();
    _test_quantum_add();
    _test_collapse();
    _test_Grover();
    
    return 0;
}