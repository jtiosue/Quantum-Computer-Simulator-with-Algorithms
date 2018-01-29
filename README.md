# Quantum-Computer-Simulator-with-Algorithms
C++ simulator of quantum registers and quantum algorithms

**Note: need std=c++11 tag to compile. In src, run g++ std=c++11 -o main \*.cpp**

*Note: I have not implemented the quantum period finding algorithm yet, currently Shor's factorization algorithm uses a classical period finding algorithm.*

*Note: My modular addition method that uses the quantum and inverse quantum fourier transform is not working.*

To see examples of how to use this code, look at src/test.cpp. To see explainations and proofs of the various algorithms, see explaination/math.pdf.

## Example

```C++
#include "quantum.h"
#include "rand.h"
#include "qalgorithms.h"
#include <iostream>

int main() {
  set_srand(); // Always need this at the start.
  
  Register reg(4); // initialize register with 4 qubits
  for(int i=0; i<4; i++) reg.Hadamard(i); // apply Hadamard gate to each qubit
  reg.ControlledNot(1, 3); // apply Controlled-Not gate to the target qubit 3 with control qubit 1
  
  QFT(&reg); // Apply quantum fourier transform to all the qubits.
  IQFT(&reg, 1, 4); // Apply inverse quantum fourier transform to qubits 1 through 3 (4 is non inclusive)
  
  reg.print_states();
  
  char c = reg.measure(0); // measure just qubit 0;
  
  std::cout << "Measured qubit 0 to be " << c << std::endl;
  
  // reg will now be collapsed into a superposition of states with the zeroth qubit being equal to c.
  std::cout << reg << std::endl; // equivalent to reg.print_states();
  
  string state = reg.measure(); // Collapse the whole system.
  
  std::cout << "Measured the register, collapsed into state " << state << std::endl;
  
  // reg will now be in state `state` with probability one.
  std::cout << reg << std::endl;
  
  return 0;
}
```

## QAlgorithms

I have implemented various algorithms that apply a series of quantum logic gates to a register in order to achieve some goal. There is quantum ripple carry addition, quantum modular arithitic using QFT and IQFT (coming soon), Grover's search algorithm, and Shor's factorization algorithm (not currently quantum; uses a classical period finding algorithm - I will implement the quantum period finding algorithm soon).

## Applying your own gates

See quantum.h/quantum.cpp and unitary.h/unitary.cpp. The Register class has a method to apply_gate(Unitary u, vector<int> v); this applies the gate given by u to the qubits given in v.
