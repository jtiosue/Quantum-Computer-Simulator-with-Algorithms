# Quantum-Computer-Simulator-with-Algorithms
*C++11 simulator of quantum registers and quantum algorithms*

*If you use this anywhere, please cite me and email (joe.iosue@yahoo.com) me so I can see it too!*

**To see explanations and proofs of the various algorithms, see explanation/math.pdf**. To see examples of how to use this code, look at src/test.cpp. To get started, compile with `g++ -std=c++11 -o test *.cpp`.

## QAlgorithms

I have implemented various algorithms that apply a series of quantum logic gates to a register in order to achieve some goal. There is quantum ripple carry addition, quantum modular arithitic using the quantum and inverse quantum Fourier transforms, Grover's search algorithm, and Shor's factorization algorithm with quantum period finding. See qalgorithms.cpp.

I am able to simulate many qubits fairly well. Finding the factors of 2813 takes on average under ten seconds (`Shor(2813)`), and 3901 (`Shor(3901)`) took under a minute on my laptop. From this, we see that simulating quantum systems on a classical computer is indeed exponentially complex.

## Example usage

```C++
#include "QuantumComputerSimulator.h"
#include <iostream>

int main() {
  
  Register reg(4); // initialize register with 4 qubits
  for(int i=0; i<4; i++) reg.Hadamard(i); // apply Hadamard gate to each qubit
  reg.ControlledNot(1, 3); // apply Controlled-Not gate to the target qubit 3 with control qubit 1
  
  QFT(&reg); // Apply quantum fourier transform to all the qubits.
  IQFT(&reg, 1, 4); // Apply inverse quantum fourier transform to qubits 1 through 3 (4 is non inclusive)
  
  reg.print_states();

  // Apply a function to the states. This particular function sends |x> to |x+1 mod 16>
  // or, in terms of qubits, |0000> goes to |0001>, |0001> goes to |0010>, ...
  // |1111> goes to |0000>.
  auto f = [](string state) {
    string r = base10_to_binary((binary_to_base10(state) + 1) % (int)pow(2, 4));
    while (r.size() < 4) r = "0" + r; // ensure resulting number is represented with 4 qubits
    return r;
  };
  // reg.apply_function will check to make sure that the function f represents a unitary transformation.
  reg.apply_function(f);

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

## Applying your own gates

See quantum.h/quantum.cpp and unitary.h/unitary.cpp. The Register class has a method to apply_gate(Unitary u, vector<int> v); this applies the gate given by u to the qubits given in v.
  
## A few things that cause weird behavior

* If you get probabilities that sum to more than one, you probably applied a controlled gate with repeated qubits; i.e. apply a controlled-not gate with control and target both being the same qubit.
* Be very careful with operations involving unsigned ints and ints. I decided to make most types in the code unsigned, but this often causes strange things to happen, and caused me many minutes of confusion. If in doubt, be very explicit; convert everything to an int when doing loops and calculations, and then convert back to unsigned int when sending values into functions.
* When computing modulo powers, use my `mod_power(a, x, N)` function from methods.h instead of `(int)pow(a, x) % N`, because the latter often causes an overflow if a or x are large enough.

## Note

A lot of the algorithms could probably be implemented more efficiently in terms of the code, maybe by combining gates or qubit operations. But for me the point was to understand what one would physically do in order to program a quantum computer. Thus, whenever I can, I only use one- or two-qubit logic gates.
