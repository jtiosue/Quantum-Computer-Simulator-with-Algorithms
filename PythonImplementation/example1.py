import QuantumSimulator.register as register
import QuantumSimulator.gates as gates


#### Initialize register of 4 qubits
reg = register.Register(4)
print("State of register", reg, "\n")


#### Create algorithm. Syntax is
####    G(angle[if necessary], qubit1, qubit2[if necessary], ...)
alg = [
    "H(0)", # Hadamard to qubit 0
    "H(1)", # Hadamard to qubit 1
    "RZ(pi/2, 2)", # Rotate qubit 2 around z axis by pi/2 radians
    "RY(.23*pi, 3)", # Rotate qubit 3 around y axis by .23*pi radians
    "T(0, 2, 1)", # Tofolli gate to qubits 0, 2, and 1
    "CN(2, 3)", # Controlled-not with control qubit 2 and target qubit 3
    "Z(1)", # Pauli-Z gate to qubit 1
    "X(2)", # Pauli-X gate to qubit 2
    "Y(0)", # Pauli-Y gate to qubit 0
    "S(2, 3)" # Swap qubits 2 and 3
]


#### Run through algorithm
for gate in alg:
    gates.apply_gate(gate, reg) # apply the next gate to the register
    # print("State of register", reg, "\n")

print("State of register", reg, "\n")


#### Find probability to measure 1001 state
print("Probability to measure 1001:", reg.probability("1001"), "\n")


#### Measure and collapse register
print("Measured register to be in", reg.measure(), "state")
print("State of register", reg)