import QuantumSimulator.register as register
import QuantumSimulator.gates as gates


#### create gates
H0 = gates.Hadamard(0) # initialize hadamard gate on qubit 0
H1 = gates.Hadamard(1) # initialize hadamard gate on qubit 1
CN24 = gates.ControlledNot(2, 4) # control qubit = 2, target qubit = 4
R = gates.ZRotation(2.3, 5) # rotate qubit 5 around z axis by 2.3 radians
T265 = gates.Tofolli(2, 6, 5) # tofolli on qubits 2, 6, 5


#### initalize register with 7 qubits
reg = register.Register(7)
print("State of register", reg, "\n")


#### apply gates. Note that 1, 2, and 3 are equivalent

## 1
# reg.apply_gate(H0)
# reg.apply_gate(H1)
# reg.apply_gate(CN24)
# reg.apply_gate(R)
# reg.apply_gate(T265)

## 2
# H0(reg)
# H1(reg)
# CN24(reg)
# R(reg)
# T265(reg)

## 3
T265(R(CN24(H1(H0(reg)))))


#### print the state of the register
print("State of register", reg, "\n")


#### measure only qubit 3 and collapse
print("Measured qubit 3 to be", reg.measure_qubit(3))
print("State of register", reg, "\n")


#### measure the entire register and collapse
print("Measured register to be", reg.measure())
print("State of register", reg)
