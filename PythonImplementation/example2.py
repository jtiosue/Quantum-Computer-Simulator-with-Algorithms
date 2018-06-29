import QuantumSimulator.register as register
import QuantumSimulator.gates as gates
import QuantumSimulator.algorithm as algorithm


num_qubits = 3

alg0 = algorithm.canonical_ordering(
    ["H(1)", "CN(0, 1)", "Y(0)", "X(1)", "T(1, 0, 2)", "S(2, 1)", "Y(2)", "P(pi/2, 0)"]
)
alg1 = algorithm.canonical_ordering(
    ["H(0)", "CN(2, 1)", "Z(0)", "Y(1)", "T(1, 2, 0)", "S(0, 1)", "X(2)", "CP(2.3*pi, 1, 2)"]
)

num_gates = len(alg0)
assert num_gates == len(alg1), "algorithms should be the same depth"


#### initialize two registers
reg0 = register.Register(num_qubits)
reg1 = register.Register(num_qubits)

D = 0.0 # D will become the average overlap of the two registers

for t in range(num_gates):
    gates.apply_gate(alg0[t], reg0) # apply t^th gate of alg0 to reg0
    gates.apply_gate(alg1[t], reg1) # apply t^th gate of alg1 to reg1

    overlap = register.overlap(reg0, reg1) # = <reg0 | reg1>
    print("<Reg0 | Reg1> after applying %d gates:" % t, overlap, "\n")
    D += abs(register.overlap(reg0, reg1))**2

print("State of register 0", reg0, "\n")
print("State of register 1", reg1, "\n")
print("Average overlap of the two registers", D / num_gates)