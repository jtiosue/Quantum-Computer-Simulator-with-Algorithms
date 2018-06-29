import QuantumSimulator.register as register
import QuantumSimulator.gates as gates

# alg = [
#     "H(1)", "CN(0, 1)", "Y(0)", "X(1)", "T(1, 0, 2)", "S(2, 1)", "Y(2)", "CN(2, 0)",
#     "H(0)", "CN(2, 1)", "Z(0)", "Y(1)", "T(1, 2, 0)", "S(0, 1)", "X(2)", "CN(1, 2)",
#     "RZ(pi/4, 2)", "RY(.21*pi, 0)", "RX(1, 1)", "CN(0, 1)", "CN(2, 0)", "CN(1, 2)"
# ]
#
# reg = register.Register(3)
#
# for g in alg:
#     gates.apply_gate(g, reg) # apply t^th gate of alg to reg
#
# print(reg)


reg = register.Register(6)

alg = ["H(0)", "CN(0, 1)", "H(2)", "CN(2, 3)"]
for g in alg:
    gates.apply_gate(g, reg)

print(reg)

alg = ["CN(2, 4)", "CN(3, 5)"]
for g in alg:
    gates.apply_gate(g, reg)

print(reg)