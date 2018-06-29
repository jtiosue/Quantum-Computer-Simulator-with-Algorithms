from numpy import exp, pi, cos, sin


class Gate:

    def __init__(self, unitary, qubits):
        """
        unitary is list of list representing unitary matrix
        qubits is tuple in order of qubits that unitary acts on
        """
        self.unitary, self.qubits = unitary, qubits
        self.dimension, self.num_qubits = len(unitary), len(qubits)

    def __getitem__(self, item):
        """ Gate[i][j] gets the (i, j) element of the unitary matrix """
        return self.unitary[item]

    def __call__(self, register):
        """
        Apply gate to register.

        :param register: Register object to apply the gate to
        :return: the register, so that we can string gate function calls.
            ie, gate1(gate2(gate3(register)))
        """
        register.apply_gate(self)
        return register


class Hadamard(Gate):
    def __init__(self, qubit):
        c = 1.0/2.0**0.5 + 0.0j
        super().__init__([
            [c, c],
            [c, -c]
        ], (qubit,))

    def __str__(self):
        return "H(%d)" % self.qubits[0]

class ControlledNot(Gate):
    def __init__(self, control_qubit, target_qubit):
        # qubits should be tuple (control, target)
        super().__init__([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0]
        ], (control_qubit, target_qubit))

    def __str__(self):
        return "CN" + str(self.qubits)

class PauliX(Gate):
    def __init__(self, qubit):
        super().__init__([
            [0.0, 1.0],
            [1.0, 0.0],
        ], (qubit,))

    def __str__(self):
        return "X%d" % self.qubits[0]

class PauliY(Gate):
    def __init__(self, qubit):
        super().__init__([
            [0.0, -1.0j],
            [1.0j, 0.0],
        ], (qubit,))

    def __str__(self):
        return "Y%d" % self.qubits[0]

class PauliZ(Gate):
    def __init__(self, qubit):
        super().__init__([
            [1.0, 0.0],
            [0.0, -1.0],
        ], (qubit,))

    def __str__(self):
        return "Z%d" % self.qubits[0]

class Tofolli(Gate):
    def __init__(self, *qubits):
        """ qubits should be a tuple of length 3 """
        unitary = [[0.0]*8 for _ in range(8)]
        for i in range(6): unitary[i][i] = 1.0
        unitary[6][7] = 1.0
        unitary[7][6] = 1.0
        super().__init__(unitary, qubits)

    def __str__(self):
        return "T" + str(self.qubits)

class Swap(Gate):
    def __init__(self, *qubits):
        """ swap two qubits. qubits should be tuple of length 2 """
        super().__init__([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0]
        ], qubits)

    def __str__(self):
        return "S" + str(self.qubits)

class PhaseShift(Gate):
    """ Phase shift P = |0><0| + exp(i theta) |1><1| """
    def __init__(self, angle, qubit):
        super().__init__([
            [1.0, 0.0],
            [0.0, exp(1.0j*angle)]
        ], (qubit,))

        self.angle = angle

    def __str__(self):
        return "P" + str((self.angle,) + self.qubits)

class ControlledPhaseShift(Gate):
    def __init__(self, angle, control_qubit, target_qubit):
        qubits = (control_qubit, target_qubit)
        super().__init__([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, exp(1.0j*angle)]
        ], qubits)

        self.angle = angle

    def __str__(self):
        return "CP" + str((self.angle,) + self.qubits)

class XRotation(Gate):
    def __init__(self, angle, qubit):
        """ rotate the qubit around the x axis by an angle """
        super().__init__([
            [cos(angle/2), -1.0j*sin(angle/2)],
            [-1.0j*sin(angle/2), cos(angle/2)]
        ], (qubit,))

        self.angle = angle

    def __str__(self):
        return "RX" + str((self.angle,) + self.qubits)

class YRotation(Gate):
    def __init__(self, angle, qubit):
        """ rotate the qubit around the y axis by an angle """
        super().__init__([
            [cos(angle/2), -sin(angle/2)],
            [sin(angle/2), cos(angle/2)]
        ], (qubit,))

        self.angle = angle

    def __str__(self):
        return "RY" + str((self.angle,) + self.qubits)

class ZRotation(Gate):
    def __init__(self, angle, qubit):
        """ rotate the qubit around the z axis by an angle """
        super().__init__([
            [exp(-1.0j*angle/2), 0.0],
            [0.0, exp(1.0j*angle/2)]
        ], (qubit,))

        self.angle = angle

    def __str__(self):
        return "RZ" + str((self.angle,) + self.qubits)


gate_map = {
    "H": Hadamard, "CN": ControlledNot, "X": PauliX, "Y": PauliY,
    "Z": PauliZ, "T": Tofolli, "S": Swap, "P": PhaseShift,
    "CP": ControlledPhaseShift, "RX": XRotation, "RY": YRotation,
    "RZ": ZRotation
}


def string_to_gate(string):
    """ string is a key of gate_map. return correct Gate object."""
    i = string.index("(")
    arg = eval(string[i:])
    if not isinstance(arg, tuple): arg = (arg,)
    return gate_map[string[:i]](*arg)


def apply_gate(string, register):
    """ apply the gate represented by string to the register """
    string_to_gate(string)(register)
