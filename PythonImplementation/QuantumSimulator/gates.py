import numpy as np

exp, PI, cos, sin = np.exp, np.pi, np.cos, np.sin

sigma_x = [[0, 1], [1, 0]]
sigma_y = [[0, -1j], [1j, 0]]
sigma_z = [[1, 0], [0, -1]]


class Gate:

    def __init__(self, unitary, qubits):
        """
        unitary is list of list representing unitary matrix
        qubits is tuple in order of qubits that unitary acts on
        """
        self.unitary, self.qubits = np.array(unitary), qubits
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

    def full_unitary(self, num_qubits):
        """
        Find the full unitary matrix of the gate on the full Hilbert space of
        dimension 2^num_qubits. ONLY WORKS FOR SINGLE QUBIT GATES RIGHT NOW
        """
        unitary = np.kron(np.eye(1 << self.qubits[0]), self.unitary)
        unitary = np.kron(unitary, np.eye(1 << (num_qubits-self.qubits[0]-1)))
        return unitary

    def __pow__(self, power):
        return Gate([list(x) for x in np.array(self.unitary)**power],
                    self.qubits)

    # def __mul__(self, other):
    #     """ self * other """




class H(Gate):
    c = 1.0/2.0**0.5 + 0.0j
    unitary = np.array([
        [c, c],
        [c, -c]
    ])
    def __init__(self, qubit):
        super().__init__(H.unitary, (qubit,))

    def __str__(self):
        return "H(%d)" % self.qubits[0]

class CX(Gate):
    unitary = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0]
    ])
    def __init__(self, control_qubit, target_qubit):
        # qubits should be tuple (control, target)
        super().__init__(CX.unitary, (control_qubit, target_qubit))

    def __str__(self):
        return "CX" + str(self.qubits)

class X(Gate):
    unitary = np.array(sigma_x)
    def __init__(self, qubit):
        super().__init__(X.unitary, (qubit,))

    def __str__(self):
        return "X%d" % self.qubits[0]

class Y(Gate):
    unitary = np.array(sigma_y)
    def __init__(self, qubit):
        super().__init__(Y.unitary, (qubit,))

    def __str__(self):
        return "Y%d" % self.qubits[0]

class Z(Gate):
    unitary = np.array(sigma_z)
    def __init__(self, qubit):
        super().__init__(Z.unitary, (qubit,))

    def __str__(self):
        return "Z%d" % self.qubits[0]

class T(Gate):
    unitary = [[0.0]*8 for _ in range(8)]
    for i in range(6): unitary[i][i] = 1.0
    unitary[6][7] = 1.0
    unitary[7][6] = 1.0
    unitary = np.array(unitary)
    def __init__(self, *qubits):
        """ qubits should be a tuple of length 3 """
        super().__init__(T.unitary, qubits)

    def __str__(self):
        return "T" + str(self.qubits)

class S(Gate):
    unitary = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ])
    def __init__(self, *qubits):
        """ swap two qubits. qubits should be tuple of length 2 """
        super().__init__(S.unitary, qubits)

    def __str__(self):
        return "S" + str(self.qubits)

class P(Gate):
    """ Phase shift P = |0><0| + exp(i theta) |1><1| """
    unitary = lambda angle: np.array([
        [1.0, 0.0],
        [0.0, exp(1.0j*angle)]
    ])
    def __init__(self, angle, qubit):
        super().__init__(P.unitary(angle), (qubit,))
        self.angle = angle

    def __str__(self):
        return "P" + str((self.angle,) + self.qubits)

class CP(Gate):
    unitary = lambda angle: np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, exp(1.0j*angle)]
    ])
    def __init__(self, angle, control_qubit, target_qubit):
        qubits = (control_qubit, target_qubit)
        super().__init__(CP.unitary(angle), qubits)
        self.angle = angle

    def __str__(self):
        return "CP" + str((self.angle,) + self.qubits)

class RX(Gate):
    unitary = lambda angle: np.array([
        [cos(angle/2), -1.0j*sin(angle/2)],
        [-1.0j*sin(angle/2), cos(angle/2)]
    ])
    def __init__(self, angle, qubit):
        """ rotate the qubit around the x axis by an angle """
        super().__init__(RX.unitary(angle), (qubit,))
        self.angle = angle

    def __str__(self):
        return "RX" + str((self.angle,) + self.qubits)

class RY(Gate):
    unitary = lambda angle: np.array([
        [cos(angle/2), -sin(angle/2)],
        [sin(angle/2), cos(angle/2)]
    ])
    def __init__(self, angle, qubit):
        """ rotate the qubit around the y axis by an angle """
        super().__init__(RY.unitary(angle), (qubit,))
        self.angle = angle

    def __str__(self):
        return "RY" + str((self.angle,) + self.qubits)

class RZ(Gate):
    unitary = lambda angle: np.array([
        [exp(-1.0j*angle/2), 0.0],
        [0.0, exp(1.0j*angle/2)]
    ])
    def __init__(self, angle, qubit):
        """ rotate the qubit around the z axis by an angle """
        super().__init__(RZ.unitary(angle), (qubit,))
        self.angle = angle

    def __str__(self):
        return "RZ" + str((self.angle,) + self.qubits)
    
class U3(Gate):
    """ u3(th, phi, lam) = Rz(phi)Ry(th)Rz(lam), see arxiv:1707.03429 """
    unitary = lambda theta, phi, lam: np.array([
        [exp(-1j*(phi+lam)/2)*cos(theta/2), 
         -exp(-1j*(phi-lam)/2)*sin(theta/2)],
        [exp(1j*(phi-lam)/2)*sin(theta/2), 
         exp(1j*(phi+lam)/2)*cos(theta/2)]
    ])
    def __init__(self, theta, phi, lam, qubit):
        super().__init__(U3.unitary(theta, phi, lam), (qubit,))
        self.params = theta, phi, lam
        
    def __str__(self):
        return "U3" + str(self.params + self.qubits)


def string_to_gate(string):
    return eval(string.upper())


def apply_gate(string, register):
    """ apply the gate represented by string to the register """
    string_to_gate(string)(register)


def apply_algorithm(algorithm, register):
    for gate in algorithm: apply_gate(gate, register)
