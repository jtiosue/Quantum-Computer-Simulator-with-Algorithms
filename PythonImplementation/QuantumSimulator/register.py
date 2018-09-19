import random
import numpy as np
import QuantumSimulator.gates as gates


def all_states(num_qubits):
    """
    Generator, yields '00', '01', '10', '11'
    for two qubits, and similar for more
    """
    if num_qubits == 1: yield from ("0", "1")
    elif num_qubits > 1:
        for s in all_states(num_qubits-1):
            yield from (s+"0", s+"1")


class Register(dict):

    def __init__(self, num_qubits):
        """ initialize register to |"0"*num_qubits> """
        super().__init__()
        self.num_qubits = num_qubits
        self["0"*num_qubits] = 1.0+0.0j

    def __getitem__(self, item):
        """
        return the amplitude of the state. setitem still works
        by inheriting from dict class
        """
        return super().__getitem__(item) if item in self else 0.0+0.0j

    def amplitude(self, state):
        return self[state]

    def probability(self, state):
        return abs(self[state])**2

    def apply_gate(self, gate):
        """ apply Gate object to the register """
        assert gate.dimension == 1 << gate.num_qubits, ( # 1 << a = 2**a
           "Unitary dimension is not correct to be applied to the input qubits"
        )

        old, temp_states = self.copy(), [x for x in all_states(gate.num_qubits)]
        for state in old:
            s = ""
            for q in gate.qubits: s += state[q]
            r = int(s, base=2)
            self[state] -= (1.0 - gate[r][r]) * old[state]
            if self.probability(state) < 1e-16: self.pop(state) # zero beyond machine precision

            j = 0
            for k in temp_states:
                if j != r:
                    s = list(state)
                    for l in range(len(k)): s[gate.qubits[l]] = k[l]
                    s = "".join(s)
                    c = gate[j][r] * old[state]
                    if s in self:
                        self[s] += c
                        if self.probability(s) < 1e-16: self.pop(s)
                    elif c != 0.0: self[s] = c
                j += 1

    def measure(self):
        """ Measure the system and collapse it into a state """
        r = random.random()
        total = 0.0
        for state in self:
            total += self.probability(state)
            if r <= total: # collapse the state
                self.clear()
                self[state] = 1.0+0.0j
                return state

    def measure_qubit(self, qubit):
        """ measure the qubit and collapse the system """
        zero_prob = 0.0
        for state in self:
            if state[qubit] == "0": zero_prob += self.probability(state)

        if random.random() < zero_prob: v, p = "0", zero_prob**0.5
        else: v, p = "1", (1-zero_prob)**0.5

        old_register = self.copy()
        self.clear()

        for (state, amp) in old_register.items():
            if state[qubit] == v: self[state] = amp / p

        return v

    def duplicate(self):
        """ return a copy of the register """
        reg = Register(self.num_qubits)
        for (key, item) in self.items(): reg[key] = item
        return reg

    def ket(self):
        """
        Returns the ket vector of the state in the computational basis.
        Returns in the form of a numpy array.
        """
        ket = np.zeros((1 << self.num_qubits, 1))*0j
        for state in self: ket[int(state, base=2)][0] = self[state]
        return ket

    def bra(self):
        """
        Returns the bra vector of the state in the computational basis.
        Returns in the form of a numpy array.
        """
        return np.conjugate(np.transpose(self.ket()))

    def density_matrix(self):
        """
        Returns the density matrix representing the state of ther resgister.
        Returns in the form of a numpy array. Note that the Register class
        only supports pure states, so the density matrix will be that of a
        pure state
        """
        ket = self.ket()
        return ket @ np.conjugate(np.transpose(ket))
        
    def apply_algorithm(self, algorithm):
        gates.apply_algorithm(algorithm, self)



class MixedRegister:

    def __init__(self, registers):
        """
        registers is a list of tuples of the form
            registers = [(.2, reg1), (.8, reg2)]
        where the 0th index of the tuple corresponds to the statistical prob
        of being in the pure state, where the 1st index of the tuple is that
        pure state (a Register object).
        """
        self.registers = registers
        self.num_qubits = registers[0][1].num_qubits
        assert(all(r[1].num_qubits == self.num_qubits for r in registers))

    def apply_gate(self, gate):
        for _, reg in self.registers: reg.apply_gate(gate)

    def density_matrix(self):
        return sum(p*reg.density_matrix() for p, reg in self.registers)

    def coefficient(self, state1, state2):
        """ Return the coefficient in front of |state1><state2| in the density matrix """
        return sum(reg.amplitude(state1)*reg.amplitude(state2).conjugate()*p
                   for p, reg in self.registers)

    def __str__(self):
        return str(self.registers)

    def dephase(self):
        rho = np.diag(np.diag(self.density_matrix()))
        self.registers = rho_to_mixedregister(rho).registers

    def duplicate(self):
        return MixedRegister(self.registers.copy())
        
    def apply_algorithm(self, algorithm):
        gates.apply_algorithm(algorithm, self)



def register_overlap(reg1:Register, reg2:Register):
    """ compute <reg1 | reg2> """
    # assert reg1.num_qubits == reg2.num_qubits, \
    #     "registers must have the same number of qubits in order to compute overlap"

    return sum((reg1[state]).conjugate() * reg2[state] for state in reg1)

def mixedregister_overlap(reg1:MixedRegister, reg2:MixedRegister):
    """ returns Tr(reg1.density_matrix()*reg2.density_matrix()) """
    # return np.trace(reg1.density_matrix()@reg2.density_matrix())
    s = 0
    for k in range(1 << reg1.num_qubits):
        for j in range(1 << reg1.num_qubits):
            state_k, state_j = bin(k)[2:], bin(j)[2:]
            state_k = "0" * (reg1.num_qubits - len(state_k)) + state_k
            state_j = "0" * (reg1.num_qubits - len(state_j)) + state_j
            s += reg1.coefficient(state_k, state_j)*reg2.coefficient(state_j, state_k)
    return s

def ket_to_register(ket):
    """

    :param ket: two-dim row vector: the ket vector of a pure state
    :return: Register object representing the state
    """
    assert abs(sum(abs(x[0])**2 for x in ket)-1) < 1e-10, "Ket not normalized"
    num_qubits = int(np.log2(len(ket)))
    reg = Register(num_qubits)
    reg.pop("0"*num_qubits)
    for i in filter(lambda j: ket[j][0], range(1 << num_qubits)):
        n = bin(i)[2:]
        n = "0"*(num_qubits-len(n)) + n
        reg[n] = complex(ket[i][0])
    return reg

def rho_to_mixedregister(rho):
    """
    Take a density matrix `rho` and return the corresponding MixedRegister
    object
    """
    assert abs(np.trace(rho) - 1) < 1e-10, "Density matrix not normalized "
    D, U = np.linalg.eig(rho)
    n = len(D)
    registers = [
        (D[i], ket_to_register(np.array([[U[j][i]] for j in range(n)])))
        for i in range(n) if D[i] > 1e-10
    ]
    return MixedRegister(registers)