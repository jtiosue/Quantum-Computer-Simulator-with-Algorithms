import random

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


def overlap(reg1:Register, reg2:Register):
    """ compute <reg1 | reg2> """
    # assert reg1.num_qubits == reg2.num_qubits, \
    #     "registers must have the same number of qubits in order to compute overlap"

    return sum((reg1[state]).conjugate() * reg2[state] for state in reg1)