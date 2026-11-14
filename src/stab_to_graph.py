"""Stabilizer state -> graph state + Clifford decomposition.

Implements the Section II construction of arXiv:2603.17743
("Fast stabilizer state preparation via AI-optimized graph decimation"):
any N-qubit stabilizer state |psi> is locally-Clifford equivalent to a graph
state |G>, i.e. |psi> = C|G> for some graph G and product of single-qubit
Clifford gates C. This module derives (G, C) from a set of stabilizer
generators, and verifies that the resulting circuit reproduces |psi>.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Sequence, Union

import networkx as nx
import numpy as np
import stim

if TYPE_CHECKING:
    from qiskit import QuantumCircuit

PauliInput = Union[str, stim.PauliString]


@dataclass
class GraphStateDecomposition:
    """|psi> = C|G>, where |G> is the graph state of `graph` and `local_clifford_circuit`
    implements C. `preparation_circuit` prepares |psi> starting from |0...0>."""

    graph: nx.Graph
    local_clifford_circuit: stim.Circuit
    preparation_circuit: stim.Circuit
    tableau: stim.Tableau


def _to_pauli_strings(stabilizers: Iterable[PauliInput]) -> list[stim.PauliString]:
    return [s if isinstance(s, stim.PauliString) else stim.PauliString(s) for s in stabilizers]


def stabilizers_to_graph_state(stabilizers: Sequence[PauliInput]) -> GraphStateDecomposition:
    """Derive a graph state and local Clifford gates equivalent to a stabilizer state.

    Uses Stim's graph-state synthesis (`stim.Tableau.to_circuit(method="graph_state")`),
    which Gaussian-eliminates the stabilizer tableau into an RX / CZ / single-qubit-Clifford
    circuit, matching the construction in arXiv:2603.17743 Section II.

    Args:
        stabilizers: N independent, commuting Pauli strings (generators) on N qubits that
            define the stabilizer state.

    Returns:
        GraphStateDecomposition with the graph G (edges = CZ gates), the local Clifford
        circuit C, and the full preparation circuit (RX -> CZ layer -> C).
    """
    pauli_strings = _to_pauli_strings(stabilizers)
    num_qubits = len(pauli_strings[0])
    if len(pauli_strings) != num_qubits:
        raise ValueError(
            f"Expected {num_qubits} independent stabilizer generators for {num_qubits} "
            f"qubits, got {len(pauli_strings)}."
        )

    tableau = stim.Tableau.from_stabilizers(pauli_strings, allow_underconstrained=False)
    circuit = tableau.to_circuit(method="graph_state")

    graph = nx.Graph()
    graph.add_nodes_from(range(num_qubits))
    local_clifford_circuit = stim.Circuit()
    for instruction in circuit:
        name = instruction.name
        if name == "CZ":
            targets = [t.value for t in instruction.targets_copy()]
            for a, b in zip(targets[0::2], targets[1::2]):
                graph.add_edge(a, b)
        elif name in ("RX", "TICK"):
            continue
        else:
            # Single-qubit Clifford correction layer (H, S, X, Y, Z, ...) = C.
            local_clifford_circuit.append(instruction)

    return GraphStateDecomposition(
        graph=graph,
        local_clifford_circuit=local_clifford_circuit,
        preparation_circuit=circuit,
        tableau=tableau,
    )


def independent_rows_gf2(rows: "np.ndarray") -> list[int]:
    """Return the indices of a maximal linearly independent subset of GF(2) matrix rows.

    Earlier rows are preferred over later, linearly-dependent ones (deterministic GF(2)
    forward elimination, no pivoting). Useful for CSS parity-check matrices that carry more
    (redundant) checks than the code's true stabilizer rank, e.g. quantum LDPC codes such as
    bivariate bicycle codes.
    """
    matrix = np.asarray(rows, dtype=np.uint8) % 2
    pivots: dict[int, "np.ndarray"] = {}
    independent_indices = []
    for i, row in enumerate(matrix):
        row = row.copy()
        for col, pivot_row in pivots.items():
            if row[col]:
                row ^= pivot_row
        nonzero = np.flatnonzero(row)
        if nonzero.size:
            pivots[int(nonzero[0])] = row
            independent_indices.append(i)
    return independent_indices


def pauli_strings_from_binary_rows(rows: "np.ndarray", letter: str, num_qubits: int) -> list[str]:
    """Turn each row of a 0/1 support matrix into a Pauli string with `letter` on the
    support and 'I' elsewhere, e.g. for building X- or Z-type CSS stabilizer generators."""
    strings = []
    for row in np.asarray(rows, dtype=np.uint8) % 2:
        chars = ["I"] * num_qubits
        for q in np.flatnonzero(row):
            chars[q] = letter
        strings.append("".join(chars))
    return strings


def graph_state_circuit(graph: nx.Graph, num_qubits: int) -> stim.Circuit:
    """Naive graph-state preparation circuit: H on all qubits, then CZ per edge (Eq. 1)."""
    circuit = stim.Circuit()
    circuit.append("H", range(num_qubits))
    for a, b in graph.edges():
        circuit.append("CZ", [a, b])
    return circuit


def verify_equivalence(stabilizers: Sequence[PauliInput], circuit: stim.Circuit) -> bool:
    """Check that `circuit`, applied to |0...0>, prepares a state stabilized by `stabilizers`.

    Direct tableau equality (`stim.Tableau.__eq__`) is not used because two tableaus can
    represent the same stabilizer state while disagreeing on the (non-unique) choice of
    destabilizers. Instead we simulate the circuit and check that every given stabilizer
    generator has expectation value +1 on the resulting state, which is exactly the
    condition for the circuit to prepare the same stabilizer group.
    """
    pauli_strings = _to_pauli_strings(stabilizers)
    sim = stim.TableauSimulator()
    sim.do_circuit(circuit)
    return all(sim.peek_observable_expectation(p) == 1 for p in pauli_strings)


def stim_circuit_to_qiskit(circuit: stim.Circuit) -> "QuantumCircuit":
    """Translate a Stim circuit into an equivalent Qiskit `QuantumCircuit` acting on |0...0>.

    Only supports the gate set emitted by `to_circuit(method="graph_state")`: RX, CZ, H, S,
    X, Y, Z. RX is Stim's reset-into-|+> operation; since qubits already start in |0>, it is
    translated as a plain H rather than a reset (Qiskit's Clifford/StabilizerState machinery
    does not support reset instructions).
    """
    from qiskit import QuantumCircuit

    num_qubits = circuit.num_qubits
    qc = QuantumCircuit(num_qubits)
    single_qubit_gates = {"H": qc.h, "S": qc.s, "X": qc.x, "Y": qc.y, "Z": qc.z}
    for instruction in circuit:
        name = instruction.name
        if name == "TICK":
            continue
        targets = [t.value for t in instruction.targets_copy()]
        if name == "RX":
            for q in targets:
                qc.h(q)
        elif name == "CZ":
            for a, b in zip(targets[0::2], targets[1::2]):
                qc.cz(a, b)
        elif name in single_qubit_gates:
            for q in targets:
                single_qubit_gates[name](q)
        else:
            raise NotImplementedError(f"stim_circuit_to_qiskit: unsupported gate {name!r}")
    return qc


def verify_equivalence_qiskit(stabilizers: Sequence[PauliInput], circuit: stim.Circuit) -> bool:
    """Cross-check of `verify_equivalence` using Qiskit instead of Stim, as an independent
    implementation to guard against bugs specific to either library.

    Stim's qubit 0 is the leftmost character of a Pauli string label, while Qiskit follows
    the little-endian convention where qubit 0 is the rightmost character, so the letters
    are reversed before constructing the `qiskit.quantum_info.Pauli` (the overall +/- sign
    is a global attribute of the operator, not tied to qubit position, and is not reversed).
    """
    from qiskit.quantum_info import PauliList, StabilizerState

    state = StabilizerState(stim_circuit_to_qiskit(circuit))
    labels = []
    for pauli_string in _to_pauli_strings(stabilizers):
        label = str(pauli_string)
        sign, letters = label[0], label[1:].replace("_", "I")
        labels.append(sign + letters[::-1])
    if not labels:
        return True
    # Evaluate every observable on U|0...0> in one Qiskit batch: U† P U
    # has expectation +1 on |0...0> iff it has no X support and phase +1.
    # This is the same check as scalar expectation_value, without its nested
    # Python loops for each of thousands of high-weight logical operators.
    evolved = PauliList(labels).evolve(state.clifford, frame="h")
    return bool(not np.any(evolved.x) and np.all(evolved.phase == 0))


if __name__ == "__main__":
    # Demo: |0>_L of the 5-qubit perfect code [[5,1,3]].
    stabilizers = ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "ZZZZZ"]

    decomposition = stabilizers_to_graph_state(stabilizers)
    print("Graph edges:", sorted(decomposition.graph.edges()))
    print("\nLocal Clifford circuit C:")
    print(decomposition.local_clifford_circuit)
    print("\nFull preparation circuit (RX -> CZ layer -> C):")
    print(decomposition.preparation_circuit)

    ok = verify_equivalence(stabilizers, decomposition.preparation_circuit)
    print("\nEquivalence check passed (Stim):", ok)

    ok_qiskit = verify_equivalence_qiskit(stabilizers, decomposition.preparation_circuit)
    print("Equivalence check passed (Qiskit):", ok_qiskit)
