"""Registry of the quantum error-correcting codes listed in Table I of arXiv:2603.17743.

Codes are built via `mqt.qecc` (github.com/munich-quantum-toolkit/qecc), which is also the
codebase behind the SAT-based state-preparation circuits (Peham et al., PRX Quantum 6, 020330
(2025), arXiv:2603.17743's ref. [45]) used in the paper as the literature baseline for most of
the smaller codes in Table I. Each entry here returns N independent |0>_L stabilizer generators
(Pauli strings), ready to pass to `stab_to_graph.stabilizers_to_graph_state`.

Only codes with an identifiable, verifiable public construction are included. Six Table I
entries ([[36,8,6]], [[40,20,6]] (non-CSS), [[42,12,5]], [[48,6,8]], [[64,12,8]], [[80,8,10]]) have no
citation in the paper and no known public source for their exact check matrices, so they are
NOT included (see README).
"""

from typing import Callable

import numpy as np

from stab_to_graph import independent_rows_gf2, pauli_strings_from_binary_rows


def _css_code_stabilizers(css_code) -> list[str]:
    """|0>_L stabilizer generators for a CSS code: independent X-checks + independent
    Z-checks + logical Z operators. Check matrix rows may be redundant (more rows than the
    code's true stabilizer rank, as for LDPC-style codes like bivariate bicycle codes), so
    only a maximal independent subset of each is kept via GF(2) elimination.
    """
    hx, hz = np.asarray(css_code.Hx), np.asarray(css_code.Hz)
    num_qubits = hx.shape[1]

    x_checks = pauli_strings_from_binary_rows(hx[independent_rows_gf2(hx)], "X", num_qubits)
    z_checks = pauli_strings_from_binary_rows(hz[independent_rows_gf2(hz)], "Z", num_qubits)
    logicals = css_code.z_logicals_as_pauli_strings()

    stabilizers = x_checks + z_checks + logicals
    if len(stabilizers) != num_qubits:
        raise ValueError(
            f"Got {len(stabilizers)} independent generators for {num_qubits} qubits "
            f"(x_checks={len(x_checks)}, z_checks={len(z_checks)}, logicals={len(logicals)})."
        )
    return stabilizers


def _mqt_codes():
    import mqt.qecc.codes as mqtcodes

    return mqtcodes


def _surface_9_1_3() -> list[str]:
    return _css_code_stabilizers(_mqt_codes().CSSCode.from_code_name("Surface", 3))


def _carbon_12_2_4() -> list[str]:
    return _css_code_stabilizers(_mqt_codes().CSSCode.from_code_name("Carbon"))


def _reed_muller_15_1_3() -> list[str]:
    # mqt.qecc names the [[15,1,3]] code "Tetrahedral"; same parameters as the paper's
    # "Reed-Muller code" entry (the punctured quantum Reed-Muller / tetrahedral code family).
    return _css_code_stabilizers(_mqt_codes().CSSCode.from_code_name("Tetrahedral"))


def _hamming_15_7_3() -> list[str]:
    return _css_code_stabilizers(_mqt_codes().CSSCode.from_code_name("Hamming"))


def _color_17_1_5() -> list[str]:
    return _css_code_stabilizers(_mqt_codes().SquareOctagonColorCode(5))


def _color_19_1_5() -> list[str]:
    return _css_code_stabilizers(_mqt_codes().HexagonalColorCode(5))


def _golay_23_1_7() -> list[str]:
    return _css_code_stabilizers(_mqt_codes().CSSCode.from_code_name("Golay"))


def _gross_144_12_12() -> list[str]:
    return _css_code_stabilizers(_mqt_codes().construct_bb_code(144))


CODES: dict[str, dict] = {
    "surface_9_1_3": {
        "table_row": "[[9,1,3]]",
        "name": "Surface code",
        "n": 9,
        "k": 1,
        "d": 3,
        "factory": _surface_9_1_3,
    },
    "carbon_12_2_4": {
        "table_row": "[[12,2,4]]",
        "name": "Carbon code",
        "n": 12,
        "k": 2,
        "d": 4,
        "factory": _carbon_12_2_4,
    },
    "reed_muller_15_1_3": {
        "table_row": "[[15,1,3]]",
        "name": "Reed-Muller code",
        "n": 15,
        "k": 1,
        "d": 3,
        "factory": _reed_muller_15_1_3,
    },
    "hamming_15_7_3": {
        "table_row": "[[15,7,3]]",
        "name": "Hamming code",
        "n": 15,
        "k": 7,
        "d": 3,
        "factory": _hamming_15_7_3,
    },
    "color_17_1_5": {
        "table_row": "[[17,1,5]]",
        "name": "Color code (4.8.8, square-octagon)",
        "n": 17,
        "k": 1,
        "d": 5,
        "factory": _color_17_1_5,
    },
    "color_19_1_5": {
        "table_row": "[[19,1,5]]",
        "name": "Color code (hexagonal)",
        "n": 19,
        "k": 1,
        "d": 5,
        "factory": _color_19_1_5,
    },
    "golay_23_1_7": {
        "table_row": "[[23,1,7]]",
        "name": "Golay code",
        "n": 23,
        "k": 1,
        "d": 7,
        "factory": _golay_23_1_7,
    },
    "gross_144_12_12": {
        "table_row": "[[144,12,12]]",
        "name": "Gross code",
        "n": 144,
        "k": 12,
        "d": 12,
        "factory": _gross_144_12_12,
    },
}

# Table I rows with no citation and no identifiable public source for their exact check
# matrices (not constructible without guessing): [[36,8,6]], [[40,20,6]] (non-CSS), [[42,12,5]],
# [[48,6,8]], [[64,12,8]], [[80,8,10]].
UNAVAILABLE_TABLE_ROWS = [
    "[[36,8,6]]",
    "[[40,20,6]] (non-CSS)",
    "[[42,12,5]]",
    "[[48,6,8]]",
    "[[64,12,8]]",
    "[[80,8,10]]",
]


def get_stabilizers(key: str) -> list[str]:
    """Return the |0>_L stabilizer generators for one of the codes in `CODES`."""
    return CODES[key]["factory"]()


if __name__ == "__main__":
    import warnings

    from stab_to_graph import (
        stabilizers_to_graph_state,
        verify_equivalence,
        verify_equivalence_qiskit,
    )

    warnings.filterwarnings("ignore")
    for key, meta in CODES.items():
        stabilizers = get_stabilizers(key)
        decomposition = stabilizers_to_graph_state(stabilizers)
        ok_stim = verify_equivalence(stabilizers, decomposition.preparation_circuit)
        ok_qiskit = verify_equivalence_qiskit(stabilizers, decomposition.preparation_circuit)
        print(
            f"{meta['table_row']:14s} {meta['name']:35s} n={meta['n']:>4} "
            f"edges={decomposition.graph.number_of_edges():>4} "
            f"stim_ok={ok_stim} qiskit_ok={ok_qiskit}"
        )
