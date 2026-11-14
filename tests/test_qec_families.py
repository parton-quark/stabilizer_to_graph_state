"""Physics/algebra regressions for the tutorial2 generators and exported samples."""

from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from qec_families import (
    export_graph_sample, load_css_code, load_mitten_code, make_bb_code,
    make_bb_preset, make_css_code, make_graph_sample, make_hgp_code,
    make_planar_surface_code, make_rotated_surface_code,
    make_color_code, make_kasai_code, make_kasai_preset, make_regular_hgp_code,
    make_regular_ldpc, make_hgp_preset, classical_seed_metrics,
)
from stab_to_graph import (independent_rows_gf2, verify_equivalence,
                           verify_equivalence_qiskit, stim_circuit_to_qiskit)


def rowspace(matrix):
    basis = matrix.toarray()
    basis = basis[independent_rows_gf2(basis)]
    words = {bytes(np.zeros(matrix.shape[1], dtype=np.uint8))}
    for row in basis:
        words |= {bytes(np.frombuffer(word, dtype=np.uint8) ^ row) for word in list(words)}
    return words


def brute_distance(code, limit):
    """Independent low-weight logical search, without trusting distance metadata."""
    sectors = [(code.hx, rowspace(code.hz)), (code.hz, rowspace(code.hx))]
    for weight in range(1, limit + 1):
        for support in combinations(range(code.n), weight):
            word = np.zeros(code.n, dtype=np.uint8)
            word[list(support)] = 1
            for checks, words in sectors:
                if not (checks @ word % 2).any() and bytes(word) not in words:
                    return weight
    return None


@pytest.mark.parametrize("factory,n", [(make_rotated_surface_code, 9),
                                      (make_planar_surface_code, 13)])
def test_surface_distance_independent_search(factory, n):
    code = factory(3)
    assert code.n == n and code.k == 1
    assert brute_distance(code, 3) == 3
    assert make_graph_sample(code).checks == {"stim": True, "qiskit": True}


@pytest.mark.parametrize("d", [25, 27, 29])
def test_target_rotated_sizes(d):
    code = make_rotated_surface_code(d)
    assert code.n == d*d and code.k == 1
    assert code.hx.shape[0] + code.hz.shape[0] == d*d - 1
    assert code.hx.sum(axis=1).max() == 4
    assert code.hz.sum(axis=1).max() == 4
    sample = make_graph_sample(code)
    assert all(sample.checks.values())


@pytest.mark.parametrize("n,k", [(72,12),(90,8),(108,8),(144,12),(288,12)])
def test_public_bb_dimensions(n, k):
    code = make_bb_preset(n)
    assert (code.n, code.k) == (n, k)
    assert make_graph_sample(code).checks["qiskit"]


def test_bb_modular_cancellation_and_hgp_sparse():
    code = make_bb_code(l=3, m=2, polynomial_a=[(0,0),(3,0)], polynomial_b=[(1,1)])
    assert code.hx[:, :6].nnz == 0
    seed = sparse.csr_matrix([[1,1,0],[0,1,1]], dtype=np.uint8)
    hgp = make_hgp_code(seed, seed)
    assert hgp.n == 13 and hgp.k == 1 and hgp.distance is None
    assert make_graph_sample(hgp).checks["stim"]


@pytest.mark.parametrize("n,k", [(150,30),(200,40),(300,60)])
def test_upstream_mitten_data(n, k):
    code = load_mitten_code(n)
    assert (code.n, code.k) == (n,k)
    assert all(make_graph_sample(code).checks.values())


def test_export_roundtrip_and_metrics(tmp_path):
    import networkx as nx
    import stim
    sample = make_graph_sample(make_bb_preset(72))
    target = export_graph_sample(sample, tmp_path / "bb72")
    payload = json.loads((target / "sample.json").read_text())
    reconstructed = nx.Graph()
    reconstructed.add_nodes_from(payload["nodes"])
    reconstructed.add_edges_from(payload["edges"])
    assert nx.utils.graphs_equal(reconstructed, sample.decomposition.graph)
    restored = load_css_code(target / "code.npz")
    assert (restored.n, restored.k) == (72, 12)
    circuit = stim.Circuit((target / "preparation.stim").read_text())
    assert verify_equivalence(sample.stabilizers, circuit)
    graph_circuit = stim.Circuit()
    graph_circuit.append("H", range(sample.code.n))
    for u, v in reconstructed.edges():
        graph_circuit.append("CZ", [u, v])
    graph_circuit += stim.Circuit((target / "local_clifford.stim").read_text())
    assert verify_equivalence(sample.stabilizers, graph_circuit)
    improvements = 0
    for vertex in reconstructed:
        changed = reconstructed.copy()
        for u, v in combinations(list(reconstructed[vertex]), 2):
            if changed.has_edge(u, v):
                changed.remove_edge(u, v)
            else:
                changed.add_edge(u, v)
        improvements += changed.number_of_edges() < reconstructed.number_of_edges()
    assert improvements == sample.metrics["improving_lc_vertices"]
    with pytest.raises(FileExistsError):
        export_graph_sample(sample, target)


@pytest.mark.parametrize("bad", [2,4,True,3.5,0])
def test_bad_rotated_distance(bad):
    with pytest.raises(ValueError):
        make_rotated_surface_code(bad)


def test_invalid_checks_and_logicals():
    with pytest.raises(ValueError, match="commute"):
        make_css_code([[1,0]], [[1,0]])
    with pytest.raises(ValueError, match="binary"):
        make_css_code([[0.5,0]], [[0,1]])
    with pytest.raises(ValueError, match="independent"):
        make_css_code([[1,1,0]], [[1,1,0]], logical_z=[[1,1,0]])
    # Zero-check / zero-logical sectors are legitimate CSS boundary cases.
    empty = sparse.csr_matrix((0,2), dtype=np.uint8)
    code = make_css_code(empty, sparse.eye(2, dtype=np.uint8))
    assert code.k == 0 and all(make_graph_sample(code).checks.values())


@pytest.mark.parametrize("lattice,n", [("4.8.8",17),("6.6.6",19)])
def test_color_distance_five(lattice, n):
    code = make_color_code(5, lattice=lattice)
    assert (code.n, code.k) == (n,1)
    assert brute_distance(code, 5) == 5
    assert all(make_graph_sample(code).checks.values())


@pytest.mark.parametrize("lattice", ["4.8.8","6.6.6"])
@pytest.mark.parametrize("d", [3,7,13,25,27,29])
def test_color_family_parameters(lattice,d):
    code = make_color_code(d,lattice=lattice)
    expected = (d*d+2*d-1)//2 if lattice=="4.8.8" else (3*d*d+1)//4
    assert code.n == expected and code.k == 1
    assert code.hx.sum(axis=1).max() <= (8 if lattice=="4.8.8" else 6)
    assert np.array_equal(code.hx.toarray(),code.hz.toarray())
    assert all(make_graph_sample(code).checks.values())


@pytest.mark.parametrize("n,k,d", [(1152,580,12),(2304,1156,14)])
def test_kasai_paper_parameters(n,k,d):
    code=make_kasai_preset(n)
    assert (code.n,code.k,code.distance,code.distance_status)==(n,k,d,"upper_bound")
    for checks in [code.hx,code.hz]:
        assert np.all(np.asarray(checks.sum(axis=1))==12)
        assert np.all(np.asarray(checks.sum(axis=0))==3)
        overlap=(checks.astype(np.int64)@checks.T).toarray()
        np.fill_diagonal(overlap,0)
        assert overlap.max()<=1  # Tanner girth >=6, independently of metadata.
    # Large graph synthesis is also run in the notebook; here check the complete
    # CSS state with Stim without paying that synthesis cost on every test run.
    import stim
    from qec_families import logical_zero_stabilizers
    t=stim.Tableau.from_stabilizers([stim.PauliString(s) for s in logical_zero_stabilizers(code)])
    assert len(t)==n


def test_kasai_inputs_and_parameter_reconstruction():
    original=make_kasai_preset(1152)
    m=original.metadata
    regenerated=make_kasai_code(m['modulus'],m['f_maps'],m['g_maps'],active_rows=m['active_rows'])
    assert (original.hx!=regenerated.hx).nnz==0
    assert (original.hz!=regenerated.hz).nnz==0
    assert regenerated.distance is None  # Coefficients do not imply a distance certificate.
    with pytest.raises(ValueError,match="coprime"):
        make_kasai_code(6,[(2,0)],[(1,0)],active_rows=1)
    with pytest.raises(ValueError,match="active_rows"):
        make_kasai_code(6,[(1,0)],[(1,0)],active_rows=2)


@pytest.mark.parametrize("method,d", [("PEG",6),("SA",8)])
def test_published_hgp_seeds(method,d):
    code=make_hgp_preset(method)
    assert (code.n,code.k,code.distance)==(625,25,d)
    assert code.distance_status=='computed'
    assert all(make_graph_sample(code).checks.values())


def test_regular_hgp_reproducible_biregular():
    a=make_regular_ldpc(16,seed=3,candidates=4)
    b=make_regular_ldpc(16,seed=3,candidates=4)
    assert (a!=b).nnz==0
    assert np.all(np.asarray(a.sum(axis=0))==3)
    assert np.all(np.asarray(a.sum(axis=1))==4)
    metrics=classical_seed_metrics(a)
    assert metrics['rank']==12 and metrics['connected']
    c=make_regular_hgp_code(16,seed=3,candidates=4)
    assert (c.n,c.k)==(400,16)
    assert c.distance==metrics['classical_distance']
    assert all(make_graph_sample(c).checks.values())
    with pytest.raises(ValueError,match='integral'):
        make_regular_ldpc(15)
    with pytest.raises(RuntimeError):
        make_regular_ldpc(16,candidates=2,max_attempts=1)


def test_qiskit_batch_matches_scalar_expectations():
    import stim
    from qiskit.quantum_info import StabilizerState,Pauli
    circuit=stim.Circuit('H 0 2\nS 0\nH 1 3\nCZ 0 1 2 3\nH 1 3\nCZ 1 3\nX 2')
    state=StabilizerState(stim_circuit_to_qiskit(circuit))
    rng=np.random.default_rng(7)
    probes=[]
    for _ in range(60):
        label=('+' if rng.integers(2) else '-')+''.join(rng.choice(list('IXYZ'),4))
        expected=state.expectation_value(Pauli(label[0]+label[1:][::-1]))==1
        assert verify_equivalence_qiskit([label],circuit)==expected
        probes.append(label)
    expected_all=all(state.expectation_value(Pauli(s[0]+s[1:][::-1]))==1 for s in probes)
    assert verify_equivalence_qiskit(probes,circuit)==expected_all
    table=stim.TableauSimulator();table.do_circuit(circuit)
    good=table.canonical_stabilizers()
    assert verify_equivalence_qiskit(good,circuit)
    wrong=list(good);wrong[0]=-wrong[0]
    assert not verify_equivalence_qiskit(wrong,circuit)
