"""Parameterized CSS code samples for tutorial2 (ldpc v2 / BP+OSD).

All states are |0>_L on every logical qubit. Code distances are construction
facts, exact small-kernel computations, or literature metadata; they are never
inferred from successful state verification.
Qubit order is the column order of hx/hz, and the left end of a Pauli string.
"""

from dataclasses import dataclass, field
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
from math import gcd
from pathlib import Path
from time import perf_counter

import networkx as nx
import numpy as np
from scipy import sparse
from ldpc import mod2
from color_lattices import color_lattice_checks, MQT_COMMIT
from algebraic_checks import reed_solomon_binary_checks, shortened_rm_checks
from semi_hyperbolic import semi_hyperbolic_color_checks

from stab_to_graph import (
    GraphStateDecomposition,
    independent_rows_gf2,
    pauli_strings_from_binary_rows,
    stabilizers_to_graph_state,
    verify_equivalence,
    verify_equivalence_qiskit,
)


def _binary_matrix(matrix, name: str) -> sparse.csr_matrix:
    """Validate before casting, so nonbinary and fractional data cannot slip in."""
    if sparse.issparse(matrix):
        values = matrix.data
        if matrix.ndim != 2 or not np.isin(values, [0, 1]).all():
            raise ValueError(f"{name} must be a two-dimensional binary matrix")
        result = matrix.astype(np.uint8).tocsr(copy=True)
        result.sum_duplicates()
        if not np.isin(result.data, [0, 1]).all():
            raise ValueError(f"{name} contains duplicate nonzero coordinates")
    else:
        values = np.asarray(matrix)
        if values.ndim != 2 or not np.isin(values, [0, 1]).all():
            raise ValueError(f"{name} must be a two-dimensional binary matrix")
        result = sparse.csr_matrix(values, dtype=np.uint8)
    result.eliminate_zeros()
    result.sort_indices()
    return result


def _nonzero_mod2(matrix) -> bool:
    return bool(np.any(matrix.data % 2))


def _integer(value, name: str, minimum: int = 1) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    if value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


@dataclass
class CodeInstance:
    name: str
    family: str
    hx: sparse.csr_matrix
    hz: sparse.csr_matrix
    logical_z: sparse.csr_matrix
    n: int
    k: int
    distance: int | None = None
    distance_status: str = "unknown"
    metadata: dict = field(default_factory=dict)

    @property
    def parameters(self) -> str:
        d = "?" if self.distance is None else str(self.distance)
        if self.distance_status == "upper_bound":
            d = "<=" + d
        elif self.distance_status == "lower_bound":
            d = ">=" + d
        return f"[[{self.n},{self.k},{d}]]"


def make_css_code(hx, hz, *, name="Custom CSS", family="custom", logical_z=None,
                  distance=None, distance_status="unknown", expected_k=None,
                  metadata=None) -> CodeInstance:
    """Validate a binary CSS code and its complete logical-Z quotient basis.

    If no logical basis is supplied, BP+OSD computes it. Supplied bases (e.g.
    upstream mitten data) are checked rather than silently replaced.
    """
    hx, hz = _binary_matrix(hx, "hx"), _binary_matrix(hz, "hz")
    if hx.shape[1] == 0 or hx.shape[1] != hz.shape[1]:
        raise ValueError("hx and hz must have the same positive column count")
    if _nonzero_mod2(hx @ hz.T):
        raise ValueError("CSS checks do not commute: hx @ hz.T != 0 mod 2")
    n = hx.shape[1]
    rank_x, rank_z = int(mod2.rank(hx)), int(mod2.rank(hz))
    k = n - rank_x - rank_z
    if expected_k is not None and k != expected_k:
        raise ValueError(f"Expected k={expected_k}, got k={k}")
    if logical_z is None:
        if k == 0:
            logical_z = sparse.csr_matrix((0, n), dtype=np.uint8)
        else:
            from bposd.css import css_code
            logical_z = css_code(hx=hx, hz=hz).lz
    logical_z = _binary_matrix(logical_z, "logical_z")
    if logical_z.shape != (k, n):
        raise ValueError(f"Expected logical_z shape {(k, n)}, got {logical_z.shape}")
    if _nonzero_mod2(hx @ logical_z.T):
        raise ValueError("Logical Z does not commute with X checks")
    if int(mod2.rank(sparse.vstack([hz, logical_z], format="csr"))) != rank_z + k:
        raise ValueError("Logical Z rows are not independent modulo Z checks")
    statuses = {"unknown", "construction", "computed", "literature", "upper_bound", "lower_bound"}
    if distance_status not in statuses:
        raise ValueError(f"distance_status must be one of {sorted(statuses)}")
    if distance is None:
        if distance_status != "unknown":
            raise ValueError("A known distance status requires a distance value")
    else:
        distance = _integer(distance, "distance")
        if distance > n or distance_status == "unknown":
            raise ValueError("Distance requires a non-unknown status and must be <= n")
    info = dict(metadata or {})
    info["rank_x"], info["rank_z"] = rank_x, rank_z
    return CodeInstance(name, family, hx, hz, logical_z, n, k,
                        distance, distance_status, info)


def make_rotated_surface_code(distance: int) -> CodeInstance:
    """Square rotated planar [[d²,1,d]] CSS code, for every odd d >= 3.

    Row-major data qubits; even (row+col) plaquettes are X. Truncated
    plaquettes are X on top/bottom and Z on left/right boundaries.
    """
    d = _integer(distance, "distance", 3)
    if d % 2 == 0:
        raise ValueError("Rotated surface distance must be odd")
    checks = {"X": [], "Z": []}
    for r in range(-1, d):
        for c in range(-1, d):
            support = [(rr, cc) for rr in (r, r + 1) for cc in (c, c + 1)
                       if 0 <= rr < d and 0 <= cc < d]
            pauli = "X" if (r + c) % 2 == 0 else "Z"
            interior = 0 <= r < d - 1 and 0 <= c < d - 1
            horizontal = r in (-1, d - 1) and 0 <= c < d - 1
            vertical = c in (-1, d - 1) and 0 <= r < d - 1
            if interior or (horizontal and pauli == "X") or (vertical and pauli == "Z"):
                checks[pauli].append([rr * d + cc for rr, cc in support])

    def matrix(rows):
        h = sparse.lil_matrix((len(rows), d * d), dtype=np.uint8)
        for i, support in enumerate(rows):
            h[i, support] = 1
        return h.tocsr()

    # Every horizontal X check overlaps this top-row Z string on 0 or 2 qubits.
    lz = sparse.lil_matrix((1, d * d), dtype=np.uint8)
    lz[0, :d] = 1
    return make_css_code(matrix(checks["X"]), matrix(checks["Z"]),
                         name=f"Rotated surface d={d}", family="rotated_surface",
                         logical_z=lz.tocsr(), expected_k=1,
                         distance=d, distance_status="construction",
                         metadata={"distance_input": d, "layout": "row-major d by d",
                                   "source": "https://arxiv.org/abs/quant-ph/0110143"})


def make_hgp_code(check_a, check_b, *, name="HGP", metadata=None) -> CodeInstance:
    """Hypergraph product of two binary classical check matrices; distance unknown."""
    from bposd.hgp import hgp
    a, b = _binary_matrix(check_a, "check_a"), _binary_matrix(check_b, "check_b")
    if min(a.shape + b.shape) < 1:
        raise ValueError("HGP seeds must be nonempty")
    code = hgp(h1=a, h2=b, compute_distance=False)
    info = dict(metadata or {})
    # Save the actual seeds, not only their sizes/seed labels.
    info.update(check_a=a.toarray().tolist(), check_b=b.toarray().tolist(),
                source="https://arxiv.org/abs/0903.0566")
    return make_css_code(code.hx, code.hz, logical_z=code.lz, name=name,
                         family="hgp", expected_k=int(code.K), metadata=info)


def classical_code_distance(checks, *, max_dimension=16) -> int | None:
    """Exact exhaustive distance of a SMALL classical kernel; None if not computed.

    Exponential in kernel dimension, not block length. Also returns None for a
    zero-dimensional kernel (no nonzero codeword). No stochastic distance label.
    """
    h = _binary_matrix(checks, "checks")
    dimension = h.shape[1] - int(mod2.rank(h))
    if dimension == 0 or dimension > max_dimension:
        return None
    basis = mod2.nullspace(h).toarray().astype(np.uint8)
    current = np.zeros(h.shape[1], dtype=np.uint8)
    best, previous_gray = h.shape[1], 0
    for i in range(1, 1 << dimension):
        gray = i ^ (i >> 1)
        changed = gray ^ previous_gray
        current ^= basis[changed.bit_length() - 1]
        best = min(best, int(current.sum()))
        previous_gray = gray
    return best


def classical_seed_metrics(checks) -> dict:
    """Finite-size diagnostics; a spectral gap is not an expansion certificate."""
    h = _binary_matrix(checks, "checks")
    overlap = (h.astype(np.int64) @ h.T).toarray()
    upper = overlap[np.triu_indices(h.shape[0], 1)]
    singular = np.linalg.svd(h.toarray().astype(float), compute_uv=False)
    tanner = nx.algorithms.bipartite.from_biadjacency_matrix(h)
    return {"rows": h.shape[0], "cols": h.shape[1], "rank": int(mod2.rank(h)),
            "row_weights": sorted(set(np.asarray(h.sum(axis=1)).ravel().tolist())),
            "column_weights": sorted(set(np.asarray(h.sum(axis=0)).ravel().tolist())),
            "four_cycles": int(np.sum(upper * (upper - 1) // 2)),
            "connected": nx.is_connected(tanner),
            "spectral_gap": float(singular[0] - singular[1]) if len(singular) > 1 else 0.0,
            "classical_distance": classical_code_distance(h)}


def make_regular_ldpc(num_bits: int, *, bit_degree=3, check_degree=4,
                      seed=0, candidates=16, max_attempts=20000) -> sparse.csr_matrix:
    """Select a simple connected biregular full-row-rank classical seed.

    Bounded configuration-model rejection sampling. Select candidates by exact
    small classical distance (when computable), then fewer 4-cycles, then larger
    spectral gap. This is NOT a reproduction of Xu et al.'s selected matrices and
    does not guarantee girth >=6 or a vertex-expansion bound.
    """
    n = _integer(num_bits, "num_bits", 2)
    dv, dc = _integer(bit_degree, "bit_degree"), _integer(check_degree, "check_degree", 2)
    candidates, max_attempts = (_integer(candidates, "candidates"),
                                _integer(max_attempts, "max_attempts"))
    seed = _integer(seed, "seed", 0)
    if dv >= dc or n * dv % dc or n < dc:
        raise ValueError("Need bit_degree < check_degree <= num_bits and integral check count")
    m = n * dv // dc
    if dv > m:
        raise ValueError("Too few check nodes for the requested bit degree")
    rng = np.random.default_rng(seed)
    bit_stubs = np.repeat(np.arange(n), dv)
    check_stubs = np.repeat(np.arange(m), dc)
    best, best_score, accepted = None, None, 0
    for _ in range(max_attempts):
        cols = rng.permutation(bit_stubs)
        if len(np.unique(check_stubs * n + cols)) != n * dv:
            continue  # Reject parallel edges; never cancel them modulo two.
        h = sparse.csr_matrix((np.ones(n * dv, dtype=np.uint8), (check_stubs, cols)),
                               shape=(m, n))
        if int(mod2.rank(h)) != m:
            continue
        metrics = classical_seed_metrics(h)
        if not metrics["connected"]:
            continue
        d = metrics["classical_distance"]
        score = (d if d is not None else 0, -metrics["four_cycles"], metrics["spectral_gap"])
        if best_score is None or score > best_score:
            best, best_score = h, score
        accepted += 1
        if accepted == candidates:
            return best
    raise RuntimeError(f"Only {accepted}/{candidates} acceptable regular seeds in {max_attempts} attempts; "
                       "increase max_attempts or change the parameters/seed")


def make_regular_hgp_code(num_bits: int, *, bit_degree=3, check_degree=4,
                          seed=0, candidates=16, max_attempts=20000) -> CodeInstance:
    """Symmetric HGP from a selected biregular seed H (H x H).

    Full row rank makes the transpose kernel trivial; the HGP distance equals
    the classical seed distance if the latter is exhaustively computed.
    """
    h = make_regular_ldpc(num_bits, bit_degree=bit_degree, check_degree=check_degree,
                          seed=seed, candidates=candidates, max_attempts=max_attempts)
    metrics = classical_seed_metrics(h)
    code = make_hgp_code(h, h, name=f"Regular HGP ({bit_degree},{check_degree}) n_seed={num_bits}",
                         metadata={"ensemble_reference": "https://arxiv.org/abs/2308.08648",
                                   "classical_seed": metrics, "selection_seed": int(seed),
                                   "candidates": int(candidates), "max_attempts": int(max_attempts),
                                   "selection": "exact small distance, fewer 4-cycles, spectral gap"})
    if metrics["classical_distance"] is not None:
        code.distance = metrics["classical_distance"]
        code.distance_status = "computed"
        code.metadata["distance_method"] = "exhaustive classical kernel + full-row-rank HGP theorem"
    return code


def make_hgp_preset(method="PEG") -> CodeInstance:
    """Published [[625,25]] HGP seed, Freire et al. arXiv:2501.09622v3.

    PEG and erasure-optimized SA are different seeds, not claims of equal
    circuit-level performance. The text files contain the classical seed H.
    """
    if method not in ("PEG", "SA"):
        raise ValueError("HGP method must be 'PEG' or 'SA'")
    root = Path(__file__).resolve().parent / "data" / "hgp"
    path = root / f"625_{method}.txt"
    manifest = json.loads((root / "manifest.json").read_text())
    if sha256(path.read_bytes()).hexdigest() != manifest["files"][path.name]["sha256"]:
        raise ValueError("Published HGP seed checksum mismatch")
    lines = path.read_text().splitlines()
    m, n = map(int, lines[0].split())
    if len(lines) != m + 1:
        raise ValueError("Unexpected HGP seed row count")
    h = sparse.lil_matrix((m, n), dtype=np.uint8)
    for row, line in enumerate(lines[1:]):
        cols = list(map(int, line.split()))
        if len(set(cols)) != len(cols) or any(col < 0 or col >= n for col in cols):
            raise ValueError("Invalid HGP seed support")
        h[row, cols] = 1
    h = h.tocsr()
    metrics = classical_seed_metrics(h)
    code = make_hgp_code(h, h, name=f"Published HGP625 {method}",
                         metadata={"source": manifest["paper"], "method": method,
                                   "input_file": manifest["files"][path.name],
                                   "classical_seed": metrics})
    # make_hgp_code sets the HGP construction source; keep the dataset separately.
    code.metadata["dataset_source"] = manifest["paper"]
    if (code.n, code.k) != (625, 25):
        raise ValueError("Published HGP parameters disagree")
    if metrics["rank"] == m and metrics["classical_distance"] is not None:
        code.distance = metrics["classical_distance"]
        code.distance_status = "computed"
        code.metadata["distance_method"] = "exhaustive classical kernel + full-row-rank HGP theorem"
    return code


def make_planar_surface_code(distance: int) -> CodeInstance:
    """Unrotated planar [[2d²-2d+1,1,d]], including even distances."""
    from ldpc.codes import rep_code
    d = _integer(distance, "distance", 2)
    seed = rep_code(d)
    code = make_hgp_code(seed, seed, name=f"Planar surface d={d}")
    code.family, code.distance, code.distance_status = "planar_surface", d, "construction"
    code.metadata["distance_input"] = d
    return code


def make_bb_code(*, l: int, m: int, polynomial_a, polynomial_b) -> CodeInstance:
    """BB CSS checks from binary monomial exponent pairs (x_power, y_power).

    Repeated monomials cancel over GF(2), including powers equal modulo l,m.
    No distance or k is promised for arbitrary polynomial choices.
    """
    l, m = _integer(l, "l"), _integer(m, "m")
    size = l * m

    def polynomial(terms):
        terms = tuple(tuple(t) for t in terms)
        mat = sparse.csr_matrix((size, size), dtype=np.uint8)
        rows = np.arange(size)
        for term in terms:
            if len(term) != 2 or any(isinstance(x, (bool, np.bool_)) or
                    not isinstance(x, (int, np.integer)) for x in term):
                raise ValueError("Monomials must be integer exponent pairs")
            a, b = map(int, term)
            cols = ((rows // m + a) % l) * m + (rows % m + b) % m
            mat = mat + sparse.csr_matrix((np.ones(size, dtype=np.uint8),
                                           (rows, cols)), shape=(size, size))
            mat.data %= 2
            mat.eliminate_zeros()
        return mat, [[int(a), int(b)] for a, b in terms]

    a, terms_a = polynomial(polynomial_a)
    b, terms_b = polynomial(polynomial_b)
    return make_css_code(sparse.hstack([a, b]), sparse.hstack([b.T, a.T]),
                         name=f"BB l={l},m={m}", family="bb",
                         metadata={"l": l, "m": m, "polynomial_a": terms_a,
                                   "polynomial_b": terms_b,
                                   "source": "https://arxiv.org/abs/2308.07915"})


# Bravyi et al. public examples; exponents agree with MQT's published table.
# k and distance are checked/annotated separately from the general constructor.
BB_PRESETS = {
    72: (6, 6, ((3, 0), (0, 1), (0, 2)), ((1, 0), (2, 0), (0, 3)), 12, 6),
    90: (15, 3, ((9, 0), (0, 1), (0, 2)), ((2, 0), (7, 0), (0, 0)), 8, 10),
    108: (9, 6, ((3, 0), (0, 1), (0, 2)), ((1, 0), (2, 0), (0, 3)), 8, 10),
    144: (12, 6, ((3, 0), (0, 1), (0, 2)), ((1, 0), (2, 0), (0, 3)), 12, 12),
    288: (12, 12, ((3, 0), (0, 2), (0, 7)), ((1, 0), (2, 0), (0, 3)), 12, 18),
}


def make_bb_preset(n: int) -> CodeInstance:
    if n not in BB_PRESETS:
        raise ValueError(f"BB preset n must be one of {list(BB_PRESETS)}")
    l, m, a, b, k, d = BB_PRESETS[n]
    code = make_bb_code(l=l, m=m, polynomial_a=a, polynomial_b=b)
    if (code.n, code.k) != (n, k):
        raise ValueError(f"BB preset parameters disagree: {code.parameters}")
    code.name = f"BB [[{n},{k},{d}]]"
    code.distance, code.distance_status = d, "literature"
    code.metadata["preset_n"] = n
    return code


def load_css_code(path, *, name=None, family="imported", metadata=None,
                  distance=None, distance_status="unknown") -> CodeInstance:
    """Read an npz (hx,hz, optional logical_z) or a yarn-style matrix directory."""
    path = Path(path)
    if path.is_dir():
        hx, hz = (np.load(path / f, allow_pickle=False) for f in ("Hx.npy", "Hz.npy"))
        lz_path = path / "Lz.npy"
        lz = np.load(lz_path, allow_pickle=False) if lz_path.exists() else None
        files = [path / "Hx.npy", path / "Hz.npy"] + ([lz_path] if lz_path.exists() else [])
    else:
        with np.load(path, allow_pickle=False) as data:
            hx, hz = data["hx"], data["hz"]
            lz = data["logical_z"] if "logical_z" in data else None
        files = [path]
    info = dict(metadata or {})
    info["input_sha256"] = {f.name: sha256(f.read_bytes()).hexdigest() for f in files}
    return make_css_code(hx, hz, logical_z=lz, name=name or path.stem,
                         family=family, metadata=info, distance=distance,
                         distance_status=distance_status)


def load_mitten_code(n: int = 150) -> CodeInstance:
    """Load vendored upstream examples; verify their pinned manifest hashes."""
    params = {150: (30, 10), 200: (40, 12), 300: (60, 14)}
    if n not in params:
        raise ValueError(f"Bundled mitten sizes: {list(params)}")
    root = Path(__file__).resolve().parent / "data" / "mitten"
    manifest = json.loads((root / "manifest.json").read_text())
    for filename in ("Hx.npy", "Hz.npy", "Lx.npy", "Lz.npy"):
        key = f"{n}/{filename}"
        if sha256((root / key).read_bytes()).hexdigest() != manifest["files"][key]["sha256"]:
            raise ValueError(f"Upstream checksum mismatch: {key}")
    k, d = params[n]
    code = load_css_code(root / str(n), name=f"Mitten [[{n},{k},{d}]]", family="mitten",
                         distance=d, distance_status="literature",
                         metadata={"source": manifest["paper"],
                                   "repository": manifest["repository"],
                                   "commit": manifest["commit"], "retrieved": manifest["retrieved"]})
    lx = _binary_matrix(np.load(root / str(n) / "Lx.npy", allow_pickle=False), "logical_x")
    pairing = (lx @ code.logical_z.T).toarray() % 2
    if code.n != n or code.k != k or lx.shape != (k, n) or _nonzero_mod2(code.hz @ lx.T):
        raise ValueError("Mitten X logicals or code parameters are inconsistent")
    if not np.array_equal(pairing, np.eye(k, dtype=np.uint8)):
        raise ValueError("Upstream logical bases are not canonically paired")
    return code


def make_color_code(distance: int, *, lattice="4.8.8") -> CodeInstance:
    """Triangular 4.8.8 or 6.6.6 [[n,1,d]] for every odd d >=3.

    Deterministic local lattice generator, with all-qubit logical Z. Requires
    no MQT installation; adapted lattice algorithm is attributed in color_lattices.
    """
    d = _integer(distance, "distance", 3)
    if d % 2 == 0 or lattice not in ("4.8.8", "6.6.6"):
        raise ValueError("Use an odd distance >=3 and lattice '4.8.8' or '6.6.6'")
    checks, coordinates = color_lattice_checks(d, lattice)
    expected_n = (d*d + 2*d - 1) // 2 if lattice == "4.8.8" else (3*d*d + 1) // 4
    if checks.shape[1] != expected_n:
        raise ValueError("Color lattice size disagrees with its construction")
    return make_css_code(checks, checks, name=f"Color {lattice} d={d}",
                         family="color", distance=d, distance_status="construction",
                         expected_k=1, logical_z=np.ones((1, expected_n), dtype=np.uint8),
                         metadata={"lattice": lattice, "distance_input": d,
                             "data_coordinates": coordinates,
                             "source": "https://github.com/munich-quantum-toolkit/qecc",
                             "lattice_algorithm_commit": MQT_COMMIT,
                             "qubit_order": "coordinates sorted by (y,x)"})


def make_semi_hyperbolic_color_code(refinement: int, *, base_tiling="bolza16") -> CodeInstance:
    """Static color code on the dual of a triangle-refined closed surface.

    L is a genuine construction parameter, not a lookup of distance instances.
    The Bolza 8.8.8-16 seed gives n=16L^2, k=8. Only L=1..4 receive the
    published distances 2/6/8/12; other refinements/custom bases stay unknown.
    A custom base is a dict with neighbors[N,3], one involution per edge color.
    """
    length = _integer(refinement, "refinement")
    checks, info = semi_hyperbolic_color_checks(base_tiling, length)
    info.update(source="https://arxiv.org/abs/2610.03682v1",
                construction="static CSS on dual of triangular Euclidean fine-graining",
                module_partition="one base triangle per module (K=1)")
    distance = {1: 2, 2: 6, 3: 8, 4: 12}.get(length) if info["base_preset"] == "bolza16" else None
    if distance is not None:
        info["distance_source"] = "arXiv:2610.03682v1, S2 and candidate tables S7-S9 (L2-L4)"
    return make_css_code(checks, checks,
                         name=f"Semi-hyperbolic color {info['base_tiling']} L={length}",
                         family="semi_hyperbolic_color", expected_k=4-2*info["euler_characteristic"],
                         distance=distance, distance_status="literature" if distance is not None else "unknown",
                         metadata=info)


def make_kasai_code(modulus: int, f_maps, g_maps, *, active_rows=3,
                    name="Kasai Code", metadata=None) -> CodeInstance:
    """Kasai Code from active-row APM CSS construction; no distance search.

    Maps are (a,b) representing x -> a*x+b mod P. Permutation convention:
    matrix[image,input]=1. Block c in X row r is F[(c-r) mod h];
    block c in Z row r is G[(r-c) mod h].T, with F/G swapped on the right.
    Arbitrary input maps are rejected when their active checks do not commute.
    """
    p, j = _integer(modulus, "modulus", 2), _integer(active_rows, "active_rows")
    f_maps, g_maps = tuple(map(tuple, f_maps)), tuple(map(tuple, g_maps))
    if not f_maps or len(f_maps) != len(g_maps) or j > len(f_maps):
        raise ValueError("Need equally sized nonempty f/g map lists and active_rows <= their length")

    def apms(maps):
        matrices, normalized = [], []
        for pair in maps:
            if len(pair) != 2 or any(isinstance(x, (bool, np.bool_)) or
                    not isinstance(x, (int, np.integer)) for x in pair):
                raise ValueError("APMs require integer (multiplier,offset) pairs")
            a, b = (int(x) % p for x in pair)
            if gcd(a, p) != 1:
                raise ValueError("APM multiplier must be coprime to modulus")
            domain = np.arange(p)
            matrices.append(sparse.csr_matrix((np.ones(p, dtype=np.uint8),
                            ((a * domain + b) % p, domain)), shape=(p,p)))
            normalized.append([a, b])
        return matrices, normalized

    f, fp = apms(f_maps)
    g, gp = apms(g_maps)
    half = len(f)
    hx = sparse.bmat([[f[(c-r) % half] for c in range(half)] +
                       [g[(c-r) % half] for c in range(half)] for r in range(j)], format="csr")
    hz = sparse.bmat([[g[(r-c) % half].T for c in range(half)] +
                       [f[(r-c) % half].T for c in range(half)] for r in range(j)], format="csr")
    info = dict(metadata or {})
    info.update(modulus=p, active_rows=j, block_columns=2*half, f_maps=fp, g_maps=gp,
                permutation_convention="matrix[image,input]=1")
    info.setdefault("construction_source", "https://arxiv.org/abs/2601.08824")
    return make_css_code(hx, hz, name=name, family="kasai", metadata=info)


def make_kasai_preset(n: int = 1152) -> CodeInstance:
    """Reconstruct the two memory examples in Zhao et al. Table A1 (v2).

    Their reported distances are UPPER bounds, never independently certified here.
    Coefficients are fixed in a provenance-bearing JSON, not a matrix lookup.
    """
    path = Path(__file__).resolve().parent / "data" / "kasai" / "zhao_table_a1.json"
    document = json.loads(path.read_text())
    if str(n) not in document["instances"]:
        raise ValueError(f"Kasai preset n must be one of {list(document['instances'])}")
    params = document["instances"][str(n)]
    code = make_kasai_code(params["modulus"], params["f_maps"], params["g_maps"],
                           active_rows=params["active_rows"], name=f"Kasai Code n={n}",
                           metadata={"source": document["source"], "source_table": "A1",
                                     "instance_variant": "Zhao et al. neutral-atom co-designed derivative",
                                     "parameters_sha256": sha256(path.read_bytes()).hexdigest()})
    if (code.n, code.k) != (n, params["k"]):
        raise ValueError(f"Kasai parameters disagree: {code.parameters}")
    code.distance, code.distance_status = params["distance_upper_bound"], "upper_bound"
    return code


def make_quantum_reed_solomon_code(field_degree: int, distance: int, *,
                                   length: int | None = None) -> CodeInstance:
    """Binary image of a quantum-MDS GRS code over GF(2^m), m=2..8.

    The parent has [[n,n-2d+2,d]]_q with n<=q and positive logical dimension.
    Binary expansion has N=m*n, K=m*(n-2d+2), and d_binary >= d. It is NOT
    labelled quantum MDS in the qubit metric. Symbol checks and basis are saved.
    """
    m = _integer(field_degree, "field_degree", 2)
    if m > 8:
        raise ValueError("Supported field degrees are 2..8")
    d = _integer(distance, "distance", 2)
    q = 1 << m
    n = q if length is None else _integer(length, "length", 3)
    k = n - 2*(d - 1)
    if n > q or k < 1:
        raise ValueError("Need 3 <= length <= 2**field_degree and length-2*(distance-1) >= 1")
    hx, hz, metadata = reed_solomon_binary_checks(m, n, d)
    metadata.update(source="https://arxiv.org/abs/quant-ph/0312164",
                    construction="Euclidean self-orthogonal GRS CSS + binary trace expansion",
                    parent_parameters={"n": n, "k": k, "distance": d, "q": q},
                    parent_quantum_mds=True, binary_quantum_mds_claimed=False,
                    distance_metric="qubit Pauli weight; lower bound from parent symbol weight",
                    binary_distance_upper_bound=m*d)
    return make_css_code(hx, hz, name=f"Quantum RS GF({q}) n={n} d={d} binary image",
                         family="quantum_reed_solomon", expected_k=m*k,
                         distance=d, distance_status="lower_bound", metadata=metadata)


def make_quantum_reed_muller_code(num_variables: int = 4, *, order: int = 1) -> CodeInstance:
    """Punctured/shortened quantum RM [[2^m-1,1,2^(r+1)-1]], m>2r.

    X checks: shortened RM(r,m); Z checks: shortened RM(m-r-1,m).
    This is the w=0 subfamily of Hastings--Haah, arXiv:1709.03543.
    r=1,m=4 gives the 15-to-1 distillation code. This generator constructs
    code states, not noisy magic states or a distillation simulation.
    """
    m, r = _integer(num_variables, "num_variables", 3), _integer(order, "order")
    if m <= 2*r:
        raise ValueError("Need num_variables > 2*order")
    hx = shortened_rm_checks(r, m)
    hz = shortened_rm_checks(m-r-1, m)
    n, d = (1 << m) - 1, (1 << (r+1)) - 1
    level = (m-1) // r
    return make_css_code(hx, hz, name=f"Quantum RM(r={r},m={m})", family="quantum_reed_muller",
                         distance=d, distance_status="construction", expected_k=1,
                         logical_z=np.ones((1, n), dtype=np.uint8),
                         metadata={"source": "https://arxiv.org/abs/1709.03543",
                                   "num_variables": m, "order": r, "puncture_weight": 0,
                                   "x_distance": (1 << (m-r)) - 1, "z_distance": d,
                                   "transversal_phase_level": level,
                                   "transversal_t": level >= 3,
                                   "uniform_phase_action": "R_level on physical qubits gives logical R_level inverse",
                                   "qubit_order": "nonzero binary points 1..2**m-1, variable 0 is LSB",
                                   "state_scope": "logical zero; distillation protocol not simulated"})


def logical_zero_stabilizers(code: CodeInstance) -> list[str]:
    hx, hz = code.hx.toarray(), code.hz.toarray()
    rows = pauli_strings_from_binary_rows(hx[independent_rows_gf2(hx)], "X", code.n)
    rows += pauli_strings_from_binary_rows(hz[independent_rows_gf2(hz)], "Z", code.n)
    rows += pauli_strings_from_binary_rows(code.logical_z.toarray(), "Z", code.n)
    if len(rows) != code.n:
        raise ValueError("The state must have n independent generators")
    return rows


@dataclass
class GraphSample:
    code: CodeInstance
    stabilizers: list[str]
    decomposition: GraphStateDecomposition
    checks: dict
    metrics: dict


def make_graph_sample(code: CodeInstance, *, verify_qiskit=True) -> GraphSample:
    start = perf_counter()
    stabilizers = logical_zero_stabilizers(code)
    decomposition = stabilizers_to_graph_state(stabilizers)
    graph = decomposition.graph
    checks = {"stim": verify_equivalence(stabilizers, decomposition.preparation_circuit),
              "qiskit": None}
    if verify_qiskit:
        checks["qiskit"] = verify_equivalence_qiskit(stabilizers, decomposition.preparation_circuit)
    if not checks["stim"] or checks["qiskit"] is False:
        raise ValueError(f"State verification failed: {checks}")
    # Count triangles in compiled sparse arithmetic, rather than Python neighbor
    # intersections for dense 1000+ qubit synthesized graphs. Use int64 to avoid
    # uint8 overflow when a pair of vertices has >255 common neighbors.
    adjacency = nx.to_scipy_sparse_array(graph, format="csr", dtype=np.int64)
    triangles = np.asarray((adjacency @ adjacency).multiply(adjacency).sum(axis=1)).ravel() // 2
    improving = sum(2 * int(triangles[v]) > degree * (degree - 1) // 2
                    for v, degree in graph.degree())
    metrics = {"n": code.n, "k": code.k, "edges": graph.number_of_edges(),
               "max_degree": max(dict(graph.degree()).values(), default=0),
               "bipartite": nx.is_bipartite(graph),
               "triangles": int(triangles.sum()) // 3,
               "improving_lc_vertices": improving,
               "max_check_weight": int(max(np.max(np.asarray(code.hx.sum(axis=1)), initial=0),
                                           np.max(np.asarray(code.hz.sum(axis=1)), initial=0))),
               "max_check_participation": int(np.max(np.asarray(code.hx.sum(axis=0) +
                                                                code.hz.sum(axis=0)), initial=0)),
               "generation_and_verification_seconds": perf_counter() - start}
    return GraphSample(code, stabilizers, decomposition, checks, metrics)


def export_graph_sample(sample: GraphSample, output_dir) -> Path:
    """Write a complete sample to a NEW directory (never overwrite another run)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    code, decomp = sample.code, sample.decomposition
    matrix_path = output_dir / "code.npz"
    np.savez_compressed(matrix_path, hx=code.hx.toarray(), hz=code.hz.toarray(),
                        logical_z=code.logical_z.toarray())
    (output_dir / "preparation.stim").write_text(str(decomp.preparation_circuit) + "\n")
    (output_dir / "local_clifford.stim").write_text(str(decomp.local_clifford_circuit) + "\n")
    (output_dir / "stabilizers.txt").write_text("\n".join(sample.stabilizers) + "\n")
    versions = {}
    for package in ("numpy", "scipy", "ldpc", "bposd", "stim", "qiskit", "networkx"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    payload = {"schema_version": 1, "name": code.name, "family": code.family,
               "n": code.n, "k": code.k, "distance": code.distance,
               "distance_status": code.distance_status,
               "state": "logical_zero_all", "pauli_order": "leftmost is qubit 0",
               "metadata": code.metadata, "metrics": sample.metrics,
               "verification": sample.checks, "versions": versions,
               "nodes": list(range(code.n)),
               "edges": sorted([min(u, v), max(u, v)] for u, v in decomp.graph.edges()),
               "files_sha256": {p.name: sha256(p.read_bytes()).hexdigest()
                                for p in output_dir.iterdir() if p.is_file()}}
    (output_dir / "sample.json").write_text(json.dumps(payload, indent=2) + "\n")
    return output_dir
