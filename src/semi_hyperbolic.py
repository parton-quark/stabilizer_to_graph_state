"""Static semi-hyperbolic color checks from a colored triangle adjacency map.

Implements the triangular fine-graining construction in arXiv:2610.03682, S2.
The base combinatorial map, not a Floquet circuit, is the input. Parallel base
edges and repeated corner triples are retained as distinct cells: on the Bolza
surface a bare list of corner triples would lose the topology.
"""

from hashlib import sha256
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse


def load_triangle_base(base_tiling):
    """Load the bundled Bolza map or validate a custom neighbors[N,3] map.

    neighbors[t][i] is the triangle across the edge opposite corner color i.
    Each column is a fixed-point-free involution. Two-color cycles are the
    vertices of the dual triangle map, including distinct parallel edges.
    """
    preset = isinstance(base_tiling, str)
    if preset:
        if base_tiling != "bolza16":
            raise ValueError("The available base_tiling preset is 'bolza16'")
        root = Path(__file__).resolve().parent / "data" / "semi_hyperbolic"
        p = root / "bolza16.json"
        expected = json.loads((root / "manifest.json").read_text())["files"][p.name]["sha256"]
        if sha256(p.read_bytes()).hexdigest() != expected:
            raise ValueError("Base tiling checksum mismatch")
        base = json.loads(p.read_text())
        digest = expected
    elif isinstance(base_tiling, dict) and "neighbors" in base_tiling:
        base = dict(base_tiling)
        digest = None
    else:
        raise ValueError("Use 'bolza16' or a dictionary containing neighbors[N,3]")
    raw = np.asarray(base["neighbors"])
    if (raw.ndim != 2 or raw.shape[1] != 3 or raw.shape[0] < 4 or
            not np.issubdtype(raw.dtype, np.integer) or raw.dtype == np.bool_):
        raise ValueError("neighbors must be an integer N-by-3 array")
    n = raw.shape[0]
    if np.any(raw < 0) or np.any(raw >= n):
        raise ValueError("Triangle neighbor index out of range")
    neighbors = raw.astype(int)
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    for face in range(n):
        if len(set(neighbors[face])) != 3:
            raise ValueError("Each triangle must have three distinct neighbors")
        for color, other in enumerate(neighbors[face]):
            if other == face or neighbors[other, color] != face:
                raise ValueError("Each neighbor color must be a fixed-point-free involution")
            graph.add_edge(face, int(other))
    if not nx.is_connected(graph):
        raise ValueError("The base tiling must be connected")
    corner_ids = np.zeros((n, 3), dtype=int)
    cycles = []
    for corner in range(3):
        other_colors = [i for i in range(3) if i != corner]
        subgraph = nx.Graph()
        subgraph.add_nodes_from(range(n))
        for t in range(n):
            for color in other_colors:
                subgraph.add_edge(t, int(neighbors[t, color]))
        for component in sorted(nx.connected_components(subgraph), key=min):
            if len(component) < 4 or len(component) % 2:
                raise ValueError("Base corner links must be even cycles of length >=4")
            for t in component:
                corner_ids[t, corner] = min(component)
            cycles.append({"corner_color": corner, "triangles": sorted(component)})
    euler = len(cycles) - n // 2
    if n % 2 or euler >= 0:
        raise ValueError("The base must describe a closed hyperbolic surface (Euler characteristic <0)")
    orientable = nx.is_bipartite(graph)
    info = {"base_tiling": base.get("name", "custom"),
            "base_preset": "bolza16" if preset else None,
            "base_neighbors": neighbors.tolist(), "base_triangles": n,
            "base_corner_cycles": cycles, "euler_characteristic": euler,
            "orientable": orientable,
            "genus": (2-euler)//2 if orientable else None,
            "base_sha256": digest, "base_provenance": base.get("provenance", {})}
    return neighbors, corner_ids, info


def semi_hyperbolic_color_checks(base_tiling, refinement):
    """Subdivide each triangle into L^2 triangles; return dual face checks.

    Barycentric points on a shared edge are glued by the edge's original cell
    identity and corner colors. No geometric tolerances or coordinate rounding.
    Qubits correspond to small triangles; checks to vertices of the refined map.
    """
    neighbors, corner_ids, info = load_triangle_base(base_tiling)
    length = refinement  # Public wrapper validates integer L>=1.
    fine_triangles, module_ids, colors = [], [], {}

    def vertex(face, weights):
        nonzero = [i for i, w in enumerate(weights) if w]
        if len(nonzero) == 1:
            corner = nonzero[0]
            key = ("v", corner, int(corner_ids[face, corner]))
        elif len(nonzero) == 2:
            opposite = next(i for i, w in enumerate(weights) if not w)
            key = ("e", opposite, min(face, int(neighbors[face, opposite])), weights[nonzero[0]])
        else:
            key = ("f", face, *weights)
        color = sum(i*w for i, w in enumerate(weights)) % 3
        if key in colors and colors[key] != color:
            raise ValueError("Inconsistent corner colors across a glued edge")
        colors[key] = color
        return key

    def add_triangle(face, weights):
        fine_triangles.append(tuple(vertex(face, w) for w in weights))
        module_ids.append(face)

    for face in range(len(neighbors)):
        for i in range(length):
            for j in range(length-i):
                a = (length-i-j, i, j)
                b = (length-i-j-1, i+1, j)
                c = (length-i-j-1, i, j+1)
                add_triangle(face, (a, b, c))
                if i+j < length-1:
                    d = (length-i-j-2, i+1, j+1)
                    add_triangle(face, (b, d, c))
    keys = sorted(colors)
    indices = {key: i for i, key in enumerate(keys)}
    rows, cols = [], []
    for qubit, triangle in enumerate(fine_triangles):
        if len(set(triangle)) != 3 or len({colors[v] for v in triangle}) != 3:
            raise ValueError("Refined triangle has invalid corner identifications/coloring")
        for v in triangle:
            rows.append(indices[v]); cols.append(qubit)
    checks = sparse.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)),
                               shape=(len(keys), len(fine_triangles)))
    expected_vertices = (len(fine_triangles) - len(neighbors)) // 2 + len(info["base_corner_cycles"])
    if len(fine_triangles) != len(neighbors)*length**2 or len(keys) != expected_vertices:
        raise ValueError("Fine-graining cell counts disagree with Euler characteristic")
    info.update(refinement=length, check_colors=[colors[k] for k in keys],
                corner_check_indices=[indices[k] for k in keys if k[0] == "v"],
                qubit_base_triangles=module_ids,
                qubit_order="base triangle, then (i,j), upward before downward small triangle",
                check_order="sorted cell keys: shared edge, face interior, original corner")
    return checks, info
