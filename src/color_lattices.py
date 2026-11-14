# Copyright (c) 2023 - 2026 Chair for Design Automation, TUM
# SPDX-License-Identifier: MIT
# Lattice coordinate algorithms adapted from Peter-Jan Derks' MQT QECC
# HexagonalColorCode and SquareOctagonColorCode, commit
# bd4b4b075e7daf5095498eae9d05632434581730.
# See third_party/MQT_LICENSE.txt for the upstream license.
"""Deterministic triangular Color lattice checks, without importing MQT.

Only the lattice coordinate construction is adapted. Binary CSS processing and
logical-state construction are handled by qec_families via ldpc/bposd.
Data coordinates and faces are sorted for reproducible generator/qubit order.
"""

import numpy as np
from scipy import sparse

MQT_COMMIT = "bd4b4b075e7daf5095498eae9d05632434581730"


def color_lattice_checks(distance, lattice):
    d = distance  # Public wrapper validates distance and lattice.
    data, faces = set(), []
    if lattice == "6.6.6":
        centers = set()
        for y in range(d + d // 2):
            for i in range(d + d // 2 - y):
                coord = (y + 2 * i, y)
                # Red, blue, green rows have their face centers at 1,2,0 mod 3.
                if i % 3 == (y + 1) % 3:
                    centers.add(coord)
                else:
                    data.add(coord)
        offsets = [(-1, 1), (-2, 0), (-1, -1), (1, -1), (2, 0), (1, 1)]
        faces = [(center, offsets) for center in sorted(centers)]
    else:
        squares = set()
        octagons = {(x, 0) for x in range(4, d // 2 * 6, 6)}

        def data_row(count, y, start):
            for i in range(count):
                data.add((start + 6 * (i // 2) + 2 * (i % 2), y))

        y, count = 1, d
        while y <= d + d // 2:
            data_row(count, y, y - 1)
            count -= 1 if y == 1 else 2
            y += 1
            if y <= d + d // 2:
                data_row(count, y, y + 1)
                y += 1
            if y <= d + d // 2:
                start = y + 1 if y % 2 == 0 else y - 2
                for i in range(count):
                    is_square = (i % 2 == 0) == (y % 2 == 0)
                    (squares if is_square else octagons).add((start + 3 * i, y))
                y += 1
        square_offsets = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        octagon_offsets = [(-2, -1), (-1, -2), (1, -2), (2, -1),
                           (2, 1), (1, 2), (-1, 2), (-2, 1)]
        faces = ([(center, square_offsets) for center in sorted(squares)] +
                 [(center, octagon_offsets) for center in sorted(octagons)])
    coords = sorted(data, key=lambda xy: (xy[1], xy[0]))
    indices = {coord: i for i, coord in enumerate(coords)}
    rows, cols = [], []
    for row, ((x, y), offsets) in enumerate(faces):
        for dx, dy in offsets:
            if (x + dx, y + dy) in indices:
                rows.append(row)
                cols.append(indices[x + dx, y + dy])
    matrix = sparse.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)),
                               shape=(len(faces), len(coords)))
    return matrix, coords
