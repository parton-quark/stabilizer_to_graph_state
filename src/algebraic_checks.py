"""Small exact algebraic check constructions; all field elements are integers.

GF(2^m) elements use polynomial coefficients, least significant coefficient
first. The RS routines deliberately keep symbol and binary metrics distinct.
"""

from itertools import combinations

import numpy as np


# Irreducible binary polynomials, including the leading x^m coefficient.
FIELD_POLYNOMIALS = {2: 0x7, 3: 0xB, 4: 0x13, 5: 0x25,
                     6: 0x43, 7: 0x83, 8: 0x11B}


class BinaryExtensionField:
    """Polynomial-basis GF(2^m), for the explicitly supported degrees 2..8."""

    def __init__(self, degree):
        if degree not in FIELD_POLYNOMIALS:
            raise ValueError("Supported field degrees are 2..8")
        self.degree = degree
        self.order = 1 << degree
        self.polynomial = FIELD_POLYNOMIALS[degree]

    def multiply(self, a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.order:
                a ^= self.polynomial
        return result

    def power(self, a, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, a)
            a = self.multiply(a, a)
            exponent >>= 1
        return result

    def trace(self, a):
        result = 0
        for _ in range(self.degree):
            result ^= a
            a = self.multiply(a, a)
        if result not in (0, 1):
            raise ValueError("Field trace did not land in GF(2)")
        return result


def reed_solomon_binary_checks(field_degree, length, distance):
    """Dual-containing GRS CSS [[n,n-2(d-1),d]]_q, then trace expansion.

    For evaluation points a_i and Lagrange weights w_i=1/prod_{j!=i}(a_i-a_j),
    v_i=sqrt(w_i) makes rows v_i*a_i^j (0<=j<d-1) self-orthogonal:
    sum_i w_i*a_i^t=0 for t<=n-2. The classical dual is MDS as well.

    X(a) is expanded in polynomial coordinates; Z(b) has binary coordinate
    j equal to Tr(b*x^j). This preserves the Pauli trace commutator even when
    the polynomial basis is NOT self-dual. Using identical bit expansions
    for X and Z would generally be incorrect.
    """
    field = BinaryExtensionField(field_degree)
    m, n, r = field_degree, length, distance - 1
    points = list(range(n))
    multipliers = []
    for a in points:
        derivative = 1
        for b in points:
            if a != b:
                derivative = field.multiply(derivative, a ^ b)
        inverse = field.power(derivative, field.order - 2)
        multipliers.append(field.power(inverse, field.order // 2))
    symbol_checks = [[field.multiply(v, field.power(a, j))
                      for a, v in zip(points, multipliers)] for j in range(r)]
    hx = np.zeros((m*r, m*n), dtype=np.uint8)
    hz = np.zeros_like(hx)
    for row, coefficients in enumerate(symbol_checks):
        for t in range(m):
            for i, coefficient in enumerate(coefficients):
                scaled = field.multiply(1 << t, coefficient)
                for j in range(m):
                    hx[m*row+t, m*i+j] = (scaled >> j) & 1
                    hz[m*row+t, m*i+j] = field.trace(field.multiply(scaled, 1 << j))
    return hx, hz, {"field_degree": m, "field_order": field.order,
                    "irreducible_polynomial": field.polynomial,
                    "evaluation_points": points, "column_multipliers": multipliers,
                    "symbol_checks": symbol_checks,
                    "binary_basis": "polynomial basis 1,x,...,x^(m-1)",
                    "binary_z_rule": "coordinate j = Tr(b*x^j)",
                    "qubit_order": "symbol-major, polynomial coefficient low degree first"}


def shortened_rm_checks(order, num_variables):
    """Shortened RM(order,m): nonconstant monomials, punctured at zero."""
    points = np.arange(1, 1 << num_variables, dtype=np.uint64)
    rows = []
    for degree in range(1, order + 1):
        for indices in combinations(range(num_variables), degree):
            mask = sum(1 << i for i in indices)
            rows.append((points & mask) == mask)
    return np.asarray(rows, dtype=np.uint8)
