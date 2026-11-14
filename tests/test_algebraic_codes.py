"""Field, distance, and transversal-phase checks for algebraic code samples."""

from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from algebraic_checks import BinaryExtensionField, FIELD_POLYNOMIALS
from qec_families import (make_quantum_reed_solomon_code,
                          make_quantum_reed_muller_code, make_graph_sample,
                          export_graph_sample, load_css_code)
from ldpc import mod2


def polynomial_remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


@pytest.mark.parametrize("m", range(2, 9))
def test_binary_fields_and_trace_pairing(m):
    f = BinaryExtensionField(m)
    # Independent trial division by every possible small monic factor.
    for degree in range(1, m//2 + 1):
        for factor in range(1 << degree, 1 << (degree+1)):
            assert polynomial_remainder(FIELD_POLYNOMIALS[m], factor) != 0
    for a in range(1, f.order):
        assert f.multiply(a, f.power(a, f.order-2)) == 1
    rng = np.random.default_rng(m)
    for _ in range(100):
        a, b = map(int, rng.integers(f.order, size=2))
        x = np.array([(a >> j) & 1 for j in range(m)])
        z = np.array([f.trace(f.multiply(b, 1 << j)) for j in range(m)])
        assert int(x @ z % 2) == f.trace(f.multiply(a, b))


@pytest.mark.parametrize("m,n,d", [(2,3,2), (2,4,2), (3,7,3), (3,8,4),
                                   (4,15,4), (4,16,5), (5,32,9),
                                   (6,9,4), (7,7,3), (8,9,4)])
def test_rs_mds_parent_and_binary_image(m,n,d):
    c = make_quantum_reed_solomon_code(m,d,length=n)
    assert (c.n,c.k)==(m*n,m*(n-2*d+2))
    assert c.distance==d and c.distance_status=='lower_bound'
    p=c.metadata['parent_parameters']
    assert p['k']+2*p['distance']==p['n']+2
    f=BinaryExtensionField(m)
    h=c.metadata['symbol_checks']
    for a in h:
        for b in h:
            overlap=0
            for x,y in zip(a,b):overlap ^= f.multiply(x,y)
            assert overlap==0
    assert all(make_graph_sample(c).checks.values())


def test_rs_distance_in_symbol_and_binary_metrics():
    c=make_quantum_reed_solomon_code(3,3,length=7)
    h=c.metadata['symbol_checks']
    f=BinaryExtensionField(3)
    # Exhaust the q-ary centralizer up to the claimed symbol distance.
    # Stabilizer words have distance n-r+1=6, so any kernel word <=3 is logical.
    first=None
    for weight in range(1,4):
        for support in combinations(range(7),weight):
            for values in product(range(1,8),repeat=weight):
                syndrome=[]
                for row in h:
                    s=0
                    for i,v in zip(support,values):s ^= f.multiply(row[i],v)
                    syndrome.append(s)
                if not any(syndrome):first=weight;break
            if first is not None:break
        if first is not None:break
    assert first==3
    # At the smallest binary image, independently find a weight-two logical.
    c=make_quantum_reed_solomon_code(2,2,length=3)
    found=None
    for weight in [1,2]:
        for support in combinations(range(c.n),weight):
            word=np.zeros(c.n,dtype=np.uint8);word[list(support)]=1
            if not (c.hx@word%2).any() and mod2.rank(np.vstack([c.hz.toarray(),word]))>mod2.rank(c.hz):
                found=weight;break
        if found is not None:break
    assert found==2
    assert not c.metadata['binary_quantum_mds_claimed']


@pytest.mark.parametrize("m,r,d", [(3,1,3),(4,1,3),(5,1,3),(6,1,3),
                                   (5,2,7),(7,2,7)])
def test_qrm_distances_and_logical_witness(m,r,d):
    c=make_quantum_reed_muller_code(m,order=r)
    assert (c.n,c.k,c.distance)==((1<<m)-1,1,d)
    # Nonzero points of an (r+1)-dimensional subspace give a weight-d logical Z.
    witness=np.zeros(c.n,dtype=np.uint8);witness[:d]=1
    assert not (c.hx@witness%2).any()
    assert int(witness.sum())==d
    assert mod2.rank(np.vstack([c.hz.toarray(),witness]))>mod2.rank(c.hz)
    # Detect every weight-one or weight-two X/Z error, independently of d metadata.
    for checks in [c.hx,c.hz]:
        columns=checks.toarray().T
        assert columns.any(axis=1).all()
        assert len({bytes(col) for col in columns})==c.n
    assert all(make_graph_sample(c).checks.values())


@pytest.mark.parametrize("m", [3,4,5,6])
def test_qrm_exact_transversal_phase_on_both_logical_states(m):
    c=make_quantum_reed_muller_code(m)
    basis=c.hx.toarray()
    modulus=1<<c.metadata['transversal_phase_level']
    # All computational basis components of |0_L> acquire phase 1; all of
    # |1_L> acquire exp(-2*pi*i/modulus). In particular m=4 implements T†.
    for coefficients in product([0,1],repeat=m):
        word=np.asarray(coefficients,dtype=np.uint8)@basis%2
        assert int(word.sum())%modulus==0
        assert int((1-word).sum())%modulus==modulus-1
    assert c.metadata['transversal_t']==(m>=4)


@pytest.mark.parametrize("args,kwargs", [((1,2),{}),((9,2),{}),((3,5),{}),
                                         ((3,2),{'length':9}),((True,2),{})])
def test_rs_invalid_parameters(args,kwargs):
    with pytest.raises(ValueError):make_quantum_reed_solomon_code(*args,**kwargs)


def test_qrm_invalid_parameters_and_rs_export(tmp_path):
    with pytest.raises(ValueError):make_quantum_reed_muller_code(4,order=2)
    c=make_quantum_reed_solomon_code(3,3,length=7)
    target=export_graph_sample(make_graph_sample(c),tmp_path/'rs')
    restored=load_css_code(target/'code.npz')
    assert (restored.hx!=c.hx).nnz==0 and (restored.hz!=c.hz).nnz==0
    import json
    data=json.loads((target/'sample.json').read_text())
    assert data['metadata']['parent_parameters']['q']==8
    assert data['distance_status']=='lower_bound'
