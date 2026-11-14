"""Topology, refinement, and independent distance checks for static color codes."""

from itertools import combinations
from pathlib import Path
import json
import sys

import networkx as nx
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from qec_families import make_semi_hyperbolic_color_code, make_graph_sample, export_graph_sample
from semi_hyperbolic import load_triangle_base


def test_bolza_seed_independent_group_model():
    neighbors, _, info = load_triangle_base("bolza16")
    published, algebraic = nx.Graph(), nx.Graph()
    # a,b,c involutory generators; reordering distinct generators changes sign.
    for word in range(16):
        sign, bits = word >> 3, word & 7
        for color in range(3):
            published.add_edge(word, int(neighbors[word,color]), color=color)
            new_sign = sign ^ ((bits >> (color+1)).bit_count() % 2)
            other = (new_sign << 3) | (bits ^ (1 << color))
            algebraic.add_edge(word, other, color=color)
    assert nx.is_isomorphic(published, algebraic,
                            edge_match=lambda a,b: a['color']==b['color'])
    assert info['genus']==2 and info['euler_characteristic']==-2
    assert len(info['base_corner_cycles'])==6
    assert all(len(c['triangles'])==8 for c in info['base_corner_cycles'])


@pytest.mark.parametrize("length", [1,2,3,4,5,6])
def test_refinement_topology_color_and_state(length):
    c=make_semi_hyperbolic_color_code(length)
    assert (c.n,c.k)==(16*length**2,8)
    assert c.hx.shape[0]==8*length**2-2
    assert (c.hx!=c.hz).nnz==0
    assert np.all(np.asarray(c.hx.sum(axis=0))==3)
    weights=np.asarray(c.hx.sum(axis=1)).ravel()
    assert (weights==8).sum()==6
    assert set(weights)<= {6,8}
    colors=np.asarray(c.metadata['check_colors'])
    for color in range(3):
        assert np.all(np.asarray(c.hx[colors==color].sum(axis=0))==1)
    assert np.array_equal(np.bincount(c.metadata['qubit_base_triangles']),np.full(16,length**2))
    if length<=4:
        assert c.distance=={1:2,2:6,3:8,4:12}[length]
        assert c.distance_status=='literature'
    else:
        assert c.distance is None and c.distance_status=='unknown'
    assert all(make_graph_sample(c).checks.values())


def test_static_distance_six_by_meet_in_middle():
    c=make_semi_hyperbolic_color_code(2)
    h=c.hx.toarray();logical=c.logical_z.toarray()
    def packed_col(matrix,i):
        return sum(int(bit)<<j for j,bit in enumerate(matrix[:,i]))
    syndrome=[packed_col(h,i) for i in range(c.n)]
    logical_syndrome=[packed_col(logical,i) for i in range(c.n)]
    small={}
    for w in range(3):
        for support in combinations(range(c.n),w):
            s=l=0
            for i in support:s^=syndrome[i];l^=logical_syndrome[i]
            small.setdefault(s,set()).add(l)
    large={};found=None
    for w in range(4):
        for support in combinations(range(c.n),w):
            s=l=mask=0
            for i in support:
                s^=syndrome[i];l^=logical_syndrome[i];mask^=1<<i
            # Any word of weight <=5 splits into pieces of weight <=2 and <=3.
            # Equal check syndrome / different logical syndrome means a logical.
            assert all(l==other for other in small.get(s,()))
            for other_l,other_mask in large.get(s,()):
                if l!=other_l:found=mask^other_mask
            large.setdefault(s,[]).append((l,mask))
    assert found is not None and found.bit_count()==6
    # Hx=Hz, so checking the pure sector suffices for the full CSS distance.


def test_custom_base_relabelled_without_borrowing_distance():
    neighbors,_,_=load_triangle_base('bolza16')
    permutation=np.random.default_rng(5).permutation(16)
    relabelled=np.empty_like(neighbors)
    for t in range(16):
        relabelled[permutation[t]]=permutation[neighbors[t,[1,0,2]]]
    c=make_semi_hyperbolic_color_code(2,base_tiling={'name':'bolza16','neighbors':relabelled})
    assert (c.n,c.k)==(64,8) and c.distance is None
    assert all(make_graph_sample(c).checks.values())


@pytest.mark.parametrize('length',[0,-1,True,2.5])
def test_bad_refinement(length):
    with pytest.raises(ValueError):make_semi_hyperbolic_color_code(length)


def test_invalid_bases_and_export(tmp_path):
    with pytest.raises(ValueError):make_semi_hyperbolic_color_code(2,base_tiling='unknown')
    neighbors,_,_=load_triangle_base('bolza16')
    broken=neighbors.copy();broken[0,0]=0
    with pytest.raises(ValueError,match='involution'):
        make_semi_hyperbolic_color_code(2,base_tiling={'neighbors':broken})
    disconnected=np.vstack([neighbors,neighbors+16])
    with pytest.raises(ValueError,match='connected'):
        make_semi_hyperbolic_color_code(2,base_tiling={'neighbors':disconnected})
    cube=[[t^1,t^2,t^4] for t in range(8)]
    with pytest.raises(ValueError,match='hyperbolic'):
        make_semi_hyperbolic_color_code(2,base_tiling={'neighbors':cube})
    c=make_semi_hyperbolic_color_code(3)
    out=export_graph_sample(make_graph_sample(c),tmp_path/'semi')
    payload=json.loads((out/'sample.json').read_text())
    assert payload['metadata']['refinement']==3
    assert payload['metadata']['base_provenance']['commit']
    assert len(payload['metadata']['qubit_base_triangles'])==144
