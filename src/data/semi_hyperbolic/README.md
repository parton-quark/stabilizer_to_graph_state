# Bolza base triangle map

`bolza16.json` stores the **base geometry**, not distance-specific CSS matrices.
`neighbors[t][i]` is the triangle across the edge opposite corner color i.
Subdividing each triangle into L^2 small triangles and taking the dual gives
static semi-hyperbolic color checks for any positive integer L.

The construction is from Higgott et al., arXiv:2610.03682v1, Supplement S2.
The Bolza 8.8.8-16 seed is identified explicitly there, and the paper states that
the triangular fine-graining geometry was also used in the earlier Floquet work.
We recover only this seed's three colored edge matchings from the first X/Y/Z
weight-two MPP layers of the author's public 16-data-qubit circuit. Commit,
URL, source-file hash and extraction description are in the JSON. The Apache-2.0
upstream license is retained. No upstream implementation code is copied.

The numerical adjacency record was derived from that source; it is not a copy
of the Floquet circuit. In particular **no Floquet stabilizers, logical operators,
or Floquet distances are used as static color-code data**. Static X and Z checks
are newly constructed as vertex-to-small-triangle incidence on the refined map.

The colored seed graph is independently isomorphic to the Cayley graph of the
16-element group with three pairwise anticommuting involutions a,b,c and central
minus one. Pair-generated corner links all have length 8, the Euler characteristic
is -2, and the surface is orientable of genus 2. The resulting color code has k=8.

Literature color-code distances for L=1,2,3,4 are 2,6,8,12, from the 2026 static
color-code paper, not the Floquet file's distance label. Other L values have no
distance assigned. `manifest.json` fixes the local numerical-record hash.
