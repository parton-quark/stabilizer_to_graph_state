# Mitten sample matrices

Upstream: https://github.com/a7b/yarn (MIT license; see `LICENSE`).
Paper: https://arxiv.org/abs/2607.28795.

The 150-, 200-, and 300-data-qubit examples contain the upstream files
`Hx.npy`, `Hz.npy`, `Lx.npy`, and `Lz.npy`, unchanged. `manifest.json` records
the exact commit, source URLs, retrieval date, and SHA-256 of every matrix.
`load_mitten_code` checks these hashes and the CSS and paired-logical algebra.
Distances are literature metadata, not independently certified by this project.

The samples are bundled so tutorial2 runs without downloading at notebook runtime.
Larger examples and construction/search tools remain in the upstream repository.
