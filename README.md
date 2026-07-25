# Stabilizer to graph state with Clifford gates
グラフ状態とClifford等価なGraph状態を導出するプログラム。
グラフ状態最適化手法をstabilizer符号の状態準備などに活用することを目指す。
- 参考: https://arxiv.org/abs/2603.17743 

- 回路として表現して等価性検証も行う。
  - stimだけだと不安だったのでqiskitも

- src/requirements.txtに必要な環境が書いてあります。

## 実装方針

参考文献 (arXiv:2603.17743, "Fast stabilizer state preparation via AI-optimized graph decimation") の Section II で示されている以下の事実をベースにする。

- 任意の stabilizer 状態 |ψ⟩ は、あるグラフ G の graph 状態 |G⟩ と局所 Clifford 演算 C を用いて |ψ⟩ = C|G⟩ と表せる (Van den Nest et al. 2004)。
  - |G⟩ = ∏_{{i,j}∈E} CZ_{i,j} |+⟩^{⊗N}
  - stabilizer generator: g_i = X_i ∏_{j∈n(i)} Z_j
- 論文では、この分解 (stabilizer tableau → graph 状態 + 局所 Clifford デコレータ) を Stim のガウス消去法の実装を使って求めている。

今回のタスクはこの Section II 相当 (stabilizer → graph 状態 + Clifford の導出、および導出結果の等価性検証) に対応する。

## スコープ外
- graph state optimizationへの接続。
    - 局所補グラフ
    - CX / CY を使った edge 削減 (論文 Eq. 7, 8) や局所補完 (LC) 最適化。
- QuSynth (MCTS + RL) による探索。

### 手法

- stabilizer tableau の表現・操作、および回路構築には [Stim](https://github.com/quantumlib/Stim) を使う (論文と同じ基盤)。
  - `stim.Tableau.from_stabilizers(...)` で stabilizer generator の集合から tableau を構築。
  - `stim.Tableau.to_circuit(method="graph_state")` で `RX (=|+> 初期化) → CZ 層 (グラフの辺) → 単一量子ビット Clifford 層` という形の回路を合成する。これはまさに Fig.1 で説明されている素朴な graph 状態準備回路 (H → CZ → 局所 Clifford) に対応する。
  - この機能は Stim 1.15 以降で利用可能 (手元環境は 1.11 だったため 1.15.0 に更新した。1.16 系はソースビルドが必要で Xcode ライセンス同意が必要なため見送り)。
- 回路から CZ 命令を集めて隣接行列 (グラフ) を、それ以外の単一量子ビットゲート列を局所 Clifford C として抽出する。
- 等価性検証: `stim.Tableau` の `==` 比較は destabilizer の取り方に依存し必ずしも一致しないため使わない。代わりに `stim.TableauSimulator` で回路を |0…0⟩ から実行し、元の stabilizer generator 全てについて `peek_observable_expectation(g) == 1` となることを確認する (= 状態が同じ stabilizer 群を持つことを直接検証する)。
- 上記はStimの実装のみに依存した検証のため、独立したライブラリによるクロスチェックとして Qiskit (`qiskit.quantum_info.StabilizerState`) でも同じ検証を行う。
  - Stim の回路 (RX, CZ, H, S, X, Y, Z のみで構成) を Qiskit の `QuantumCircuit` に変換する (`stim_circuit_to_qiskit`)。RX は |0⟩ からの |+⟩ リセットだが、量子ビットは既に |0⟩ から始まっているため単なる H に読み替える (Qiskit の Clifford/StabilizerState は reset命令をサポートしないため)。
  - Pauli文字列の量子ビット順序規約が Stim (左端が qubit 0) と Qiskit (little-endian、右端が qubit 0) で逆になる点に注意し、`Pauli` に渡す前に文字列を反転する (符号は反転しない)。非対称なランダム stabilizer 状態で実際に検証し、この反転が必要であることを確認済み。
  - `state.expectation_value(pauli) == 1` を全生成元について確認する (`verify_equivalence_qiskit`)。

## Table I の符号への対応 (`src/qec_codes.py`)

論文 Table I に登場する 14 個の符号 [[n,k,d]] のうち、公開されている情報から stabilizer generator を確実に再構成できるものを [MQT QECC](https://github.com/munich-quantum-toolkit/qecc) (`mqt.qecc`) を使って実装している。MQT QECC は論文が比較対象としている SAT ベースの手法 (Peham et al., PRX Quantum 6, 020330 (2025), 論文の ref. [45]) と同じグループが公開しているツールであり、Steane/Golay/Carbon/Hamming/Tetrahedral の厳密な parity-check 行列や、rotated surface code・color code・bivariate bicycle code の構成的生成器を提供している。

### 出典一覧 (Table I の14符号すべて)

| Table I | 名前 | 状態 | 出典・実装方法 |
|---|---|---|---|
| [[9,1,3]] | Surface code | ✅ 対応済み | `mqt.qecc` `CSSCode.from_code_name("Surface", 3)` (rotated surface code) |
| [[12,2,4]] | Carbon code | ✅ 対応済み | `mqt.qecc` `CSSCode.from_code_name("Carbon")` ([[4,2,2]] と C₆ の concatenation, [Error Correction Zoo](https://errorcorrectionzoo.org/c/carbon)) |
| [[15,1,3]] | Reed-Muller code | ✅ 対応済み | `mqt.qecc` `CSSCode.from_code_name("Tetrahedral")` (同じ [[15,1,3]] パラメータ。mqt.qecc 内での呼称は Tetrahedral) |
| [[15,7,3]] | Hamming code | ✅ 対応済み | `mqt.qecc` `CSSCode.from_code_name("Hamming")` |
| [[17,1,5]] | Color code | ✅ 対応済み | `mqt.qecc` `SquareOctagonColorCode(5)` (4.8.8 格子) |
| [[19,1,5]] | Color code | ✅ 対応済み | `mqt.qecc` `HexagonalColorCode(5)` |
| [[23,1,7]] | Golay code | ✅ 対応済み | 論文 ref. [43] (Steane, *Simple quantum error-correcting codes*, PRA 54, 4741 (1996))。実装は `mqt.qecc` `CSSCode.from_code_name("Golay")` |
| [[144,12,12]] | Gross code | ✅ 対応済み | 論文 ref. [44] (Bravyi et al., *High-threshold and low-overhead fault-tolerant quantum memory*, Nature 627, 778 (2024))。bivariate bicycle code, 多項式 A=x³+y+y², B=y³+x+x², (l,m)=(12,6)。実装は `mqt.qecc` `construct_bb_code(144)` |
| [[36,8,6]] | (名称なし) | ❌ 未対応 | 論文中に引用なし。[codetables.de](https://www.codetables.de/QECC/QECCTable.php) の best-known [[36,8]] は d=8〜11 で d=6 と不一致 (=最適符号ではない別の符号)。[quantumgizmos/ldpc](https://github.com/quantumgizmos/ldpc) にも該当インスタンスなし |
| [[40,20,6]] (non-CSS) | (名称なし) | ❌ 未対応 | 論文中に引用なし。codetables.de の best-known [[40,20]] は d=6〜7 で distance は一致するが、そちらは (通常) MAGMA 生成の代数的符号で non-CSS 指定もなく、論文の符号と同一とは確認できない |
| [[42,12,5]] | (名称なし) | ❌ 未対応 | 論文中に引用なし。codetables.de の best-known [[42,12]] は d=8〜11 で d=5 と不一致 |
| [[48,6,8]] | (名称なし) | ❌ 未対応 | 論文中に引用なし。codetables.de の best-known [[48,6]] は d=11〜15 で d=8 と不一致 |
| [[64,12,8]] | (名称なし) | ❌ 未対応 | 論文中に引用なし。codetables.de の best-known [[64,12]] は d=13〜19 で d=8 と不一致 |
| [[80,8,10]] | (名称なし) | ❌ 未対応 | 論文中に引用なし。codetables.de の best-known [[80,8]] は d=15〜25 で d=10 と不一致 |

未対応の6符号について: 論文 Table I ではこれらに引用がなく (脚注3 = 著者ら独自の LC-optimization ベースラインというだけで、符号自体の出典ではない)。上記の通り codetables.de の best-known 符号は同じ [[n,k]] でも distance が (ほぼ全て) 一致せず、論文の符号は distance を最適化した符号ではなく QLDPC 系 (hypergraph product / bivariate bicycle など) の何らかの具体的なパラメータを用いていそう。`quantumgizmos/ldpc` にも一致するインスタンスは入っていなかった。exact match が確認できない以上、、現時点では対応を保留。

**検証方法**: 各符号について CSS の check 行列 (Hx, Hz) が GF(2) 上で直交する (`Hx Hz^T = 0 mod 2`) ことと、rank から求めた k が符号のパラメータと一致することを確認した上で、`_css_code_stabilizers` で独立な stabilizer generator 集合 (X-check の独立行 + Z-check の独立行 + 論理 Z 演算子、合計 n 個) を構成している。LDPC 系の符号 (Gross code など) は check 行列に冗長な行が含まれるため、`independent_rows_gf2` で GF(2) 上の最大独立行を選んでいる。得られた stabilizer は全て `verify_equivalence` / `verify_equivalence_qiskit` の両方で検証済み (144量子ビットの Gross code でも Stim 側 0.2秒未満、Qiskit 側 0.2秒程度で完了)。

なお、論文 Table I の "Graph edges" 列とは一致しない (例: Golay code は 77 で一致するが、他は異なる)。これは stabilizer generator の選び方や順序、符号の具体的な生成子基底の違いによって Gaussian 消去が導くグラフが変わるためで (同じ状態でも LC 軌道上の異なるグラフになりうる)、状態そのものの正しさとは無関係。

## 構成

- `src/stab_to_graph.py`: stabilizer → graph 状態 + Clifford の導出、および等価性検証を行う (`independent_rows_gf2`, `pauli_strings_from_binary_rows` など CSS 符号向けの補助関数も含む)。
- `src/qec_codes.py`: 論文 Table I の符号レジストリ (`mqt.qecc` ベース)。
- `src/requirements.txt`: 実行に必要な Python パッケージ。
- `src/tutorial.ipynb`: 使い方のデモ (QEC符号の例など)。