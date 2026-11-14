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

## 手法

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
- `src/tutorial2.ipynb`: 距離指定の回転表面符号、通常の平面表面符号、HGP、一般BB/公開BB例、公開mitten例、距離指定Color、Semi-hyperbolic Color、Kasai Code、Quantum Reed–Muller／Reed–Solomonからのサンプル生成・検証・保存。
- `src/qec_families.py`: 上記の再利用可能な生成・読込み関数。共通CSS処理はquantumgizmosのldpc/bposdを使用する。
- `src/data/mitten/`: 著者公開の150・200・300量子ビット例。commit・SHA-256・元のライセンスを同梱。

## 新しい符号族のサンプル (`tutorial2.ipynb`)

旧tutorialのldpc v1環境と分け、Python 3.11以上の新しい仮想環境を用意する。

```bash
python3 -m venv .venv-tutorial2
.venv-tutorial2/bin/python -m pip install -r src/requirements-tutorial2.txt
.venv-tutorial2/bin/python -m ipykernel install --user --name stabilizer-tutorial2 --display-name "QEC tutorial2"
```

`src/tutorial2.ipynb` を開き、上記カーネルで順に実行する。回転表面符号は奇数距離3以上を引数に取り、notebookの標準系列には目標の距離25・27・29も含む。
各サンプルをStim/Qiskitで検証し、`output/tutorial2/<run>/<sample>/` にチェック行列・論理Z・状態生成元・グラフ辺リスト・局所Clifford・準備回路・構成メタデータを保存する。実行時のネット接続は不要。

Colorは `make_color_code(d, lattice="4.8.8")` / `lattice="6.6.6"` で奇数距離3以上を生成する。MQT由来の格子生成部分をMITライセンス付きでローカルに取り込み、追加のMQT依存は不要。

HGPは任意行列の `make_hgp_code(h1, h2)` に加え、(3,4)-biregularの `make_regular_hgp_code(num_bits, seed=...)` と、論文公開のPEG/SA成分による `make_hgp_preset("PEG")` / `"SA"` を用意した。公開625量子ビット例の距離6/8は小規模古典kernelの全探索と対称HGPの距離公式で確認する。SAの選択対象は消失ノイズ。

Kasai Codeは、affine permutation matrices（APM）を使うQLDPC符号で、 `make_kasai_code(modulus, f_maps, g_maps, active_rows=3)` でAPM係数から生成する。`make_kasai_preset(1152)` / `(2304)` は [Zhao et al. Table A1](https://arxiv.org/html/2604.16209v2) の係数を固定し、それぞれ `[[1152,580,≤12]]` / `[[2304,1156,≤14]]` を再構成する。距離の上界を確定値として扱わない。これらはKasaiのAPM原構成に基づく**Zhaoらの中性原子向け派生インスタンス**であり、Kasai原論文の9216量子ビット例は未収録。構成原典とインスタンスの出典は別々に記録する。論文中の名称は **Kasai Code** に統一する。初出では[Kasaiの原構成](https://arxiv.org/abs/2601.08824)を引用し、使用する有限サイズインスタンスの出典としてZhaoらも引用する。notebookの標準実行に両方を含み、大きい例は数分かかる場合がある。

代数符号は次の生成関数を使う。追加パッケージは不要。

- `make_quantum_reed_muller_code(m, order=r)`：m>2rのpunctured RMから `[[2**m-1,1,2**(r+1)-1]]` を生成する。標準例は `[[7,1,3]]`、15-to-1蒸留で使う `[[15,1,3]]`、`[[31,1,3]]`、`[[63,1,3]]`、`[[31,1,7]]`、`[[127,1,7]]`。[Hastings–Haah](https://arxiv.org/abs/1709.03543) のw=0族。横断Tの可否も記録する。生成対象は論理ゼロ状態で、魔法状態や蒸留プロトコルのシミュレーションではない。
- `make_quantum_reed_solomon_code(field_degree=m, distance=d, length=n)`：GF(2**m)上のGRSからQuantum MDS `[[n,n-2d+2,d]]_q` を作り、トレース双対のX/Z座標でqubitへ展開する。m=2〜8、3≤n≤2**m、n−2(d−1)>0に対応。length省略時はn=2**m。展開後はN=mn、K=m(n−2d+2)、**qubit距離はd以上の下界**として保存する。MDS性と確定距離は元のqudit符号に対する主張で、qubit符号のMDS性は主張しない。[Grassl–Beth–Rötteler](https://arxiv.org/abs/quant-ph/0312164)

Semi-hyperbolic Colorは `make_semi_hyperbolic_color_code(L, base_tiling="bolza16")` で任意の整数L≥1を生成する。[Higgott et al.](https://arxiv.org/abs/2610.03682) の三角形細分化と双対化を実装し、n=16L²、k=8。L=2/3/4の文献距離は6/8/12、その他の未確認Lでは距離を推測しない。notebookにはL=1〜5を含める。カスタム基礎タイリングは色付き三角形隣接 `neighbors[N,3]` を渡せる。
基礎データは著者の旧Floquetデータに含まれるBolza隣接関係のみを参照し、静的ColorのHx=Hzを新たに生成する。測定回路・論理基底・Floquet距離は流用しない。出典commit・ハッシュ・ライセンスを `src/data/semi_hyperbolic/` に保存した。

一般BB多項式構成のセルはコメントアウトしてあり、公開BBプリセットを標準で使う。GALAは追加予定から外した。Semi-hyperbolic Surfaceと新規APM係数の探索は未実装。詳細は [add_codes.md](add_codes.md)。

生成コードの回帰確認は、上記環境で `python -m pytest tests -q` を実行する。
