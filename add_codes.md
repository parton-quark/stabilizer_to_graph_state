# QEC符号データセットの実装状況と構成方針

更新日: 2026-11-14。公開実装の調査日: 2026-11-13。
現在利用できる符号族、生成API、構成・距離の根拠、検証方法をまとめる。実行例は [src/tutorial2.ipynb](src/tutorial2.ipynb)、生成APIは [src/qec_families.py](src/qec_families.py) を参照。

## 実装の進捗 (`tutorial2.ipynb`)

- 距離指定の回転表面符号、通常の平面表面符号、行列指定のHGP、周期サイズ・多項式指定のBBを実装した。
- BBの公開5例とmittenの150・200・300量子ビット例を利用できる。Mittenの行列・ライセンス・commit・ハッシュを同梱した。
- 共通CSS検証、論理ゼロ状態の生成、Stim/Qiskit検証、初期グラフ指標の集計、サンプルの保存を実装した。
- Colorは4.8.8 / 6.6.6を奇数距離からローカル生成する。MQTのMITライセンス格子アルゴリズムを取り込み、追加のMQT依存は不要。標準notebookでd=3,5,7,9,11,13を実行する。
- Kasai Codeは一般係数からの生成とZhao Table A1の1152・2304量子ビット例を実装。HGPには(3,4)-biregular系列と公開PEG/SA成分を追加。
- Quantum Reed–Mullerはm,r指定のpunctured RM族、Quantum Reed–SolomonはGF(2^m)上のQuantum MDS親符号とqubit展開を実装。RSのqudit距離とqubit距離の下界を区別して保存する。
- Semi-hyperbolic ColorはBolza基礎タイリングの任意L細分化とカスタム基礎隣接入力を実装。L=2/3/4の64・144・256量子ビット例を含む。
- Semi-hyperbolic Surface、新規APM係数の探索とmittenの一般構成は未実装。任意CSS行列を読み込む共通入口は利用できる。
- 旧コードと依存を混在させないよう、ldpc v2 / bposd v2用の `src/requirements-tutorial2.txt` を分離した。旧 `tutorial.ipynb` と `qec_codes.py` は維持する。

## 目的と対象

グラフ状態最適化に使う、QEC符号由来のグラフ状態データセットを生成する。
最適化実験は共著者が担当し、このリポジトリでは符号の構成、状態・グラフの生成、等価性検証と再現可能なデータ出力を扱う。

対象は、チェック演算子に全論理Zを加えて一意に定めた論理ゼロ状態
`|0>_L^{⊗k}` の準備。任意の論理入力に対するエンコーダや、耐故障な準備回路の評価は別の課題とする。
各符号の応用上の意義、サイズに伴う変化、局所的な符号と高レートQLDPCの違いを、量子情報分野外の読者にも説明できる構成を目指す。

## 利用できる符号族

| 符号族 | 生成・読込みAPI | 標準サンプル・対応範囲 |
|---|---|---|
| 回転表面符号 | `make_rotated_surface_code(d)` | 任意の奇数d≥3。標準系列にd=25,27,29を含む |
| 通常の平面表面符号 | `make_planar_surface_code(d)` | 整数d≥2。反復符号同士のHGP |
| BB | `make_bb_preset(n)` / `make_bb_code(...)` | 公開n=72,90,108,144,288。一般多項式例のセルはコメントアウト |
| Color | `make_color_code(d, lattice=...)` | 4.8.8 / 6.6.6、任意の奇数d≥3 |
| Kasai Code | `make_kasai_preset(n)` / `make_kasai_code(...)` | Zhaoらの1152・2304量子ビット例、一般APM係数入力 |
| HGP | `make_hgp_code(h1,h2)` / `make_regular_hgp_code(...)` / `make_hgp_preset(...)` | 任意の古典行列、biregular系列、公開PEG/SA例 |
| Mitten | `load_mitten_code(n)` | 公開n=150,200,300の行列・論理基底 |
| Semi-hyperbolic Color | `make_semi_hyperbolic_color_code(L, base_tiling=...)` | Bolza基礎の任意L≥1、カスタム基礎隣接入力 |
| Quantum Reed–Muller | `make_quantum_reed_muller_code(m, order=r)` | m>2r≥2のpunctured RM族 |
| Quantum Reed–Solomon | `make_quantum_reed_solomon_code(m,d,length=n)` | GF(2^m)、m=2〜8のQuantum MDS親符号とqubit展開 |

共通CSS処理は `ldpc==2.4.1` / `bposd==2.1`。依存は `src/requirements-tutorial2.txt` にまとめ、旧tutorialのldpc v1環境と分離する。Color格子はローカル生成し、MQTの実行時依存は不要。

旧 `src/tutorial.ipynb` / `src/qec_codes.py` は、MQTを使った既存Table Iの8例を保持する。Carbon・Hamming・Golayなどの固定例はそちらから利用できる。同一のチェック行列を特定できていない残り6例は未対応。

直近の確認では回帰テスト90件が成功し、標準notebookの64サンプルすべてでStim/Qiskitによる状態等価性検証、保存、ファイルハッシュ確認を完了した。

## 距離・サイズの目標と引用根拠

**回転表面符号の主目標は d=25,27,29。** 小規模からの系列も用意する。

推奨する初期系列: `d = 3, 5, 7, 9, 11, 13, 15, 19, 21, 25, 27, 29`。
これは生成対象のプリセットであり、関数自体はこのリスト以外の有効な距離も受け付ける。

| 距離 | データ量子ビット数 / グラフ頂点数 `n=d²` | 標準的な測定補助込み `2d²−1` |
|---|---:|---:|
| 13 | 169 | 337 |
| 25 | 625 | 1249 |
| 27 | 729 | 1457 |
| 29 | 841 | 1681 |

主な引用根拠:

1. Majaniemi & Matekole, *Reducing quantum error correction overhead using soft information*, Quantum Sci. Technol. 11, 025024 (2026), [arXiv:2504.03504v3, §II.2・Table 2](https://arxiv.org/html/2504.03504v3)。超伝導向けSI1000系回路レベルノイズ、2量子ビットのノイズ確率 `p=0.2%`、LCDを仮定すると、MegaQuop目標に対して測定分類エラー `p_S=p` ではd=25、`p_S=5p` ではhard decodingでd=29、soft decodingでd=27と推定される。同論文の基準はdラウンドの論理メモリエラー確率 `P_L≤10^-6` で、格子手術を想定した指標。距離9を超える推定は外挿であり、一般のアプリ全体の成功率の保証ではない。
2. Gidney, *How to factor 2048 bit RSA integers with less than a million noisy qubits* (2025), [arXiv:2505.15917v1, §3.2](https://arxiv.org/html/2505.15917v1#S3.SS2)。一様脱分極ノイズ0.1%、1 µsの訂正サイクル、10 µsの制御反応時間、格子手術等を仮定し、通常のhot-storageパッチにd=25を採用する。これは大規模アルゴリズム資源推定との接点として引用する。

他の符号族すべてに距離25〜29を要求するわけではない。Colorには距離系列、semi-hyperbolicには基礎タイリングと細分化の系列、BB・HGP・mitten・Kasai Codeには出典と構成が明確な複数サイズを用意する。
同じ距離でも符号族・ノイズ・操作が異なれば論理性能は異なるため、同距離を同等の耐故障性能とは扱わない。

## 途中の最適化結果を踏まえた選定理由

共有された最適化原稿のTable 3では、Grossの辺数比がLSで0.55、SA-ILPで0.73、GBFSで0.45、最大次数比がそれぞれ0.79、1.00、0.57となっている。BBを複数例へ拡張して、この改善が単一例に依存しないか確認する価値が高い。
一方、通常の平面表面符号では距離6でGBFSが優位でも、距離7ではSA-ILPが優位。したがって、大距離への有効性は仮説として扱い、生成可能なサイズと探索が有効なサイズを区別する。原稿の表面符号系列は通常版であり、既存コードの回転版とは別に扱う。

符号名だけでなく、変換後の初期グラフの構造に基づいて選ぶ。
次数qの頂点vの近傍にt_v本の辺があると、一手のLCによる辺数変化は `Δ|E|=q(q−1)/2−2t_v`。三角形のないグラフでは一手で辺数を減らせない。
対象のCSS論理ゼロ状態には、独立なXチェックを系統形 `[I_r|P]` にしたとき、隣接行列 `[[0,P],[P.T,0]]` の二部グラフ表現がある。符号と二部グラフ・pivotの対応は [Danielsen & Parker](https://arxiv.org/abs/0710.2243) を参照。
ただし、現在のStim変換が返したグラフが実際に二部グラフかは個別に確認する。一手のLSが改善しないことだけでLC軌道上の最適性を結論しない。

選定時に比較したい構造:

- BB、mitten、Kasai Code: 異なる代数的構成で、消去による密化とLCでの削減がどう変わるか。
- 通常・回転表面、Color、semi-hyperbolic: 局所格子、大域的な接続、細分化の違いによる変化。
- HGP: 成分の古典LDPC行列・seedを変えた、条件を制御しやすい比較群。
- Golayなど: 改善しない例・小規模で最適性を確認できる例も残す。

高レート、大きなTanner girth、符号の自己同型が多いことだけから、状態グラフのLC最適化が効くとは判断しない。符号の対称性が選択した状態・初期グラフにも残るかを確認する。
共著者側の実験候補として、時間に対する最良辺数、メモリ使用量、pivot（隣接頂点でのLCを3回まとめた操作）を用いる局所探索との比較も挙げる。最適化アルゴリズム自体はこのリポジトリで実装しない。

## 符号族ごとの構成と出典

### 共通基盤: quantumgizmos

[quantumgizmos/ldpc](https://github.com/quantumgizmos/ldpc) のGF(2)線形代数と、[quantumgizmos/bp_osd](https://github.com/quantumgizmos/bp_osd) のCSS・HGP処理を使う。

- `ldpc.mod2`: rank・nullspaceなどの行列処理。
- `bposd.css.css_code`: Hx,Hzからの論理基底計算。
- `bposd.hgp.hgp`: 古典チェック行列2個からのHGP構成。
- 格子、多項式、有限体、APM、三角形細分化は、本リポジトリの生成関数でHx,Hzへ変換する。

処理の流れは `符号族ごとの生成器 → Hx,Hzと論理Z → 論理ゼロ状態の生成元 → Stimのグラフ変換 → Stim/Qiskit検証 → 保存`。共通表現には疎行列と構成・距離のメタデータを保持する。

### 1. 回転表面符号と通常の平面表面符号

回転版は資源推定との接続とスケーリングの基準にする。`make_rotated_surface_code(d)` は任意の奇数d≥3から格子の内側・境界のチェックを生成し、`[[d²,1,d]]` を構成する。d=25,27,29を含めてStim/Qiskit検証済み。長方形のd_X,d_Z指定は未対応。

通常版は `make_planar_surface_code(d)`。`ldpc.codes.rep_code(d)` 同士を `bposd.hgp.hgp` に渡して `[[2d²−2d+1,1,d]]` を構成し、偶数距離にも対応する（[BP+OSD公式例](https://github.com/quantumgizmos/bp_osd#hypergraph-product-codes)）。d=3なら通常版は13量子ビット、回転版は9量子ビットとなる。

比較用の外部実装として [MQT QECC](https://github.com/munich-quantum-toolkit/qecc) と [qLDPC SurfaceCode](https://qldpc.readthedocs.io/autoapi/qldpc/codes/quantum/index.html#qldpc.codes.quantum.SurfaceCode) がある。

### 2. BB符号

既存Grossを基準に、具体的な低オーバーヘッドQLDPCの応用例を増やす。
応用上の根拠は [Bravyi et al., Nature 627, 778–782 (2024)](https://arxiv.org/abs/2308.07915)。

- MQTの [`construct_bb_code(n)`](https://mqt.readthedocs.io/projects/qecc/en/latest/api/mqt/qecc/codes/constructions/bb_codes/index.html) は `n=72,90,108,144,288` の既知例を提供する。これは既知例の選択APIであり、任意n・任意距離の生成器ではない。
- 一般構成にはqLDPCの [`BBCode`](https://qldpc.readthedocs.io/autoapi/qldpc/codes/quantum/index.html#qldpc.codes.quantum.BBCode) と[公式リポジトリの例](https://github.com/qLDPCOrg/qLDPC/tree/main/examples)が候補。周期サイズと2変数多項式から生成する。

実装済み: `make_bb_preset(n)` は公開5例の構成パラメータを選び、`make_bb_code(l=...,m=...,polynomial_a=...,polynomial_b=...)` は周期シフト多項式から `Hx=[A|B]`、`Hz=[B.T|A.T]` を生成する。標準notebookは公開プリセットを使い、一般多項式の実行例はコメントアウトしてある。
`distance` だけで多項式を決めることはしない。多項式と周期サイズを変えた場合、kと距離は別途確認する。
複数の既知例をサイズ順に並べても、同一の構成規則による距離系列とは限らないことを明記する。

### 3. Color符号

4.8.8と6.6.6の距離系列で、格子構造の違いを比較する。

- MQT [`SquareOctagonColorCode(distance)`](https://mqt.readthedocs.io/projects/qecc/en/latest/api/mqt/qecc/codes/constructions/color_codes/square_octagon_color_code/index.html)。
- MQT [`HexagonalColorCode(distance)`](https://mqt.readthedocs.io/projects/qecc/en/latest/api/mqt/qecc/codes/constructions/color_codes/hexagonal_color_code/index.html)。

実装済み: `make_color_code(d, lattice="4.8.8" / "6.6.6")`。MQTの格子生成アルゴリズムを `src/color_lattices.py` に取り込み、出典commit・MITライセンスを保存した。MQTのインストールは不要。任意の奇数d≥3からHx=Hzを作り、k=1と可換性を検証する。4.8.8はn=(d²+2d−1)/2、6.6.6はn=(3d²+1)/4。notebookではd=3,5,7,9,11,13を生成し、テストでは25/27/29までの構成と状態等価性、d=5の独立した低重み論理探索も確認する。
その後は実測した生成時間・メモリ量に基づいて拡張する。表面符号の資源推定からColorの必要距離を直接決めない。

### 4. Kasai Code（Zhaoインスタンス）

高レート・多数論理量子ビットを持つQLDPCとして採用する。原構成と中性原子向け派生を区別して記録する。

- 原論文: [Kasai, Breaking the Orthogonality Barrier in Quantum LDPC Codes](https://arxiv.org/abs/2601.08824)。
- 著者公開の構成プログラム: [`kasaikenta/construct_apm_css_code`](https://github.com/kasaikenta/construct_apm_css_code)。C++17の探索プログラムで、`P,J,L` と非可換条件などを入力する。READMEによれば `--save-mtx` で `HX.mtx` と `HZ.mtx` を出力できる。距離を直接入力するPython APIではない。
- 著者公開のデモ: [`kasaikenta/joint_BP_plus_PP`](https://github.com/kasaikenta/joint_BP_plus_PP)。原論文の `[[9216,4612,≤48]]` 例のパラメータを含む。まず構成・入出力形式の参考にし、9216頂点例は大規模側の追加候補とする。
- 中性原子向け派生: [Zhao et al., Towards Ultra-High-Rate Quantum Error Correction with Reconfigurable Atom Arrays](https://arxiv.org/abs/2604.16209)。収録例は `[[1152,580,≤12]]` と `[[2304,1156,≤14]]`。v2のAppendix Table A1にAPM係数が明記されており、専用GitHubを必要とせず再構成できる。本実装は表の最初のP=96とP=192の係数を使用する。原構成の別デモや同じn,kの別カタログ例を代用しない。
- [著者の符号カタログ](https://kasai.ict.eng.isct.ac.jp/kasai-code-catalog.html)にも個別パラメータと距離の記録がある。使用する場合は配布ファイル・seed・取得日を固定する。同じn,kでも別インスタンスの距離記録を流用しない。

名称は **Kasai Code**。原構成としてKasai、使用した有限サイズインスタンスの出典としてZhaoらを引用する。APM（affine permutation matrices）によるCSS構成を `make_kasai_*` APIで生成する。

現在のプリセットは**Kasai原構成を基にしたZhaoらの中性原子向け派生例**。Kasai原論文の9216量子ビット例や、Kasai著者の探索プログラムをそのまま移植したものではない。構成原典はKasai、係数・インスタンスの出典はZhao Table A1として別々に記録する。

実装済み: `make_kasai_code(modulus, f_maps, g_maps, active_rows=3)` はアフィン写像x→ax+b mod Pから疎APMを生成し、Hx,Hzの可換性とrank、論理Zを検証する。`make_kasai_preset(1152)` / `(2304)` は `src/data/kasai/zhao_table_a1.json` の係数と出典から再構成する。新規係数の探索・デコーダ性能の再現は含まない。notebookでは両例を共通のグラフ変換・Stim/Qiskit検証・保存へ接続する。
`d≤12` を `d=12` として登録しない。確定距離・下限・上限・未確認を分ける。

### 5. HGP符号

古典符号の選択・seedを変え、構造とサイズの多様性を増やす。実用性を意識した主系列は、全行rankの(3,4)-biregular古典LDPC同士の対称HGPとする。
原構成は [Tillich–Zémor](https://arxiv.org/abs/0903.0566)。高レートQLDPCの比較群にする場合は、成分にも適切な古典LDPC符号を使う。

- [`quantumgizmos/bp_osd`](https://github.com/quantumgizmos/bp_osd) の `bposd.hgp.hgp(h1=..., h2=..., compute_distance=False)` を使用する。tutorial2はbposd 2.1で検証済み。
- qLDPCの `HGPCode(code_a, code_b)` も候補（[公式チュートリアル](https://qldpc.readthedocs.io/examples/basics.html)）。

実装済み: 任意の古典チェック行列から `make_hgp_code(h1,h2)` で構成する。`make_regular_ldpc` / `make_regular_hgp_code` はbit degree 3 / check degree 4を標準とし、有限回の候補から小規模で全探索できる古典距離、4-cycle数、スペクトルギャップの順で選ぶ。生成seedと実際の古典行列を保存する。大きいkernelでは距離を未確認とする。

実用性とベンチマークの両面からの選定:

- **主系列: (3,4)-biregular HGP**。[Xu et al., Nat. Phys. 20, 1084 (2024)](https://www.nature.com/articles/s41567-024-02479-z) / [arXiv](https://arxiv.org/abs/2308.08648) の中性原子・回路レベル評価と対応する構成クラス。全行rankならN=n_c²+m_c²、K=(n_c−m_c)²、m_c=3n_c/4、レート1/25、チェック重み7。原論文は距離・girth≥6・スペクトルギャップで成分を選ぶ。本実装のランダム候補は同一行列ではなく、girth≥6や展開性・論文の論理誤り率を保証しない。
- **サイズ系列**: 古典長12/16/20/24からN=225/400/625/900、K=9/16/25/36。これはグラフ最適化を試しやすい系列であり、量子距離25〜29を満たす系列ではない。複数seedの分布も調べる。
- **公開行列の基準点: PEG / SA**。[Freire–Delfosse–Leverrier, ISIT 2025](https://arxiv.org/abs/2501.09622) の付属15×20古典行列を同梱。ともに[[625,25]]で、古典kernelの全探索と全行rankの対称HGPの距離公式によりPEGはd=6、SAはd=8と確認する。SAは消失ノイズの評価で選択されており、一般の回路ノイズでの優位性とは分ける。SAには4-cycleがあるので、girthだけで実用性能を決めない。
- **高レート比較: (3,6)**。同じ生成関数で指定可能。全行rankでレート1/5、チェック重み9となる。高レートだけを理由に主系列より良いとはしない。回路実装・距離・デコーダ性能も評価する。

長期的には既存のPEGや距離・展開性を改善した古典成分を一般行列APIへ入力する。新規探索の目的関数（消失・回路ノイズ・グラフCZ数）を明記し、QEC性能とグラフ最適化効果を別々に記録する。
反復符号同士のHGPは通常の平面表面符号の確認例になるが、それだけで高レートHGPの代表としない。

### 6. Mitten codes

非可換群に基づくlifted-product系のCSS符号。BBとは異なる代数構造を持つ比較群として採用する。
原論文は [Bhardwaj et al., High-rate qLDPC processors](https://arxiv.org/abs/2607.28795)。符号化率20%、チェック重み9の具体的な例を、論理操作・回路レベルノイズも含めて評価している。

- 著者公開実装: [a7b/yarn](https://github.com/a7b/yarn)。
- [processor_codesの形式説明](https://github.com/a7b/yarn/blob/main/processor_codes/README.md)に、`Hx.npy`、`Hz.npy`、`Lx.npy`、`Lz.npy` と論理測定用gadgetsの配布が明記されている。今回まず使うのは元の符号のチェックと論理基底。
- 収録例は `[[150,30,10]]`、`[[200,40,12]]`、`[[300,60,14]]`。
- 配布ディレクトリの距離ラベルだけで認証済みと扱わない。例えば論文では975量子ビット例を `d≤24` と記載しているため、採用版ごとに距離の根拠・上下限を確認する。

実装済み: `load_mitten_code(n)` は `src/data/mitten/` の公開行列と付属論理Zを読み込む。元ファイルのcommit・ハッシュを固定し、チェックの可換性、論理基底の独立性、付属論理X/Zのペアリングを検証する。一般の群・base matrixからの生成は未実装。

### 7. Semi-hyperbolic Color

[Higgott et al., Low-Overhead Quantum Error Correction with Boundary-Connected Planar Modules](https://arxiv.org/abs/2610.03682) のsemi-hyperbolic Color構成を使う。平面モジュールと境界の非局所接続を持ち、通常の表面・Colorと代数的QLDPCの間を比較する題材にする。
同論文はモジュール化したメモリの回路レベル評価と論理操作向けの構成を含む。そこでの物理資源削減と、今回の状態グラフの辺数削減は別の評価として扱う。

理論・構成上の利点は、基礎タイリングを固定して細分化パラメータLを増やすと、トポロジーとkを保ちつつnを増やせること。Color構成では `n=L²×基礎三角形数` で、距離はおおむねLとともに増える（§S2）。同じ大域構造を保った系列と、同じLで基礎構造を変える系列を分ける。

標準サンプルは、同じ8.8.8の基礎構造から得られるColor系列（同論文Tables S7–S9）:

| 細分化L | パラメータ | 用途 |
|---|---|---|
| 2 | `[[64,8,6]]` | 小規模で導入・検証 |
| 3 | `[[144,8,8]]` | Grossと同じ頂点数で比較 |
| 4 | `[[256,8,12]]` | 同一基礎構造のサイズ拡張 |

実装済み: `make_semi_hyperbolic_color_code(refinement=L, base_tiling="bolza16")`。L≥1の各基礎三角形をL²個へ細分化し、双対ColorのチェックHx=Hzを新規生成する。隣接・辺・角の同一視を整数のcell keyで処理し、Bolzaの並行辺・重複する角頂点組を潰さない。n=16L²、k=8。標準notebookはL=1〜5を含み、テストはL=6まで生成する。L=1/2/3/4のみ文献距離2/6/8/12を付け、それ以外の距離は未確認とする。L=2の距離6は、独立したmeet-in-the-middle探索でも確認する。

基礎タイリングは、同論文S2で共用すると記される2024年Floquet論文の [著者公開データ](https://github.com/oscarhiggott/hyperbolic-floquet-data) にある16量子ビットBolza例から、X/Y/Zの色付き辺の隣接関係だけを復元した。Floquetの測定チェック・論理基底・距離を静的CSSへ代用していない。commit・原ファイルhash・数値データhash・Apacheライセンスを `src/data/semi_hyperbolic/` に保存する。基礎の色付きグラフは、3つの反可換なinvolutionで生成される16元群のCayleyグラフとも独立に同型確認済み。

カスタム基礎は `base_tiling={"name": ..., "neighbors": ...}` で指定する。neighbors[t][i]を三角形tの色iの角の対辺を越えた隣と定義し、各色がfixed-point-free involution、全体が連結、角リンクが偶数cycle、Euler characteristic<0であることを検証する。genus・チェック色・qubitの所属する基礎三角形（K=1モジュール）も保存する。カスタム入力にはプリセットの距離を流用しない。
初期グラフにモジュール構造や対称性が残るか、それがLC削減率や対称性削減の効果に関係するかは実験仮説とする。

### 8. Quantum Reed–Muller：小規模な魔法状態蒸留用符号

旧 `qec_codes.py` はMQTのTetrahedral `[[15,1,3]]` を既に含む。新APIはこれを固定例から、非ゼロ2進評価点上のshortened RMによる生成関数に拡張する。共通CSS処理はldpc/bposdで、MQT依存は不要。

`make_quantum_reed_muller_code(m, order=r)` はm>2r≥2から `[[2^m−1,1,2^(r+1)−1]]` を構成する。[Hastings–Haah, Distillation with sublogarithmic overhead](https://arxiv.org/abs/1709.03543) のw=0族に対応する。Xチェックはshortened RM(r,m)、Zチェックはshortened RM(m−r−1,m)。

標準notebookの例:

| m | r | パラメータ | 横断位相の最大保証レベルν=floor((m−1)/r) |
|---|---|---|---|
| 3 | 1 | [[7,1,3]] | 2（Clifford） |
| 4 | 1 | [[15,1,3]] | 3（15-to-1 T蒸留で重要） |
| 5 | 1 | [[31,1,3]] | 4 |
| 6 | 1 | [[63,1,3]] | 5 |
| 5 | 2 | [[31,1,7]] | 2（この例では横断Tを主張しない） |
| 7 | 2 | [[127,1,7]] | 3 |

物理一様Rν=diag(1,exp(2πi/2^ν))は論理Rν逆を実装する。ν≥3ならTも横断的に実装可能。テストではr=1系列の全Xチェック語に対して、論理0/1の計算基底成分が正しい位相を得ることを直接確認する。

距離・横断ゲート・サイズを制御できる代数的比較群として有用。ただし今回の最適化対象は論理ゼロ状態の準備グラフであり、非Clifford魔法状態自体や15-to-1蒸留回路のCZ数・歩留まり・出力誤り率とは同一視しない。

### 9. Quantum Reed–Solomon：Quantum MDS親符号とqubit展開

[Grassl–Beth–Rötteler, On optimal quantum codes](https://arxiv.org/abs/quant-ph/0312164) の長さn≤qのQuantum MDSクラスを、GF(q=2^m)上の自己直交GRS CSSで構成する。

`make_quantum_reed_solomon_code(field_degree=m, distance=d, length=n)`：m=2〜8、3≤n≤q、d≥2、k=n−2(d−1)>0。length省略時はn=q。評価点a_i=0..n−1とw_i=1/∏(a_i−a_j)、v_i=sqrt(w_i)から、Hの各行をv_i a_i^j（j=0..d−2）とする。H H^T=0を満たし、親符号は `[[n,n−2d+2,d]]_q`。Quantum Singleton boundを飽和する。有限体演算は小規模な多項式基底実装で、既約多項式と評価点・列係数・元のq元チェックを保存する。

現在のStim/Qiskitはqubit用なので、1 quditをm qubitへ展開する。X(a)はaの多項式座標、Z(b)はTr(b x^j)をZ座標に用いる。X/Z両方を単純に同じビット列にする写像は一般に可換性を保存しないため採用しない。

展開後のパラメータはN=mn、K=m(n−2d+2)、d_binary≥d。親のqudit距離は確定値・MDSとしてmetadataに保存し、`CodeInstance.distance` はqubit距離の下界として `lower_bound` を付ける。qubitのMDS符号と誤記しない。

標準系列はGF(4/8/16/32/64)上の複数長・距離。GF(64)の `[[64,16,25]]_64` と `[[64,8,29]]_64` は、それぞれqubit展開でN=384、K=96/48、距離下界25/29となる。ただし表面符号の距離目標や論理誤り率・耐故障資源をそのまま対応させるものではない。密な代数構成とLDPC構成のグラフ最適化を比較する用途に使う。

## 生成API

以下は実装済みの公開API。詳細な引数・検証条件は [src/qec_families.py](src/qec_families.py) を参照。

```python
# 距離・構成パラメータから生成する。
make_rotated_surface_code(distance: int) -> CodeInstance
make_planar_surface_code(distance: int) -> CodeInstance
make_color_code(distance: int, *, lattice: str) -> CodeInstance
make_quantum_reed_muller_code(num_variables=4, *, order=1) -> CodeInstance
make_quantum_reed_solomon_code(field_degree: int, distance: int, *, length=None) -> CodeInstance

# 行列・多項式・APMなどの構成パラメータを受け取る。
make_bb_code(*, l: int, m: int, polynomial_a, polynomial_b) -> CodeInstance
make_bb_preset(n: int) -> CodeInstance
make_hgp_code(check_a, check_b) -> CodeInstance
make_regular_hgp_code(num_bits: int, *, bit_degree=3, check_degree=4, seed=0) -> CodeInstance
make_hgp_preset(method="PEG") -> CodeInstance
make_kasai_code(modulus: int, f_maps, g_maps, *, active_rows=3) -> CodeInstance
make_kasai_preset(n=1152) -> CodeInstance
make_semi_hyperbolic_color_code(refinement: int, *, base_tiling="bolza16") -> CodeInstance
load_mitten_code(n=150) -> CodeInstance
load_css_code(path, *, name=None, family="imported", metadata=None,
              distance=None, distance_status="unknown") -> CodeInstance

# 構成と、特定の論理状態を選ぶ処理を分ける。
logical_zero_stabilizers(code: CodeInstance) -> list[str]
make_graph_sample(code: CodeInstance, *, verify_qiskit=True) -> GraphSample
export_graph_sample(sample: GraphSample, output_dir) -> Path
```

- `CodeInstance` はチェック行列、論理演算子、n,k、構成・距離のメタデータを保持する。生の行列を保持し、グラフ生成前にも調査できるようにする。
- 外部のMQT/qLDPC/bposdオブジェクトは共通表現へ変換する。共通CSS処理はldpc/bposdを優先する。既存の `_css_code_stabilizers` はMQTの属性・メソッドに依存しているため、他ライブラリをそのまま渡す前提にしない。
- `CODES` と `get_stabilizers(key)` は、旧tutorialのTable I再現用APIとして維持する。
- 距離ごとにfactoryや行列を追加する方式を主実装にしない。既知例のプリセットは、一般生成関数へ渡すパラメータ集合として持つ。
- Mittenは公開行列を読み込む。格子・多項式・APMからの一般生成と、公開インスタンスの読込みを区別する。
- 距離不明の族で `distance=27` と指定すると適当にパラメータを選ぶAPIは作らない。将来探索機能を作るなら、目標距離と認証済み距離を区別した別処理とする。
- 行列から独立な論理Zを求め、チェックと合わせて状態を一意に定める。多数の論理量子ビットを持つKasai Codeでも同じ定義を使う。

## 検証と共著者に渡すデータ

実装済みの確認:

- GF(2)上で `Hx @ Hz.T == 0`、`n-rank(Hx)-rank(Hz)==k` を確認する。
- 論理ZがXチェックと可換で、Zチェックを法としてk個独立であることを確認する。最終生成元がn個独立で、矛盾する符号を含まないことも確認する。
- 生成した状態をStim/Qiskitで検証する。小規模例に加え、表面符号の主目標d=25,27,29と採用した新規例で確認する。
- 状態等価性・rankの検証と、符号距離の認証は別。前者の成功を距離の証明とは扱わない。

保存する項目:

- 符号族、n,k、距離の値/下限/上限とその根拠、論文・配布元。
- 構成パラメータ、seed、元ファイルのハッシュ、依存ライブラリのバージョン/commit。
- Hx,Hz、採用した論理Zと状態生成元、その順序、量子ビットの順序。
- グラフの辺リスト、局所Clifford、準備回路、等価性検証結果。
- 辺数、次数、チェック重み、状態生成・変換・検証時間など。元のチェックが疎でも、変換後のグラフが疎とは限らない。
- 二部グラフか、三角形数、一手のLCで改善できる頂点数など、局所探索の効きやすさを説明する初期グラフの指標。
- Semi-hyperbolicでは基礎タイリングの識別子、細分化L、利用できる場合は元のモジュール分割。Mitten/Kasai Codeでは群・生成元・公開基底の識別情報。

同じ状態でも基底・順序・ライブラリ版によって初期グラフが変わりうるため、生成条件を固定する。
頂点数だけでなく、符号族別・サイズ別に集計できるデータを渡す。

## 追加評価・拡張

- 複数seedのHGPで、同じサイズでも初期グラフと最適化効果がどれだけ変わるか比較する。
- Semi-hyperbolic Colorは同一基礎タイリングのL系列を主比較にし、別の基礎タイリングは構造差を調べる追加比較に使う。
- 実行規模は状態生成・グラフ変換・最適化の実測時間とメモリ量に基づいて決める。ピークメモリと各段階の個別時間は追加の測定対象。
- 最適化後のグラフ・LC列が得られたら、元の論理ゼロ状態との等価性を検証する。
- Semi-hyperbolic Surface、新規APM係数の探索、Mittenの一般構成生成器は未実装。量子BCH・Icebergは後続候補。

実行可能な標準サンプルは [src/tutorial2.ipynb](src/tutorial2.ipynb)、回帰テストは `tests/` で管理する。
