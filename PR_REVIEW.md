# PR コードレビュー: exmc update by composer-1.5 #1

## 総評

Exchange MC（レプリカ交換法）の実装と overlap の追加は、物理的にも実装としても妥当です。全体としてよく整理されていますが、いくつか修正・改善を推奨する点があります。

---

## 良い点

### 1. Exchange MC の受理確率

`attempt_exchange` の受理確率
$$P = \min\left(1, \exp\left[\left(\frac{1}{T_i} - \frac{1}{T_j}\right)(E_i - E_j)\right]\right)$$
は、レプリカ交換法の標準的な式と一致しています。

### 2. Overlap の定義と計算

`compute_overlap` の定義
$$q(t) = \frac{1}{N}\sum_i \mathbf{S}_i(t)\cdot\mathbf{S}_i(0)$$
は、`docs/observables_and_averaging.md` と整合しています。`copy_spins_to_init` で burn-in 直後の配置を保存し、測定開始時点を t=0 としているのも妥当です。

### 3. メモリ管理

`prev_sx/prev_sy/prev_sz` と `sample_overlap/accum_overlap` の allocate/free が対応しており、メモリリークは見当たりません。

### 4. ベンチマークスクリプト

`benchmark_exchange_mc.py` は、Exchange MC の有無による比較がしやすく、`--preset strong` や `--no-plot` などオプションも整理されています。

---

## 要修正・改善点

### 1. シンボリックリンクの絶対パス（重要）

`samples/exchange_mc_benchmark/default/` と `strong/` 内の `lattice.def`、`interaction.def` が絶対パスを指しています。

```
lattice.def -> /Users/yuichiro/.../samples/square_L16_Ising/lattice.def
```

他環境や別パスにクローンするとリンクが壊れます。

**推奨**: シンボリックリンクをリポジトリに含めず、`benchmark_exchange_mc.py` 実行時に毎回作成する運用にするか、相対パスでリンクするように変更してください。既存の絶対パスリンクは削除して、`.gitignore` に `samples/exchange_mc_benchmark/*/lattice.def` などを追加するのも一案です。

### 2. `memory.h` のインデント揺れ

```c
double *accum_E;      /* Sum over samples of per-sample <E/N> */
double *accum_C;     /* Sum over samples of per-sample specific heat */
double *accum_M2;    /* Sum over samples of per-sample <M^2> */
```

`accum_C` と `accum_M2` のコメント前のスペース数が他と揃っていません。スタイルを統一すると読みやすくなります。

### 3. ベンチマーク出力のコミット可否

`samples/exchange_mc_benchmark/default/exchange_mc_benchmark.png` や `MC_simple_result.dat` などは実行結果のため、再現性の観点からは除外する方がよい場合があります。

**推奨**: `.gitignore` に `samples/exchange_mc_benchmark/*/MC_simple_result.dat` や `samples/exchange_mc_benchmark/*/exchange_mc_benchmark.png` を追加するか、サンプル出力はリポジトリに含めない方針を検討してください。

---

## 軽微な指摘・提案

### 4. `attempt_exchange` の `exp` オーバーフロー

`(1/T_i - 1/T_j) * (E_i - E_j)` が大きいと `exp` がオーバーフローする可能性があります。現状は `prob > 1.0` で `prob = 1.0` にしているため、`inf` は 1 に丸められており、受理判定には影響しません。ただし、より安全にするなら `exp` の前に `exp` の引数が閾値（例: 700）を超える場合は `prob = 1.0` に設定する処理を入れることも検討できます。

### 5. Python の `parse_result_file` の列数チェック

`len(parts) >= 4` で T, E, C, M2 を読み、`len(parts) >= 6` で overlap を読みます。現在の出力は常に 6 列なので問題ありませんが、`len(parts) >= 5` の列（acceptance）がある場合のフォーマットをコメントで明示しておくと、将来の変更時に分かりやすくなります。

### 6. `samples/square_L16_Ising/param.def.bak`

`param.def.bak` に `enable_exchange = 0` のみが含まれており、用途が不明です。不要であれば削除するか、意図があれば README に説明を追記することを推奨します。

### 7. `README.md` の `enable_exchange` の記述

`enable (default)` とある一方で、デフォルト値の説明がやや短いです。`param.def` に書かない場合のデフォルトが 1 であることを明示すると、より分かりやすくなります。

---

## 動作確認の提案

以下を確認してください。

1. `enable_exchange = 0` と `1` で、同じ乱数シードで E, C, M² が期待される挙動になるか
2. `scripts/benchmark_exchange_mc.py --n-runs 2` が正常に完了するか
3. `spin_dim = 2`（XY）や `spin_dim = 3`（Heisenberg）でも Exchange MC が正常に動作するか

---

## まとめ

| 項目 | 評価 |
|------|------|
| 物理的妥当性 | 良好 |
| 実装の一貫性 | 良好 |
| メモリ管理 | 妥当 |
| ドキュメント | 十分 |
| 移植性 | シンボリックリンクの絶対パス要修正 |

**結論**: シンボリックリンクの絶対パス問題を修正すれば、マージしてよい内容だと思います。その他の修正は任意です。
