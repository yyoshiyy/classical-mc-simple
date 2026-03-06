# scripts/

## run_ner_and_plot.py

NER モードで計算を実行し、完了後に自動でプロットする。

```bash
python scripts/run_ner_and_plot.py
```

- `--sample DIR`: サンプルディレクトリ (デフォルト: square_L16_Ising)
- `--param FILE`: param ファイル (デフォルト: param_ner.def)
- `--no-save`: プロットをファイルに保存せず表示のみ
- `--no-analyze`: 臨界指数推定をスキップ（デフォルトは自動実行）

## analyze_ner_critical.py

NER_result.dat から臨界指数を推定。m(t) ~ A*t^(-α) をフィットし、α = β/(νz) から動的臨界指数 z を推定。

```bash
python scripts/analyze_ner_critical.py samples/square_L16_Ising/NER_result.dat
python scripts/analyze_ner_critical.py NER_result.dat --short-time --plot
```

- `--t-min`, `--t-max`: フィット範囲
- `--short-time`: 短時間領域のみフィット（推奨）
- `--plot`: フィット結果をプロット

※ べき乗減衰は T≈Tc で成立。`param_ner_tc.def` (Ini_T=2.27) を使用。

## plot_ner_result.py

NER_result.dat をプロットする（単体実行用）。

```bash
cd samples/square_L16_Ising
python ../../scripts/plot_ner_result.py -o NER_plot.png
```

## benchmark_exchange_mc.py

Exchange MC の効果をエラーバー比較で検証するスクリプト。

### 使い方

```bash
# リポジトリルートで実行
python scripts/benchmark_exchange_mc.py

# 効果をはっきり見たい場合（L=32, 低温含む, 長い burn-in）
python scripts/benchmark_exchange_mc.py --preset strong
```

### オプション

- `--n-runs N`: 各モード（Exchange MC 有/無）での実行回数（デフォルト: 20）
- `--preset {default|strong}`:
  - `default`: L=16, Burn_in=2000, Total_Step=5000（短時間）
  - `strong`: L=32, T=1.2〜2.6, Burn_in=5000, Total_Step=10000（効果が顕著になりやすい）
- `-o FILE`: 出力プロットファイル名
- `--no-plot`: プロットをスキップし、数値のみ出力

### 結果の見方

- **エラーバー比 < 1**: Exchange MC で統計誤差が減少
- **平均エラーバー比**: 全温度・全観測量の幾何平均。<1 なら全体的に効果あり
- **比が 1 前後でばらつく場合**: L=16 では critical slowing down が弱く、差が出にくい。`--preset strong` を試す

### 必要なもの

- numpy
- matplotlib（プロット時のみ）

### 出力

- `samples/exchange_mc_benchmark/{preset}/exchange_mc_benchmark.png`: E, C, M², overlap の比較プロット
- 標準出力: 各温度での mean ± stderr、overlap 比較、エラーバー比

### overlap による EXMC 効果

- **overlap** = t=0 と t の間の類似度 `(1/N)Σ S_i(t)·S_i(0)`
- 高い overlap = 初期状態からあまり離れていない（critical slowing down）
- **Exchange MC ありで overlap が低い** = 温度交換で phase space をよく探索できている
