# scripts/

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
