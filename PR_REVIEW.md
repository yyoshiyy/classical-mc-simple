# Pull Request Review: Feature/ner mpi #3

## 総合評価: ✅ Approve（軽微な提案あり）

NER (Non-Equilibrium Relaxation) モードと MPI 並列化を追加する PR です。実装はよく整理されており、既存の equilibrium モードとの分離も明確です。ビルドも成功し、設計方針も妥当です。

---

## 良い点 👍

### 1. アーキテクチャ
- **NER と equilibrium の分岐が明確**: `enable_ner` による分岐が main.c で分かりやすく、責務が分離されている
- **lattice の初期化順序の整理**: `build_lattice` を main で呼び、`initial()` 内から削除したことで、CLI から interaction ファイルパスを渡せるようになっている
- **MPI の optional ビルド**: `USE_MPI` オプションで従来の単一プロセスビルドを維持しつつ、MPI を選択可能

### 2. NER モードの実装
- **initial_ner**: FM, AF, Stripe1/2, RANDOM など複数の初期状態に対応
- **Ini_sx/Ini_sy/Ini_sz**: overlap q(t) = S(t)·S(0) の計算に必要な初期配置を適切に保存
- **時系列出力**: m(t), m²(t), e_per_site(t), q(t) を NER_result.dat に出力

### 3. MPI 並列化
- **independent-run 並列化**: 各 rank が異なる RNG seed で独立実行し、MPI_Reduce で集約する設計は妥当
- **MPI_SEED_STRIDE**: ランク間で seed が重ならないよう適切にオフセット
- **MASTER のみの出力**: printf/fprintf を rank 0 に限定して重複を防止

### 4. Python スクリプト
- **run_ner_and_plot.py**: 計算→プロット→臨界指数推定を一括実行するワークフローが便利
- **analyze_ner_critical.py**: log-log 空間でのべき乗減衰フィットと z 推定が実装されている
- **scripts/README.md**: 各スクリプトの使い方が明記されている

### 5. ドキュメント
- README に NER モードと MPI ビルド手順が追記されている
- param_ner.def, param_ner_tc.def のサンプルが用意されている

---

## 軽微な指摘・提案

### 1. input_parser.c: run_mode の大文字小文字
`run_mode = ner` は `strcmp(val, "ner")` で完全一致のため、`run_mode = NER` と書くと NER モードにならない。README に「小文字の `ner` を指定」と明記するか、`strcasecmp` で大文字小文字を無視することを検討してください。

### 2. main.c: NER メモリ解放時のエラーハンドリング
```c
if (ner_m == NULL || ner_m2 == NULL || ner_e == NULL || ner_q == NULL) {
    fprintf(stderr, "Error: NER accumulator allocation failed\n");
    free(ner_m);  // NULL の可能性あり → free(NULL) は安全だが、部分的に確保済みの場合の解放順序は問題なし
    ...
}
```
`free(NULL)` は C では安全なので問題ありませんが、`calloc` 失敗時は複数が NULL になり得るため、現状の実装で問題ありません。

### 3. バイナリ・データファイルのコミット
`samples/square_L16_Ising/NER_result.dat`, `NER_plot.png`, `NER_critical_fit.png` がリポジトリに含まれています。サンプル出力として残す意図であれば問題ありませんが、再現性の観点から `.gitignore` に含め、必要に応じて `run_ner_and_plot.py` で生成する運用も検討できます。

### 4. analyze_ner_critical.py: q(t) のフィット
`mask_q = (t >= args.t_min) & (t <= t_max) & (q > 0.01)` でチェックしている一方、`fit_power_law` 内では `y > 0` でマスクしています。q が負になり得る overlap では、現状の `y > 0` で十分ですが、「q>0 の範囲のみ」というコメントと実装の対応は明確です。

---

## 動作確認

- デフォルトビルド（NO_MPI）: 成功
- 変更ファイル: 16 ファイル、+938/-40 行

---

## 結論

マージして問題ない内容です。上記の軽微な点は、必要に応じて後続の PR や issue で対応可能です。

---

*This review was conducted by **Composer 1.5** (Cursor AI).*
