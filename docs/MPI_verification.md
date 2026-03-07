# MPIモード動作確認 (Issue #4)

## 環境

- macOS (arm64), Homebrew
- Open MPI: `brew install open-mpi` (5.0.9)

## 手順

### 1. MPIビルド

```bash
cd src
cmake -S . -B build -DUSE_MPI=ON
cmake --build build -j
```

### 2. NERサンプル実行

```bash
cd samples/square_L16_Ising
mpirun -np 2 ../../src/build/MC_simple param_ner.def
mpirun -np 4 ../../src/build/MC_simple param_ner.def
```

## 確認結果

- **ビルド**: `USE_MPI=ON` でビルド成功（MPI_C 検出済み）
- **実行**: `mpirun -np 2` / `-np 4` で正常終了
- **出力**: `NER_result.dat` が生成され、列 `t, T, m, m2, e_per_site, q` が期待どおり
- **集約**: 各 rank が異なる seed で独立実行し、MASTER のみが集約後の平均をファイル出力

以上で MPI モードの動作確認を完了した。
