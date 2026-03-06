# 物理量の計算と統計平均

このドキュメントでは、MC_simple で計算される各物理量の定義・計算方法と、統計平均の取り方を説明する。

---

## 1. 各物理量の定義と計算

### 1.1 E per site（サイトあたりエネルギー）

**定義**:
$$
\frac{E}{N} = \frac{\mathcal{H}}{N}
$$

**ハミルトニアン**:
$$
\mathcal{H} = J \sum_{\langle i,j \rangle} \mathbf{S}_i \cdot \mathbf{S}_j - H \sum_i S_i^z
$$

**計算方法**:
- 全エネルギー \(E\) は隣接テーブルから計算。テーブルに i→j と j→i の両方を含むため、和を 2 で割る:
  $$
  E = \frac{1}{2} \sum_i \sum_{j \in \text{neighbors}(i)} J_{ij} \, \mathbf{S}_i \cdot \mathbf{S}_j - H \sum_i S_i^z
  $$
- MC 更新時は差分 \(\Delta E\) のみ計算し、\(E \leftarrow E + \Delta E\) で更新
- 各 MC ステップ \(t\) で \(E(t)/N\) を記録

---

### 1.2 M²（磁化の二乗）

**定義**:
$$
M^2 = \frac{1}{N^2} \left( M_x^2 + M_y^2 + M_z^2 \right)
$$

**計算方法**:
$$
M_\alpha = \sum_{i=1}^{N} S_i^\alpha \quad (\alpha = x, y, z)
$$
$$
M^2(t) = \frac{M_x^2 + M_y^2 + M_z^2}{N^2}
$$

- **Ising** (\(d=1\)): \(M_x = M_y = 0\) なので \(M^2 = M_z^2 / N^2\)
- **XY** (\(d=2\)): \(M_z = 0\)
- **Heisenberg** (\(d=3\)): 3 成分すべて使用

---

### 1.3 Overlap

**定義**:
$$
q(t) = \frac{1}{N} \sum_{i=1}^{N} \mathbf{S}_i(t) \cdot \mathbf{S}_i(0)
$$

**計算方法**:
- measurement 開始時（burn-in 直後）の配置 \(\mathbf{S}_i(0)\) を保存
- 各 MC ステップ \(t\) で:
  $$
  q(t) = \frac{1}{N} \sum_i \left( S_i^x(t) S_i^x(0) + S_i^y(t) S_i^y(0) + S_i^z(t) S_i^z(0) \right)
  $$

---

### 1.4 C per site（サイトあたり比熱）

**定義**:
$$
C = \frac{\langle E^2 \rangle - \langle E \rangle^2}{N \, T^2}
$$

**計算方法**:
- 各ステップで \(E(t)\) と \(E^2(t)\) を記録
- 時間平均 \(\langle E \rangle\), \(\langle E^2 \rangle\) を計算して上式に代入

---

### 1.5 acceptance（Metropolis 受理率）

**定義**:
$$
\text{acceptance} = \frac{\text{受理した更新数}}{N \times \text{Total\_Step}}
$$

---

## 2. 統計平均の取り方

### 2.1 二段階の平均

MC_simple では **二段階** で統計平均を取る。

| 段階 | 内容 |
|------|------|
| **時間平均** | 1 回の run 内で、Total_Step ステップにわたって平均 |
| **Sample 平均** | Sample 回の独立 run の平均 |

---

### 2.2 時間平均（1 run 内）

各 run で、measurement フェーズの各ステップ \(t = 1, \ldots, \text{Total\_Step}\) について観測量 \(O(t)\) を記録し、その平均を取る:

$$
\bar{O}^{(\text{sample})} = \frac{1}{\text{Total\_Step}} \sum_{t=1}^{\text{Total\_Step}} O(t)
$$

- \(O\): \(E/N\), \(M^2\), \(q\) など
- 各 run は 1 つの \(\bar{O}^{(\text{sample})}\) を生成

---

### 2.3 Sample 平均（独立 run の平均）

各 run は **異なる乱数シード** で独立に実行される:

$$
\text{seed}^{(s)} = \text{base\_seed} + 8945 \times s \quad (s = 0, 1, \ldots, \text{Sample}-1)
$$

最終的な熱力学的平均は:

$$
\langle O \rangle = \frac{1}{\text{Sample}} \sum_{s=1}^{\text{Sample}} \bar{O}^{(s)}
$$

---

### 2.4 Sample の役割

| パラメータ | 意味 |
|------------|------|
| **Sample** | 独立な MC run の回数 |
| **Total_Step** | 各 run の measurement ステップ数 |

- Sample を増やすと、**統計誤差**（run 間のばらつき）が減る
- Total_Step を増やすと、**時間平均の精度**が上がる（相関の影響も含む）

---

## 3. 計算フロー概要

```
for s = 1 to Sample:
    seed = base_seed + 8945 * s
    初期化（ランダムスピン）
    
    Burn-in: Total_Step ステップ（測定なし）
    
    S(0) を保存（overlap 用）
    
    for t = 1 to Total_Step:
        MC スイープ（全温度）
        Exchange MC（隣接温度間の交換試行）
        E(t), E²(t), M²(t), q(t) を蓄積
    
    sample_E[s] = (1/Total_Step) × Σ E(t)/N
    sample_M2[s] = (1/Total_Step) × Σ M²(t)
    sample_overlap[s] = (1/Total_Step) × Σ q(t)

最終出力:
    ⟨E/N⟩ = (1/Sample) × Σ sample_E[s]
    ⟨M²⟩  = (1/Sample) × Σ sample_M2[s]
    ⟨q⟩   = (1/Sample) × Σ sample_overlap[s]
```

---

## 4. strong プリセット時の例

| パラメータ | 値 |
|------------|-----|
| Sample | 5 |
| Burn_in | 5000 |
| Total_Step | 10000 |

- 1 run あたり: 5000 (burn-in) + 10000 (measurement) = 15000 ステップ
- 5 run 合計: 75,000 MC ステップ / 温度
- 最終値は 5 本の独立 run の平均
