#!/usr/bin/env python3
"""
Exchange MC 効果検証スクリプト

Exchange MC の有無でシミュレーションを複数回実行し、
エラーバー（標準誤差）を比較してプロットする。

使い方:
  cd classical-mc-simple
  python scripts/benchmark_exchange_mc.py

  # 効果をはっきり見たい場合（L=32, 低温含む, 長い burn-in）:
  python scripts/benchmark_exchange_mc.py --preset strong

  または:
  python scripts/benchmark_exchange_mc.py --n-runs 10  # 高速テスト用
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

# matplotlib はオプション（なければプロットをスキップ）
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# リポジトリルート（このスクリプトの親の親）
REPO_ROOT = Path(__file__).resolve().parent.parent
BUILD_DIR = REPO_ROOT / "src" / "build"
MC_EXEC = BUILD_DIR / "MC_simple"

# プリセット: Exchange MC の効果が顕著になりやすい条件
PRESETS = {
    "default": {
        "sample_dir": "square_L16_Ising",
        "num_temp": 5,
        "ini_t": 1.5,
        "delta_t": 0.35,
        "burn_in": 2000,
        "total_step": 5000,
        "sample": 5,
    },
    "strong": {
        # L=32: critical slowing down が顕著
        # 低温 T=1.2 含む: 通常 MC が動きにくい
        # 長い burn-in: 低温での平衡化が困難
        "sample_dir": "square_L32_Ising",
        "num_temp": 5,
        "ini_t": 1.2,
        "delta_t": 0.35,  # T = 1.2, 1.55, 1.9, 2.25, 2.6 (Tc≈2.27 付近を含む)
        "burn_in": 5000,
        "total_step": 10000,
        "sample": 5,
    },
}


def make_param_content(enable_exchange: int, num_temp: int = 5,
                       ini_t: float = 1.5, delta_t: float = 0.35,
                       burn_in: int = 2000, total_step: int = 5000,
                       sample: int = 5) -> str:
    """param.def の内容を生成"""
    return f"""# Benchmark: Exchange MC comparison
Burn_in    = {burn_in}
Total_Step = {total_step}
Sample     = {sample}
num_temp   = {num_temp}
Ini_T      = {ini_t}
Delta_T    = {delta_t}
lambda     = 1.0
H          = 0.0
spin_dim   = 1
enable_exchange = {enable_exchange}
"""


def parse_result_file(result_path: Path) -> dict:
    """MC_simple_result.dat をパースして T, E, C, M2, overlap を返す"""
    data = {"T": [], "E": [], "C": [], "M2": [], "overlap": []}
    with open(result_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                data["T"].append(float(parts[0]))
                data["E"].append(float(parts[1]))
                data["C"].append(float(parts[2]))
                data["M2"].append(float(parts[3]))
                data["overlap"].append(float(parts[5]) if len(parts) >= 6 else float("nan"))
    return data


def run_single(seed: int, enable_exchange: int, param_path: Path,
               lattice_path: Path, interaction_path: Path,
               exec_path: Path, cwd: Path) -> dict:
    """1回のシミュレーションを実行し、結果を返す"""
    env = os.environ.copy()
    env["MC_SIMPLE_SEED"] = str(seed)

    result = subprocess.run(
        [str(exec_path), str(param_path), str(lattice_path),
         str(interaction_path)],
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"MC_simple failed with seed={seed}")

    return parse_result_file(cwd / "MC_simple_result.dat")


def collect_runs(n_runs: int, enable_exchange: int, param_path: Path,
                 lattice_path: Path, interaction_path: Path,
                 exec_path: Path, cwd: Path) -> dict:
    """n_runs 回実行し、各温度での mean ± stderr を計算"""
    all_T = []
    all_E = []
    all_C = []
    all_M2 = []

    all_overlap = []
    for run in range(n_runs):
        seed = 10000 + run * 1007  # 異なるシード
        d = run_single(seed, enable_exchange, param_path, lattice_path,
                       interaction_path, exec_path, cwd)
        all_T.append(d["T"])
        all_E.append(d["E"])
        all_C.append(d["C"])
        all_M2.append(d["M2"])
        all_overlap.append(d["overlap"])

    # 温度は全 run で同じ
    T = np.array(all_T[0])
    E = np.array(all_E)
    C = np.array(all_C)
    M2 = np.array(all_M2)
    overlap = np.array(all_overlap)

    return {
        "T": T,
        "E_mean": E.mean(axis=0),
        "E_err": E.std(axis=0) / np.sqrt(n_runs) if n_runs > 1 else np.zeros_like(T),
        "C_mean": C.mean(axis=0),
        "C_err": C.std(axis=0) / np.sqrt(n_runs) if n_runs > 1 else np.zeros_like(T),
        "M2_mean": M2.mean(axis=0),
        "M2_err": M2.std(axis=0) / np.sqrt(n_runs) if n_runs > 1 else np.zeros_like(T),
        "overlap_mean": overlap.mean(axis=0),
        "overlap_err": overlap.std(axis=0) / np.sqrt(n_runs) if n_runs > 1 else np.zeros_like(T),
    }


def plot_results(data_no_ex, data_with_ex, out_path: Path):
    """E, C, M2, overlap をプロット"""
    if not HAS_MATPLOTLIB:
        print("matplotlib がインストールされていません。プロットをスキップします。")
        return

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    # E_per_site
    ax = axes[0, 0]
    ax.errorbar(data_no_ex["T"], data_no_ex["E_mean"], yerr=data_no_ex["E_err"],
                fmt="o-", capsize=3, label="no Exchange MC")
    ax.errorbar(data_with_ex["T"], data_with_ex["E_mean"],
                yerr=data_with_ex["E_err"],
                fmt="s--", capsize=3, label="with Exchange MC")
    ax.set_xlabel("T")
    ax.set_ylabel("E per site")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # C_per_site
    ax = axes[0, 1]
    ax.errorbar(data_no_ex["T"], data_no_ex["C_mean"], yerr=data_no_ex["C_err"],
                fmt="o-", capsize=3, label="no Exchange MC")
    ax.errorbar(data_with_ex["T"], data_with_ex["C_mean"],
                yerr=data_with_ex["C_err"],
                fmt="s--", capsize=3, label="with Exchange MC")
    ax.set_xlabel("T")
    ax.set_ylabel("C per site")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # M2
    ax = axes[1, 0]
    ax.errorbar(data_no_ex["T"], data_no_ex["M2_mean"],
                yerr=data_no_ex["M2_err"],
                fmt="o-", capsize=3, label="no Exchange MC")
    ax.errorbar(data_with_ex["T"], data_with_ex["M2_mean"],
                yerr=data_with_ex["M2_err"],
                fmt="s--", capsize=3, label="with Exchange MC")
    ax.set_xlabel("T")
    ax.set_ylabel("M²")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # overlap（低い = よく mix、EXMC の効果）
    ax = axes[1, 1]
    ax.errorbar(data_no_ex["T"], data_no_ex["overlap_mean"],
                yerr=data_no_ex["overlap_err"],
                fmt="o-", capsize=3, label="no Exchange MC")
    ax.errorbar(data_with_ex["T"], data_with_ex["overlap_mean"],
                yerr=data_with_ex["overlap_err"],
                fmt="s--", capsize=3, label="with Exchange MC")
    ax.set_xlabel("T")
    ax.set_ylabel("overlap (consecutive steps)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"プロットを保存: {out_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Exchange MC の効果をエラーバー比較で検証"
    )
    parser.add_argument("--n-runs", type=int, default=20,
                       help="各モードでの実行回数（デフォルト: 20）")
    parser.add_argument("--preset", choices=list(PRESETS), default="default",
                       help="default: L=16, 短い burn-in. strong: L=32, 低温含む, 長い burn-in（効果が顕著）")
    parser.add_argument("--output", "-o", type=str,
                       default="exchange_mc_benchmark.png",
                       help="出力プロットファイル名")
    parser.add_argument("--no-plot", action="store_true",
                       help="プロットをスキップし、数値のみ出力")
    args = parser.parse_args()

    preset = PRESETS[args.preset]
    sample_dir = REPO_ROOT / "samples" / preset["sample_dir"]

    if not MC_EXEC.exists():
        print(f"エラー: MC_simple が見つかりません。先にビルドしてください:")
        print("  cd src && cmake -S . -B build && cmake --build build -j")
        sys.exit(1)

    if not sample_dir.exists():
        print(f"エラー: サンプルディレクトリが見つかりません: {sample_dir}")
        sys.exit(1)

    print(f"プリセット: {args.preset} (sample_dir={preset['sample_dir']})")
    print(f"  温度: T={preset['ini_t']}..{preset['ini_t']+(preset['num_temp']-1)*preset['delta_t']:.1f}, "
          f"Burn_in={preset['burn_in']}, Total_Step={preset['total_step']}")

    # benchmark 用ディレクトリを作成（プリセットごとに分離）
    bench_dir = REPO_ROOT / "samples" / "exchange_mc_benchmark" / args.preset
    bench_dir.mkdir(parents=True, exist_ok=True)

    # lattice, interaction をコピー（またはシンボリックリンク）
    import shutil
    for name in ["lattice.def", "interaction.def"]:
        dst = bench_dir / name
        src = sample_dir / name
        if dst.exists():
            dst.unlink()
        try:
            dst.symlink_to(src.resolve())
        except OSError:
            shutil.copy(src, dst)

    lattice_path = bench_dir / "lattice.def"
    interaction_path = bench_dir / "interaction.def"
    if not lattice_path.exists() or not interaction_path.exists():
        print("エラー: lattice.def または interaction.def が見つかりません")
        sys.exit(1)

    param_no_ex = bench_dir / "param_no_exchange.def"
    param_with_ex = bench_dir / "param_with_exchange.def"

    param_no_ex.write_text(make_param_content(
        enable_exchange=0, **{k: preset[k] for k in ["num_temp", "ini_t", "delta_t",
                                                     "burn_in", "total_step", "sample"]}
    ))
    param_with_ex.write_text(make_param_content(
        enable_exchange=1, **{k: preset[k] for k in ["num_temp", "ini_t", "delta_t",
                                                     "burn_in", "total_step", "sample"]}
    ))

    print(f"Exchange MC なし: {args.n_runs} 回実行...")
    data_no_ex = collect_runs(
        args.n_runs, 0, param_no_ex, lattice_path, interaction_path,
        MC_EXEC, bench_dir
    )

    print(f"Exchange MC あり: {args.n_runs} 回実行...")
    data_with_ex = collect_runs(
        args.n_runs, 1, param_with_ex, lattice_path, interaction_path,
        MC_EXEC, bench_dir
    )

    # 数値結果を表示
    print("\n=== 結果 (mean ± stderr) ===")
    print("\n[Exchange MC なし]")
    for i in range(len(data_no_ex["T"])):
        ov = data_no_ex["overlap_mean"][i]
        print(f"  T={data_no_ex['T'][i]:.2f}: E={data_no_ex['E_mean'][i]:.4f}±{data_no_ex['E_err'][i]:.4f}, "
              f"C={data_no_ex['C_mean'][i]:.4f}±{data_no_ex['C_err'][i]:.4f}, "
              f"M2={data_no_ex['M2_mean'][i]:.4f}±{data_no_ex['M2_err'][i]:.4f}, overlap={ov:.4f}")

    print("\n[Exchange MC あり]")
    for i in range(len(data_with_ex["T"])):
        ov = data_with_ex["overlap_mean"][i]
        print(f"  T={data_with_ex['T'][i]:.2f}: E={data_with_ex['E_mean'][i]:.4f}±{data_with_ex['E_err'][i]:.4f}, "
              f"C={data_with_ex['C_mean'][i]:.4f}±{data_with_ex['C_err'][i]:.4f}, "
              f"M2={data_with_ex['M2_mean'][i]:.4f}±{data_with_ex['M2_err'][i]:.4f}, overlap={ov:.4f}")

    # 平均エラーバー比（全温度・全観測量の幾何平均）
    ratios = []
    for i in range(len(data_no_ex["T"])):
        for err_no, err_with in [
            (data_no_ex["E_err"][i], data_with_ex["E_err"][i]),
            (data_no_ex["C_err"][i], data_with_ex["C_err"][i]),
            (data_no_ex["M2_err"][i], data_with_ex["M2_err"][i]),
        ]:
            if err_no > 1e-12:
                ratios.append(err_with / err_no)
    if ratios:
        geo_mean = np.exp(np.mean(np.log(ratios)))
        print(f"\n=== 平均エラーバー比（幾何平均）: {geo_mean:.3f} "
              f"(<1 なら Exchange MC で全体的にエラー減少) ===")

    # overlap 比較（低い = よく mix、EXMC の効果）
    print("\n=== overlap 比較（Exchange MC ありの方が低いと mix が良い）===")
    for i in range(len(data_no_ex["T"])):
        T = data_no_ex["T"][i]
        ov_no = data_no_ex["overlap_mean"][i]
        ov_ex = data_with_ex["overlap_mean"][i]
        diff = ov_ex - ov_no
        print(f"  T={T:.2f}: no_EX={ov_no:.4f}, with_EX={ov_ex:.4f}, 差={diff:+.4f} "
              f"(負なら EXMC で overlap 減少 = 効果あり)")

    # エラーバー比較の要約
    print("\n=== エラーバー比較（Exchange MC ありの方が小さいと効果あり）===")
    for i in range(len(data_no_ex["T"])):
        T = data_no_ex["T"][i]
        e_ratio = data_with_ex["E_err"][i] / data_no_ex["E_err"][i] if data_no_ex["E_err"][i] > 0 else 0
        c_ratio = data_with_ex["C_err"][i] / data_no_ex["C_err"][i] if data_no_ex["C_err"][i] > 0 else 0
        m2_ratio = data_with_ex["M2_err"][i] / data_no_ex["M2_err"][i] if data_no_ex["M2_err"][i] > 0 else 0
        print(f"  T={T:.2f}: E_err比={e_ratio:.3f}, C_err比={c_ratio:.3f}, M2_err比={m2_ratio:.3f} "
              f"(<1 なら Exchange MC でエラーが減少)")

    # プロット
    if not args.no_plot and HAS_MATPLOTLIB:
        out_path = bench_dir / args.output
        plot_results(data_no_ex, data_with_ex, out_path)
    elif args.no_plot:
        print("\n--no-plot のためプロットをスキップしました。")


if __name__ == "__main__":
    main()
