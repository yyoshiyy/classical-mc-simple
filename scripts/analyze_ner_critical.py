#!/usr/bin/env python3
"""
NER_result.dat から臨界指数を推定するスクリプト

臨界温度付近で m(t) ~ A*t^(-α) のべき乗減衰が期待される。
α = β/(νz) より、2D Ising の既知値 β=1/8, ν=1 から動的臨界指数 z を推定。

使い方:
  python scripts/analyze_ner_critical.py samples/square_L16_Ising/NER_result.dat
  python scripts/analyze_ner_critical.py NER_result.dat --t-min 2 --t-max 50
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# 2D Ising の静的臨界指数（厳密解）
BETA_2D_ISING = 1.0 / 8.0
NU_2D_ISING = 1.0


def power_law(t, A, alpha):
    """y = A * t^(-alpha). t=0 は inf になるので 1e-10 でクリップ"""
    t_safe = np.maximum(t, 1e-10)
    return A * np.power(t_safe, -alpha)


def fit_power_law(t, y, t_min, t_max):
    """
    y(t) = A * t^(-α) を log-log 空間で線形回帰。
    log(y) = log(A) - α*log(t)
    返り値: (A, alpha, alpha_err)
    """
    mask = (t >= t_min) & (t <= t_max) & (y > 0)
    t_fit = t[mask]
    y_fit = y[mask]

    if len(t_fit) < 3:
        return None, None, None

    log_t = np.log(t_fit)
    log_y = np.log(y_fit)

    # 線形回帰: log_y = b - alpha * log_t
    coeffs, residuals, rank, s, cond = np.polyfit(log_t, log_y, 1, full=True)
    slope = coeffs[0]  # = -alpha
    intercept = coeffs[1]  # = log(A)
    alpha = -slope
    A = np.exp(intercept)

    # 標準誤差（簡易）
    n = len(t_fit)
    if n > 2 and residuals.size > 0:
        mse = residuals[0] / (n - 2)
        var_slope = mse / np.sum((log_t - np.mean(log_t)) ** 2)
        alpha_err = np.sqrt(var_slope)
    else:
        alpha_err = np.nan

    return A, alpha, alpha_err


def main():
    parser = argparse.ArgumentParser(
        description="NER 時系列から臨界指数（べき乗減衰）を推定"
    )
    parser.add_argument(
        "file",
        nargs="?",
        default="NER_result.dat",
        help="NER_result.dat のパス",
    )
    parser.add_argument(
        "--t-min",
        type=int,
        default=2,
        help="フィット開始ステップ (t=0 は除外, デフォルト: 2)",
    )
    parser.add_argument(
        "--t-max",
        type=int,
        default=None,
        help="フィット終了ステップ (未指定時は全データ)",
    )
    parser.add_argument(
        "--short-time",
        action="store_true",
        help="短時間領域のみフィット (t_max=min(50, 全データの1/3))",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="フィット結果をプロットして保存",
    )
    parser.add_argument(
        "-o", "--output",
        default="NER_critical_fit.png",
        help="プロット出力ファイル (--plot 時)",
    )
    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        print(f"Error: ファイルが見つかりません: {path}", file=sys.stderr)
        sys.exit(1)

    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)

    t = data[:, 0]
    T = data[:, 1]
    m = data[:, 2]
    q = data[:, 5]

    if args.t_max is not None:
        t_max = args.t_max
    elif args.short_time:
        t_max = min(50, int(len(t) / 3))
    else:
        t_max = int(t[-1])

    # m(t) のべき乗フィット
    A_m, alpha_m, alpha_m_err = fit_power_law(t, m, args.t_min, t_max)

    if A_m is None:
        print("Error: フィットに十分なデータがありません。--t-min / --t-max を調整してください。")
        sys.exit(1)

    print(f"=== NER 臨界指数推定 (T = {T[0]:.3f}) ===")
    print(f"フィット範囲: t = {args.t_min} .. {t_max}")
    print()
    print("【m(t) のべき乗減衰】 m(t) ~ A * t^(-α)")
    print(f"  A     = {A_m:.6f}")
    print(f"  α     = {alpha_m:.6f} ± {alpha_m_err:.6f}")
    print()
    print("【2D Ising からの z 推定】 α = β/(νz) より z = β/(ν*α)")
    print(f"  β = 1/8, ν = 1 (2D Ising 厳密値)")
    z_est = BETA_2D_ISING / (NU_2D_ISING * alpha_m)
    z_err = (
        z_est * (alpha_m_err / alpha_m)
        if not np.isnan(alpha_m_err) and alpha_m > 0
        else np.nan
    )
    print(f"  z_est = {z_est:.4f} ± {z_err:.4f}")
    print(f"  (文献値: z ≈ 2.17)")
    print()
    print("※ べき乗減衰は臨界温度 T=Tc 付近で成立。T が Tc から離れると指数減衰になる。")
    if alpha_m <= 0:
        print("※ α≤0: T が Tc から離れている可能性。T≈2.27 付近で再実行を推奨。")

    # q(t) のフィット（q>0 の範囲のみ）
    mask_q = (t >= args.t_min) & (t <= t_max) & (q > 0.01)
    if np.sum(mask_q) >= 3:
        A_q, alpha_q, alpha_q_err = fit_power_law(t, q, args.t_min, t_max)
        if A_q is not None:
            print()
            print("【q(t) のべき乗減衰】 q(t) ~ A * t^(-α) (q>0 の範囲)")
            print(f"  α_q   = {alpha_q:.6f} ± {alpha_q_err:.6f}")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("--plot には matplotlib が必要です", file=sys.stderr)
            sys.exit(1)

        m_fit = power_law(t, A_m, alpha_m)
        m_fit[t < args.t_min] = np.nan
        m_fit[t > t_max] = np.nan

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        axes[0].plot(t, m, "b.-", markersize=4, label="m(t)")
        axes[0].plot(t, m_fit, "r--", lw=2, label=rf"fit: $A t^{{-{alpha_m:.3f}}}$")
        axes[0].axvspan(args.t_min, t_max, alpha=0.2, color="gray")
        axes[0].set_xlabel("MC step t")
        axes[0].set_ylabel(r"$m(t)$")
        axes[0].set_title("Order parameter decay + power-law fit")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        axes[0].set_xlim(0, t[-1])

        axes[1].plot(t, np.log(m + 1e-10), "b.-", markersize=4, label="log m(t)")
        log_m_fit = np.log(m_fit + 1e-10)
        axes[1].plot(t, log_m_fit, "r--", lw=2, label="fit (log)")
        axes[1].axvspan(args.t_min, t_max, alpha=0.2, color="gray")
        axes[1].set_xlabel("MC step t")
        axes[1].set_ylabel(r"$\ln m(t)$")
        axes[1].set_title("Log-log (linear region = power law)")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        axes[1].set_xlim(0, t[-1])

        fig.suptitle(f"NER critical exponent fit (T={T[0]:.2f}, α={alpha_m:.4f}, z_est={z_est:.3f})")
        plt.tight_layout()

        out_path = Path(args.output)
        if not out_path.is_absolute():
            out_path = path.parent / out_path.name
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
