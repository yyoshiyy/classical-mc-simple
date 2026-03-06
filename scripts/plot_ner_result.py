#!/usr/bin/env python3
"""
NER_result.dat をプロットするスクリプト

使い方:
  cd samples/square_L16_Ising
  python ../../scripts/plot_ner_result.py

  またはファイルパスを指定:
  python ../../scripts/plot_ner_result.py NER_result.dat
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib が必要です: pip install matplotlib", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="NER_result.dat をプロット")
    parser.add_argument(
        "file",
        nargs="?",
        default="NER_result.dat",
        help="入力ファイル (デフォルト: NER_result.dat)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="出力画像ファイル (指定しない場合は表示)",
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
    m2 = data[:, 3]
    e_per_site = data[:, 4]
    q = data[:, 5]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)

    axes[0, 0].plot(t, m, "b.-", markersize=4)
    axes[0, 0].set_ylabel(r"$m(t)$")
    axes[0, 0].set_title("Order parameter magnitude")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(t, m2, "g.-", markersize=4)
    axes[0, 1].set_ylabel(r"$m^2(t)$")
    axes[0, 1].set_title("Second moment")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(t, e_per_site, "r.-", markersize=4)
    axes[1, 0].set_xlabel("MC step t")
    axes[1, 0].set_ylabel(r"$E/N$")
    axes[1, 0].set_title("Energy per site")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(t, q, "m.-", markersize=4)
    axes[1, 1].set_xlabel("MC step t")
    axes[1, 1].set_ylabel(r"$q(t)$")
    axes[1, 1].set_title("Overlap with initial state")
    axes[1, 1].grid(True, alpha=0.3)

    T_val = T[0] if len(T) > 0 else 0
    fig.suptitle(f"NER relaxation (T = {T_val})", fontsize=12)
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150, bbox_inches="tight")
        print(f"Saved: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
