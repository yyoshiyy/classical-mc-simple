#!/usr/bin/env python3
"""
NER モードで計算を実行し、完了後に自動でプロットするスクリプト

使い方:
  cd classical-mc-simple
  python scripts/run_ner_and_plot.py

  サンプルディレクトリを指定:
  python scripts/run_ner_and_plot.py --sample square_L24_Ising

  プロットを表示のみ（ファイル保存しない）:
  python scripts/run_ner_and_plot.py --no-save
"""

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BUILD_DIR = REPO_ROOT / "src" / "build"
MC_EXEC = BUILD_DIR / "MC_simple"
PLOT_SCRIPT = REPO_ROOT / "scripts" / "plot_ner_result.py"


def main():
    parser = argparse.ArgumentParser(description="NER 計算 + プロットを一括実行")
    parser.add_argument(
        "--sample",
        default="square_L16_Ising",
        help="サンプルディレクトリ名 (デフォルト: square_L16_Ising)",
    )
    parser.add_argument(
        "--param",
        default="param_ner.def",
        help="param ファイル (デフォルト: param_ner.def)",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="プロットをファイルに保存せず表示のみ",
    )
    parser.add_argument(
        "--no-analyze",
        action="store_true",
        help="臨界指数推定をスキップ（デフォルトは自動実行）",
    )
    args = parser.parse_args()

    sample_dir = REPO_ROOT / "samples" / args.sample
    if not sample_dir.exists():
        print(f"Error: サンプルディレクトリが見つかりません: {sample_dir}", file=sys.stderr)
        sys.exit(1)

    param_path = sample_dir / args.param
    if not param_path.exists():
        print(f"Error: param ファイルが見つかりません: {param_path}", file=sys.stderr)
        sys.exit(1)

    if not MC_EXEC.exists():
        print(f"Error: MC_simple がビルドされていません: {MC_EXEC}", file=sys.stderr)
        print("  cd src && cmake -S . -B build && cmake --build build -j", file=sys.stderr)
        sys.exit(1)

    # 1. NER 計算を実行
    print("Running NER simulation...", flush=True)
    rc = subprocess.run(
        [str(MC_EXEC), str(param_path), "lattice.def", "interaction.def"],
        cwd=sample_dir,
    )
    if rc.returncode != 0:
        print(f"Error: MC_simple が終了コード {rc.returncode} で終了しました", file=sys.stderr)
        sys.exit(rc.returncode)

    ner_result = sample_dir / "NER_result.dat"
    if not ner_result.exists():
        print("Error: NER_result.dat が生成されませんでした (enable_ner=1 の param を使用していますか?)", file=sys.stderr)
        sys.exit(1)

    # 2. プロットを実行
    print("Plotting NER_result.dat...", flush=True)
    plot_args = [sys.executable, str(PLOT_SCRIPT), str(ner_result)]
    if not args.no_save:
        plot_args.extend(["-o", str(sample_dir / "NER_plot.png")])

    rc = subprocess.run(plot_args, cwd=sample_dir)
    if rc.returncode != 0:
        print("Error: プロットに失敗しました", file=sys.stderr)
        sys.exit(rc.returncode)

    if not args.no_analyze:
        print("Analyzing critical exponent...", flush=True)
        analyze_script = REPO_ROOT / "scripts" / "analyze_ner_critical.py"
        rc = subprocess.run(
            [sys.executable, str(analyze_script), str(ner_result), "--short-time", "--plot", "-o", "NER_critical_fit.png"],
            cwd=sample_dir,
        )
        if rc.returncode != 0:
            print("Warning: 臨界指数推定に失敗しました", file=sys.stderr)

    print("Done.")


if __name__ == "__main__":
    main()
