#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
postprocess.py - Hong et al. (2017) ZDPlasKin 结果后处理与可视化

功能：
  1. 读取 ZDPlasKin 输出文件（时间演化的物种密度）
  2. 绘制物种密度随时间演化（对比 Figure 3）
  3. 计算电子能量分支（对比 Figure 4）
  4. 绘制 NH3 产率 vs 条件（对比 Figure 5）
  5. 验证关键结论

使用方法：
  python postprocess.py --input data/results/baseline.out
  python postprocess.py --compare --dir data/results
  python postprocess.py --energy_branching --input data/results/baseline.out

依赖：numpy, matplotlib, pandas
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse
import json
import sys

# 设置中文字体（如果需要）
plt.rcParams['font.size'] = 12

# ============================================================================
# 验证目标（从 params.json 读取）
# ============================================================================
VALIDATION_TARGETS = {
    "NH3_steady_density_cm3": (1e14, 1e15),  # 范围
    "electron_temperature_ev": (1.0, 5.0),
    "N2v1_vs_N_ratio": 100,  # N2(v=1) / N(4S) > 100
}

SPECIES_COLORS = {
    "N2": "blue",
    "H2": "green",
    "NH3": "red",
    "N": "orange",
    "H": "purple",
    "e": "black",
    "N2_v1": "cyan",
    "NH": "magenta",
    "NH2": "brown",
}


def parse_zdplaskin_output(filepath):
    """
    解析 ZDPlasKin 输出文件
    
    ZDPlasKin 输出格式（QtPlaskin 兼容格式）：
    第一行：标题/列名
    后续行：time, species1, species2, ..., Te, Tgas, EN, ...
    
    注意：实际格式取决于主程序的输出格式。需要根据你的输出格式调整。
    
    这里假设为 CSV-like 格式：
    # time(s)  N2(cm-3)  H2(cm-3)  NH3(cm-3)  N(cm-3)  H(cm-3)  e(cm-3)  Te(K)  Tgas(K)
    0.000e+00  2.50e+18  7.50e+18  0.00e+00  ...
    """
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"ERROR: File not found: {filepath}")
        return None
    
    # 尝试读取 CSV 格式
    try:
        # 跳过注释行（以 # 开头）
        with open(filepath, "r") as f:
            lines = f.readlines()
        
        # 找到列名行（第一个不以 # 开头的行）
        header_line = None
        data_lines = []
        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if header_line is None:
                # 检查是否是数据行（包含数字）还是标题行
                parts = line.split()
                try:
                    float(parts[0])  # 如果第一个是数字，说明没有标题行
                    header_line = "time " + " ".join([f"col_{i}" for i in range(len(parts))])
                    data_lines.append(line)
                except ValueError:
                    header_line = line
            else:
                data_lines.append(line)
        
        if header_line is None:
            print("ERROR: No data found in file")
            return None
        
        # 解析数据
        columns = header_line.split()
        data = []
        for line in data_lines:
            parts = line.split()
            if len(parts) >= len(columns):
                data.append([float(p) for p in parts[:len(columns)]])
        
        df = pd.DataFrame(data, columns=columns)
        return df
        
    except Exception as e:
        print(f"ERROR parsing file: {e}")
        return None


def plot_species_evolution(df, output_path=None, species_list=None):
    """
    绘制物种密度随时间演化（对比论文 Figure 3）
    
    典型物种：N2, H2, NH3, N, H, e, N2_v1
    """
    if df is None:
        return
    
    time_col = df.columns[0]  # 假设第一列是时间
    time = df[time_col].values
    
    if species_list is None:
        # 默认绘制主要物种
        species_list = ["N2", "H2", "NH3", "N", "H", "e", "N2_v1"]
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # 上图：线性尺度（短时间）
    ax1 = axes[0]
    for sp in species_list:
        if sp in df.columns:
            color = SPECIES_COLORS.get(sp, "gray")
            ax1.plot(time, df[sp].values, label=sp, color=color, linewidth=1.5)
    ax1.set_ylabel("Density (cm$^{-3}$)")
    ax1.set_title("Species Density Evolution - Linear Scale")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)
    
    # 下图：对数尺度（全时间范围）
    ax2 = axes[1]
    for sp in species_list:
        if sp in df.columns:
            color = SPECIES_COLORS.get(sp, "gray")
            y = df[sp].values
            y = np.maximum(y, 1e-20)  # 避免 log(0)
            ax2.semilogy(time, y, label=sp, color=color, linewidth=1.5)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("Density (cm$^{-3}$)")
    ax2.set_title("Species Density Evolution - Log Scale")
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"Saved: {output_path}")
    else:
        plt.show()


def calculate_energy_branching(df):
    """
    计算电子能量分支（对比论文 Figure 4）
    
    需要从 ZDPlasKin 输出中提取各反应的能量消耗。
    ZDPlasKin 的诊断输出可以给出每个反应的速率。
    然后乘以每个反应的阈值能量，得到能量分支。
    
    或者从 BOLSIG+ 输出中直接读取能量分支。
    """
    print("\n" + "=" * 60)
    print("ENERGY BRANCHING ANALYSIS")
    print("=" * 60)
    print("Note: This requires reaction rate output from ZDPlasKin")
    print("or BOLSIG+ energy partitioning data.")
    print("\nExpected energy branching (from paper):")
    print("  Vibrational excitation  > 50% (dominant)")
    print("  Electronic excitation     ~ 20-30%")
    print("  Dissociation              ~ 10-20%")
    print("  Ionization                < 5%")
    print("  Elastic scattering        < 5%")
    print("\nTODO: Implement extraction from ZDPlasKin diagnostic output")


def validate_results(df):
    """
    验证结果是否在论文报告的范围内
    """
    print("\n" + "=" * 60)
    print("VALIDATION CHECK")
    print("=" * 60)
    
    if df is None:
        print("No data to validate")
        return
    
    # 检查 NH3 稳态密度
    if "NH3" in df.columns:
        nh3_final = df["NH3"].values[-1]
        target_min, target_max = VALIDATION_TARGETS["NH3_steady_density_cm3"]
        status = "PASS" if target_min <= nh3_final <= target_max else "FAIL"
        print(f"  NH3 steady-state density: {nh3_final:.2e} cm-3  [{status}]")
        print(f"    Target: {target_min:.2e} - {target_max:.2e} cm-3")
    else:
        print("  NH3 column not found in output")
    
    # 检查电子温度
    if "Te" in df.columns or "T_e" in df.columns:
        te_col = "Te" if "Te" in df.columns else "T_e"
        te_final = df[te_col].values[-1]
        # 如果单位是 K，转换为 eV
        if te_final > 100:  # 假设是 K
            te_ev = te_final / 11604.5
        else:
            te_ev = te_final
        target_min, target_max = VALIDATION_TARGETS["electron_temperature_ev"]
        status = "PASS" if target_min <= te_ev <= target_max else "FAIL"
        print(f"  Electron temperature: {te_ev:.2f} eV  [{status}]")
        print(f"    Target: {target_min:.1f} - {target_max:.1f} eV")
    else:
        print("  Te column not found in output")
    
    # 检查 N2(v=1) vs N 密度比
    if "N2_v1" in df.columns and "N" in df.columns:
        n2v1_final = df["N2_v1"].values[-1]
        n_final = df["N"].values[-1]
        if n_final > 0:
            ratio = n2v1_final / n_final
            target = VALIDATION_TARGETS["N2v1_vs_N_ratio"]
            status = "PASS" if ratio > target else "FAIL"
            print(f"  N2(v=1)/N ratio: {ratio:.1f}  [{status}]")
            print(f"    Target: > {target}")
        else:
            print("  N density is zero, cannot calculate ratio")
    else:
        print("  N2_v1 or N column not found")
    
    print("=" * 60)


def compare_parameter_scan(results_dir):
    """
    对比参数扫描结果
    """
    results_dir = Path(results_dir)
    if not results_dir.exists():
        print(f"ERROR: Directory not found: {results_dir}")
        return
    
    print(f"\n{'='*60}")
    print("PARAMETER SCAN COMPARISON")
    print(f"{'='*60}")
    
    # 查找所有输出文件
    out_files = sorted(results_dir.glob("*.out"))
    
    if not out_files:
        print("No output files found")
        return
    
    print(f"Found {len(out_files)} result files")
    
    # 分类：按参数类型
    pressure_files = [f for f in out_files if "pressure" in f.name]
    ratio_files = [f for f in out_files if "N2H2" in f.name]
    en_files = [f for f in out_files if "EN_" in f.name]
    temp_files = [f for f in out_files if "T_" in f.name]
    
    # 绘制压力扫描结果
    if pressure_files:
        fig, ax = plt.subplots(figsize=(8, 6))
        pressures = []
        nh3_densities = []
        for f in pressure_files:
            df = parse_zdplaskin_output(f)
            if df is not None and "NH3" in df.columns:
                p = float(f.stem.split("_")[1].replace("atm", ""))
                pressures.append(p)
                nh3_densities.append(df["NH3"].values[-1])
        if pressures:
            ax.semilogy(pressures, nh3_densities, "bo-", markersize=8, linewidth=2)
            ax.set_xlabel("Pressure (atm)")
            ax.set_ylabel("NH$_3$ Steady-State Density (cm$^{-3}$)")
            ax.set_title("NH$_3$ Yield vs Pressure")
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(results_dir / "pressure_scan.png", dpi=300)
            print(f"Saved: {results_dir / 'pressure_scan.png'}")
            plt.close()
    
    # 绘制 E/N 扫描结果
    if en_files:
        fig, ax = plt.subplots(figsize=(8, 6))
        ens = []
        nh3_densities = []
        for f in en_files:
            df = parse_zdplaskin_output(f)
            if df is not None and "NH3" in df.columns:
                en = float(f.stem.split("_")[1].replace("Td", ""))
                ens.append(en)
                nh3_densities.append(df["NH3"].values[-1])
        if ens:
            ax.plot(ens, nh3_densities, "ro-", markersize=8, linewidth=2)
            ax.set_xlabel("E/N (Td)")
            ax.set_ylabel("NH$_3$ Steady-State Density (cm$^{-3}$)")
            ax.set_title("NH$_3$ Yield vs Reduced Electric Field")
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(results_dir / "EN_scan.png", dpi=300)
            print(f"Saved: {results_dir / 'EN_scan.png'}")
            plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Hong et al. (2017) ZDPlasKin Post-Processing"
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Input ZDPlasKin output file to process"
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Compare parameter scan results in directory"
    )
    parser.add_argument(
        "--dir",
        type=str,
        default="data/results",
        help="Directory containing results"
    )
    parser.add_argument(
        "--energy_branching",
        action="store_true",
        help="Calculate energy branching"
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate results against paper targets"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="Species to plot (default: major species)"
    )
    
    args = parser.parse_args()
    
    if args.compare:
        compare_parameter_scan(args.dir)
    elif args.input:
        df = parse_zdplaskin_output(args.input)
        if df is not None:
            plot_species_evolution(df, species_list=args.species)
            if args.validate:
                validate_results(df)
            if args.energy_branching:
                calculate_energy_branching(df)
    else:
        print("Usage:")
        print("  python postprocess.py --input data/results/baseline.out --validate")
        print("  python postprocess.py --compare --dir data/results")
        print("  python postprocess.py --input data/results/baseline.out --energy_branching")


if __name__ == "__main__":
    main()
