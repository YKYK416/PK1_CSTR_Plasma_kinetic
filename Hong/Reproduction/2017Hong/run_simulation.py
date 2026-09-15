#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_simulation.py - Hong et al. (2017) ZDPlasKin 复现运行脚本

功能：
  1. 读取 params.json 配置
  2. 批量运行参数扫描（压力、N2:H2比、E/N、温度）
  3. 调用编译后的 ZDPlasKin 可执行文件
  4. 收集结果并保存

使用方法：
  python run_simulation.py --mode baseline
  python run_simulation.py --mode scan --parameter pressure
  python run_simulation.py --mode scan --parameter N2_H2_ratio
  python run_simulation.py --mode scan --parameter reduced_field

依赖：
  - ZDPlasKin 已编译（test_NH3_baseline 可执行文件）
  - params.json 配置文件
"""

import json
import subprocess
import os
import sys
import argparse
from pathlib import Path
import numpy as np

# ============================================================================
# 配置路径
# ============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
PARAMS_FILE = SCRIPT_DIR / "params.json"
RESULTS_DIR = SCRIPT_DIR / "data" / "results"
CROSS_SECTIONS_DIR = SCRIPT_DIR / "data" / "cross_sections"

# ZDPlasKin 可执行文件（编译后）
EXECUTABLE = SCRIPT_DIR / "test_NH3_baseline"
# 如果编译后的文件名不同，请修改这里

# ============================================================================
# 默认运行参数（从 params.json 读取，但可以被命令行覆盖）
# ============================================================================
DEFAULT_BASELINE = {
    "pressure_pa": 101325.0,
    "gas_temperature_k": 300.0,
    "N2_fraction": 0.25,
    "H2_fraction": 0.75,
    "reduced_field_townsend": 150.0,
    "initial_electron_density_cm3": 1.0e8,
    "simulation_time_s": 0.01,
    "output_interval_s": 1.0e-6,
}


def load_params():
    """加载 params.json 配置文件"""
    with open(PARAMS_FILE, "r", encoding="utf-8") as f:
        return json.load(f)


def write_input_file(params, output_dir, run_name):
    """
    生成 ZDPlasKin 输入文件（或参数传递文件）
    
    注意：ZDPlasKin 的参数传递方式取决于主程序的实现。
    如果主程序从命令行参数读取，直接传递参数。
    如果主程序从输入文件读取，生成相应的文件。
    """
    input_file = output_dir / f"{run_name}.inp"
    with open(input_file, "w") as f:
        f.write(f"# ZDPlasKin input for {run_name}\n")
        f.write(f"pressure_pa = {params['pressure_pa']}\n")
        f.write(f"gas_temperature_k = {params['gas_temperature_k']}\n")
        f.write(f"N2_fraction = {params['N2_fraction']}\n")
        f.write(f"H2_fraction = {params['H2_fraction']}\n")
        f.write(f"reduced_field_townsend = {params['reduced_field_townsend']}\n")
        f.write(f"initial_electron_density_cm3 = {params['initial_electron_density_cm3']}\n")
        f.write(f"simulation_time_s = {params['simulation_time_s']}\n")
        f.write(f"output_interval_s = {params['output_interval_s']}\n")
    return input_file


def run_zdplaskin(input_file, output_dir, run_name):
    """
    调用 ZDPlasKin 可执行文件运行模拟
    
    返回：subprocess.CompletedProcess 对象
    """
    if not EXECUTABLE.exists():
        print(f"ERROR: ZDPlasKin executable not found: {EXECUTABLE}")
        print("Please compile the Fortran code first:")
        print("  ./preprocessor kinet.inp")
        print("  gfortran -o test_NH3_baseline test_NH3_baseline.F90 zdplaskin_m.F90 ...")
        sys.exit(1)
    
    output_file = output_dir / f"{run_name}.out"
    
    # 构建命令行参数
    # 根据你的主程序接口调整参数传递方式
    cmd = [
        str(EXECUTABLE),
        str(input_file),
        str(output_file)
    ]
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1小时超时
            cwd=str(SCRIPT_DIR)
        )
        
        if result.returncode != 0:
            print(f"ERROR: ZDPlasKin failed with return code {result.returncode}")
            print("STDERR:", result.stderr)
            return None
        
        print(f"Completed: {run_name}")
        return result
        
    except subprocess.TimeoutExpired:
        print(f"ERROR: ZDPlasKin timeout for {run_name}")
        return None
    except Exception as e:
        print(f"ERROR: {e}")
        return None


def run_baseline():
    """运行基准条件模拟"""
    params = load_params()
    baseline = params["simulation"]["baseline_conditions"]
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 生成输入文件
    input_file = write_input_file(baseline, RESULTS_DIR, "baseline")
    
    # 运行模拟
    result = run_zdplaskin(input_file, RESULTS_DIR, "baseline")
    
    if result:
        print("\n" + "=" * 60)
        print("BASELINE SIMULATION COMPLETE")
        print("=" * 60)
        print(f"Results saved to: {RESULTS_DIR / 'baseline.out'}")
        print("Next: Run postprocess.py to visualize results")
    else:
        print("\nBaseline simulation failed!")


def run_pressure_scan():
    """压力扫描：验证常压 vs 低压的效率差异"""
    params = load_params()
    baseline = params["simulation"]["baseline_conditions"]
    pressure_values = params["simulation"]["parameter_scans"]["pressure"]["values_atm"]
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print("PRESSURE SCAN")
    print(f"{'='*60}")
    print(f"Values: {pressure_values} atm")
    print(f"Expected: Atmospheric pressure > low pressure by 10-100x")
    
    for p_atm in pressure_values:
        p_pa = p_atm * 101325.0
        run_params = baseline.copy()
        run_params["pressure_pa"] = p_pa
        
        run_name = f"pressure_{p_atm:.2f}atm"
        input_file = write_input_file(run_params, RESULTS_DIR, run_name)
        result = run_zdplaskin(input_file, RESULTS_DIR, run_name)
        
        if result is None:
            print(f"Warning: Pressure scan failed at {p_atm} atm")
    
    print(f"\nPressure scan complete. Results in: {RESULTS_DIR}")


def run_composition_scan():
    """气体组成扫描：N2:H2 比例"""
    params = load_params()
    baseline = params["simulation"]["baseline_conditions"]
    ratios = params["simulation"]["parameter_scans"]["N2_H2_ratio"]["values"]
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print("N2:H2 RATIO SCAN")
    print(f"{'='*60}")
    print(f"Ratios: {ratios}")
    
    for n2, h2 in ratios:
        total = n2 + h2
        run_params = baseline.copy()
        run_params["N2_fraction"] = n2 / total
        run_params["H2_fraction"] = h2 / total
        
        run_name = f"N2H2_{n2}_{h2}"
        input_file = write_input_file(run_params, RESULTS_DIR, run_name)
        result = run_zdplaskin(input_file, RESULTS_DIR, run_name)
        
        if result is None:
            print(f"Warning: Composition scan failed at N2:H2 = {n2}:{h2}")
    
    print(f"\nComposition scan complete. Results in: {RESULTS_DIR}")


def run_EN_scan():
    """约化电场 E/N 扫描"""
    params = load_params()
    baseline = params["simulation"]["baseline_conditions"]
    en_values = params["simulation"]["parameter_scans"]["reduced_field"]["values_townsend"]
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print("REDUCED FIELD E/N SCAN")
    print(f"{'='*60}")
    print(f"Values: {en_values} Td")
    
    for en in en_values:
        run_params = baseline.copy()
        run_params["reduced_field_townsend"] = en
        
        run_name = f"EN_{en}Td"
        input_file = write_input_file(run_params, RESULTS_DIR, run_name)
        result = run_zdplaskin(input_file, RESULTS_DIR, run_name)
        
        if result is None:
            print(f"Warning: E/N scan failed at {en} Td")
    
    print(f"\nE/N scan complete. Results in: {RESULTS_DIR}")


def run_temperature_scan():
    """气体温度扫描"""
    params = load_params()
    baseline = params["simulation"]["baseline_conditions"]
    T_values = params["simulation"]["parameter_scans"]["gas_temperature"]["values_k"]
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print("GAS TEMPERATURE SCAN")
    print(f"{'='*60}")
    print(f"Values: {T_values} K")
    
    for T in T_values:
        run_params = baseline.copy()
        run_params["gas_temperature_k"] = T
        
        run_name = f"T_{T}K"
        input_file = write_input_file(run_params, RESULTS_DIR, run_name)
        result = run_zdplaskin(input_file, RESULTS_DIR, run_name)
        
        if result is None:
            print(f"Warning: Temperature scan failed at {T} K")
    
    print(f"\nTemperature scan complete. Results in: {RESULTS_DIR}")


def main():
    parser = argparse.ArgumentParser(
        description="Hong et al. (2017) ZDPlasKin Reproduction Runner"
    )
    parser.add_argument(
        "--mode",
        choices=["baseline", "scan"],
        default="baseline",
        help="Run mode: baseline or parameter scan"
    )
    parser.add_argument(
        "--parameter",
        choices=["pressure", "N2_H2_ratio", "reduced_field", "temperature", "all"],
        default="all",
        help="Parameter to scan (only used with --mode scan)"
    )
    
    args = parser.parse_args()
    
    if args.mode == "baseline":
        run_baseline()
    else:
        if args.parameter in ["pressure", "all"]:
            run_pressure_scan()
        if args.parameter in ["N2_H2_ratio", "all"]:
            run_composition_scan()
        if args.parameter in ["reduced_field", "all"]:
            run_EN_scan()
        if args.parameter in ["temperature", "all"]:
            run_temperature_scan()


if __name__ == "__main__":
    main()
