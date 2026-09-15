#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Windows-native representative reproduction for Shao & Mesbah (JACS Au 2024).

The upstream workflow is a Slurm array.  This script reproduces the four
published field/temperature benchmark conditions with the author-supplied,
pre-generated DFT ZDPlasKin module.  It never rewrites the upstream mechanism.
"""
from __future__ import annotations

import argparse
import csv
import json
import pickle
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SHAO_ROOT = HERE.parent
CASE_TEMPLATE = SHAO_ROOT / "dft_microkinetic_windows_case"
REFERENCE = HERE / "author_reference"
CASES_DIR = HERE / "cases_4000s_finalsample"
RESULTS_DIR = HERE / "results"
TIMEPOINT_S = 4000.0
ELECTRON_DENSITY = 8.27e7
N2_FRACTION = 1.0 / 3.0
CASES = (
    ("EF0p06_T350", 0.06, 350.0),
    ("EF0p06_T550", 0.06, 550.0),
    ("EF0p11_T350", 0.11, 350.0),
    ("EF0p11_T550", 0.11, 550.0),
)
CASE_FILES = (
    "Const_E.F90", "bolsig_x86_64_g.dll", "bolsig_x86_64_g.lib", "bolsigdb.dat",
    "dvode_f90_m.f90", "zdplaskin_m.F90",
)


def _write_table(path: Path, header: list[str], values: list[float]) -> None:
    path.write_text(" ".join(header) + "\n" + " ".join(f"{v:.16g}" for v in values) + "\n",
                    encoding="utf-8")


def _load_parameter_model():
    with (REFERENCE / "pred_collect.pkl").open("rb") as handle:
        pred_collect = pickle.load(handle)
    with (REFERENCE / "basis.pkl").open("rb") as handle:
        basis = np.asarray(pickle.load(handle), dtype=float).reshape(-1)
    with (REFERENCE / "entropyPara.pkl").open("rb") as handle:
        entropy = np.asarray(pickle.load(handle), dtype=float)
    with (REFERENCE / "EleFieldforEntropy.pkl").open("rb") as handle:
        entropy_field = np.asarray(pickle.load(handle), dtype=float).reshape(-1)
    return pred_collect, basis, entropy, entropy_field


def _dft_parameters(field_v_per_angstrom: float):
    """Reproduce the parameter interpolation in the authors' Slurm worker."""
    pred_collect, basis, entropy, entropy_field = _load_parameter_model()
    sample_fields = np.asarray([-0.407617857, 0.0, 0.405267196, 0.814412714, 0.081528])
    reaction_in, reaction_basis = [], []
    for index, values in enumerate(pred_collect):
        if index < 12:
            continue
        polynomial = np.poly1d(np.polyfit(sample_fields, np.asarray(values, dtype=float)[1:], 3))
        reaction_in.append(float(polynomial(-field_v_per_angstrom * 10.0)))
        reaction_basis.append(float(polynomial(basis[3])))
    entropy_in, entropy_basis = [], []
    for values in entropy:
        polynomial = np.poly1d(np.polyfit(entropy_field, values, 3))
        entropy_in.append(float(polynomial(field_v_per_angstrom * 10.0)))
        entropy_basis.append(float(polynomial(0.81528)))
    if len(reaction_in) != 4 or len(entropy_in) != 21:
        raise RuntimeError("作者参数插值维度异常，拒绝生成算例")
    return reaction_in, reaction_basis, entropy_in, entropy_basis


def prepare_case(case_name: str, field: float, temperature: float) -> Path:
    case = CASES_DIR / case_name
    if case.exists():
        raise RuntimeError(f"算例目录已存在，拒绝覆盖：{case}")
    case.mkdir(parents=True)
    for name in CASE_FILES:
        source = CASE_TEMPLATE / name
        if not source.is_file():
            raise FileNotFoundError(f"模板缺少：{source}")
        shutil.copy2(source, case / name)
    main_source = case / "Const_E.F90"
    source_text = main_source.read_text(encoding="utf-8")
    source_text, replacements = re.subn(r"time_end\s*=\s*9\.0d3", "time_end = 4.0d3", source_text)
    if replacements != 1:
        raise RuntimeError("未能唯一替换作者主程序的 time_end=9.0d3")
    marker = "  enddo\n  !\n  ! for quasineutrality"
    if source_text.count(marker) != 1:
        raise RuntimeError("未能在主积分循环后插入 4000 s 最终采样")
    source_text = source_text.replace(
        marker,
        "  enddo\n  call ZDPlasKin_write_qtplaskin(time, LFORCE_WRITE=.true.)\n  !\n  ! for quasineutrality",
        1,
    )
    main_source.write_text(source_text, encoding="utf-8")
    reaction_in, reaction_basis, entropy_in, entropy_basis = _dft_parameters(field)
    electric_field_v_cm = field * 1.0e7
    _write_table(case / "Ele.dat", ["EField", "Electrons_cm-3"],
                 [electric_field_v_cm, ELECTRON_DENSITY])
    _write_table(case / "Tgas.dat", ["K"], [temperature])
    _write_table(case / "other_para.dat", ["N2_frac"], [N2_FRACTION])
    _write_table(case / "REACTION_E_IN.DAT",
                 ["H2_NHs_ace", "Ns_Hs_ACT_E", "NHs_Hs_ACT_E", "NH2s_Hs_ACT_E"], reaction_in)
    _write_table(case / "REACTION_E_BASIS.DAT",
                 ["H2_NHs_ace", "Ns_Hs_ACT_E", "NHs_Hs_ACT_E", "NH2s_Hs_ACT_E"], reaction_basis)
    entropy_names = [
        "H2_vib_1", "N2_vib_1", "NH_vib_1", "NH2_vib_1", "NH2_vib_2", "NH2_vib_3",
        "NH3_vib_1", "NH3_vib_2", "NH3_vib_3", "NH3_vib_4", "NH3_vib_5", "NH3_vib_6",
        "H2_inertia_1", "N2_inertia_1", "NH_inertia_1", "NH2_inertia_1", "NH2_inertia_2",
        "NH2_inertia_3", "NH3_inertia_1", "NH3_inertia_2", "NH3_inertia_3",
    ]
    _write_table(case / "ENTROPY_PARA_IN.DAT", entropy_names, entropy_in)
    _write_table(case / "ENTROPY_INFO_BASIS.DAT",
                 [name.replace("_1", "_basis_1") for name in entropy_names], entropy_basis)
    (case / "case.json").write_text(json.dumps({
        "case": case_name, "field_v_per_angstrom": field, "temperature_k": temperature,
        "electron_density_cm3": ELECTRON_DENSITY, "n2_fraction": N2_FRACTION,
        "target_time_s": TIMEPOINT_S,
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tool_dir = SHAO_ROOT.parent / "工具脚本"
    sys.path.insert(0, str(tool_dir))
    import zdp_runtime  # pylint: disable=import-outside-toplevel
    if zdp_runtime.build(case, mode="full", log_cb=None) != 0:
        raise RuntimeError(f"派生的 4000 s 算例构建失败：{case}")
    return case


def run_case(case: Path, timeout_s: int) -> dict:
    started = time.monotonic()
    log_path = case / "console.log"
    with log_path.open("w", encoding="utf-8", newline="") as log:
        try:
            result = subprocess.run([str(case / "Const_E.exe")], cwd=case, stdout=log,
                                    stderr=subprocess.STDOUT, text=True, timeout=timeout_s)
            return {"case": case.name, "returncode": result.returncode,
                    "elapsed_s": time.monotonic() - started, "timeout": False}
        except subprocess.TimeoutExpired:
            return {"case": case.name, "returncode": None,
                    "elapsed_s": time.monotonic() - started, "timeout": True}


def _load_qtplaskin():
    package_root = SHAO_ROOT.parent / "qtplaskin-master"
    if not package_root.is_dir():
        raise FileNotFoundError(f"找不到已迁入的 qtplaskin：{package_root}")
    sys.path.insert(0, str(package_root))
    from qtplaskin import FastDirData  # pylint: disable=import-outside-toplevel
    return FastDirData


def _value_at(data, series, timepoint: float) -> float:
    if len(data.t) == 0 or data.t[-1] < timepoint:
        return float("nan")
    return float(np.interp(timepoint, data.t, series))


def metrics_from_dir(directory: Path) -> dict:
    FastDirData = _load_qtplaskin()
    data = FastDirData(str(directory))
    result = {"completed_time_s": float(data.t[-1]), "n_time_points": len(data.t)}
    for species in ("NH3", "N2", "H2", "E"):
        result[species] = _value_at(data, data.get_spec(species), TIMEPOINT_S)
    result["power_w_cm3"] = _value_at(data, data.get_cond("Power density [W/cm3]"), TIMEPOINT_S)
    return result


def author_metrics(field: float, temperature: float, kind: str) -> dict:
    FastDirData = _load_qtplaskin()
    name = f"{kind}_{field:.2f}_{int(temperature)}.pkl"
    with (REFERENCE / name).open("rb") as handle:
        data = pickle.load(handle)
    result = {"completed_time_s": float(data.t[-1]), "n_time_points": len(data.t)}
    for species in ("NH3", "N2", "H2", "E"):
        result[species] = _value_at(data, data.get_spec(species), TIMEPOINT_S)
    result["power_w_cm3"] = _value_at(data, data.get_cond("Power density [W/cm3]"), TIMEPOINT_S)
    return result


def write_comparison(prepared: dict[str, tuple[float, float]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fields = ("NH3", "N2", "H2", "E", "power_w_cm3")
    rows = []
    for name, (field, temperature) in prepared.items():
        base = author_metrics(field, temperature, "base")
        reference = author_metrics(field, temperature, "extract")
        rerun = metrics_from_dir(CASES_DIR / name)
        for metric in fields:
            ref_value, rerun_value = reference[metric], rerun[metric]
            rel_error = (rerun_value - ref_value) / ref_value if np.isfinite(ref_value) and ref_value else float("nan")
            rows.append({"case": name, "field_v_per_angstrom": field, "temperature_k": temperature,
                         "metric": metric, "author_base": base[metric], "author_dft": ref_value,
                         "rerun_dft": rerun_value, "relative_error_to_author_dft": rel_error,
                         "rerun_completed_time_s": rerun["completed_time_s"]})
    with (RESULTS_DIR / "Shao2024复现对比.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    try:
        import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel
        nh3 = [row for row in rows if row["metric"] == "NH3"]
        x = np.arange(len(nh3))
        fig, ax = plt.subplots(figsize=(10, 5), dpi=160)
        ax.bar(x - 0.25, [row["author_base"] for row in nh3], 0.25, label="Author base")
        ax.bar(x, [row["author_dft"] for row in nh3], 0.25, label="Author DFT")
        ax.bar(x + 0.25, [row["rerun_dft"] for row in nh3], 0.25, label="Local DFT rerun")
        ax.set_yscale("log")
        ax.set_ylabel("NH3 density [cm$^{-3}$] at 4000 s")
        ax.set_xticks(x, [row["case"] for row in nh3])
        ax.legend()
        fig.tight_layout()
        fig.savefig(RESULTS_DIR / "NH3_4000s_对比.png")
        plt.close(fig)
    except Exception as exc:  # noqa: BLE001
        (RESULTS_DIR / "绘图警告.txt").write_text(str(exc) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Shao 2024 四工况 DFT 微观动力学复现")
    parser.add_argument("--prepare-only", action="store_true", help="仅生成四个算例目录")
    parser.add_argument("--workers", type=int, default=2, help="并行运行数（默认 2）")
    parser.add_argument("--timeout-s", type=int, default=3600, help="单算例最长墙钟时间（默认 3600 s）")
    args = parser.parse_args(argv)
    if CASES_DIR.exists() and any(CASES_DIR.iterdir()):
        raise RuntimeError(f"算例目录非空，拒绝覆盖：{CASES_DIR}")
    CASES_DIR.mkdir(parents=True, exist_ok=True)
    prepared = {name: (field, temperature) for name, field, temperature in CASES}
    for name, field, temperature in CASES:
        prepare_case(name, field, temperature)
    (HERE / "prepared_cases.json").write_text(json.dumps(prepared, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if args.prepare_only:
        print("已生成 4 个 Shao 复现算例；未运行。")
        return 0
    outcomes = []
    with ThreadPoolExecutor(max_workers=max(1, min(args.workers, len(CASES)))) as executor:
        futures = [executor.submit(run_case, CASES_DIR / name, args.timeout_s) for name in prepared]
        for future in as_completed(futures):
            outcome = future.result()
            outcomes.append(outcome)
            print(json.dumps(outcome, ensure_ascii=False))
    (HERE / "run_outcomes.json").write_text(json.dumps(outcomes, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    failures = [item for item in outcomes if item["timeout"] or item["returncode"] != 0]
    if failures:
        print("存在失败或超时算例；保留日志，不生成误导性的比较结果。", file=sys.stderr)
        return 1
    write_comparison(prepared)
    print(f"复现完成：{RESULTS_DIR / 'Shao2024复现对比.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
